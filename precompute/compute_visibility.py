#!/usr/bin/env python3
"""
Eclipse 2026 Terrain Visibility Precomputation
===============================================
For each cell in the SRTM 90m DEM covering Spain, determines whether
the sun is terrain-blocked at the moment of the August 12 2026 solar eclipse.

Outputs two GeoJSON files:
  - terrain_blocked.geojson   — areas where mountains obstruct the sun
  - terrain_clear.geojson     — areas with unobstructed view

Usage:
    pip install -r requirements.txt
    python compute_visibility.py
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
import rasterio
import rasterio.warp
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.features import shapes
from rasterio.merge import merge
from rasterio.transform import from_bounds
from shapely.geometry import shape, mapping

# ── Configuration ─────────────────────────────────────────────────────────────

# Bounding box for Spain + small margin so rays don't fall off the edge
# (west, south, east, north)
BBOX = (-10.0, 34.5, 5.5, 44.5)

# Eclipse time: August 12 2026, approximately 18:30 UTC (20:30 CEST)
# Maximum totality over Spain is ~20:29:59 CEST.
ECLIPSE_UTC = (2026, 8, 12, 18, 30, 0)

# How far ahead (in the sun's direction) to scan for blocking terrain.
# At sun altitude ~25°, a 2500m peak 10km away would just barely block you.
# 30km covers all realistic cases.
MAX_RAY_DIST_M = 30_000

# Polygon simplification tolerance in degrees (~200m at Spain's latitude)
SIMPLIFY_TOL = 0.002

# Minimum polygon area to keep (degrees²). Filters out tiny isolated pixels.
MIN_AREA_DEG2 = 0.0001  # roughly 1 km²

DEM_PATH = Path("spain_dem.tif")
BLOCKED_TIF = Path("terrain_blocked.tif")
BLOCKED_GEOJSON = Path("../data/terrain_blocked.geojson")
CLEAR_GEOJSON = Path("../data/terrain_clear.geojson")

SRTM_DIR = Path("../data/srtm")   # .hgt.gz files downloaded by download_srtm.sh
RESAMPLE_FACTOR = 3               # 30m → ~90m

# ── Step 1: Build merged DEM from local SRTM tiles ───────────────────────────

def download_dem():
    if DEM_PATH.exists():
        print(f"DEM already exists at {DEM_PATH}, skipping.")
        return

    tile_files = sorted(SRTM_DIR.glob("*.hgt.gz"))
    if not tile_files:
        sys.exit(f"No .hgt.gz files found in {SRTM_DIR}. Run download_srtm.sh first.")

    print(f"Opening {len(tile_files)} local SRTM tiles from {SRTM_DIR} …")
    datasets = []
    for i, f in enumerate(tile_files, 1):
        vsi_path = f"/vsigzip/{f.resolve()}"
        try:
            ds = rasterio.open(vsi_path)
            datasets.append(ds)
        except rasterio.errors.RasterioIOError as e:
            print(f"  WARNING: could not open {f.name}: {e}")
        print(f"  {i}/{len(tile_files)} opened", end="\r", flush=True)

    print(f"\n  Merging {len(datasets)} tiles …")
    mosaic, transform = merge(datasets, method="first", nodata=-32768)
    for ds in datasets:
        ds.close()

    src_h, src_w = mosaic.shape[1], mosaic.shape[2]
    dst_h = round(src_h / RESAMPLE_FACTOR)
    dst_w = round(src_w / RESAMPLE_FACTOR)
    crs = CRS.from_epsg(4326)

    west  = transform.c
    north = transform.f
    east  = west  + transform.a * src_w
    south = north + transform.e * src_h
    dst_transform = from_bounds(west, south, east, north, dst_w, dst_h)

    print(f"  Resampling {src_h}×{src_w} → {dst_h}×{dst_w} (~90m) …")
    resampled = np.empty((1, dst_h, dst_w), dtype=np.float32)
    mem_profile = dict(driver="MEM", count=1, dtype="float32", crs=crs,
                       transform=transform, width=src_w, height=src_h, nodata=-32768)
    with rasterio.MemoryFile() as mem:
        with mem.open(**mem_profile) as src_ds:
            src_ds.write(mosaic.astype(np.float32))
            rasterio.warp.reproject(
                source=rasterio.band(src_ds, 1),
                destination=resampled,
                dst_transform=dst_transform,
                dst_crs=crs,
                resampling=Resampling.average,
            )

    profile = dict(driver="GTiff", count=1, dtype="float32", crs=crs,
                   transform=dst_transform, width=dst_w, height=dst_h,
                   compress="lzw", nodata=-32768)
    with rasterio.open(DEM_PATH, "w", **profile) as dst:
        dst.write(resampled)

    size_mb = DEM_PATH.stat().st_size / 1e6
    print(f"  Saved {DEM_PATH}  ({dst_h}×{dst_w} @ ~90m,  {size_mb:.0f} MB)")


# ── Step 2: Sun position at eclipse time ──────────────────────────────────────

def get_sun_position(lat, lon):
    """
    Returns (azimuth_deg, altitude_deg) for the centre of Spain at eclipse time.
    Pysolar convention: azimuth is degrees clockwise from North.
    """
    from datetime import datetime, timezone
    from pysolar.solar import get_altitude, get_azimuth

    dt = datetime(*ECLIPSE_UTC, tzinfo=timezone.utc)
    alt = get_altitude(lat, lon, dt)
    az  = get_azimuth(lat, lon, dt)
    return az, alt


# ── Step 3: Ray-cast every DEM cell ──────────────────────────────────────────

def _build_numba_kernel():
    """Compile and return the numba ray-casting kernel (fast path)."""
    import numba

    @numba.njit(parallel=True)
    def _raycast(dem, dr, dc, tan_sun_alt, cell_size_m, max_steps):
        rows, cols = dem.shape
        blocked = np.zeros((rows, cols), dtype=numba.boolean)

        for r in numba.prange(rows):
            for c in range(cols):
                obs_h = dem[r, c]
                if obs_h < -200:          # ocean / nodata — skip
                    continue

                for k in range(1, max_steps + 1):
                    sr = r + k * dr
                    sc = c + k * dc

                    ir = int(sr)
                    ic = int(sc)
                    if ir < 0 or ir >= rows - 1 or ic < 0 or ic >= cols - 1:
                        break           # ray left the DEM — no blocker found

                    # Bilinear interpolation
                    fr = sr - ir
                    fc = sc - ic
                    h = (dem[ir,     ic    ] * (1 - fr) * (1 - fc) +
                         dem[ir + 1, ic    ] * fr        * (1 - fc) +
                         dem[ir,     ic + 1] * (1 - fr) * fc        +
                         dem[ir + 1, ic + 1] * fr        * fc)

                    if (h - obs_h) / (k * cell_size_m) > tan_sun_alt:
                        blocked[r, c] = True
                        break

        return blocked

    return _raycast


def _raycast_numpy(dem, dr, dc, tan_sun_alt, cell_size_m, max_steps):
    """
    Pure-numpy fallback. Vectorises across all cells for each step k,
    skipping cells already determined to be blocked.
    ~10–30× slower than numba but correct.
    """
    from scipy.ndimage import map_coordinates

    rows, cols = dem.shape
    blocked  = np.zeros((rows, cols), dtype=bool)
    settled  = np.zeros((rows, cols), dtype=bool)   # blocked OR ray left DEM

    col_idx = np.arange(cols, dtype=np.float32)
    row_idx = np.arange(rows, dtype=np.float32)
    col_grid, row_grid = np.meshgrid(col_idx, row_idx)

    print(f"    numpy fallback: {max_steps} ray steps …")
    for k in range(1, max_steps + 1):
        if k % 50 == 0:
            pct = 100 * settled.mean()
            print(f"    step {k}/{max_steps}  ({pct:.0f}% cells settled)")

        active = ~settled
        if not active.any():
            break

        sr = (row_grid + k * dr)[active]
        sc = (col_grid + k * dc)[active]

        in_bounds = (sr >= 0) & (sr < rows - 1) & (sc >= 0) & (sc < cols - 1)

        # Cells whose ray has left the DEM are settled (not blocked)
        out_mask = np.zeros((rows, cols), dtype=bool)
        out_mask[active] = ~in_bounds
        settled |= out_mask

        if not in_bounds.any():
            continue

        # Only interpolate the in-bounds subset
        active_ib = active.copy()
        active_ib[active] = in_bounds

        coords = np.array([sr[in_bounds], sc[in_bounds]])
        sampled_h = map_coordinates(dem, coords, order=1, mode="nearest")

        obs_h     = dem[active_ib]
        dist_m    = k * cell_size_m
        new_block = (sampled_h - obs_h) / dist_m > tan_sun_alt

        newly_blocked = np.zeros((rows, cols), dtype=bool)
        newly_blocked[active_ib] = new_block

        blocked  |= newly_blocked
        settled  |= newly_blocked

    return blocked


def compute_blocking(dem, cell_size_m, sun_az_deg, sun_alt_deg):
    az_rad  = math.radians(sun_az_deg)
    # In raster coords: rows increase southward, cols increase eastward
    dr = -math.cos(az_rad)   # north component  → decreasing row
    dc =  math.sin(az_rad)   # east component   → increasing col

    tan_sun_alt = math.tan(math.radians(sun_alt_deg))
    max_steps   = int(MAX_RAY_DIST_M / cell_size_m)

    print(f"  Sun: azimuth={sun_az_deg:.1f}°  altitude={sun_alt_deg:.1f}°")
    print(f"  Ray: dr={dr:.4f}  dc={dc:.4f}  max_steps={max_steps}")
    print(f"  DEM: {dem.shape[0]}×{dem.shape[1]} = {dem.size:,} cells")

    try:
        kernel = _build_numba_kernel()
        print("  Using numba (parallel) …")
        blocked = kernel(
            dem.astype(np.float32),
            np.float32(dr), np.float32(dc),
            np.float32(tan_sun_alt),
            np.float32(cell_size_m),
            max_steps,
        )
    except ImportError:
        print("  numba not installed — using numpy fallback …")
        blocked = _raycast_numpy(
            dem.astype(np.float32),
            dr, dc, tan_sun_alt, cell_size_m, max_steps,
        )

    pct = 100 * blocked.mean()
    print(f"  Result: {blocked.sum():,} / {blocked.size:,} cells blocked ({pct:.1f}%)")
    return blocked


# ── Step 4: Vectorise to GeoJSON ──────────────────────────────────────────────

def vectorise(blocked_arr, transform, crs, output_path, export_blocked: bool):
    """
    Convert a boolean raster to simplified GeoJSON polygons.
    export_blocked=True  → polygons where terrain blocks the sun
    export_blocked=False → polygons where the sun is visible
    """
    from scipy.ndimage import binary_closing

    label = "blocked" if export_blocked else "clear"
    print(f"  Vectorising {label} areas → {output_path}")

    mask = blocked_arr if export_blocked else ~blocked_arr

    # Small morphological closing removes single-cell speckle
    mask = binary_closing(mask, structure=np.ones((3, 3))).astype(np.uint8)

    polys = []
    for geom_dict, val in shapes(mask, mask=mask, transform=transform):
        if val == 0:
            continue
        poly = shape(geom_dict).simplify(SIMPLIFY_TOL, preserve_topology=True)
        if poly.is_empty or poly.area < MIN_AREA_DEG2:
            continue
        polys.append(poly)

    print(f"    {len(polys):,} polygons after simplification")

    geojson = {
        "type": "FeatureCollection",
        "features": [
            {"type": "Feature",
             "geometry": mapping(p),
             "properties": {"terrain_blocked": export_blocked}}
            for p in polys
        ],
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(geojson, f)

    size_mb = output_path.stat().st_size / 1e6
    print(f"    Written {output_path}  ({size_mb:.1f} MB)")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("\n=== Step 1: DEM ===")
    download_dem()

    print("\n=== Step 2: Load DEM ===")
    with rasterio.open(DEM_PATH) as src:
        dem       = src.read(1).astype(np.float32)
        transform = src.transform
        crs       = src.crs
        # Degrees per pixel → metres per pixel (approximate, mid-Spain latitude)
        deg_per_px   = abs(transform.a)
        cell_size_m  = deg_per_px * 111_320 * math.cos(math.radians(40.0))
    print(f"  Shape: {dem.shape}   CRS: {crs}   Cell ≈ {cell_size_m:.0f} m")

    # Clamp nodata / nonsense values to sea level
    dem = np.clip(dem, 0, 8_850)

    print("\n=== Step 3: Sun position ===")
    centre_lat = (BBOX[1] + BBOX[3]) / 2
    centre_lon = (BBOX[0] + BBOX[2]) / 2
    sun_az, sun_alt = get_sun_position(centre_lat, centre_lon)
    print(f"  Centre ({centre_lat:.1f}°N, {centre_lon:.1f}°E): "
          f"azimuth={sun_az:.1f}°  altitude={sun_alt:.1f}°")

    if sun_alt <= 0:
        sys.exit("ERROR: Sun is below the horizon at the given time — check ECLIPSE_UTC.")

    print("\n=== Step 4: Ray casting ===")
    blocked = compute_blocking(dem, cell_size_m, sun_az, sun_alt)

    # Save intermediate GeoTIFF (useful for debugging in QGIS etc.)
    with rasterio.open(
        BLOCKED_TIF, "w",
        driver="GTiff", height=blocked.shape[0], width=blocked.shape[1],
        count=1, dtype="uint8", crs=crs, transform=transform,
        compress="lzw",
    ) as dst:
        dst.write(blocked.astype(np.uint8), 1)
    print(f"  Saved intermediate raster → {BLOCKED_TIF}")

    print("\n=== Step 5: Vectorise ===")
    vectorise(blocked, transform, crs, BLOCKED_GEOJSON, export_blocked=True)
    vectorise(blocked, transform, crs, CLEAR_GEOJSON,   export_blocked=False)

    print("\n=== Done ===")
    print(f"  {BLOCKED_GEOJSON}  — terrain-blocked areas (overlay these as 'no view')")
    print(f"  {CLEAR_GEOJSON}    — terrain-clear areas (overlay as 'visible')")
    print()
    print("Next: download the 2026 eclipse path GeoJSON from")
    print("  https://svs.gsfc.nasa.gov/5073  (NASA eclipse visualisation)")
    print("  or Xavier Jubier's tool: http://xjubier.free.fr/en/site_pages/solar_eclipses/TSE_2026_GoogleMapFull.html")
    print("Place the umbra/penumbra GeoJSON in data/ and update index.html.")


if __name__ == "__main__":
    main()
