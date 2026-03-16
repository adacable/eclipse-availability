#!/usr/bin/env bash
# Download SRTM1 tiles for Spain from Mapzen/AWS public bucket.
# Tiles are 1°×1° HGT files at 1 arc-second (~30m) resolution.
# Resume-safe: wget -c skips files that are already complete.

set -euo pipefail

DEST="../data/srtm"
BASE="https://s3.amazonaws.com/elevation-tiles-prod/skadi"

mkdir -p "$DEST"

# Spain bounding box: lon -10→5, lat 34→44
LATS=(34 35 36 37 38 39 40 41 42 43 44)
LONS=(-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5)

total=$(( ${#LATS[@]} * ${#LONS[@]} ))
done=0

for lat in "${LATS[@]}"; do
    ns="N"
    lat_fmt=$(printf "%02d" $lat)

    for lon in "${LONS[@]}"; do
        done=$(( done + 1 ))

        if (( lon < 0 )); then
            ew="W"
            lon_fmt=$(printf "%03d" $(( -lon )))
        else
            ew="E"
            lon_fmt=$(printf "%03d" $lon)
        fi

        name="${ns}${lat_fmt}${ew}${lon_fmt}"
        url="${BASE}/${ns}${lat_fmt}/${name}.hgt.gz"
        out="${DEST}/${name}.hgt.gz"

        echo "[${done}/${total}] ${name}  →  ${out}"

        # Skip if already fully downloaded (check size > 0)
        if [[ -s "${out}" ]]; then
            echo "  already exists, skipping"
            continue
        fi

        http_code=$(curl \
            --silent \
            --show-error \
            --retry 5 \
            --retry-delay 3 \
            --connect-timeout 10 \
            --max-time 120 \
            --write-out "%{http_code}" \
            --output "${out}" \
            "${url}")

        if [[ "${http_code}" == "404" ]]; then
            rm -f "${out}"
            echo "  404 (ocean/void tile) — skipped"
        elif [[ "${http_code}" != "200" ]]; then
            rm -f "${out}"
            echo "  ERROR: HTTP ${http_code} for ${url}"
        else
            size=$(wc -c < "${out}" | tr -d ' ')
            echo "  OK  (${size} bytes)"
        fi
    done
done

echo ""
echo "Done. $(ls "$DEST"/*.hgt.gz 2>/dev/null | wc -l | tr -d ' ') tiles in ${DEST}/"
