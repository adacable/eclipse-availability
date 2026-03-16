// Paths relative to web/ — adjust if serving from repo root
const DATA = {
  // NASA / Xavier Jubier eclipse path polygons
  // Expected properties: "type" = "umbra" | "penumbra" (or adapt STYLE_FOR below)
  eclipsePath:    "../data/eclipse_path_2026.geojson",

  // Output of precompute/compute_visibility.py
  terrainBlocked: "../data/terrain_blocked_clipped.geojson",
};

// ── Map init ─────────────────────────────────────────────────────────────────

const map = L.map("map", {
  center: [40.4, -3.7],   // Madrid
  zoom: 6,
  zoomControl: true,
});

L.tileLayer("https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png", {
  attribution: '© <a href="https://www.openstreetmap.org/copyright">OSM</a> © <a href="https://carto.com/">CARTO</a>',
  subdomains: "abcd",
  maxZoom: 19,
}).addTo(map);

// ── Styles ───────────────────────────────────────────────────────────────────

const STYLES = {
  umbra: {
    color:       "#e0c060",
    weight:      2,
    fillColor:   "#1a1a2e",
    fillOpacity: 0.75,
  },
  penumbra: {
    color:       "#c8a028",
    weight:      1,
    fillColor:   "#c8a028",
    fillOpacity: 0.20,
  },
  blocked: {
    color:       "#b43c3c",
    weight:      0,
    fillColor:   "#b43c3c",
    fillOpacity: 0.45,
  },
};

// ── Layer loading ─────────────────────────────────────────────────────────────

/**
 * Fetch a GeoJSON file and add it to the map.
 * styleFn: (feature) => L.PathOptions, or a plain options object.
 * Returns the layer, or null if the file is missing.
 */
async function addGeoJsonLayer(url, styleFn, options = {}) {
  let data;
  try {
    const res = await fetch(url);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    data = await res.json();
  } catch (err) {
    console.warn(`Could not load ${url}:`, err.message);
    return null;
  }

  const layer = L.geoJSON(data, {
    style: typeof styleFn === "function" ? styleFn : () => styleFn,
    ...options,
  }).addTo(map);

  return layer;
}

function eclipseStyle(feature) {
  const layer = feature.properties?.layer ?? "";
  switch (layer) {
    case "umbra":        return STYLES.umbra;
    case "north_limit":
    case "south_limit":
    case "central":      return { color: "#e0c060", weight: 1.5, fill: false };
    case "penumbra_south": return { color: "#c8a028", weight: 1, dashArray: "4 4", fill: false };
    case "mag_0.8":      return { color: "#c8a028", weight: 1, opacity: 0.7, fill: false };
    case "mag_0.6":      return { color: "#c8a028", weight: 1, opacity: 0.5, fill: false };
    case "mag_0.4":      return { color: "#c8a028", weight: 1, opacity: 0.3, fill: false };
    default:             return { color: "#c8a028", weight: 1, fill: false };
  }
}

// ── Boot ──────────────────────────────────────────────────────────────────────

(async () => {
  // Draw terrain-blocked areas first (bottom layer)
  await addGeoJsonLayer(DATA.terrainBlocked, STYLES.blocked);

  // Eclipse path on top
  await addGeoJsonLayer(DATA.eclipsePath, eclipseStyle);
})();
