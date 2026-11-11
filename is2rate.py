#!/usr/bin/env python3
"""is2rate: rates of vertical land motion from ICESat-2 on a hexagon grid.

ICESat-2 ground segments (SlideRule ATL06-SR) are differenced from a reference DEM. For each hexagon and pass,
the median residual of the target land cover in the hexagon is differenced from the median residual of
reference land cover within +-ref_halfwidth_m on the same pass, and a least-squares trend (with a seasonal
cycle and one constant per DEM) is fitted to these per-pass differences. Output: a 2-band GeoTIFF
(band 1 rate, band 2 SE, mm/yr), a QGIS style, hexagons (GeoPackage) and a quick-look PNG.

Usage:
    python is2rate.py --init params.txt        # write a parameter file with defaults, then edit it
    python is2rate.py params.txt               # run
    python is2rate.py params.txt --set hex_width_m=5000 --set end=2022-06-30
"""
import argparse, glob, json, math, os, sys, time

import numpy as np
import pandas as pd

__version__ = "0.1.0"

# (key, default, comment). Defaults are those used for the Waikato (Harvey, in prep.).
PARAMS = [
    ("# --- area and time", None, None),
    ("bbox", "", "lon_min, lat_min, lon_max, lat_max (WGS84 degrees), e.g. 175.2, -37.95, 175.7, -37.6. Required"),
    ("start", "2018-10-13", "first pass (YYYY-MM-DD, or YYYY = 1 January)"),
    ("end", "2026-12-31", "last pass (YYYY-MM-DD, or YYYY = 31 December)"),
    ("out_dir", "is2rate_out", "output directory"),
    ("name", "is2rate", "prefix of output files"),
    ("crs", "auto", "working and output CRS: auto (UTM zone of the bbox centre) or an EPSG code, e.g. 2193"),
    ("# --- reference surface", None, None),
    ("dem", "copernicus", "copernicus (GLO-30, read from AWS), DEM GeoTIFF/VRT path(s) comma-separated in priority order, or none (repeat pairs, no DEM)"),
    ("pair_radius_m", "10", "dem = none: segments from different passes within this distance form a repeat pair"),
    ("dem_geoid", "egm2008", "vertical datum of the DEM(s): egm2008, ellipsoid, or a geoid-height GeoTIFF (N; H = h - N)"),
    ("# --- land cover and reference", None, None),
    ("landcover", "worldcover", "worldcover (ESA WorldCover 2021 10 m, read from AWS) or a GeoTIFF with WorldCover class codes"),
    ("target_classes", "30, 40", "land-cover classes measured in each hexagon (WorldCover 30 grassland, 40 cropland)"),
    ("ref_classes", "30, 40", "land-cover classes used as the reference"),
    ("min_target_share", "0.25", "hexagons with less than this fraction of their area in target classes are not reported"),
    ("exclude", "none", "vector (gpkg, shp, geojson) or raster (nonzero) of areas not used as reference, e.g. peat soils"),
    ("exclude_buffer_m", "500", "reference segments must lie at least this far from 'exclude'"),
    ("reference", "local", "local: relative to reference classes on the same passes; none: absolute rates"),
    ("ref_halfwidth_m", "10000", "reference window, +- this distance from the hexagon centre (hexagon itself excluded)"),
    ("# --- hexagons and output grid", None, None),
    ("hex_width_m", "3000", "hexagon width, flat to flat (m)"),
    ("cell_shape", "hex", "hex (hexagons) or square (square cells of width hex_width_m)"),
    ("hex_orientation", "flat", "flat (flat top, as the published maps) or pointy (pointy top)"),
    ("grid_origin", "auto", "x, y of the first hexagon centre in the working CRS (auto: lower-left corner of the bbox)"),
    ("out_res_m", "100", "pixel size of the output GeoTIFF (m)"),
    ("# --- ICESat-2 segments (SlideRule ATL06-SR)", None, None),
    ("segments_file", "none", "existing ATL06-SR GeoParquet file(s), glob allowed; none: request from SlideRule (cached in out_dir/cache)"),
    ("sliderule_url", "slideruleearth.io", "SlideRule server"),
    ("seg_len_m", "20", "segment length (m)"),
    ("seg_res_m", "10", "segment spacing (m)"),
    ("atl08_ground", "yes", "use only photons classified as ground by ATL08"),
    ("tile_deg", "0.5", "size of SlideRule request tiles (degrees)"),
    ("# --- segment quality", None, None),
    ("min_photons", "10", "minimum photons per segment"),
    ("max_h_sigma_m", "0.15", "maximum height uncertainty (m)"),
    ("max_along_slope", "0.03", "maximum along-track slope (m/m)"),
    ("max_window_m", "4", "maximum photon window (m)"),
    ("max_terrain_slope_deg", "2", "maximum DEM slope over +-10 m (degrees)"),
    ("max_resid_dev_m", "0.5", "drop segments whose residual is further than this from the pass median (m)"),
    ("# --- trend fit", None, None),
    ("min_segs_per_pass", "10", "minimum target and reference segments on a pass"),
    ("min_passes", "8", "minimum passes for a rate"),
    ("seasonal", "yes", "fit an annual sinusoid"),
    ("se_method", "bootstrap", "band 2: bootstrap (whole passes resampled) or formula (least-squares SE)"),
    ("n_boot", "300", "bootstrap resamples"),
    ("seed", "0", "random seed"),
    ("# --- outputs", None, None),
    ("write_gpkg", "yes", "write hexagons as GeoPackage"),
    ("write_png", "yes", "write a quick-look PNG"),
]

# Colour classes of the published maps (mm/yr): upper bound, colour, label
BREAKS = [(-40, "#08306b", "≤ −40"), (-30, "#1f5aa6", "−40 to −30"), (-20, "#2a78d6", "−30 to −20"), (-10, "#6fa6e8", "−20 to −10"),
          (-5, "#b5d2f4", "−10 to −5"), (5, "#efeee9", "−5 to +5"), (10, "#f5c4c3", "+5 to +10"), (20, "#e8807f", "+10 to +20"),
          (np.inf, "#b22b2b", "> +20")]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ---------------------------------------------------------------- parameters
def write_template(path):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# is2rate {__version__} parameter file. Lines: key = value  # comment\n")
        for k, v, c in PARAMS:
            if v is None: f.write(f"\n{k}\n"); continue
            f.write(f"{k + ' = ' + v:<34} # {c}\n")
    print("wrote", path)


def read_params(path, overrides=()):
    p = {k: v for k, v, _ in PARAMS if v is not None}
    lines = open(path, encoding="utf-8").read().splitlines() + list(overrides)
    for ln in lines:
        ln = ln.split("#", 1)[0].strip()
        if not ln: continue
        if "=" not in ln: sys.exit(f"bad parameter line: {ln}")
        k, v = [s.strip() for s in ln.split("=", 1)]
        if k not in p: sys.exit(f"unknown parameter: {k}")
        p[k] = v
    yes = lambda s: s.lower() in ("yes", "y", "true", "1")
    lst = lambda s: [float(t) for t in s.replace(";", ",").split(",") if t.strip()]
    q = dict(p)
    if not p["bbox"]: sys.exit("bbox is required")
    q["bbox"] = lst(p["bbox"])
    if len(q["bbox"]) != 4 or q["bbox"][0] >= q["bbox"][2] or q["bbox"][1] >= q["bbox"][3]: sys.exit("bbox must be lon_min, lat_min, lon_max, lat_max")
    q["start"] = pd.Timestamp(p["start"] + ("-01-01" if len(p["start"]) == 4 else ""))
    q["end"] = pd.Timestamp(p["end"] + ("-12-31" if len(p["end"]) == 4 else "")) + pd.Timedelta(hours=23, minutes=59, seconds=59) \
        if len(p["end"]) <= 10 else pd.Timestamp(p["end"])
    q["target_classes"] = [int(v) for v in lst(p["target_classes"])]
    q["ref_classes"] = [int(v) for v in lst(p["ref_classes"])]
    for k in ("min_target_share", "exclude_buffer_m", "ref_halfwidth_m", "hex_width_m", "out_res_m", "seg_len_m", "seg_res_m", "tile_deg",
              "max_h_sigma_m", "max_along_slope", "max_window_m", "max_terrain_slope_deg", "max_resid_dev_m"):
        q[k] = float(p[k])
    for k in ("min_photons", "min_segs_per_pass", "min_passes", "n_boot", "seed"):
        q[k] = int(p[k])
    for k in ("atl08_ground", "seasonal", "write_gpkg", "write_png"):
        q[k] = yes(p[k])
    q["dem"] = p["dem"].lower() if p["dem"].lower() in ("copernicus", "none") else [s.strip() for s in p["dem"].split(",")]
    q["pair_radius_m"] = float(p["pair_radius_m"])
    q["grid_origin"] = None if p["grid_origin"].lower() == "auto" else lst(p["grid_origin"])
    q["exclude"] = None if p["exclude"].lower() == "none" else p["exclude"]
    q["segments_file"] = None if p["segments_file"].lower() == "none" else p["segments_file"]
    if p["reference"] not in ("local", "none"): sys.exit("reference must be local or none")
    if p["cell_shape"] not in ("hex", "square"): sys.exit("cell_shape must be hex or square")
    if p["hex_orientation"] not in ("flat", "pointy"): sys.exit("hex_orientation must be flat or pointy")
    if p["se_method"] not in ("bootstrap", "formula"): sys.exit("se_method must be bootstrap or formula")
    return q, p


def working_crs(p):
    from pyproj import CRS
    if p["crs"].lower() != "auto": return CRS.from_user_input(int(p["crs"]) if p["crs"].isdigit() else p["crs"])
    lon = (p["bbox"][0] + p["bbox"][2]) / 2; lat = (p["bbox"][1] + p["bbox"][3]) / 2
    return CRS.from_epsg((32600 if lat >= 0 else 32700) + int((lon + 180) // 6) + 1)


# ---------------------------------------------------------------- raster input
def gdal_env():
    import rasterio
    return rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif", GDAL_HTTP_MAX_RETRY="4", GDAL_HTTP_RETRY_DELAY="2")


def fetch_tiles(urls, bbox, out, nodata):
    """Read the bbox window of cloud-optimised GeoTIFF tiles (geographic) into one local GeoTIFF."""
    import rasterio
    from rasterio.merge import merge
    if os.path.exists(out): return out
    pad = 0.02; b = (bbox[0] - pad, bbox[1] - pad, bbox[2] + pad, bbox[3] + pad)
    with gdal_env():
        src = []
        for u in urls:
            try: src.append(rasterio.open(u))
            except Exception: log("  no tile", u.rsplit("/", 1)[-1], "(sea or outside coverage)")
        if not src: sys.exit("no tiles found for the bbox")
        log("  reading", len(src), "tiles")
        a, tr = merge(src, bounds=b, nodata=nodata)
        prof = dict(driver="GTiff", height=a.shape[1], width=a.shape[2], count=1, dtype=a.dtype, crs=src[0].crs, transform=tr, nodata=nodata,
                    compress="deflate", tiled=True)
        with rasterio.open(out, "w", **prof) as d: d.write(a[0], 1)
    return out


def copernicus(bbox, cache):
    urls = []
    for la in range(math.floor(bbox[1]), math.ceil(bbox[3])):
        for lo in range(math.floor(bbox[0]), math.ceil(bbox[2])):
            n = f"Copernicus_DSM_COG_10_{'N' if la >= 0 else 'S'}{abs(la):02d}_00_{'E' if lo >= 0 else 'W'}{abs(lo):03d}_00_DEM"
            urls.append(f"/vsicurl/https://copernicus-dem-30m.s3.amazonaws.com/{n}/{n}.tif")
    log("Copernicus GLO-30")
    return fetch_tiles(urls, bbox, f"{cache}/copernicus_glo30.tif", -32767.0)


def worldcover(bbox, cache):
    urls = []
    for la in range(math.floor(bbox[1] / 3) * 3, math.ceil(bbox[3] / 3) * 3, 3):
        for lo in range(math.floor(bbox[0] / 3) * 3, math.ceil(bbox[2] / 3) * 3, 3):
            t = f"{'N' if la >= 0 else 'S'}{abs(la):02d}{'E' if lo >= 0 else 'W'}{abs(lo):03d}"
            urls.append(f"/vsicurl/https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_{t}_Map.tif")
    log("ESA WorldCover 2021")
    return fetch_tiles(urls, bbox, f"{cache}/worldcover_2021.tif", 0)


def sample(path, x, y, xy_crs, order=1, offsets=((0.0, 0.0),), block=1024):
    """Sample band 1 of a raster at points (x, y) in xy_crs, plus offset points (dx, dy in xy_crs units).
    order 1: bilinear (NaN where any neighbour is nodata); order 0: nearest. Returns array (len(offsets), n)."""
    import rasterio
    from rasterio.windows import Window
    from pyproj import Transformer
    from scipy import ndimage
    out = np.full((len(offsets), len(x)), np.nan)
    with rasterio.open(path) as d:
        tr = Transformer.from_crs(xy_crs, d.crs, always_xy=True); inv = ~d.transform
        C = np.empty((len(offsets), len(x))); R = np.empty_like(C)
        for k, (dx, dy) in enumerate(offsets):
            u, v = tr.transform(x + dx, y + dy); c, r = inv * (np.asarray(u), np.asarray(v)); C[k] = c - 0.5; R[k] = r - 0.5
        ok = np.isfinite(C).all(0) & np.isfinite(R).all(0)
        pad = int(np.ceil(np.nanmax(np.abs(C - C[0])) if len(offsets) > 1 else 0)) + 2
        bc = np.full(len(x), -10**9); br = np.full(len(x), -10**9)
        bc[ok] = np.floor(C[0, ok] / block).astype(int); br[ok] = np.floor(R[0, ok] / block).astype(int)
        inside = ok & (bc >= -1) & (br >= -1) & (bc * block <= d.width) & (br * block <= d.height)
        key = pd.Series(np.arange(len(x))[inside]).groupby([br[inside], bc[inside]])
        nod = d.nodata; fill = nod if nod is not None else (np.nan if d.dtypes[0].startswith("float") else 0)
        for (rb, cb), idx in key:
            idx = idx.values; c0 = cb * block - pad; r0 = rb * block - pad
            n_ = block + 2 * pad; a = np.full((n_, n_), np.nan, np.float32)
            ca, cz = max(c0, 0), min(c0 + n_, d.width); ra, rz = max(r0, 0), min(r0 + n_, d.height)
            if cz <= ca or rz <= ra: continue
            v_ = d.read(1, window=Window(ca, ra, cz - ca, rz - ra)).astype(np.float32)   # plain read: boundless reads of VRTs are very slow
            if nod is not None: v_[v_ == nod] = np.nan
            a[ra - r0:rz - r0, ca - c0:cz - c0] = v_
            for k in range(len(offsets)):
                cc = C[k, idx] - c0; rr = R[k, idx] - r0
                if order == 0:
                    ci = np.clip(np.rint(cc).astype(int), 0, a.shape[1] - 1); ri = np.clip(np.rint(rr).astype(int), 0, a.shape[0] - 1); out[k, idx] = a[ri, ci]
                else:
                    v = ndimage.map_coordinates(np.nan_to_num(a, nan=-1e12), [rr, cc], order=1, mode="constant", cval=-1e12)
                    v[v < -1e11] = np.nan; out[k, idx] = v
    return out


def to_grid(path, G, resampling):
    """Reproject band 1 of a raster to the working grid G."""
    import rasterio
    from rasterio.warp import reproject
    a = np.zeros(G["shape"], np.float32)
    with rasterio.open(path) as d:
        reproject(rasterio.band(d, 1), a, dst_transform=G["transform"], dst_crs=G["crs"], resampling=resampling, dst_nodata=0)
    return a


# ---------------------------------------------------------------- ICESat-2
def fetch_segments(p, cache):
    from sliderule import sliderule, icesat2
    b = p["bbox"]; td = p["tile_deg"]
    nx = max(1, math.ceil((b[2] - b[0]) / td - 1e-9)); ny = max(1, math.ceil((b[3] - b[1]) / td - 1e-9))
    sliderule.init(p["sliderule_url"], verbose=False)
    files = []
    for i in range(nx):
        for j in range(ny):
            x0 = b[0] + i * (b[2] - b[0]) / nx; x1 = b[0] + (i + 1) * (b[2] - b[0]) / nx
            y0 = b[1] + j * (b[3] - b[1]) / ny; y1 = b[1] + (j + 1) * (b[3] - b[1]) / ny
            f = f"{cache}/atl06sr_{i}_{j}.parquet"; files.append(f)
            parms = {"poly": [{"lon": x0, "lat": y0}, {"lon": x1, "lat": y0}, {"lon": x1, "lat": y1}, {"lon": x0, "lat": y1}, {"lon": x0, "lat": y0}],
                     "srt": icesat2.SRT_LAND, "cnf": icesat2.CNF_SURFACE_LOW, "len": p["seg_len_m"], "res": p["seg_res_m"],
                     "ats": 5.0, "cnt": 5, "maxi": 6, "H_min_win": 3.0, "sigma_r_max": 5.0,
                     "t0": p["start"].strftime("%Y-%m-%dT%H:%M:%SZ"), "t1": p["end"].strftime("%Y-%m-%dT%H:%M:%SZ")}
            if p["atl08_ground"]: parms["atl08_class"] = ["atl08_ground"]
            pf = f.replace(".parquet", "_parms.json")
            if os.path.exists(f) and os.path.exists(pf) and json.load(open(pf)) == json.loads(json.dumps(parms)):
                log(f"  tile {i},{j}: cached"); continue
            t0 = time.time(); g = icesat2.atl06p(parms)
            log(f"  tile {i},{j}: {len(g)} segments, {time.time() - t0:.0f} s")
            if len(g): g.to_parquet(f)
            else: files.remove(f)
            json.dump(parms, open(pf, "w"), indent=1)
    return files


def load_segments(files, p):
    import geopandas as gpd
    fs = sorted(set(sum([glob.glob(f) for f in files], [])))
    if not fs: sys.exit("no segment files")
    g = pd.concat([gpd.read_parquet(f).reset_index() for f in fs], ignore_index=True)
    g["lon"] = g.geometry.x.values; g["lat"] = g.geometry.y.values; g = pd.DataFrame(g.drop(columns="geometry"))
    g = g.drop_duplicates(subset=["rgt", "cycle", "spot", "x_atc"])
    b = p["bbox"]; g = g[(g.lon >= b[0]) & (g.lon <= b[2]) & (g.lat >= b[1]) & (g.lat <= b[3])]
    t = pd.to_datetime(g.time); t = t.dt.tz_localize(None) if t.dt.tz is not None else t
    g["time"] = t; g = g[(t >= p["start"]) & (t <= p["end"])].reset_index(drop=True)
    return g


def geoid_height(p, lon, lat):
    """N (m) at lon/lat for the DEM's vertical datum, so that H = h - N."""
    if p["dem"] == "none": return np.zeros(len(lon))   # repeat pairs: the geoid cancels
    gd = p["dem_geoid"] if p["dem"] != "copernicus" else "egm2008"
    if gd.lower() == "ellipsoid": return np.zeros(len(lon))
    if gd.lower() == "egm2008":
        import pyproj
        from pyproj import Transformer
        pyproj.network.set_network_enabled(True)   # PROJ fetches the EGM2008 grid from cdn.proj.org if it is not installed
        t = Transformer.from_crs("EPSG:4979", "EPSG:4326+3855", always_xy=True)
        _, _, H = t.transform(lon, lat, np.zeros(len(lon)))
        N = -np.asarray(H)
        if not np.isfinite(N).all() or np.nanmax(np.abs(N)) < 1e-6:
            sys.exit("EGM2008 grid not available to PROJ; set dem_geoid to a geoid GeoTIFF (e.g. us_nga_egm08_25.tif from cdn.proj.org)")
        return N
    from pyproj import CRS
    return sample(gd, lon, lat, CRS.from_epsg(4326))[0]


# ---------------------------------------------------------------- hexagons
def hexgrid(x0, y0, x1, y1, width, flat=True, square=False):
    """Cell centres (width = flat-to-flat distance); a point belongs to its nearest centre, so square centres give square cells."""
    if square:
        X, Y = np.meshgrid(np.arange(x0, x1 + width, width), np.arange(y0, y1 + width, width))
        return np.column_stack([X.ravel(), Y.ravel()]), width / 2
    r = width / np.sqrt(3); d = 1.5 * r; pts = []
    if flat:
        for i, x in enumerate(np.arange(x0, x1 + d, d)):
            for y in np.arange(y0 - (i % 2) * width / 2, y1 + width, width): pts.append((x, y))
    else:
        for j, y in enumerate(np.arange(y0, y1 + d, d)):
            for x in np.arange(x0 - (j % 2) * width / 2, x1 + width, width): pts.append((x, y))
    return np.array(pts), r


def cell_vertices(xc, yc, r, p):
    """Vertices of a cell: r is the hexagon circumradius, or half the square width."""
    if p["cell_shape"] == "square":
        a = np.radians(45 + 90 * np.arange(4)); r = r * np.sqrt(2)
    else:
        a = np.radians((0 if p["hex_orientation"] == "flat" else 30) + 60 * np.arange(6))
    return np.column_stack([xc + r * np.cos(a), yc + r * np.sin(a)])


def hex_polygon(xc, yc, r, p):
    from shapely.geometry import Polygon
    return Polygon(cell_vertices(xc, yc, r, p))


# ---------------------------------------------------------------- trend fit
def design(t, dem, seasonal):
    cols = [t - 2022.0]
    if seasonal: cols += [np.sin(2 * np.pi * t), np.cos(2 * np.pi * t)]
    for dm in np.unique(dem): cols.append((dem == dm).astype(float))
    return np.column_stack(cols)


def fit(t, dem, d_mm, seasonal, min_passes):
    n = len(t)
    if n < min_passes: return np.nan, np.nan, 0
    X = design(t, dem, seasonal); dof = n - X.shape[1]
    if dof < 2: return np.nan, np.nan, dof
    beta, *_ = np.linalg.lstsq(X, d_mm, rcond=None)
    s2 = np.sum((d_mm - X @ beta) ** 2) / dof
    return beta[0], float(np.sqrt(s2 * np.linalg.pinv(X.T @ X)[0, 0])), dof


def boot_se(t, dem, d_mm, seasonal, n_boot, rng):
    v = []; k = design(t, dem, seasonal).shape[1]
    for _ in range(n_boot):
        i = rng.integers(0, len(t), len(t))
        if len(np.unique(i)) < k + 2: continue   # too few distinct passes: near-singular fit, absurd rates
        r, _, _ = fit(t[i], dem[i], d_mm[i], seasonal, 1)
        if np.isfinite(r): v.append(r)
    return float(np.std(v, ddof=1)) if len(v) > 10 else np.nan


# ---------------------------------------------------------------- main
def run(p, raw):
    import rasterio
    from rasterio.enums import Resampling
    from pyproj import CRS, Transformer
    from scipy.spatial import cKDTree
    from scipy.ndimage import distance_transform_edt
    t_start = time.time()
    os.makedirs(p["out_dir"], exist_ok=True); cache = os.path.join(p["out_dir"], "cache"); os.makedirs(cache, exist_ok=True)
    W = working_crs(p); WGS = CRS.from_epsg(4326); log("working CRS", W.to_string())
    b = p["bbox"]; tw = Transformer.from_crs(WGS, W, always_xy=True)
    ex, ey = tw.transform(np.r_[np.linspace(b[0], b[2], 50), np.full(50, b[2]), np.linspace(b[2], b[0], 50), np.full(50, b[0])],
                          np.r_[np.full(50, b[1]), np.linspace(b[1], b[3], 50), np.full(50, b[3]), np.linspace(b[3], b[1], 50)])
    GR = 30.0
    X0 = math.floor(min(ex) / GR) * GR; X1 = math.ceil(max(ex) / GR) * GR; Y0 = math.floor(min(ey) / GR) * GR; Y1 = math.ceil(max(ey) / GR) * GR
    G = dict(crs=W, transform=rasterio.transform.from_origin(X0, Y1, GR, GR), shape=(int(round((Y1 - Y0) / GR)), int(round((X1 - X0) / GR))))

    # inputs
    PAIR = p["dem"] == "none"
    dems = [copernicus(b, cache)] if p["dem"] == "copernicus" else ([] if PAIR else p["dem"])
    lc = worldcover(b, cache) if p["landcover"].lower() == "worldcover" else p["landcover"]
    files = [p["segments_file"]] if p["segments_file"] else fetch_segments(p, cache)
    # heights, DEM residuals, land cover (cached: sampling a 1 m DEM takes minutes)
    import hashlib
    fs_ = sorted(set(sum([glob.glob(f) for f in files], [])))
    key = json.dumps(dict(f=[(f, os.path.getsize(f), os.path.getmtime(f)) for f in fs_], d=[(f, os.path.getmtime(f)) for f in dems], g=raw["dem_geoid"],
                          l=(lc, os.path.getmtime(lc)), c=W.to_string(), b=b, t=(str(p["start"]), str(p["end"]))), default=str)
    sc = f"{cache}/segments_sampled_{hashlib.md5(key.encode()).hexdigest()[:12]}.parquet"
    if os.path.exists(sc):
        s = pd.read_parquet(sc); n = {"segments": len(s)}; log("sampled segments from cache", len(s))
    else:
        s = load_segments(files, p); n = {"segments": len(s)}; log("segments in bbox and period", len(s))
        if not len(s): sys.exit("no segments")
        s["H"] = s.h_mean.values - geoid_height(p, s.lon.values, s.lat.values)
        s["x"], s["y"] = tw.transform(s.lon.values, s.lat.values)
        s["dem"] = -1; s["zdem"] = np.nan; s["slope"] = np.nan
        OFF = ((0, 0), (10, 0), (-10, 0), (0, 10), (0, -10))
        for k, dp in enumerate(dems):
            m = (s.dem.values < 0)
            if not m.any(): break
            log("sampling DEM", os.path.basename(dp))
            z = sample(dp, s.x.values[m], s.y.values[m], W, 1, OFF)
            sl = np.degrees(np.arctan(np.hypot((z[1] - z[2]) / 20, (z[3] - z[4]) / 20)))
            ok = np.isfinite(z[0]); idx = np.where(m)[0][ok]
            s.loc[idx, "dem"] = k; s.loc[idx, "zdem"] = z[0][ok]; s.loc[idx, "slope"] = sl[ok]
        if PAIR: s["dem"] = 0; s["zdem"] = 0.0; s["slope"] = 0.0   # no DEM: heights compared between passes at the same place
        s = s[s.dem >= 0].reset_index(drop=True)
        log("sampling land cover"); s["lc"] = sample(lc, s.x.values, s.y.values, W, 0)[0]
        keep = ["time", "lon", "lat", "x", "y", "H", "dem", "zdem", "slope", "lc", "rgt", "cycle", "spot", "h_mean", "n_fit_photons", "h_sigma",
                "dh_fit_dx", "w_surface_window_final"]
        s = s[keep]; s.to_parquet(sc)
    n["on_dem"] = len(s)
    s["xdist"] = np.inf
    if p["exclude"]:
        log("exclusion zone", os.path.basename(p["exclude"]))
        if p["exclude"].lower().endswith((".tif", ".tiff", ".vrt")):
            em = to_grid(p["exclude"], G, Resampling.max) > 0
        else:
            import geopandas as gpd
            from rasterio.features import rasterize
            eg = gpd.read_file(p["exclude"]).to_crs(W)
            em = rasterize(((g, 1) for g in eg.geometry if g is not None and not g.is_empty), out_shape=G["shape"], transform=G["transform"], fill=0, dtype="uint8", all_touched=True) > 0
        dist = distance_transform_edt(~em) * GR
        c = np.clip(((s.x.values - X0) // GR).astype(int), 0, G["shape"][1] - 1); r = np.clip(((Y1 - s.y.values) // GR).astype(int), 0, G["shape"][0] - 1)
        s["xdist"] = dist[r, c]

    # quality and residuals
    q = (s.n_fit_photons >= p["min_photons"]) & (s.h_sigma < p["max_h_sigma_m"]) & (s.dh_fit_dx.abs() < p["max_along_slope"]) & \
        (s.w_surface_window_final < p["max_window_m"]) & (s.slope < p["max_terrain_slope_deg"])
    s = s[q].copy(); n["quality"] = len(s)
    s["r"] = s.H - s.zdem
    s["pid"] = s.rgt.astype(np.int64) * 100 + s.cycle.astype(np.int64)
    if not PAIR:   # repeat pairs are filtered on pair differences instead
        med = s.groupby(["pid", "dem"]).r.transform("median")
        s = s[(s.r - med).abs() <= p["max_resid_dev_m"]]
    s = s.reset_index(drop=True); n["residual_filter"] = len(s)
    ts = s.time; y0 = pd.to_datetime(ts.dt.year.astype(str) + "-01-01"); y1 = pd.to_datetime((ts.dt.year + 1).astype(str) + "-01-01")
    s["t"] = ts.dt.year + (ts - y0).dt.total_seconds() / (y1 - y0).dt.total_seconds()
    log("segments after quality filters", len(s), "passes", s.pid.nunique())

    # hexagons and target share
    ox, oy = p["grid_origin"] if p["grid_origin"] else (X0, Y0)
    C, rr = hexgrid(ox, oy, X1, Y1, p["hex_width_m"], p["hex_orientation"] == "flat", p["cell_shape"] == "square"); tree = cKDTree(C)
    log("hexagons", len(C), "land-cover share")
    lcg = to_grid(lc, G, Resampling.nearest)
    gx, gy = np.meshgrid(X0 + GR * (np.arange(G["shape"][1]) + 0.5), Y1 - GR * (np.arange(G["shape"][0]) + 0.5))
    tm = np.isin(lcg, p["target_classes"])
    hg = tree.query(np.column_stack([gx[tm], gy[tm]]))[1]
    share = np.bincount(hg, minlength=len(C)) * GR * GR / (p["hex_width_m"] ** 2 * (1 if p["cell_shape"] == "square" else np.sqrt(3) / 2))
    del gx, gy, lcg, tm, hg
    s["hex"] = tree.query(np.column_stack([s.x.values, s.y.values]))[1]
    is_t = s.lc.isin(p["target_classes"]).values
    is_r = s.lc.isin(p["ref_classes"]).values & (s.xdist.values >= p["exclude_buffer_m"])

    # per-hexagon series: one value per pass
    rng = np.random.default_rng(p["seed"]); R = p["ref_halfwidth_m"]; MN = p["min_segs_per_pass"]

    def dem_series(h):
        ti = np.asarray(T[h]); tp = pd.DataFrame({k: v[ti] for k, v in cols.items()})
        major = tp.groupby("pid").dem.agg(lambda v: v.value_counts().index[0])
        tp = tp[tp.dem.values == major.reindex(tp.pid).values]
        tg = tp.groupby("pid").agg(t=("t", "mean"), dem=("dem", "first"), m=("r", "median"), n=("r", "size")); tg = tg[tg.n >= MN]
        if len(tg) < p["min_passes"]: return None, 0
        if p["reference"] == "local":
            xc, yc = C[h]; a, z = np.searchsorted(xs, [xc - R, xc + R]); w = xs_o[a:z]
            w = w[(np.abs(Y[w] - yc) <= R) & is_r[w] & (cols["hex"][w] != h) & np.isin(cols["pid"][w], tg.index.values)]
            rg = pd.DataFrame({"pid": cols["pid"][w], "dem": cols["dem"][w], "r": cols["r"][w]}).groupby(["pid", "dem"]).r.agg(["median", "size"])
            key = pd.MultiIndex.from_arrays([tg.index.values, tg.dem.values])
            rm = rg.reindex(key); ok = (rm["size"].values >= MN)
            tg = tg[ok].copy(); tg["d"] = tg.m.values - rm["median"].values[ok]
        else:
            tg["d"] = tg.m.values
        return tg, len(tp)

    def pair_series(h):
        """Repeat pairs: median height change between each two passes (target minus reference pairs), inverted to one height per pass.
        Each connected network of passes has its own constant (passed to the fit as 'dem')."""
        from scipy.sparse.csgraph import connected_components
        pi = np.asarray(T[h])
        g = pd.DataFrame(dict(pa=PR["pa"][pi], pb=PR["pb"][pi], dh=PR["dh"][pi])).groupby(["pa", "pb"]).dh.agg(["median", "size"])
        g = g[g["size"] >= MN]
        if p["reference"] == "local" and len(g):
            xc, yc = C[h]; a, z = np.searchsorted(pxs, [xc - R, xc + R]); w = px_o[a:z]
            w = w[(np.abs(PR["my"][w] - yc) <= R) & PR["is_r"][w] & (PR["ha"][w] != h) & (PR["hb"][w] != h)]
            rg = pd.DataFrame(dict(pa=PR["pa"][w], pb=PR["pb"][w], dh=PR["dh"][w])).groupby(["pa", "pb"]).dh.agg(["median", "size"]).reindex(g.index)
            ok = rg["size"].values >= MN; g = g[ok]; D = g["median"].values - rg["median"].values[ok]
        else:
            D = g["median"].values
        if not len(g): return None, 0
        pa = g.index.get_level_values(0).values; pb = g.index.get_level_values(1).values
        ps = np.unique(np.r_[pa, pb])
        if len(ps) < p["min_passes"]: return None, 0
        ia = np.searchsorted(ps, pa); ib = np.searchsorted(ps, pb)
        A = np.zeros((len(g), len(ps))); A[np.arange(len(g)), ia] = -1; A[np.arange(len(g)), ib] = 1
        hgt = np.linalg.lstsq(A, D, rcond=None)[0]   # minimum norm: zero mean within each connected network
        _, comp = connected_components((np.abs(A.T) @ np.abs(A)) > 0, directed=False)
        return pd.DataFrame(dict(t=tpass.reindex(ps).values, dem=comp, d=hgt), index=ps), int(g["size"].sum())

    if PAIR:
        P = cKDTree(np.column_stack([s.x.values, s.y.values])).query_pairs(p["pair_radius_m"], output_type="ndarray")
        P = P[s.pid.values[P[:, 0]] != s.pid.values[P[:, 1]]]
        sw = s.t.values[P[:, 0]] > s.t.values[P[:, 1]]; ia = np.where(sw, P[:, 1], P[:, 0]); ib = np.where(sw, P[:, 0], P[:, 1])   # a earlier than b
        pr = pd.DataFrame(dict(pa=s.pid.values[ia], pb=s.pid.values[ib], dh=s.H.values[ib] - s.H.values[ia], mx=(s.x.values[ia] + s.x.values[ib]) / 2,
                               my=(s.y.values[ia] + s.y.values[ib]) / 2, ha=s.hex.values[ia], hb=s.hex.values[ib],
                               is_t=is_t[ia] & is_t[ib], is_r=is_r[ia] & is_r[ib]))
        dev = pr.dh - pr.groupby(["pa", "pb"]).dh.transform("median"); pr = pr[dev.abs() <= p["max_resid_dev_m"]].reset_index(drop=True)
        n["pairs"] = len(pr); log("repeat pairs within", p["pair_radius_m"], "m:", len(pr), "between", pr.groupby(["pa", "pb"]).ngroups, "pass pairs")
        PR = {k: pr[k].values for k in pr.columns}; px_o = np.argsort(PR["mx"]); pxs = PR["mx"][px_o]
        tpass = s.groupby("pid").t.mean()
        T = pr[pr.is_t & (pr.ha == pr.hb)].groupby("ha").groups; series = pair_series
    else:
        xs_o = np.argsort(s.x.values); xs = s.x.values[xs_o]; Y = s.y.values
        cols = dict(pid=s.pid.values, dem=s.dem.values, r=s.r.values, t=s.t.values, hex=s.hex.values)
        T = s[is_t].groupby("hex").groups; series = dem_series
    rows = []; cand = [h for h in T if share[h] >= p["min_target_share"]]
    log("fitting", len(cand), "hexagons with target data and share >=", p["min_target_share"])
    for i, h in enumerate(cand):
        tg, nseg = series(h)
        if tg is None or len(tg) < p["min_passes"]: continue
        d = tg.d.values
        rate, se_f, dof = fit(tg.t.values, tg.dem.values, d * 1000.0, p["seasonal"], p["min_passes"])
        if not np.isfinite(rate): continue
        se_b = boot_se(tg.t.values, tg.dem.values, d * 1000.0, p["seasonal"], p["n_boot"], rng) if p["se_method"] == "bootstrap" else np.nan
        rows.append(dict(hid=int(h), x=C[h][0], y=C[h][1], rate=rate, se=se_b if p["se_method"] == "bootstrap" else se_f, se_formula=se_f,
                         se_boot=se_b, n_pass=len(tg), dof=dof, target_share=round(float(share[h]), 3), n_target_seg=nseg))
        if (i + 1) % 100 == 0: log(f"  {i + 1}/{len(cand)}")
    Hx = pd.DataFrame(rows); log("hexagons with a rate", len(Hx))
    if not len(Hx): sys.exit("no hexagon has enough passes")

    # outputs
    stem = os.path.join(p["out_dir"], p["name"])
    res = p["out_res_m"]; gx0 = math.floor(X0 / res) * res; gy1 = math.ceil(Y1 / res) * res
    nx_ = int(math.ceil((X1 - gx0) / res)); ny_ = int(math.ceil((gy1 - Y0) / res))
    px, py = np.meshgrid(gx0 + res * (np.arange(nx_) + 0.5), gy1 - res * (np.arange(ny_) + 0.5))
    hp = tree.query(np.column_stack([px.ravel(), py.ravel()]))[1].reshape(px.shape)
    rate = np.full(len(C), np.nan, np.float32); se = np.full(len(C), np.nan, np.float32)
    rate[Hx.hid.values] = Hx.rate.values; se[Hx.hid.values] = Hx.se.values
    prof = dict(driver="GTiff", height=ny_, width=nx_, count=2, dtype="float32", crs=W, transform=rasterio.transform.from_origin(gx0, gy1, res, res),
                nodata=np.nan, compress="deflate", tiled=True, blockxsize=256, blockysize=256)
    with rasterio.open(stem + "_rate.tif", "w", **prof) as d:
        d.write(rate[hp], 1); d.write(se[hp], 2)
        d.set_band_description(1, "rate (mm/yr)"); d.set_band_description(2, f"SE (mm/yr, {p['se_method']})")
        d.update_tags(software=f"is2rate {__version__}", units="mm/yr", reference=("relative to " + ",".join(map(str, p["ref_classes"])) +
                      f" within +-{R:.0f} m on the same passes") if p["reference"] == "local" else "absolute", **{f"param_{k}": v for k, v in raw.items()})
    write_qml(stem + "_rate.qml"); log("wrote", stem + "_rate.tif")
    if p["write_gpkg"]:
        import geopandas as gpd
        gpd.GeoDataFrame(Hx, geometry=[hex_polygon(a, c, rr, p) for a, c in zip(Hx.x, Hx.y)], crs=W).to_file(stem + "_hexagons.gpkg"); log("wrote", stem + "_hexagons.gpkg")
    if p["write_png"]: quicklook(Hx, rr, W, p, stem + "_rate.png", (X0, Y0, X1, Y1))
    summ = dict(version=__version__, params=raw, crs=W.to_string(), counts=n, passes=int(s.pid.nunique()), hexagons=len(Hx),
                median_passes=float(Hx.n_pass.median()), median_rate=float(Hx.rate.median()), median_se=float(Hx.se.median()),
                runtime_s=round(time.time() - t_start))
    json.dump(summ, open(stem + "_summary.json", "w"), indent=1, default=str); log("done", json.dumps({k: summ[k] for k in ("hexagons", "median_passes", "runtime_s")}))


def write_qml(path):
    items = "\n".join(f'          <item alpha="255" value="{"inf" if np.isinf(v) else v}" label="{l.replace(">", "&gt;")}" color="{c}"/>' for v, c, l in BREAKS)
    open(path, "w", encoding="utf-8").write(f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.28" styleCategories="AllStyleCategories">
  <pipe>
    <rasterrenderer type="singlebandpseudocolor" band="1" opacity="1" alphaBand="-1" classificationMin="-40" classificationMax="20" nodataColor="">
      <rastershader>
        <colorrampshader colorRampType="DISCRETE" classificationMode="3" clip="0" labelPrecision="0">
{items}
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
  </pipe>
</qgis>
""")


def quicklook(Hx, rr, W, p, path, ext):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Patch
    polys = [cell_vertices(x, y, rr, p) for x, y in zip(Hx.x, Hx.y)]
    up = np.array([v for v, _, _ in BREAKS]); cl = [BREAKS[int(np.searchsorted(up, v))][1] for v in Hx.rate.values]
    w = ext[2] - ext[0]; h = ext[3] - ext[1]
    fig, ax = plt.subplots(figsize=(7, 7 * h / w + 1.2))
    ax.add_collection(PolyCollection(polys, facecolors=cl, edgecolors="white", linewidths=0.4))
    ax.set_xlim(ext[0], ext[2]); ax.set_ylim(ext[1], ext[3]); ax.set_aspect("equal"); ax.set_facecolor("#d4d4d4")
    ax.ticklabel_format(useOffset=False, style="plain"); ax.tick_params(labelsize=7)
    ref = "relative to reference land on the same passes" if p["reference"] == "local" else "absolute"
    ax.set_title(f"ICESat-2, {p['hex_width_m'] / 1000:g} km {'squares' if p['cell_shape'] == 'square' else 'hexagons'} ({ref}), {p['start']:%Y-%m-%d} to {p['end']:%Y-%m-%d}", loc="left", fontsize=8)
    ax.legend(handles=[Patch(facecolor=c, edgecolor="#999999", lw=0.4, label=l) for _, c, l in BREAKS], title="Rate (mm year⁻¹)", fontsize=6.5,
              title_fontsize=7, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.06), frameon=False)
    ax.set_xlabel(W.name, fontsize=7)
    fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig); log("wrote", path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("params", help="parameter file")
    ap.add_argument("--init", action="store_true", help="write a parameter file with default values and exit")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="override a parameter")
    a = ap.parse_args()
    if a.init:
        if os.path.exists(a.params): sys.exit(f"{a.params} exists")
        write_template(a.params); return
    p, raw = read_params(a.params, a.set)
    run(p, raw)


if __name__ == "__main__":
    main()
