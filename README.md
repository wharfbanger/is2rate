# is2rate

Rates of vertical land motion from ICESat-2 on a hexagon grid, written as a GeoTIFF. Band 1 is the rate and band 2 its standard error (SE), both in mm year⁻¹.

The method is the one used for the Waikato, New Zealand (Harvey, in prep., *ICESat-2 confirms peat subsidence in the Waikato, New Zealand*). The defaults are the Waikato settings.

## Method

1. ICESat-2 ground segments come from SlideRule as ATL06-SR: 20 m segments every 10 m, ATL08 ground photons only, for the bounding box and period.
2. Segment heights are converted to the DEM's vertical datum and differenced from a reference DEM to give a residual. Ground that lowers after the DEM epoch gives increasingly negative residuals on later passes.
3. A segment is kept if it meets all of these thresholds:
   - at least 10 photons;
   - height uncertainty below 0.15 m;
   - along-track slope below 0.03;
   - photon window below 4 m;
   - DEM slope below 2° (measured over ±10 m).

   Segments more than 0.5 m from the median residual of their pass are then removed. These fall on buildings, trees, water or misclassified photons.
4. On each pass, the median residual of the target land cover in a hexagon (grassland and cropland by default) is differenced from the median residual of reference land cover. The reference is taken within ±10 km of the hexagon centre on the same pass, excluding the hexagon itself and any `exclude` zone plus its buffer. This removes orbit, atmospheric and other errors common to a pass.

   Where DEMs overlap, each pass uses the DEM that holds most of the hexagon's segments.
5. The rate is the least-squares trend of these per-pass differences. The fit also includes an annual sinusoid and a constant for each DEM. A rate needs at least 8 passes, each with at least 10 target and 10 reference segments.
6. The SE is computed in one of two ways:
   - by resampling whole passes, 300 times (default);
   - from the least-squares formula.

   Passes, not segments, are the independent measurements. Rates are not reported for hexagons with less than 25% of their area in target land cover.

Rates are relative to the reference land. Choose a reference that is stable, or check its stability separately (for example with `reference = none`, which gives absolute rates).

## Install

Requires Python 3.10 or later.

```
pip install -r requirements.txt
```

No account is needed. The Copernicus DEM and ESA WorldCover are read from the AWS open-data buckets, and ICESat-2 segments come from the public SlideRule server. Only the window covering the bounding box is read.

## Use

```
python is2rate.py --init params.txt     # write a parameter file with defaults
# edit params.txt: at least bbox, start, end
python is2rate.py params.txt
python is2rate.py params.txt --set hex_width_m=5000 --set end=2022-06-30   # override a parameter
```

SlideRule requests are cached in `out_dir/cache` and reused when the parameters are unchanged.

## Parameters

See `params_example.txt` for every parameter and its default. The main ones:

| Parameter | Default | Meaning |
|---|---|---|
| `bbox` | — | lon_min, lat_min, lon_max, lat_max (WGS84) |
| `start`, `end` | 2018-10-13, 2026-12-31 | period (YYYY-MM-DD or YYYY) |
| `crs` | auto | output CRS (UTM zone of the bbox centre, or an EPSG code) |
| `dem` | copernicus | Copernicus GLO-30, or your own DEM(s) in priority order |
| `dem_geoid` | egm2008 | DEM vertical datum: egm2008, ellipsoid, or a geoid-height GeoTIFF |
| `landcover` | worldcover | ESA WorldCover 2021, or a GeoTIFF with WorldCover codes |
| `target_classes`, `ref_classes` | 30, 40 | grassland and cropland |
| `exclude`, `exclude_buffer_m` | none, 500 | areas not used as reference (e.g. peat soils) |
| `reference` | local | local (relative to reference land) or none (absolute) |
| `ref_halfwidth_m` | 10000 | reference window, ±m from the hexagon centre |
| `hex_width_m` | 3000 | hexagon width, flat to flat (m) |
| `hex_orientation` | flat | flat-topped (as the published maps) or pointy-topped |
| `out_res_m` | 100 | GeoTIFF pixel size (m) |
| `se_method` | bootstrap | band 2: bootstrap or formula |

The DEM is only a fixed reference surface. It is also used for the slope filter. Static DEM errors are constant in time and are absorbed by the DEM constant, so they do not bias the trend. They do add scatter, because tracks fall in different places on each pass. Use LiDAR where it exists; Copernicus works anywhere but is noisier (see Validation).

## Outputs (in `out_dir`)

| File | Content |
|---|---|
| `NAME_rate.tif` | band 1 rate, band 2 SE (mm year⁻¹), float32, NaN where there is no rate |
| `NAME_rate.qml` | QGIS style (loaded automatically) with the colour classes of the published maps |
| `NAME_hexagons.gpkg` | hexagons with rate, SE (bootstrap and formula), number of passes and target share |
| `NAME_rate.png` | quick-look map |
| `NAME_summary.json` | parameters, segment counts and summary statistics |

## Validation (Waikato, 2018-10-13 to 2026-10-21)

These runs used the published hexagon grid, 2024 and 2021 LINZ LiDAR (NZGeoid2016), WorldCover, and S-map Organic Soils as `exclude`.

- **Against the published map.** 429 of the 440 published 3 km hexagons are reproduced, with rate r = 0.82, a median difference of +0.2 mm year⁻¹ and a median |difference| of 3.6 mm year⁻¹. The published hexagon SEs are about 6–9 mm year⁻¹.
- **Against an independent re-implementation with `reference = none`.** 551 of 551 hexagons agree, with r = 0.999 and a median difference of 0.00 mm year⁻¹.
- **Copernicus GLO-30 in place of LiDAR.** The bias is small (median +0.5 mm year⁻¹), but the hexagon rates are noisy (r = 0.34 against the LiDAR run) and SEs are larger (median 13 against 8 mm year⁻¹). Treat single hexagons from Copernicus with caution, and consider `hex_width_m = 5000`.

The runtime was 11 min for a 100 × 100 km area with 1 m LiDAR, most of it spent sampling the DEM. Sampled segments are cached, so a re-run with different hexagon or reference settings takes under 2 min.

An online test also ran successfully, on a 0.2° × 0.15° area near Kopuatai with all defaults. It used Copernicus and WorldCover read from AWS, EGM2008 from PROJ, and segments from SlideRule, and took 4 min.

## Caveats

- **SlideRule read failures.** SlideRule can report occasional granule read failures, logged as `Alert <-7>: H5Coro::Future read failure`. Those passes are then partly missing. Delete `out_dir/cache/atl06sr_*` and re-run to retry them.

- **Few passes.** Hexagons have few passes (about 13 at 3 km in the Waikato), so their SEs are approximate. Map patterns are more reliable than single hexagons.
- **Correlated passes.** Passes close in time can share weather and surface conditions. SEs that treat passes as independent may be too small.
- **Off-pointing.** ICESat-2 tracks are displaced between passes (off-pointing), so a hexagon is sampled at different places on each pass. Within-hexagon variability adds noise, and can bias the result if it is correlated with time.
- **Unsuitable ground.** Steep, forested, built or wet ground is removed by the quality and land-cover filters. Results are for the target land cover only.

## Licence

MIT. See `LICENSE`.

## Experimental: repeat pairs (`dem = none`)

This mode needs no DEM. Segments from different passes within `pair_radius_m` (default 10 m) of each other form repeat pairs. The pairs are differenced against the reference land and inverted to give one height per pass.

It was tested in the Waikato against the LiDAR run:

- **Bias** is small (median difference +0.1 to +1.2 mm year⁻¹).
- **Noise** is high: r = 0.3–0.44 at 3 km and 0.27 at 10 km.
- **Coverage** is lower, because only 30–45% of segments have a repeat. Rates are reported for 205–461 of 562 hexagons, depending on pair radius and minimum pairs.

This is no better than Copernicus, so Copernicus remains the default where no LiDAR exists.
