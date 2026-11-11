# Parameter files for Harvey (in prep.), Waikato

These are the parameter files for the paper *ICESat-2 confirms peat subsidence in the Waikato, New Zealand*. Paths under `data/` stand for the user's local copies of these inputs:

- LINZ Waikato LiDAR 1 m DEMs, 2024 and 2021, mosaicked as VRTs;
- NZGeoid2016;
- S-map Organic Soils;
- target polygons.

ESA WorldCover is read automatically.

## Hexagon maps

These give Figures 4d, 5 and 6a and Table 4.

| File | Grid |
|---|---|
| `waikato_paper_3km.txt` | 3 km hexagons on the published grid, `seasonal = regional` |
| `waikato_paper_5km.txt` | 5 km hexagons on the published grid, `seasonal = regional` |

## Maximum rates

These give Figures 7 and 8 and Table 5. They use `waikato_paper_targets.txt`, with one run per target polygon (`--set target_area=... --set name=...`) for each of three periods:

- `--set start=2018-10-13 --set end=2022-06-30`
- `--set start=2022-07-01 --set end=2026-10-20`
- the whole record, 2018-10-13 to 2026-10-20 (the file's default)

`ref_area` is the study-area rectangle (NZTM 1750000–1846000 E, 5785000–5886000 N). The reference is farmland more than 500 m from S-map Organic Soils, on the same passes.

### Target polygons

- **Distance bands inside S-map Organic Soils:** inward buffers of the dissolved Organic Soils polygons. The bands are 0–250, 250–500, 500–1000, 1000–1500, 1500–2000 and 2000–3000 m. For example, the 0–250 m band is `organic − buffer(organic, −250 m)`.
- **Farmland selected by map rate:** 240 m cells of the paper's maps at or below a threshold, polygonised and dissolved.
  - InSAR map, relative to mineral-soil farmland (Figure 4c): thresholds −5 to −70 mm/yr.
  - LiDAR − Copernicus map (Figure 4a): thresholds −5 to −40 mm/yr.

Each run writes `NAME_target.json` with the rate, its formula and bootstrap SEs, the 95% CI and the number of passes. Table 5 reports the formula SE.
