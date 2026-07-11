<!--
<img src="https://drive.google.com/uc?id=1fgYuG7jpnekZrkoL_PdVUnSiUFBFX-vI" alt="Logo" width="150" style="float: left; margin-right: 10px;">
-->

<img src="https://drive.google.com/uc?id=1szqLViKqTX5C1XF8uV7HbIst0i6Xvv7g" alt="Logo" width="300">

# GEE Python Toolkit

![In Development](https://img.shields.io/badge/Status-In%20Development-yellow)
![Languages](https://img.shields.io/badge/Languages-Python%20%7C%20GEE-blue)

A collection of Python scripts for driving [Google Earth Engine](https://earthengine.google.com/) (GEE) from the **Earth Engine Python API** locally, rather than the JavaScript Code Editor. The scripts preprocess and extract remote sensing and geospatial covariates — digital elevation models, terrain and hydrological indices, soil properties, and multi-sensor spectral time series — and export the results (GeoTIFFs and CSVs) to Google Drive.


---

## Repository structure

```
gee_python_toolkit/
├── python/                  # GEE Python scripts (run these)
│   ├── _gee_config.py       # Shared config: EE project + Drive folder
│   ├── *.py                 # Dataset/workflow scripts
│   └── utils/               # Shared helper modules
├── gee_compute_reports/     # EECU compute-usage reports (txt)
├── AGENTS.md                # Coding conventions
└── README.md
```

---

## GEE via VS Code (Python API)

Python API workflow (not the JS Code Editor). `geemap` handles interactive maps, since VS Code doesn't render `ee` objects natively.

### Prerequisites (Cloud side, one-time)

- **Cloud project + EE API enabled.** Every request routes through a Google Cloud project with the Earth Engine API turned on.
- **Registration + tier.** Noncommercial projects default to the Community Tier; the tier can be changed anytime.

### Environment setup

- **Extensions:** Python (Microsoft) + Jupyter. Pylance optional (typing).
- **Isolated environment** (keeps EE deps off the ArcPy/base env):

```bash
python -m venv .venv
# or: conda create -n gee python=3.11
```

- **Install packages:**

```bash
pip install earthengine-api geemap
# geemap pulls in folium/ipyleaflet, ipywidgets, etc.
```

- **Select the interpreter:** Command Palette → *Python: Select Interpreter* → `.venv`/conda env.
- **Optional:** install the `gcloud` CLI. GEE uses Google Cloud for auth; `gcloud` gives the cleanest flow and lets you manage/switch projects from the terminal. See [Spatial Thoughts](https://courses.spatialthoughts.com/install-gee-python-api.html).

### Authentication & initialization

One-time auth (opens a browser, stores a token at `~/.config/earthengine/credentials`):

```python
import ee
ee.Authenticate()            # or run `earthengine authenticate` in the terminal
ee.Initialize(project='your-project-id')
```

Set a default project to drop the argument later:

```bash
earthengine set_project your-project-id
```

Sanity check:

```python
print(ee.String('Hello from the Earth Engine servers!').getInfo())
```

The scripts in this repository handle initialization for you via `utils.gee_utils.initialize_ee()`, which reads the project ID from [`_gee_config.py`](python/_gee_config.py).

### Interactive workflow
How to see Earth Engine data on a map when working from the Python API using Jupyter Notebooks (`.ipynb` files)


- **`.ipynb` / interactive window** — best for exploration; `geemap` renders inline:

```python
import geemap
m = geemap.Map(center=[54.5, -113.5], zoom=6)   # north-central AB
m.add_layer(ee.Image('USGS/SRTMGL1_003'), {'min': 0, 'max': 3000}, 'DEM')
m
```

---

## Configuration

All scripts read shared settings from [`python/_gee_config.py`](python/_gee_config.py):

| Setting | Description |
|---------|-------------|
| `EE_PROJECT` | Google Cloud project ID registered for Earth Engine. Find it in the [Code Editor](https://code.earthengine.google.com) (profile icon, top right) or register at [code.earthengine.google.com/register](https://code.earthengine.google.com/register). |
| `DRIVE_FOLDER` | Default Google Drive folder for exports (e.g. `gee_exports`). |

Set `EE_PROJECT` to your own project before running any script. Most scripts also expose user parameters near the top (export scale, CRS, test vs. full AOI, whether to write a compute report).

---

## Running a script

1. Activate the environment and select the interpreter (above).
2. Set `EE_PROJECT` in [`python/_gee_config.py`](python/_gee_config.py).
3. Open a script in `python/`, review its user parameters, and run it.
4. Exports appear as tasks in the [Earth Engine Tasks tab](https://code.earthengine.google.com/tasks) and land in your `DRIVE_FOLDER` on Google Drive.

---

## Contents

### Scripts (`python/`)

| File | Description |
|------|-------------|
| [_gee_config.py](python/_gee_config.py) | Shared configuration (Earth Engine project ID and default Drive export folder) used by all scripts. |
| [fabdem.py](python/fabdem.py) | Mosaics the FABDEM DEM, clips to the US + Canada, and exports a GeoTIFF. |
| [fabdem_twi_alberta.py](python/fabdem_twi_alberta.py) | Computes the Topographic Wetness Index (TWI = ln(α/tanβ)) for Alberta using FABDEM slope and MERIT Hydro upslope area. |
| [global_geomorphometric_layers.py](python/global_geomorphometric_layers.py) | Loads Geomorpho90m geomorphometric variables, mosaics and clips them, and exports a multiband GeoTIFF for Alberta. |
| [hydrologically_adjusted_elevation.py](python/hydrologically_adjusted_elevation.py) | Extracts Height Above Nearest Drainage (HAND) from MERIT Hydro and exports it for Alberta. |
| [nrcan_topographic_indices.py](python/nrcan_topographic_indices.py) | Derives terrain metrics (elevation, slope, aspect, northness, eastness) from the NRCan/CDEM DEM. |
| [hihydrosoil_v2.py](python/hihydrosoil_v2.py) | Exports HiHydroSoil v2.0 soil hydraulic properties (native ~250 m and 1000 m) and extracts point-level values to CSV. |
| [soil_grids_250.py](python/soil_grids_250.py) | Exports ISRIC SoilGrids 250m v2.0 soil properties with unit rescaling and extracts point-level values to CSV. |
| [landsat_time_series.py](python/landsat_time_series.py) | Builds an annual Landsat 5/7/8/9 spectral-index time series and exports multiband GeoTIFFs (native and focal scales). |
| [landsat_time_series_to_poly.py](python/landsat_time_series_to_poly.py) | Summarizes the Landsat spectral-index time series to polygons and exports a per-polygon-per-date CSV. |
| [sentinel2_time_series.py](python/sentinel2_time_series.py) | Builds an annual Sentinel-2 spectral-index time series (with NDRS forest bands) and exports multiband GeoTIFFs. |
| [modis_land_cover_dynamics.py](python/modis_land_cover_dynamics.py) | Extracts annual MODIS MCD12Q2 phenology bands and exports multiband GeoTIFFs (native and focal scales). |
| [land_cover_time_series.py](python/land_cover_time_series.py) | Exports annual forest land cover (Canada VLCE2) GeoTIFFs, one per year. |

### Helper modules (`python/utils/`)

| File | Description |
|------|-------------|
| [gee_utils.py](python/utils/gee_utils.py) | Authentication/initialization and a Drive export wrapper with optional task monitoring. |
| [gee_helpers.py](python/utils/gee_helpers.py) | General imagery helpers: date lists, AOI tiling, band filtering, point/focal reductions, image statistics, and Drive export (adapted in part from geeTools). |
| [compute_report.py](python/utils/compute_report.py) | Collects Earth Engine EECU compute-usage information and writes plain-text reports. |
| [calculate_twi.py](python/utils/calculate_twi.py) | Calculates TWI from FABDEM. |
| [geomorpho90m.py](python/utils/geomorpho90m.py) | Extracts Geomorpho90m layers. |
| [canopy_height.py](python/utils/canopy_height.py) | Retrieves canopy height. |
| [gap_filling.py](python/utils/gap_filling.py) | Gap-filling functions for image collections. |
| [annual_forest_land_cover.py](python/utils/annual_forest_land_cover.py) | Retrieves annual forest land cover (Canada VLCE2). |
| [get_annual_forest_tree_species.py](python/utils/get_annual_forest_tree_species.py) | Retrieves annual forest tree species. |
| [proportion_forested_land_cover.py](python/utils/proportion_forested_land_cover.py) | Computes land cover proportions. |
| [proportion_of_leading_tree_species.py](python/utils/proportion_of_leading_tree_species.py) | Computes leading tree species proportions. |
| [landsat_time_series.py](python/utils/landsat_time_series.py) | Builds a Landsat image time series. |
| [landsat_indices_and_masks.py](python/utils/landsat_indices_and_masks.py) | Landsat spectral indices and mask functions. |
| [sentinel_time_series.py](python/utils/sentinel_time_series.py) | Builds a Sentinel-2 image time series. |
| [sentinel_indices_and_masks.py](python/utils/sentinel_indices_and_masks.py) | Sentinel-2 spectral indices and mask functions. |
| [masks.py](python/utils/masks.py) | General image masking functions. |
| [image_collection_to_features.py](python/utils/image_collection_to_features.py) | Extracts image-collection values to features/polygons. |
| [image_collection_to_points.py](python/utils/image_collection_to_points.py) | Extracts image-collection values to points. |
| [image_to_points.py](python/utils/image_to_points.py) | Extracts single-image values to points. |

---

## Compute reports

Scripts with `COMPUTE_REPORT` enabled write EECU usage summaries to `gee_compute_reports/` via [`compute_report.py`](python/utils/compute_report.py). These capture per-algorithm EECU profiles for computations and total batch EECU-seconds for export tasks — useful for finding compute choke points. 

---

## Data sources

| Dataset | Used by |
|---------|---------|
| FABDEM (`projects/sat-io/open-datasets/FABDEM`) | `fabdem.py`, `fabdem_twi_alberta.py` |
| MERIT Hydro (`MERIT/Hydro/v1_0_1`) | `fabdem_twi_alberta.py`, `hydrologically_adjusted_elevation.py` |
| Geomorpho90m (`projects/sat-io/open-datasets/Geomorpho90m`) | `global_geomorphometric_layers.py` |
| NRCan/CDEM | `nrcan_topographic_indices.py` |
| HiHydroSoil v2.0 (FutureWater / sat-io) | `hihydrosoil_v2.py` |
| SoilGrids 250m v2.0 (`projects/soilgrids-isric/*_mean`) | `soil_grids_250.py` |
| Landsat 5/7/8/9 SR (`LANDSAT/*/C02/T1_L2`) | `landsat_time_series.py`, `landsat_time_series_to_poly.py` |
| Sentinel-2 SR (`COPERNICUS/S2_SR_HARMONIZED`) | `sentinel2_time_series.py` |
| MODIS MCD12Q2 (`MODIS/061/MCD12Q2`) | `modis_land_cover_dynamics.py` |
| Canada Forest LC VLCE2 (`projects/sat-io/open-datasets/CA_FOREST_LC_VLCE2`) | `land_cover_time_series.py`, `sentinel2_time_series.py` |
| FAO GAUL boundaries (`FAO/GAUL/2015/*`) | AOI clipping (most scripts) |

---

