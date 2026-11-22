"""Convert the SPARTACUS JSON download into a gridded NetCDF file.

Input:  data/raw/spartacus-v2-1y-1km_waldviertel_1996-2026.json
        (GeoSphere Austria Data Hub, SPARTACUS v2 yearly, 1 km; GeoJSON point features)
Output: data/processed/spartacus_yearly_waldviertel.nc  – dims (time, y, x)

The grid points are the cell centres of the SPARTACUS 1 km grid in
ETRS89 / Austria Lambert (EPSG:3416), so they map one-to-one onto a regular raster.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "raw" / "spartacus-v2-1y-1km_waldviertel_1996-2026.json"
OUT = ROOT / "data" / "processed" / "spartacus_yearly_waldviertel.nc"
CRS = "EPSG:3416"


def main():
    raw = json.loads(SRC.read_text())
    time = pd.to_datetime(raw["timestamps"]).tz_localize(None)
    feats = raw["features"]
    lon = np.array([f["geometry"]["coordinates"][0] for f in feats])
    lat = np.array([f["geometry"]["coordinates"][1] for f in feats])
    x, y = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform(lon, lat)
    x, y = np.round(x, -1), np.round(y, -1)            # snap to the 1 km grid (cell centres)

    xs, ys = np.unique(x), np.unique(y)[::-1]          # y from north to south
    ix, iy = np.searchsorted(xs, x), len(ys) - 1 - np.searchsorted(ys[::-1], y)

    def cube(param, scale=1.0):
        arr = np.full((len(time), len(ys), len(xs)), np.nan, dtype="float32")
        for f, i, j in zip(feats, iy, ix):
            arr[:, i, j] = np.asarray(f["properties"]["parameters"][param]["data"], float) * scale
        return arr

    grid_lon = np.full((len(ys), len(xs)), np.nan)
    grid_lat = np.full((len(ys), len(xs)), np.nan)
    grid_lon[iy, ix], grid_lat[iy, ix] = lon, lat

    ds = xr.Dataset(
        {
            "TM": (("time", "y", "x"), cube("TM"), {"long_name": "yearly mean air temperature", "units": "degC"}),
            "RR": (("time", "y", "x"), cube("RR"), {"long_name": "yearly precipitation sum", "units": "mm"}),
            "SA": (("time", "y", "x"), cube("SA", 1 / 3600), {"long_name": "yearly sunshine duration", "units": "h"}),
        },
        coords={"time": time, "y": ys, "x": xs,
                "lon": (("y", "x"), grid_lon), "lat": (("y", "x"), grid_lat)},
        attrs={"title": "SPARTACUS v2 yearly 1 km – Waldviertel / Horn basin subset",
               "source": "GeoSphere Austria, SPARTACUS v2 (spartacus-v2-1y-1km), CC BY 4.0",
               "crs": CRS},
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(OUT)
    print(ds)
    print(f"\n{int(ds.TM.isel(time=0).notnull().sum())} grid cells, "
          f"{len(xs)} x {len(ys)} at 1 km -> {OUT.name}")


if __name__ == "__main__":
    main()
