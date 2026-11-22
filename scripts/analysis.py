"""Gridded climate analysis (SPARTACUS, 1 km) for the eastern Waldviertel / Horn basin, 1996–2025.

- 30-year means of temperature, precipitation and sunshine (maps)
- per-grid-cell trends: Theil-Sen slope + Mann-Kendall significance (maps)
- regional mean time series and anomalies
- relation between local climate (cold plateau vs. warm basin) and warming rate
- comparison with GeoSphere station data (Horn, Gars am Kamp)

Run scripts/prepare_data.py first.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import TwoSlopeNorm
from pyproj import Transformer
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
NC = ROOT / "data" / "processed" / "spartacus_yearly_waldviertel.nc"
STATIONS_CSV = ROOT / "data" / "external" / "geosphere_stations_annual_horn_gars.csv"
PROC = ROOT / "data" / "processed"
FIG = ROOT / "figures"

FIRST, LAST = 1996, 2025                       # 30 complete years
# Approximate station positions (town centres); used to pick the nearest grid cell
STATIONS = {"Horn": (15.657, 48.663), "Gars am Kamp": (15.653, 48.590)}
TOWNS = {"Zwettl": (15.169, 48.604), "Krems": (15.609, 48.410), "Eggenburg": (15.818, 48.639),
         **STATIONS}


# ------------------------------------------------------------------ helpers
def load():
    ds = xr.open_dataset(NC)
    return ds.sel(time=slice(f"{FIRST}-01-01", f"{LAST}-12-31"))


def to_xy(lon, lat):
    return Transformer.from_crs("EPSG:4326", "EPSG:3416", always_xy=True).transform(lon, lat)


def extent(ds):
    h = 500
    return [ds.x.min() - h, ds.x.max() + h, ds.y.min() - h, ds.y.max() + h]


def decorate(ax, ds, towns=True):
    if towns:
        for name, (lo, la) in TOWNS.items():
            x, y = to_xy(lo, la)
            ax.plot(x, y, "o", ms=4, mfc="white", mec="black", mew=0.8)
            right = x > ds.x.max() - 12000
            ax.annotate(name, (x, y), xytext=(-4 if right else 4, 3), textcoords="offset points",
                        fontsize=7.5, ha="right" if right else "left",
                        bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.6, lw=0))
    x0, x1, y0, y1 = extent(ds)
    ax.plot([x1 - 12000, x1 - 2000], [y0 + 2500] * 2, color="black", lw=2.5)
    ax.text(x1 - 7000, y0 + 3300, "10 km", ha="center", fontsize=7.5)
    ax.set_xticks([])
    ax.set_yticks([])


def trend_maps(da):
    """Theil-Sen slope (per decade) and Mann-Kendall p-value for every grid cell."""
    years = da.time.dt.year.values.astype(float)
    arr = da.values.reshape(len(years), -1)
    slope = np.full(arr.shape[1], np.nan)
    pval = np.full(arr.shape[1], np.nan)
    for k in range(arr.shape[1]):
        y = arr[:, k]
        if np.isnan(y).any():
            continue
        slope[k] = stats.theilslopes(y, years)[0] * 10
        pval[k] = stats.kendalltau(years, y)[1]
    shape = da.shape[1:]
    return (xr.DataArray(slope.reshape(shape), coords={"y": da.y, "x": da.x}),
            xr.DataArray(pval.reshape(shape), coords={"y": da.y, "x": da.x}))


def nearest_cell(ds, lon, lat):
    x, y = to_xy(lon, lat)
    return ds.sel(x=x, y=y, method="nearest")


# ------------------------------------------------------------------ figures
def fig_climatology(ds, path):
    mean = ds.mean("time")
    panels = [("TM", "Mean temperature (°C)", "RdYlBu_r"),
              ("RR", "Annual precipitation (mm)", "YlGnBu"),
              ("SA", "Sunshine duration (h/yr)", "YlOrBr")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, (var, label, cmap) in zip(axes, panels):
        im = ax.imshow(mean[var], extent=extent(ds), cmap=cmap, origin="upper")
        cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cb.set_label(label)
        decorate(ax, ds)
    fig.suptitle(f"Climate {FIRST}–{LAST} (SPARTACUS 1 km)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_trends(ds, t_slope, t_p, r_slope_pct, r_p, path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    im = axes[0].imshow(t_slope, extent=extent(ds), cmap="Reds", origin="upper",
                        vmin=float(np.nanmin(t_slope)), vmax=float(np.nanmax(t_slope)))
    fig.colorbar(im, ax=axes[0], shrink=0.8, pad=0.02).set_label("Temperature trend (°C per decade)")
    lim = float(np.nanmax(np.abs(r_slope_pct)))
    im = axes[1].imshow(r_slope_pct, extent=extent(ds), cmap="BrBG",
                        norm=TwoSlopeNorm(0, -lim, lim), origin="upper")
    fig.colorbar(im, ax=axes[1], shrink=0.8, pad=0.02).set_label("Precipitation trend (% per decade)")
    for ax, p, title in [(axes[0], t_p, "Temperature"), (axes[1], r_p, "Precipitation")]:
        # stippling: cells where the trend is NOT significant
        yy, xx = np.meshgrid(ds.y, ds.x, indexing="ij")
        ns = (p.values >= 0.05)
        ax.scatter(xx[ns][::3], yy[ns][::3], s=1, color="black", alpha=0.35, lw=0)
        share = 100 * np.mean(p.values < 0.05)
        ax.set_title(f"{title} trend {FIRST}–{LAST}  ({share:.0f} % of cells significant)",
                     fontsize=10, fontweight="bold")
        decorate(ax, ds)
    fig.text(0.5, 0.01, "Theil-Sen slope; dots = not significant (Mann-Kendall, α = 0.05)",
             ha="center", fontsize=8, color="grey")
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_regional(reg, path):
    years = reg.index.values
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True,
                             gridspec_kw={"height_ratios": [0.45, 1, 1]})
    anom = reg.TM - reg.TM.mean()
    axes[0].bar(years, 1, width=1, color=plt.cm.RdBu_r(TwoSlopeNorm(0, -1.5, 1.5)(anom.values)))
    axes[0].set_yticks([])
    axes[0].set_title("Regional temperature anomaly (warming stripes)", fontsize=10, loc="left")
    for ax, col, color, unit in [(axes[1], "TM", "#c53030", "°C"), (axes[2], "RR", "#2b6cb0", "mm")]:
        a = reg[col] - reg[col].mean()
        ax.bar(years, a, color=np.where(a > 0, color, "#a0aec0") if col == "TM" else
               np.where(a > 0, color, "#dd6b20"), alpha=0.85)
        sl, ic, *_ = stats.theilslopes(reg[col], years)
        p = stats.kendalltau(years, reg[col])[1]
        ax.plot(years, ic + sl * years - reg[col].mean(), color="black", lw=1.5, ls="--")
        ax.axhline(0, color="black", lw=0.6)
        ax.set_ylabel(f"Anomaly ({unit})")
        ax.set_title(f"{'Mean temperature' if col == 'TM' else 'Precipitation'}: "
                     f"trend {10 * sl:+.2f} {unit}/decade ({'p < 0.001' if p < 0.001 else f'p = {p:.3f}'}), "
                     f"mean {reg[col].mean():.1f} {unit}",
                     fontsize=10, loc="left")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Regional mean {FIRST}–{LAST} (relative to 30-year mean)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_warming_vs_climate(mean_tm, t_slope, mean_rr, path):
    x, y, c = mean_tm.values.ravel(), t_slope.values.ravel(), mean_rr.values.ravel()
    rho, p = stats.spearmanr(x, y)
    fig, ax = plt.subplots(figsize=(7, 4.6))
    sc = ax.scatter(x, y, c=c, cmap="YlGnBu", s=8, alpha=0.8)
    fig.colorbar(sc, ax=ax).set_label("Mean precipitation (mm)")
    ax.set_xlabel(f"Mean temperature {FIRST}–{LAST} (°C)")
    ax.set_ylabel("Temperature trend (°C per decade)")
    ax.set_title(f"Warming rate vs. local climate (Spearman ρ = {rho:.2f})", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return rho, p


def fig_validation(ds, path):
    st = pd.read_csv(STATIONS_CSV)
    rows = []
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for col_i, (name, (lo, la)) in enumerate(STATIONS.items()):
        cell = nearest_cell(ds, lo, la).to_dataframe()[["TM", "RR"]]
        cell.index = cell.index.year
        s = st[st.station == name].set_index("year")
        both = cell.join(s[["tmean_c", "precip_mm"]], how="inner").dropna()
        for row_i, (g, o, unit) in enumerate([("TM", "tmean_c", "°C"), ("RR", "precip_mm", "mm")]):
            ax = axes[row_i, col_i]
            ax.plot(cell.index, cell[g], color="#2b6cb0", lw=1.8, label="SPARTACUS grid cell")
            ax.plot(s.index, s[o], "o-", color="#c53030", ms=3.5, lw=1, label="Station")
            bias = (both[g] - both[o]).mean()
            r = both[g].corr(both[o])
            ax.set_title(f"{name}: bias {bias:+.2f} {unit}, r = {r:.2f}", fontsize=10, loc="left")
            ax.set_ylabel("Mean temperature (°C)" if g == "TM" else "Precipitation (mm)")
            ax.spines[["top", "right"]].set_visible(False)
            rows.append({"station": name, "variable": g, "n_years": len(both),
                         "bias_grid_minus_station": round(bias, 2), "pearson_r": round(r, 3)})
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("SPARTACUS vs. GeoSphere station data (overlapping years)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------ main
def main():
    FIG.mkdir(exist_ok=True)
    ds = load()
    mean = ds.mean("time")

    t_slope, t_p = trend_maps(ds.TM)
    r_slope, r_p = trend_maps(ds.RR)
    r_slope_pct = 100 * r_slope / mean.RR

    reg = ds.mean(["x", "y"]).to_dataframe()[["TM", "RR", "SA"]]
    reg.index = reg.index.year
    reg.index.name = "year"

    fig_climatology(ds, FIG / "climatology_maps.png")
    fig_trends(ds, t_slope, t_p, r_slope_pct, r_p, FIG / "trend_maps.png")
    fig_regional(reg, FIG / "regional_timeseries.png")
    rho, p_rho = fig_warming_vs_climate(mean.TM, t_slope, mean.RR, FIG / "warming_vs_climate.png")
    val = fig_validation(ds, FIG / "station_validation.png")

    reg.round(2).to_csv(PROC / "regional_annual_means.csv")
    val.to_csv(PROC / "station_validation.csv", index=False)
    xr.Dataset({"TM_mean": mean.TM, "RR_mean": mean.RR, "SA_mean": mean.SA,
                "TM_trend_per_decade": t_slope, "TM_trend_p": t_p,
                "RR_trend_pct_per_decade": r_slope_pct, "RR_trend_p": r_p},
               attrs={"crs": "EPSG:3416", "period": f"{FIRST}-{LAST}"}
               ).to_netcdf(PROC / "climatology_and_trends.nc")

    yrs = reg.index.values
    summary = {}
    for col in ["TM", "RR", "SA"]:
        sl = stats.theilslopes(reg[col], yrs)[0] * 10
        summary[col] = {"mean": round(reg[col].mean(), 1),
                        "trend_per_decade": round(sl, 2),
                        "p_value": round(stats.kendalltau(yrs, reg[col])[1], 4),
                        "min_cell_mean": round(float(mean[col].min()), 1),
                        "max_cell_mean": round(float(mean[col].max()), 1)}
    summary["TM"]["cells_significant_pct"] = round(100 * float((t_p < 0.05).mean()), 1)
    summary["TM"]["cell_trend_range"] = f"{float(t_slope.min()):.2f} – {float(t_slope.max()):.2f}"
    summary["RR"]["cells_significant_pct"] = round(100 * float((r_p < 0.05).mean()), 1)
    summary["RR"]["cell_trend_range_pct"] = f"{float(r_slope_pct.min()):.1f} – {float(r_slope_pct.max()):.1f}"
    summ = pd.DataFrame(summary).T
    summ.to_csv(PROC / "summary.csv")

    pd.set_option("display.width", 160)
    print("=== Regional summary", FIRST, "–", LAST, "===\n", summ.to_string())
    print(f"\nWarming rate vs. mean temperature: Spearman ρ = {rho:.2f} (p = {p_rho:.3g})")
    print("\n=== Station validation ===\n", val.to_string(index=False))
    print("\nWarmest years:", reg.TM.nlargest(4).round(2).to_dict())
    print("Driest years:", reg.RR.nsmallest(3).round(0).to_dict(), "| wettest:", reg.RR.nlargest(3).round(0).to_dict())
    first, last = reg.TM.loc[FIRST:FIRST + 9].mean(), reg.TM.loc[LAST - 9:LAST].mean()
    print(f"Decade means: {FIRST}–{FIRST + 9}: {first:.2f} °C, {LAST - 9}–{LAST}: {last:.2f} °C (+{last - first:.2f})")


if __name__ == "__main__":
    main()
