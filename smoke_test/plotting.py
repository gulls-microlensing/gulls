"""Plotting helpers for smoke test lightcurve products."""
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from astropy.coordinates import SkyCoord
import astropy.units as u

from .constants import REPO_ROOT
from .errors import SmokeTestError

candidate = (REPO_ROOT.parent / "VBMicrolensing").resolve()
if candidate.is_dir():
    sys.path.append(str(candidate))

from VBMicrolensing import VBMicrolensing as VBMicrolensingClass  # type: ignore[attr-defined]

VBM_CLASS = VBMicrolensingClass  # type: ignore


def _derive_event_key(lc_file: Path) -> Tuple[int, int, int] | None:
    stem = lc_file.stem.split(".", 1)[0]
    parts = stem.rsplit("_", 3)
    if len(parts) < 4:
        return None
    return tuple(int(part) for part in parts[-3:])


def _format_metric(value: float | None, precision: int = 3) -> str:
    if value is None or math.isnan(value):
        return "n/a"
    return f"{value:.{precision}f}"


def _galactic_pm_to_icrs(l_deg: float, b_deg: float, mu_l: float, mu_b: float) -> Tuple[float, float]:
    coord = SkyCoord(
        l=l_deg * u.deg,
        b=b_deg * u.deg,
        pm_l_cosb=mu_l * u.mas / u.yr,
        pm_b=mu_b * u.mas / u.yr,
        frame="galactic",
    )
    icrs = coord.icrs
    return (
        icrs.pm_ra_cosdec.to_value(u.mas / u.yr),
        icrs.pm_dec.to_value(u.mas / u.yr),
    )


def _parse_header(lc_file: Path) -> Tuple[List[float] | None, List[float] | None]:
    planet_vals: List[float] | None = None
    event_vals: List[float] | None = None
    with lc_file.open(encoding="utf-8") as header_reader:
        for raw in header_reader:
            if not raw.startswith("#"):
                break
            stripped = raw.strip()
            if stripped.startswith("#Planet:"):
                planet_vals = [float(x) for x in stripped.split()[1:]]
            elif stripped.startswith("#Event:"):
                event_vals = [float(x) for x in stripped.split()[1:]]
    return planet_vals, event_vals


def _plot_photometry_only(
    lc_file: Path,
    output_dir: Path,
    title: str,
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray,
    true_flux: np.ndarray | None,
) -> Path:
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    fig.suptitle(title, fontsize=14)
    ax.errorbar(
        time,
        flux,
        yerr=flux_err,
        fmt="o",
        markersize=2,
        alpha=0.5,
        color="C0",
        label="Measured",
        zorder=1,
    )
    if true_flux is not None:
        ax.plot(
            time,
            true_flux,
            "-",
            linewidth=1.5,
            color="red",
            label="True",
            zorder=2,
            alpha=0.8,
        )
    ax.axhline(1.0, color="k", linestyle="--", linewidth=1.5, label="Baseline", zorder=3)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Relative Flux")
    ax.set_title("Light Curve")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plot_file = output_dir / f"{lc_file.stem}_plot.png"
    fig.tight_layout()
    fig.savefig(plot_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Generated plot: {plot_file.name}")
    return plot_file


def _compute_vbm_model(
    summary: Dict[str, float] | None,
    planet_vals: List[float] | None,
    event_vals: List[float] | None,
    source_pm_icrs: Tuple[float, float] | None,
    lens_pm_icrs: Tuple[float, float] | None,
    theta_e_float: float | None,
    source_dist_float: float | None,
    event_ra_float: float | None,
    event_dec_float: float | None,
    alpha_deg_float: float,
    sim_zero_offset: float,
    time: np.ndarray,
    true_x_vals: np.ndarray | None,
    true_y_vals: np.ndarray | None,
) -> Tuple[Dict[str, np.ndarray | str] | None, str | None]:
    if summary is None:
        return None, "summary metrics missing for this lightcurve"
    if not (
        planet_vals
        and len(planet_vals) >= 6
        and event_vals
        and len(event_vals) >= 8
    ):
        return None, "header missing #Planet/#Event metadata"

    q_val = float(planet_vals[4])
    s_val = float(planet_vals[5])
    if q_val <= 0 or s_val <= 0:
        return None, f"unphysical planet parameters (q={q_val}, s={s_val})"

    rho_val = summary.get("rho")
    if rho_val is None or math.isnan(rho_val):
        rho_val = float(event_vals[7])
    tE_val = summary.get("tE_ref")
    if tE_val is None or math.isnan(tE_val):
        tE_val = float(event_vals[6])
    u0_val = summary.get("u0")
    if u0_val is None or math.isnan(u0_val):
        u0_val = float(event_vals[0])
    alpha_deg = summary.get("alpha_event")
    if alpha_deg is None or math.isnan(alpha_deg):
        alpha_deg = float(event_vals[1])
    t0_val = summary.get("t0")
    if t0_val is None or math.isnan(t0_val):
        t0_val = float(event_vals[2])
    pi_n_val = summary.get("pi_n")
    pi_e_val = summary.get("pi_e")
    if pi_n_val is None or math.isnan(pi_n_val):
        pi_n_val = 0.0
    if pi_e_val is None or math.isnan(pi_e_val):
        pi_e_val = 0.0

    if (
        theta_e_float is None
        or theta_e_float <= 0
        or source_dist_float is None
        or source_dist_float <= 0
        or source_pm_icrs is None
        or event_ra_float is None
        or event_dec_float is None
    ):
        return None, "insufficient astrometric metadata (theta_E, source distance, or PM)"

    pi_s_val = 1.0 / source_dist_float
    if (
        pi_s_val <= 0
        or rho_val is None
        or float(rho_val) <= 0
        or tE_val is None
        or float(tE_val) <= 0
    ):
        return None, "missing positive rho/tE/source distance for VBM evaluation"

    vbm = VBM_CLASS()  # type: ignore[operator]
    skycoord = SkyCoord(
        ra=float(event_ra_float) * u.deg,
        dec=float(event_dec_float) * u.deg,
    )
    coord_str = (
        f"{skycoord.ra.to_string(unit=u.hour, sep=':', pad=True)} "
        f"{skycoord.dec.to_string(unit=u.deg, sep=':', pad=True, alwayssign=True)}"
    )
    vbm.SetObjectCoordinates(coord_str)
    params_vbm = [
        math.log(float(s_val)),
        math.log(float(q_val)),
        float(u0_val),
        math.radians(float(alpha_deg)),
        math.log(float(rho_val)),
        math.log(float(tE_val)),
        float(t0_val) + sim_zero_offset,
        float(pi_n_val),
        float(pi_e_val),
        float(source_pm_icrs[1]),
        float(source_pm_icrs[0]),
        float(pi_s_val),
        float(theta_e_float),
    ]
    results = vbm.BinaryAstroLightCurve(params_vbm, time + sim_zero_offset)

    lens_dec_deg = np.array(results[3], dtype=float)
    lens_ra_deg = np.array(results[4], dtype=float)
    y1 = np.array(results[5], dtype=float)
    y2 = np.array(results[6], dtype=float)
    lensframe_label = "Source Trajectory BinaryAstroLightCurve"

    vbm_x = y1
    vbm_y = y2
    if (
        true_x_vals is not None
        and true_y_vals is not None
        and len(true_x_vals) == len(y1)
    ):
        combos = [
            ("x= y1, y= y2", y1, y2),
            ("x=-y1, y= y2", -y1, y2),
            ("x= y1, y=-y2", y1, -y2),
            ("x=-y1, y=-y2", -y1, -y2),
            ("x= y2, y= y1", y2, y1),
            ("x=-y2, y= y1", -y2, y1),
            ("x= y2, y=-y1", y2, -y1),
            ("x=-y2, y=-y1", -y2, -y1),
        ]
        best = None
        best_err = None
        for label, cand_x, cand_y in combos:
            err = np.nanmean((cand_x - true_x_vals) ** 2 + (cand_y - true_y_vals) ** 2)
            if best_err is None or err < best_err:
                best_err = err
                best = (label, cand_x, cand_y)
        if best:
            lensframe_label = f"Source trajectory BinaryAstroLightCurve ({best[0]})"
            vbm_x, vbm_y = best[1], best[2]

    if len(vbm_x) != len(time):
        return None, "VBM returned mismatched array lengths"

    model: Dict[str, np.ndarray | str] = {
        "lens_x": vbm_x,
        "lens_y": vbm_y,
        "lens_label": lensframe_label,
        "sky_ra": lens_ra_deg,
        "sky_dec": lens_dec_deg,
        "sky_label": "VBM BinaryAstroLightCurve (sky)",
    }
    return model, None


def _render_lensframe(
    lc_file: Path,
    output_dir: Path,
    time: np.ndarray,
    cmap: matplotlib.colors.Colormap,
    norm: Normalize,
    vbm_model: Dict[str, np.ndarray | str],
    true_x_vals: np.ndarray | None,
    true_y_vals: np.ndarray | None,
    meas_x: np.ndarray | None,
    meas_y: np.ndarray | None,
) -> Path:
    vbm_x = np.asarray(vbm_model["lens_x"])  # type: ignore[index]
    vbm_y = np.asarray(vbm_model["lens_y"])  # type: ignore[index]
    vbm_label = str(vbm_model["lens_label"])

    fig2, ax2 = plt.subplots(figsize=(6, 6))

    if meas_x is not None and meas_y is not None:
        mask_meas = np.isfinite(meas_x) & np.isfinite(meas_y)
        if np.any(mask_meas):
            ax2.scatter(
                meas_x[mask_meas],
                meas_y[mask_meas],
                c=time[mask_meas],
                cmap=cmap,
                norm=norm,
                s=25,
                alpha=0.6,
                label="Measured centroid",
                zorder=1,
            )
    else:
        ax2.text(
            0.02,
            0.98,
            "Measured centroid not available",
            ha="left",
            va="top",
            transform=ax2.transAxes,
            fontsize=8,
        )

    if true_x_vals is not None and true_y_vals is not None:
        mask_true = np.isfinite(true_x_vals) & np.isfinite(true_y_vals)
        if np.any(mask_true):
            ax2.scatter(
                true_x_vals[mask_true],
                true_y_vals[mask_true],
                c=time[mask_true],
                cmap=cmap,
                norm=norm,
                s=20,
                marker="x",
                linewidths=0.9,
                alpha=1.0,
                label="True centroid",
                zorder=2,
            )

    ax2.plot(
        vbm_x,
        vbm_y,
        color="black",
        linewidth=1.5,
        label=vbm_label,
        zorder=3,
    )
    ax2.set_xlabel("x_centroid (Einstein radii)")
    ax2.set_ylabel("y_centroid (Einstein radii)")
    ax2.set_title(f"Lens-frame Centroid: {lc_file.stem}")
    ax2.grid(True, alpha=0.3)

    segments_x = [vbm_x]
    segments_y = [vbm_y]
    if meas_x is not None and meas_y is not None:
        mask_meas = np.isfinite(meas_x) & np.isfinite(meas_y)
        if np.any(mask_meas):
            segments_x.append(np.asarray(meas_x)[mask_meas])
            segments_y.append(np.asarray(meas_y)[mask_meas])
    if true_x_vals is not None and true_y_vals is not None:
        mask_true = np.isfinite(true_x_vals) & np.isfinite(true_y_vals)
        if np.any(mask_true):
            segments_x.append(true_x_vals[mask_true])
            segments_y.append(true_y_vals[mask_true])
    all_x = np.concatenate(segments_x) if segments_x else np.array([0.0])
    all_y = np.concatenate(segments_y) if segments_y else np.array([0.0])
    x_min = float(np.nanmin(all_x))
    x_max = float(np.nanmax(all_x))
    y_min = float(np.nanmin(all_y))
    y_max = float(np.nanmax(all_y))
    x_c = 0.5 * (x_min + x_max)
    y_c = 0.5 * (y_min + y_max)
    half_span = max(x_max - x_min, y_max - y_min) * 0.5
    half_span = max(half_span, 1e-6)
    half_span *= 1.15
    ax2.set_xlim(x_c - half_span, x_c + half_span)
    ax2.set_ylim(y_c - half_span, y_c + half_span)
    ax2.set_aspect("equal", adjustable="box")
    ax2.legend(loc="upper left")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar2 = fig2.colorbar(sm, ax=ax2, fraction=0.046, pad=0.04)
    cbar2.set_label("Time (days)")

    lensframe_file = output_dir / f"{lc_file.stem}_lensframe_plot.png"
    fig2.tight_layout()
    fig2.savefig(lensframe_file, dpi=150, bbox_inches="tight")
    plt.close(fig2)
    return lensframe_file


def _render_astrometric_figure(
    lc_file: Path,
    output_dir: Path,
    title: str,
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray,
    true_flux: np.ndarray | None,
    true_N_mas: np.ndarray,
    true_E_mas: np.ndarray,
    meas_N_mas: np.ndarray,
    meas_E_mas: np.ndarray,
    meas_N_err_mas: np.ndarray,
    meas_E_err_mas: np.ndarray,
    true_ra_deg: np.ndarray,
    true_dec_deg: np.ndarray,
    meas_ra_deg: np.ndarray,
    meas_dec_deg: np.ndarray,
    meas_ra_err_deg: np.ndarray,
    meas_dec_err_deg: np.ndarray,
    vector_specs: List[Dict[str, float | str]],
    vbm_model: Dict[str, np.ndarray | str] | None,
    true_x_vals: np.ndarray | None,
    true_y_vals: np.ndarray | None,
    meas_x: np.ndarray | None,
    meas_y: np.ndarray | None,
) -> Tuple[Path, Path | None]:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(title, fontsize=14)
    ax_light = axes[0, 0]
    ax_time = axes[0, 1]
    ax_radec = axes[1, 0]
    ax_ne = axes[1, 1]

    ax_light.errorbar(
        time,
        flux,
        yerr=flux_err,
        fmt="o",
        markersize=2,
        alpha=0.5,
        color="C0",
        label="Measured",
        zorder=1,
    )
    if true_flux is not None:
        ax_light.plot(
            time,
            true_flux,
            "-",
            linewidth=1.5,
            color="red",
            label="True",
            zorder=2,
            alpha=0.8,
        )
    ax_light.axhline(
        1.0,
        color="k",
        linestyle="--",
        linewidth=1.5,
        label="Baseline",
        zorder=3,
    )
    ax_light.set_xlabel("Time (days)")
    ax_light.set_ylabel("Relative Flux")
    ax_light.set_title("Light Curve")
    ax_light.legend()
    ax_light.grid(True, alpha=0.3)

    norm = Normalize(vmin=np.min(time), vmax=np.max(time)) if len(time) else Normalize(0, 1)
    cmap = plt.get_cmap("plasma")

    span_years: float | None = None
    if len(time):
        span_days = float(time.max() - time.min())
        if span_days > 0:
            span_years = span_days / 365.25

    ax_radec.errorbar(
        meas_ra_deg,
        meas_dec_deg,
        xerr=meas_ra_err_deg,
        yerr=meas_dec_err_deg,
        fmt="none",
        ecolor="lightgray",
        alpha=0.5,
        capsize=2,
        zorder=0,
    )
    sc_ra = ax_radec.scatter(
        meas_ra_deg,
        meas_dec_deg,
        c=time,
        cmap=cmap,
        norm=norm,
        s=25,
        alpha=0.5,
        label="Measured",
        zorder=1,
    )
    ax_radec.plot(
        true_ra_deg,
        true_dec_deg,
        color="black",
        linewidth=1.2,
        alpha=0.8,
        label="True track",
        zorder=4,
    )
    ax_radec.scatter(
        true_ra_deg,
        true_dec_deg,
        c=time,
        cmap=cmap,
        norm=norm,
        s=18,
        marker="x",
        linewidths=0.8,
        alpha=1.0,
        label="True samples",
        zorder=3,
    )
    if vbm_model is not None:
        deg_to_mas = 3600.0 * 1000.0
        baseline_ra = true_ra_deg[0]
        baseline_dec = true_dec_deg[0]
        cos_dec0 = math.cos(math.radians(baseline_dec))
        if abs(cos_dec0) < 1e-6:
            cos_dec0 = 1e-6 if cos_dec0 >= 0 else -1e-6

        vbm_ra_offset_mas = np.asarray(vbm_model["sky_ra"], dtype=float)  # type: ignore[index]
        vbm_dec_offset_mas = np.asarray(vbm_model["sky_dec"], dtype=float)  # type: ignore[index]
        vbm_ra_abs = baseline_ra + (vbm_ra_offset_mas / (deg_to_mas * cos_dec0))
        vbm_dec_abs = baseline_dec + (vbm_dec_offset_mas / deg_to_mas)

        ax_radec.plot(
            vbm_ra_abs,
            vbm_dec_abs,
            color="tab:purple",
            linewidth=1.2,
            alpha=0.9,
            label=vbm_model.get("sky_label", "VBM BinaryAstroLightCurve (sky)"),
        )

    ax_radec.set_xlabel("RA (degrees)")
    ax_radec.set_ylabel("Dec (degrees)")
    ax_radec.set_title("Absolute Astrometric Position")
    ax_radec.grid(True, alpha=0.3)
    ax_radec.axis("equal")

    if span_years and vector_specs:
        start_ra = true_ra_deg[0]
        start_dec = true_dec_deg[0]
        cos_dec = math.cos(math.radians(start_dec))
        if abs(cos_dec) < 1e-6:
            cos_dec = 1e-6 if cos_dec >= 0 else -1e-6
        for spec in vector_specs:
            pm_ra = spec["pm_ra"]
            pm_dec = spec["pm_dec"]
            delta_ra_deg = (pm_ra * span_years) / (3600000.0 * cos_dec)
            delta_dec_deg = (pm_dec * span_years) / 3600000.0
            end_ra = start_ra + delta_ra_deg
            end_dec = start_dec + delta_dec_deg
            ax_radec.annotate(
                "",
                xy=(end_ra, end_dec),
                xytext=(start_ra, start_dec),
                arrowprops=dict(color=spec["color"], arrowstyle="->", linewidth=1),
                zorder=5,
            )
            ax_radec.plot([], [], color=spec["color"], linewidth=2, label=spec["label"])
    ax_radec.legend(ncol=2, fontsize=8)

    ax_ne.plot(
        true_E_mas,
        true_N_mas,
        color="black",
        linewidth=1.2,
        alpha=0.8,
        label="True track",
        zorder=4,
    )
    ax_ne.scatter(
        true_E_mas,
        true_N_mas,
        c=time,
        cmap=cmap,
        norm=norm,
        s=18,
        marker="x",
        linewidths=0.8,
        alpha=1.0,
        label="True samples",
        zorder=3,
    )
    ax_ne.errorbar(
        meas_E_mas,
        meas_N_mas,
        xerr=meas_E_err_mas,
        yerr=meas_N_err_mas,
        fmt="none",
        ecolor="lightgray",
        alpha=0.5,
        capsize=2,
        zorder=0,
    )
    ax_ne.scatter(
        meas_E_mas,
        meas_N_mas,
        c=time,
        cmap=cmap,
        norm=norm,
        s=25,
        alpha=0.5,
        label="Measured",
        zorder=1,
    )
    if vbm_model is not None:
        vbm_E_mas = np.asarray(vbm_model["sky_ra"], dtype=float)  # type: ignore[index]
        vbm_N_mas = np.asarray(vbm_model["sky_dec"], dtype=float)  # type: ignore[index]
        ax_ne.plot(
            vbm_E_mas,
            vbm_N_mas,
            color="tab:purple",
            linewidth=1.2,
            alpha=0.9,
            label=vbm_model.get("sky_label", "VBM BinaryAstroLightCurve (sky)"),
        )
    ax_ne.set_xlabel("ΔEast (mas)")
    ax_ne.set_ylabel("ΔNorth (mas)")
    ax_ne.set_title("Astrometric Centroid (N/E), Relative to the Lens")
    ax_ne.grid(True, alpha=0.3)
    ax_ne.axis("equal")

    if span_years and vector_specs:
        start_E = true_E_mas[0]
        start_N = true_N_mas[0]
        for spec in vector_specs:
            pm_ra = spec["pm_ra"]
            pm_dec = spec["pm_dec"]
            delta_E_mas = pm_ra * span_years
            delta_N_mas = pm_dec * span_years
            end_E = start_E + delta_E_mas
            end_N = start_N + delta_N_mas
            ax_ne.annotate(
                "",
                xy=(end_E, end_N),
                xytext=(start_E, start_N),
                arrowprops=dict(color=spec["color"], arrowstyle="->", linewidth=1),
                zorder=5,
            )
            ax_ne.plot([], [], color=spec["color"], linewidth=2, label=spec["label"])
    ax_ne.legend(ncol=2, fontsize=8)

    fig.tight_layout(rect=[0, 0.12, 1, 1])
    cbar_ax = fig.add_axes([0.25, 0.06, 0.5, 0.025])
    cbar = fig.colorbar(sc_ra, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Time (days)")

    ax_time.plot(time, true_N_mas, "b-", label="True N", linewidth=2, alpha=0.7)
    ax_time.plot(time, meas_N_mas, "r.", label="Meas N", markersize=2, alpha=0.5)
    ax_time.plot(time, true_E_mas, "g-", label="True E", linewidth=2, alpha=0.7)
    ax_time.plot(time, meas_E_mas, "m.", label="Meas E", markersize=2, alpha=0.5)
    ax_time.set_xlabel("Time (days)")
    ax_time.set_ylabel("Centroid Shift (mas)")
    ax_time.set_title("Astrometric Timeseries (N and E)")
    ax_time.legend(fontsize=8)
    ax_time.grid(True, alpha=0.3)

    plot_file = output_dir / f"{lc_file.stem}_plot.png"
    fig.savefig(plot_file, dpi=150, bbox_inches="tight")
    plt.close(fig)

    lensframe_path: Path | None = None
    if vbm_model is not None:
        lensframe_path = _render_lensframe(
            lc_file,
            output_dir,
            time,
            cmap,
            norm,
            vbm_model,
            true_x_vals,
            true_y_vals,
            meas_x,
            meas_y,
        )

    return plot_file, lensframe_path


def plot_lightcurves(
    output_dir: Path,
    summaries: Dict[Tuple[int, int, int], Dict[str, float]] | None = None,
    params: Dict[str, str] | None = None,
) -> None:
    lc_files = sorted(output_dir.rglob("*.lc"))
    if not lc_files:
        return

    sim_zero_offset = 0.0
    if params:
        sz = params.get("SIMULATION_ZERO_TIME")
        if sz:
            sim_zero_offset = float(sz) - 2450000.0

    astrometry_expected = False
    if params:
        val = params.get("ASTROMETRY_ON")
        if val is not None:
            astrometry_expected = str(val).strip().lower() not in {"0", "false", "off"}

    for lc_file in lc_files:
        planet_vals, event_vals = _parse_header(lc_file)
        df = pd.read_csv(lc_file, sep=r"\s+", comment="#")
        if df.empty:
            continue

        column_names = list(df.columns)

        def _require_column(name: str) -> np.ndarray:
            if name not in column_names:
                raise SmokeTestError(f"Smoke test failed: column '{name}' missing in {lc_file.name}")
            return df[name].to_numpy(dtype=float, copy=False)

        def _optional_column(name: str) -> np.ndarray | None:
            if name not in column_names:
                return None
            return df[name].to_numpy(dtype=float, copy=False)

        summary = None
        if summaries:
            event_key = _derive_event_key(lc_file)
            if event_key is not None:
                summary = summaries.get(event_key)
        if summary is None:
            raise SmokeTestError(f"Smoke test failed: summary metrics missing for {lc_file.name}")

        title = f"Smoke Test: {lc_file.stem}"
        theta_e_float: float | None = None
        alpha_deg_float: float | None = None
        pi_n_float: float | None = None
        pi_e_float: float | None = None
        source_dist_float: float | None = None
        event_ra_float: float | None = None
        event_dec_float: float | None = None
        if summary:
            lens_mass = _format_metric(summary.get("lens_mass"))
            lens_dist = _format_metric(summary.get("lens_dist"))
            source_dist = _format_metric(summary.get("source_dist"))
            theta_e = _format_metric(summary.get("theta_e"))
            pm_alpha = _format_metric(summary.get("pm_helio_alpha"))
            pm_delta = _format_metric(summary.get("pm_helio_delta"))
            subtitle = (
                f"Lens M={lens_mass} Msun, Lens D={lens_dist} pc, "
                f"Source D={source_dist} pc, theta_E={theta_e}, "
                f"mu_rel=({pm_alpha}, {pm_delta}) mas/yr"
            )
            title = f"{title}\n{subtitle}"
            val = summary.get("theta_e")
            if val is not None and not math.isnan(val):
                theta_e_float = float(val)
            val = summary.get("alpha_event")
            if val is not None and not math.isnan(val):
                alpha_deg_float = float(val)
            val = summary.get("pi_n")
            if val is not None and not math.isnan(val):
                pi_n_float = float(val)
            val = summary.get("pi_e")
            if val is not None and not math.isnan(val):
                pi_e_float = float(val)
            val = summary.get("source_dist")
            if val is not None and not math.isnan(val) and val > 0:
                source_dist_float = float(val)
            val = summary.get("event_ra")
            if val is not None and not math.isnan(val):
                event_ra_float = float(val)
            val = summary.get("event_dec")
            if val is not None and not math.isnan(val):
                event_dec_float = float(val)
        if alpha_deg_float is None and event_vals and len(event_vals) >= 2:
            alpha_deg_float = float(event_vals[1])
        if alpha_deg_float is None:
            raise SmokeTestError(f"Smoke test failed: missing alpha_event for {lc_file.name}")

        time = _require_column("Simulation_time")
        flux = _require_column("measured_relative_flux")
        flux_err = _require_column("measured_relative_flux_error")
        true_flux = _require_column("true_relative_flux")

        astrom_cols = [
            "true_N_centroid_mas",
            "true_E_centroid_mas",
            "measured_N_centroid_mas",
            "measured_E_centroid_mas",
            "measured_N_centroid_error_mas",
            "measured_E_centroid_error_mas",
            "true_centroid_ra_deg",
            "true_centroid_dec_deg",
            "measured_centroid_ra_deg",
            "measured_centroid_dec_deg",
            "measured_centroid_ra_error_deg",
            "measured_centroid_dec_error_deg",
        ]
        missing_astrom_cols = [col for col in astrom_cols if col not in column_names]
        has_astrom = not missing_astrom_cols
        if astrometry_expected and missing_astrom_cols:
            print(
                f"  Debug: {lc_file.name} missing astrometry columns: "
                + ", ".join(missing_astrom_cols)
            )
        if astrometry_expected and not has_astrom:
            raise SmokeTestError(
                f"Smoke test failed: astrometric columns missing in {lc_file.name}"
            )

        if not has_astrom:
            _plot_photometry_only(lc_file, output_dir, title, time, flux, flux_err, true_flux)
            continue

        true_N_mas = _require_column("true_N_centroid_mas")
        true_E_mas = _require_column("true_E_centroid_mas")
        meas_N_mas = _require_column("measured_N_centroid_mas")
        meas_E_mas = _require_column("measured_E_centroid_mas")
        meas_N_err_mas = _require_column("measured_N_centroid_error_mas")
        meas_E_err_mas = _require_column("measured_E_centroid_error_mas")
        true_ra_deg = _require_column("true_centroid_ra_deg")
        true_dec_deg = _require_column("true_centroid_dec_deg")
        meas_ra_deg = _require_column("measured_centroid_ra_deg")
        meas_dec_deg = _require_column("measured_centroid_dec_deg")
        meas_ra_err_deg = _require_column("measured_centroid_ra_error_deg")
        meas_dec_err_deg = _require_column("measured_centroid_dec_error_deg")
        true_x_vals = _optional_column("true_x_centroid")
        true_y_vals = _optional_column("true_y_centroid")
        meas_x = _optional_column("x_centroid")
        meas_y = _optional_column("y_centroid")

        pm_ref_alpha_float = None
        pm_ref_delta_float = None
        pm_ref_alpha_val = summary.get("pm_ref_alpha")
        pm_ref_delta_val = summary.get("pm_ref_delta")
        if pm_ref_alpha_val is not None and pm_ref_delta_val is not None:
            pm_ref_alpha_float = float(pm_ref_alpha_val)
            pm_ref_delta_float = float(pm_ref_delta_val)
            if math.isnan(pm_ref_alpha_float) or math.isnan(pm_ref_delta_float):
                pm_ref_alpha_float = None
                pm_ref_delta_float = None

        source_pm_icrs: Tuple[float, float] | None = None
        lens_pm_icrs: Tuple[float, float] | None = None
        if summary:
            src_vals = (
                summary.get("source_mul"),
                summary.get("source_mub"),
                summary.get("source_l"),
                summary.get("source_b"),
            )
            if all(v is not None and not math.isnan(v) for v in src_vals):
                source_pm_icrs = _galactic_pm_to_icrs(
                    float(src_vals[2]),
                    float(src_vals[3]),
                    float(src_vals[0]),
                    float(src_vals[1]),
                )
            lens_vals = (
                summary.get("lens_mul"),
                summary.get("lens_mub"),
                summary.get("lens_l"),
                summary.get("lens_b"),
            )
            if all(v is not None and not math.isnan(v) for v in lens_vals):
                lens_pm_icrs = _galactic_pm_to_icrs(
                    float(lens_vals[2]),
                    float(lens_vals[3]),
                    float(lens_vals[0]),
                    float(lens_vals[1]),
                )

        span_years = None
        if len(time):
            span_days = float(time.max() - time.min())
            if span_days > 0:
                span_years = span_days / 365.25

        vector_specs: List[Dict[str, float | str]] = []
        if span_years:
            if pm_ref_alpha_float is not None and pm_ref_delta_float is not None:
                vector_specs.append(
                    {
                        "label": "Relative proper motion (geocentric)",
                        "color": "black",
                        "pm_ra": pm_ref_alpha_float,
                        "pm_dec": pm_ref_delta_float,
                    }
                )
            if source_pm_icrs:
                vector_specs.append(
                    {
                        "label": "Source proper motion (heliocentric)",
                        "color": "tab:blue",
                        "pm_ra": source_pm_icrs[0],
                        "pm_dec": source_pm_icrs[1],
                    }
                )
            if lens_pm_icrs:
                vector_specs.append(
                    {
                        "label": "Lens proper motion (heliocentric)",
                        "color": "tab:red",
                        "pm_ra": lens_pm_icrs[0],
                        "pm_dec": lens_pm_icrs[1],
                    }
                )

        vbm_model, vbm_reason = _compute_vbm_model(
            summary,
            planet_vals,
            event_vals,
            source_pm_icrs,
            lens_pm_icrs,
            theta_e_float,
            source_dist_float,
            event_ra_float,
            event_dec_float,
            alpha_deg_float,
            sim_zero_offset,
            time,
            true_x_vals,
            true_y_vals,
        )
        if astrometry_expected and vbm_model is None:
            detail = f" ({vbm_reason})" if vbm_reason else ""
            raise SmokeTestError(
                f"Smoke test failed: missing VBM lens-frame plot for {lc_file.name}{detail}"
            )

        plot_file, lensframe_path = _render_astrometric_figure(
            lc_file,
            output_dir,
            title,
            time,
            flux,
            flux_err,
            true_flux,
            true_N_mas,
            true_E_mas,
            meas_N_mas,
            meas_E_mas,
            meas_N_err_mas,
            meas_E_err_mas,
            true_ra_deg,
            true_dec_deg,
            meas_ra_deg,
            meas_dec_deg,
            meas_ra_err_deg,
            meas_dec_err_deg,
            vector_specs,
            vbm_model,
            true_x_vals,
            true_y_vals,
            meas_x,
            meas_y,
        )

        if astrometry_expected:
            if lensframe_path is None or not lensframe_path.exists():
                raise SmokeTestError(
                    f"Smoke test failed: missing lens-frame plot for {lc_file.name}"
                )
        print(f"  Generated plot: {plot_file.name}")
        if lensframe_path is not None:
            print(f"  Generated plot: {lensframe_path.name}")

__all__ = ["plot_lightcurves"]
