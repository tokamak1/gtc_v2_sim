#!/usr/bin/env python3
"""Generic plotting utilities for GTC history.out and snap*.out files.

Examples:
  python3 gtc_plot.py case --case-dir run_dir --out-dir run_dir/plots

  python3 gtc_plot.py case \
    --input run_dir/gtc.input \
    --history run_dir/history.out \
    --snap-glob 'run_dir/snap*.out' \
    --out-dir run_dir/plots

  python3 gtc_plot.py compare-history \
    --case-dir run_no_zonal \
    --case-dir run_zonal \
    --label "no zonal" \
    --label "zonal" \
    --out-dir compare_plots

  python3 gtc_plot.py compare-history \
    --input case1/gtc.input --history case1/history.out \
    --input case2/gtc.input --history case2/history.out \
    --label case1 --label case2 \
    --out-dir compare_plots
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

mpl_config = Path(os.environ.get("MPLCONFIGDIR", Path.cwd() / ".matplotlib")).expanduser()
if not mpl_config.exists() or not os.access(mpl_config, os.W_OK):
    mpl_config = Path.cwd() / ".matplotlib"
    mpl_config.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(mpl_config)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


TRACER_NAMES = [
    "tracer_x_minus_1",
    "tracer_z",
    "tracer_theta",
    "tracer_zeta",
    "tracer_weight",
    "tracer_v_parallel_norm",
    "tracer_energy_error",
    "tracer_momentum_error",
]

PHYS_NAMES = [
    "ddeni",
    "eradial",
    "efield_rms",
    "entropy_i",
    "dflowi_over_vthi",
    "pfluxi_over_vthi",
    "efluxi_norm",
    "heat_flux_bin_1",
    "heat_flux_bin_2",
    "heat_flux_bin_3",
    "heat_flux_bin_4",
    "heat_flux_bin_5",
]

PHYS_LABELS = {
    "ddeni": r"$\delta n_i$",
    "eradial": r"$E_r$",
    "efield_rms": r"$|\delta E|_{\mathrm{rms}}$",
    "entropy_i": r"$S_i$",
    "dflowi_over_vthi": r"$\delta u_{\parallel i}/v_{thi}$",
    "pfluxi_over_vthi": r"$\Gamma_i/v_{thi}$",
    "efluxi_norm": r"$Q_i$",
    "heat_flux_bin_1": r"$\chi_i(r_1)$",
    "heat_flux_bin_2": r"$\chi_i(r_2)$",
    "heat_flux_bin_3": r"$\chi_i(r_3)$",
    "heat_flux_bin_4": r"$\chi_i(r_4)$",
    "heat_flux_bin_5": r"$\chi_i(r_5)$",
}

TRACER_LABELS = {
    "tracer_x_minus_1": r"$X_{\mathrm{tr}}-1$",
    "tracer_z": r"$Z_{\mathrm{tr}}$",
    "tracer_theta": r"$\theta_{\mathrm{tr}}$",
    "tracer_zeta": r"$\zeta_{\mathrm{tr}}$",
    "tracer_weight": r"$w_{\mathrm{tr}}$",
    "tracer_v_parallel_norm": r"$v_{\parallel,\mathrm{tr}}/v_{thi}$",
    "tracer_energy_error": r"$\Delta E_{\mathrm{tr}}/E_0$",
    "tracer_momentum_error": r"$\Delta P_{\zeta,\mathrm{tr}}$",
}

RADIAL_NAMES = [
    "zonali",
    "phip00_over_rho",
    "marker",
    "fflows",
    "dflows",
    "ftem",
    "dtem",
]

RADIAL_LABELS = {
    "zonali": r"$\delta n_{i,00}$",
    "phip00_over_rho": r"$E_{r,00}/\rho_i$",
    "marker": r"$N_{\mathrm{marker}}$",
    "fflows": r"$f_{\mathrm{flow}}$",
    "dflows": r"$\delta u_{\parallel i}$",
    "ftem": r"$T_i$",
    "dtem": r"$\delta T_i$",
}


@dataclass(frozen=True)
class GtcInput:
    path: Path | None
    text: str

    @classmethod
    def from_path(cls, path: Path | None) -> "GtcInput":
        if path is None or not path.exists():
            return cls(path, "")
        return cls(path.resolve(), path.read_text())

    def real(self, name: str, default: float) -> float:
        match = re.search(
            rf"\b{name}\s*=\s*([-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?)",
            self.text,
            flags=re.IGNORECASE,
        )
        if match is None:
            return default
        return float(match.group(1).replace("D", "E").replace("d", "E"))

    def integer(self, name: str, default: int) -> int:
        match = re.search(rf"\b{name}\s*=\s*([-+]?\d+)", self.text, flags=re.IGNORECASE)
        return int(match.group(1)) if match else default

    def int_list(self, name: str) -> list[int]:
        match = re.search(rf"\b{name}\s*=\s*([^/\n]+)", self.text, flags=re.IGNORECASE)
        if match is None:
            return []
        return [int(value) for value in re.findall(r"[-+]?\d+", match.group(1))]


@dataclass(frozen=True)
class CaseSpec:
    label: str
    case_dir: Path | None
    input_path: Path | None
    history_path: Path | None
    snapshot_paths: tuple[Path, ...]

    @property
    def gtc_input(self) -> GtcInput:
        return GtcInput.from_path(self.input_path)


def resolve_path(path: str | Path | None) -> Path | None:
    if path is None:
        return None
    return Path(path).expanduser().resolve()


def snapshot_step(path: Path) -> int:
    match = re.search(r"(\d+)(?!.*\d)", path.stem)
    return int(match.group(1)) if match else 0


def collect_snapshot_paths(case_dir: Path | None, globs: list[str], explicit: list[str]) -> tuple[Path, ...]:
    paths: list[Path] = []
    if explicit:
        paths.extend(resolve_path(path) for path in explicit if path)
    if globs:
        for pattern in globs:
            paths.extend(resolve_path(path) for path in glob.glob(pattern))
    elif case_dir is not None:
        paths.extend(case_dir.glob("snap*.out"))
    existing = {path.resolve() for path in paths if path is not None and path.exists()}
    return tuple(sorted(existing, key=snapshot_step))


def case_from_dir(case_dir: Path, label: str | None = None) -> CaseSpec:
    case_dir = case_dir.resolve()
    return CaseSpec(
        label=label or case_dir.name,
        case_dir=case_dir,
        input_path=case_dir / "gtc.input",
        history_path=case_dir / "history.out",
        snapshot_paths=tuple(sorted(case_dir.glob("snap*.out"), key=snapshot_step)),
    )


def case_from_paths(
    *,
    label: str | None,
    case_dir: Path | None,
    input_path: Path | None,
    history_path: Path | None,
    snapshot_paths: tuple[Path, ...],
) -> CaseSpec:
    inferred_dir = case_dir
    if inferred_dir is None:
        for path in (history_path, input_path):
            if path is not None:
                inferred_dir = path.parent
                break
    inferred_label = label or (inferred_dir.name if inferred_dir is not None else "case")
    return CaseSpec(
        label=inferred_label,
        case_dir=inferred_dir,
        input_path=input_path,
        history_path=history_path,
        snapshot_paths=snapshot_paths,
    )


def parse_case_arg(spec: str) -> CaseSpec:
    if "=" in spec:
        label, value = spec.split("=", 1)
        label = label.strip() or None
    else:
        label, value = None, spec
    value = value.strip()
    parts = value.split(":")
    if len(parts) == 1:
        return case_from_dir(Path(parts[0]).expanduser(), label)
    if len(parts) == 2:
        input_path = resolve_path(parts[0])
        history_path = resolve_path(parts[1])
        return case_from_paths(
            label=label,
            case_dir=history_path.parent if history_path is not None else None,
            input_path=input_path,
            history_path=history_path,
            snapshot_paths=(),
        )
    raise ValueError(
        "--case expects LABEL=CASE_DIR, CASE_DIR, LABEL=GTC_INPUT:HISTORY, or GTC_INPUT:HISTORY"
    )


def mode_numbers(gtc_input: GtcInput) -> tuple[list[int], list[int]]:
    return gtc_input.int_list("nmode"), gtc_input.int_list("mmode")


def mode_title(gtc_input: GtcInput, slot: int, field_symbol: str) -> str:
    nmode, mmode = mode_numbers(gtc_input)
    n_value = nmode[slot - 1] if slot <= len(nmode) else slot
    if slot <= len(mmode) and mmode[slot - 1] != 0:
        m_value = mmode[slot - 1]
        return rf"${field_symbol}_{{m={m_value},\,n={n_value}}}$"
    return rf"${field_symbol}_{{n={n_value}}}$"


def display_label(column: str, gtc_input: GtcInput | None = None) -> str:
    if column in PHYS_LABELS:
        return PHYS_LABELS[column]
    if column in TRACER_LABELS:
        return TRACER_LABELS[column]
    if column in RADIAL_LABELS:
        return RADIAL_LABELS[column]
    match = re.match(r"phi_mode_(\d+)_(amp|real|imag)$", column)
    if match:
        suffix = {"amp": r"|\delta\phi|", "real": r"\Re\,\delta\phi", "imag": r"\Im\,\delta\phi"}[
            match.group(2)
        ]
        return mode_title(gtc_input or GtcInput(None, ""), int(match.group(1)), suffix)
    match = re.match(r"density_mode_(\d+)_(amp|real|imag)$", column)
    if match:
        suffix = {"amp": r"|\delta n|", "real": r"\Re\,\delta n", "imag": r"\Im\,\delta n"}[
            match.group(2)
        ]
        return mode_title(gtc_input or GtcInput(None, ""), int(match.group(1)), suffix)
    return column.replace("_", " ")


def history_mode_order(mquantity: int, requested: str) -> str:
    if requested != "auto":
        return requested
    return "interleaved" if mquantity == len(PHYS_NAMES) else "field_blocks"


def history_columns(mquantity: int, num_mode: int, order: str) -> list[str]:
    if mquantity >= len(TRACER_NAMES) + len(PHYS_NAMES):
        columns = TRACER_NAMES + PHYS_NAMES
        if mquantity > len(columns):
            columns += [f"quantity_{i + 1}" for i in range(len(columns), mquantity)]
    elif mquantity == len(PHYS_NAMES):
        columns = PHYS_NAMES.copy()
    else:
        columns = [f"quantity_{i + 1}" for i in range(mquantity)]

    if order == "interleaved":
        for slot in range(1, num_mode + 1):
            for field_kind in ("phi", "density"):
                columns.append(f"{field_kind}_mode_{slot}_real")
                columns.append(f"{field_kind}_mode_{slot}_imag")
    elif order == "field_blocks":
        for field_kind in ("phi", "density"):
            for slot in range(1, num_mode + 1):
                columns.append(f"{field_kind}_mode_{slot}_real")
                columns.append(f"{field_kind}_mode_{slot}_imag")
    else:
        raise ValueError(f"unknown history mode order: {order}")
    return columns


def read_history(
    history_path: Path,
    gtc_input: GtcInput,
    *,
    mode_order: str = "auto",
) -> tuple[pd.DataFrame, dict[str, float | int | str]]:
    values = [float(line.strip()) for line in history_path.read_text().splitlines() if line.strip()]
    if len(values) < 6:
        raise ValueError(f"{history_path}: history file has fewer than 6 header values")
    irun = int(values[0])
    mquantity = int(values[1])
    mflux = int(values[2])
    num_mode = int(values[3])
    n_outputs = int(values[4])
    diagnostic_dt = values[5]
    record_len = mquantity + 4 * num_mode
    payload = np.asarray(values[6:], dtype=float)
    if payload.size % record_len != 0:
        raise ValueError(f"{history_path}: payload size {payload.size} is not divisible by {record_len}")

    order = history_mode_order(mquantity, mode_order)
    columns = history_columns(mquantity, num_mode, order)
    df = pd.DataFrame(payload.reshape((-1, record_len)), columns=columns)
    ndiag = gtc_input.integer("ndiag", 5)
    tstep = gtc_input.real("tstep", 0.2)
    df.insert(0, "step", np.arange(1, len(df) + 1) * ndiag)
    df.insert(1, "time", df["step"].to_numpy(dtype=float) * tstep)
    df.insert(2, "diagnostic_time", np.arange(1, len(df) + 1) * diagnostic_dt)

    for field_kind in ("phi", "density"):
        for slot in range(1, num_mode + 1):
            real = df.get(f"{field_kind}_mode_{slot}_real")
            imag = df.get(f"{field_kind}_mode_{slot}_imag")
            if real is not None and imag is not None:
                df[f"{field_kind}_mode_{slot}_amp"] = np.hypot(real, imag)

    return df, {
        "irun": irun,
        "mquantity": mquantity,
        "mflux": mflux,
        "num_mode": num_mode,
        "n_outputs_header": n_outputs,
        "n_outputs_parsed": len(df),
        "diagnostic_dt": diagnostic_dt,
        "record_len": record_len,
        "history_mode_order": order,
    }


def time_values(df: pd.DataFrame, axis: str) -> np.ndarray:
    if axis == "step":
        return df["step"].to_numpy(dtype=float)
    if axis == "diagnostic_time":
        return df["diagnostic_time"].to_numpy(dtype=float)
    return df["time"].to_numpy(dtype=float)


def x_axis_label(axis: str) -> str:
    if axis == "step":
        return "time step"
    if axis == "diagnostic_time":
        return "diagnostic time"
    return r"$t\,[(C_s/L_n)^{-1}]$"


def plot_panel_lines(
    df: pd.DataFrame,
    columns: list[str],
    title: str,
    path: Path,
    gtc_input: GtcInput,
    *,
    log_abs: bool = False,
    x_axis: str = "time",
) -> Path:
    if not columns:
        return path
    ncols = 2
    nrows = int(np.ceil(len(columns) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 2.55 * nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    x = time_values(df, x_axis)
    for ax, column in zip(axes, columns):
        y = df[column].to_numpy(dtype=float)
        if log_abs:
            ax.semilogy(x, np.maximum(np.abs(y), np.finfo(float).tiny), linewidth=1.15)
            ax.set_ylabel(r"$|\mathrm{value}|$")
        else:
            ax.plot(x, y, linewidth=1.15)
            ax.set_ylabel("value")
        ax.set_title(display_label(column, gtc_input), fontsize=10)
        ax.grid(True, alpha=0.3)
    for ax in axes[len(columns) :]:
        ax.axis("off")
    fig.supxlabel(x_axis_label(x_axis))
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    if x.size < 2 or y.size < 2:
        return float("nan"), float("nan"), float("nan")
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept
    ss_res = float(np.sum((y - predicted) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")
    return float(slope), float(intercept), float(r2)


def mode_fit_mask(time_axis: np.ndarray, amplitude: np.ndarray) -> np.ndarray:
    y = np.asarray(amplitude, dtype=float)
    finite = np.isfinite(time_axis) & np.isfinite(y) & (y > 0.0)
    if not np.any(finite):
        return finite
    ymax = float(np.nanmax(y[finite]))
    if ymax <= 0.0:
        return finite
    for lower, upper in ((1.0e-6, 0.98), (1.0e-8, 0.995), (1.0e-10, 1.0), (0.0, 1.0)):
        mask = finite & (y >= ymax * lower) & (y <= ymax * upper)
        if np.count_nonzero(mask) >= 20:
            return mask
    return finite


def fit_complex_mode(time_axis: np.ndarray, complex_amplitude: np.ndarray, tstep: float) -> dict[str, float]:
    z = np.asarray(complex_amplitude, dtype=np.complex128)
    amplitude = np.abs(z)
    mask = mode_fit_mask(time_axis, amplitude)
    if np.count_nonzero(mask) < 2:
        return {
            "gamma": float("nan"),
            "omega": float("nan"),
            "frequency": float("nan"),
            "fit_time_start": float("nan"),
            "fit_time_end": float("nan"),
            "fit_step_start": float("nan"),
            "fit_step_end": float("nan"),
            "r2_growth": float("nan"),
            "r2_phase": float("nan"),
        }
    fit_time = time_axis[mask]
    gamma, _, r2_growth = linear_fit(fit_time, np.log(amplitude[mask]))
    phase = np.unwrap(np.angle(z[mask]))
    omega, _, r2_phase = linear_fit(fit_time, phase)
    return {
        "gamma": gamma,
        "omega": omega,
        "frequency": omega / (2.0 * np.pi),
        "fit_time_start": float(fit_time[0]),
        "fit_time_end": float(fit_time[-1]),
        "fit_step_start": float(fit_time[0] / tstep) if tstep else float("nan"),
        "fit_step_end": float(fit_time[-1] / tstep) if tstep else float("nan"),
        "r2_growth": r2_growth,
        "r2_phase": r2_phase,
    }


def density_mode_fits(
    df: pd.DataFrame, meta: dict[str, float | int | str], gtc_input: GtcInput
) -> list[dict[str, float | int]]:
    num_mode = int(meta["num_mode"])
    nmode, mmode = mode_numbers(gtc_input)
    t_axis = time_values(df, "time")
    tstep = gtc_input.real("tstep", 0.2)
    rows: list[dict[str, float | int]] = []
    for slot in range(1, num_mode + 1):
        real_col = f"density_mode_{slot}_real"
        imag_col = f"density_mode_{slot}_imag"
        if real_col not in df.columns or imag_col not in df.columns:
            continue
        z = df[real_col].to_numpy(dtype=float) + 1j * df[imag_col].to_numpy(dtype=float)
        fit = fit_complex_mode(t_axis, z, tstep)
        fit.update(
            {
                "slot": slot,
                "n": nmode[slot - 1] if slot <= len(nmode) else slot,
                "m": mmode[slot - 1] if slot <= len(mmode) else slot,
            }
        )
        rows.append(fit)
    return rows


def write_density_mode_fit_summary(rows: list[dict[str, float | int]], path: Path) -> Path:
    columns = [
        "slot",
        "n",
        "m",
        "gamma",
        "omega",
        "frequency",
        "fit_step_start",
        "fit_step_end",
        "fit_time_start",
        "fit_time_end",
        "r2_growth",
        "r2_phase",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fp:
        fp.write(",".join(columns) + "\n")
        for row in rows:
            values = []
            for column in columns:
                value = row[column]
                if isinstance(value, (int, np.integer)):
                    values.append(str(int(value)))
                else:
                    values.append(f"{float(value):.9e}")
            fp.write(",".join(values) + "\n")
    return path


def plot_mode_amplitudes(
    df: pd.DataFrame,
    meta: dict[str, float | int | str],
    gtc_input: GtcInput,
    out_dir: Path,
) -> list[Path]:
    num_mode = int(meta["num_mode"])
    output_paths: list[Path] = []
    if num_mode <= 0:
        return output_paths

    def plot_group(field_kind: str, field_symbol: str, out_name: str, title: str) -> None:
        ncols = 2
        nrows = int(np.ceil(num_mode / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(12, 2.55 * nrows), sharex=True)
        axes = np.atleast_1d(axes).ravel()
        t_axis = time_values(df, "time")
        for slot, ax in enumerate(axes[:num_mode], start=1):
            column = f"{field_kind}_mode_{slot}_amp"
            if column not in df.columns:
                ax.axis("off")
                continue
            y = df[column].to_numpy(dtype=float)
            ax.semilogy(t_axis, np.maximum(y, np.finfo(float).tiny), linewidth=1.1)
            ax.set_title(mode_title(gtc_input, slot, field_symbol), fontsize=10)
            ax.set_ylabel(rf"$|{field_symbol}|$")
            ax.grid(True, alpha=0.3)
        for ax in axes[num_mode:]:
            ax.axis("off")
        for ax in axes[-ncols:]:
            ax.set_xlabel(x_axis_label("time"))
        fig.suptitle(title, fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        out_path = out_dir / out_name
        fig.savefig(out_path, dpi=180)
        plt.close(fig)
        output_paths.append(out_path)

    plot_group("phi", r"\delta\phi", "delta_phi_mode_amplitudes.png", r"$\delta\phi$ mode amplitudes")
    plot_group("density", r"\delta n", "delta_n_mode_amplitudes.png", r"$\delta n$ mode amplitudes")
    return output_paths


def plot_growth_normalized_density_modes(
    df: pd.DataFrame,
    meta: dict[str, float | int | str],
    fit_rows: list[dict[str, float | int]],
    out_dir: Path,
) -> Path:
    num_mode = int(meta["num_mode"])
    if num_mode <= 0:
        return out_dir / "delta_n_growth_normalized_oscillations.png"
    ncols = 2
    nrows = int(np.ceil(num_mode / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 2.7 * nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    t_axis = time_values(df, "time")
    fits_by_slot = {int(row["slot"]): row for row in fit_rows}
    for slot, ax in enumerate(axes[:num_mode], start=1):
        fit = fits_by_slot.get(slot)
        real_col = f"density_mode_{slot}_real"
        imag_col = f"density_mode_{slot}_imag"
        if fit is None or real_col not in df.columns or imag_col not in df.columns:
            ax.axis("off")
            continue
        z = df[real_col].to_numpy(dtype=float) + 1j * df[imag_col].to_numpy(dtype=float)
        gamma = float(fit["gamma"])
        t0 = float(fit["fit_time_start"])
        t1 = float(fit["fit_time_end"])
        window = (t_axis >= t0) & (t_axis <= t1) if np.isfinite(t0) and np.isfinite(t1) else np.ones_like(t_axis, dtype=bool)
        normalized = z[window] * np.exp(-gamma * (t_axis[window] - t0)) if np.isfinite(gamma) and np.isfinite(t0) else z[window]
        scale = np.nanmax(np.abs(normalized)) if normalized.size else 1.0
        if not np.isfinite(scale) or scale <= 0.0:
            scale = 1.0
        normalized = normalized / scale
        ax.plot(t_axis[window], normalized.real, linewidth=1.0, label=r"$\Re$")
        ax.plot(t_axis[window], normalized.imag, linewidth=1.0, linestyle="--", label=r"$\Im$")
        ax.set_title(
            rf"$\delta n_{{m={int(fit['m'])},\,n={int(fit['n'])}}}$, "
            rf"$\gamma={float(fit['gamma']):.3g}$, $\omega={float(fit['omega']):.3g}$",
            fontsize=10,
        )
        ax.set_ylabel("normalized")
        ax.grid(True, alpha=0.3)
    for ax in axes[num_mode:]:
        ax.axis("off")
    for ax in axes[-ncols:]:
        ax.set_xlabel(x_axis_label("time"))
    if num_mode > 0:
        axes[0].legend(loc="best")
    fig.suptitle(r"Growth-normalized $\delta n$ oscillations and fitted frequencies", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_path = out_dir / "delta_n_growth_normalized_oscillations.png"
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def read_int(tokens: list[str], idx: int) -> tuple[int, int]:
    return int(float(tokens[idx])), idx + 1


def read_floats(tokens: list[str], idx: int, count: int) -> tuple[np.ndarray, int]:
    arr = np.asarray([float(x) for x in tokens[idx : idx + count]], dtype=float)
    return arr, idx + count


def parse_snapshot(path: Path) -> dict[str, object]:
    tokens = path.read_text().split()
    idx = 0
    time_value = float(tokens[idx])
    idx += 1
    q_mid = float(tokens[idx])
    idx += 1

    mbin_u, idx = read_int(tokens, idx)
    mbin_psi, idx = read_int(tokens, idx)
    mpsi, idx = read_int(tokens, idx)
    jm, idx = read_int(tokens, idx)
    mzetamax, idx = read_int(tokens, idx)
    midiag, idx = read_int(tokens, idx)

    n_u_quant, idx = read_int(tokens, idx)
    idx += 3 * mbin_u
    idx += mbin_psi * 6 * mbin_u

    n_r_quant, idx = read_int(tokens, idx)
    radial, idx = read_floats(tokens, idx, mpsi)
    radial_quantities = {}
    for name in RADIAL_NAMES:
        radial_quantities[name], idx = read_floats(tokens, idx, mpsi)

    n_pol_quant, idx = read_int(tokens, idx)
    shape_pol = (jm + 1, mpsi)
    x, idx = read_floats(tokens, idx, (jm + 1) * mpsi)
    z, idx = read_floats(tokens, idx, (jm + 1) * mpsi)
    phi_xz, idx = read_floats(tokens, idx, (jm + 1) * mpsi)

    n_flux_quant, idx = read_int(tokens, idx)
    flux_phi, idx = read_floats(tokens, idx, jm * mzetamax)

    num_mode, idx = read_int(tokens, idx)
    m_poloidal, idx = read_int(tokens, idx)
    nmode, idx = read_floats(tokens, idx, num_mode)
    eigenmode, idx = read_floats(tokens, idx, mpsi * num_mode * m_poloidal)

    return {
        "time": time_value,
        "q_mid": q_mid,
        "mbin_u": mbin_u,
        "mbin_psi": mbin_psi,
        "mpsi": mpsi,
        "jm": jm,
        "mzetamax": mzetamax,
        "midiag": midiag,
        "n_u_quant": n_u_quant,
        "n_r_quant": n_r_quant,
        "n_pol_quant": n_pol_quant,
        "n_flux_quant": n_flux_quant,
        "radial": radial,
        "radial_quantities": radial_quantities,
        "x": x.reshape(shape_pol),
        "z": z.reshape(shape_pol),
        "phi_xz": phi_xz.reshape(shape_pol),
        "flux_phi": flux_phi.reshape((jm, mzetamax)),
        "num_mode": num_mode,
        "m_poloidal": m_poloidal,
        "nmode": nmode,
        "eigenmode": eigenmode.reshape((mpsi, num_mode, m_poloidal)),
    }


def symmetric_limits(arr: np.ndarray) -> tuple[float, float]:
    vmax = float(np.nanmax(np.abs(arr)))
    if not np.isfinite(vmax) or vmax == 0.0:
        vmax = 1.0
    return -vmax, vmax


def snapshot_time_label(path: Path, gtc_input: GtcInput) -> str:
    step = snapshot_step(path)
    if step > 0:
        return f"{gtc_input.real('tstep', 0.2) * step:g}"
    return ""


def plot_snapshot_eigenmode(data: dict[str, object], step: int, out_dir: Path, suffix: str) -> dict[str, Path]:
    num_mode = int(data["num_mode"])
    m_poloidal = int(data["m_poloidal"])
    if num_mode <= 0 or m_poloidal <= 0:
        return {}

    radial = np.asarray(data["radial"], dtype=float)
    nmode = np.asarray(data["nmode"], dtype=int)
    eigenmode = np.abs(np.asarray(data["eigenmode"], dtype=float))
    safe = np.maximum(eigenmode, np.finfo(float).tiny)
    log_amp = np.log10(safe)
    vmin = float(np.nanpercentile(log_amp, 2.0))
    vmax = float(np.nanpercentile(log_amp, 99.5))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
        vmin, vmax = -16.0, 0.0

    ncols = 2
    nrows = int(np.ceil(num_mode / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(13, 2.75 * nrows), sharex=True, sharey=True, constrained_layout=True
    )
    axes = np.atleast_1d(axes).ravel()
    for ax, mode_idx in zip(axes, range(num_mode)):
        image = log_amp[:, mode_idx, :].T
        mesh = ax.imshow(
            image,
            origin="lower",
            aspect="auto",
            interpolation="nearest",
            cmap="magma",
            vmin=vmin,
            vmax=vmax,
            extent=[radial[0], radial[-1], 0, m_poloidal - 1],
        )
        ax.set_title(rf"$n={int(nmode[mode_idx])}$", fontsize=10)
        ax.set_xlabel(r"$r/\rho_i$")
        ax.set_ylabel(r"poloidal mode index")
    for ax in axes[num_mode:]:
        ax.axis("off")
    cbar = fig.colorbar(mesh, ax=axes[:num_mode].tolist(), shrink=0.92)
    cbar.set_label(r"$\log_{10}|\delta\phi_{n,m}(r)|$")
    fig.suptitle(rf"Snapshot selected-$n$ eigenmode spectrum, step {step}", fontsize=14)
    heatmap_path = out_dir / f"eigenmode_selected_n_{suffix}.png"
    fig.savefig(heatmap_path, dpi=180)
    plt.close(fig)

    rows = []
    for mode_idx, n_value in enumerate(nmode):
        plane = eigenmode[:, mode_idx, :]
        if plane.size == 0 or not np.isfinite(plane).any():
            continue
        flat = int(np.nanargmax(plane))
        radial_idx, poloidal_idx = np.unravel_index(flat, plane.shape)
        rows.append(
            {
                "step": step,
                "n": int(n_value),
                "peak_radial_index": int(radial_idx + 1),
                "peak_radial": float(radial[radial_idx]),
                "peak_poloidal_index": int(poloidal_idx),
                "peak_amplitude": float(plane[radial_idx, poloidal_idx]),
            }
        )
    summary_path = out_dir / f"eigenmode_selected_n_summary_{suffix}.csv"
    pd.DataFrame(rows).to_csv(summary_path, index=False)
    return {"eigenmode_selected_n": heatmap_path, "eigenmode_selected_n_summary": summary_path}


def plot_snapshot(path: Path, out_dir: Path, gtc_input: GtcInput) -> dict[str, Path]:
    data = parse_snapshot(path)
    step = snapshot_step(path)
    time_label = snapshot_time_label(path, gtc_input)
    suffix = f"{step:05d}" if step > 0 else path.stem
    outputs: dict[str, Path] = {}

    x = data["x"]
    z = data["z"]
    phi = data["phi_xz"]
    fig, ax = plt.subplots(figsize=(7, 6))
    vmin, vmax = symmetric_limits(phi)
    triang = mtri.Triangulation(x.ravel(), z.ravel())
    mesh = ax.tricontourf(triang, phi.ravel(), levels=80, cmap="RdBu_r", vmin=vmin, vmax=vmax)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X")
    ax.set_ylabel("Z")
    ax.set_title(rf"$\delta\phi(X,Z)$, $t={time_label}\,(C_s/L_n)^{{-1}}$")
    fig.colorbar(mesh, ax=ax, label=r"$\delta\phi/\rho_i^2$")
    fig.tight_layout()
    outputs["poloidal_phi"] = out_dir / f"poloidal_phi_{suffix}.png"
    fig.savefig(outputs["poloidal_phi"], dpi=180)
    plt.close(fig)

    flux_phi = data["flux_phi"]
    fig, ax = plt.subplots(figsize=(8, 5))
    vmin, vmax = symmetric_limits(flux_phi)
    mesh = ax.imshow(
        flux_phi,
        origin="lower",
        aspect="auto",
        cmap="RdBu_r",
        vmin=vmin,
        vmax=vmax,
        extent=[0, data["mzetamax"], 1, data["jm"]],
    )
    ax.set_xlabel(r"$\zeta$ grid index")
    ax.set_ylabel(r"$\theta$ grid index")
    ax.set_title(rf"$\delta\phi(\theta,\zeta)$, step {step}")
    fig.colorbar(mesh, ax=ax, label=r"$\delta\phi/\rho_i^2$")
    fig.tight_layout()
    outputs["flux_phi"] = out_dir / f"flux_phi_{suffix}.png"
    fig.savefig(outputs["flux_phi"], dpi=180)
    plt.close(fig)

    outputs.update(plot_snapshot_eigenmode(data, step, out_dir, suffix))

    radial = data["radial"]
    radial_quantities = data["radial_quantities"]
    ncols = 2
    nrows = int(np.ceil(len(RADIAL_NAMES) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12.5, 2.65 * nrows), sharex=False)
    axes = np.atleast_1d(axes).ravel()
    for ax, name in zip(axes, RADIAL_NAMES):
        ax.plot(radial, radial_quantities[name], linewidth=1.1)
        ax.set_title(RADIAL_LABELS[name], fontsize=10)
        ax.set_xlabel(r"$r/\rho_i$")
        ax.set_ylabel("value")
        ax.grid(True, alpha=0.3)
    for ax in axes[len(RADIAL_NAMES) :]:
        ax.axis("off")
    fig.suptitle(rf"Snapshot radial 1-D profiles, step {step}, $t={time_label}$", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    outputs["radial_1d"] = out_dir / f"radial_1d_{suffix}.png"
    fig.savefig(outputs["radial_1d"], dpi=180)
    plt.close(fig)

    return outputs


def write_history_summary(df: pd.DataFrame, columns: list[str], path: Path) -> Path:
    rows = []
    for column in columns:
        values = df[column].to_numpy(dtype=float)
        rows.append(
            {
                "quantity": column,
                "final": values[-1] if values.size else np.nan,
                "abs_max": np.nanmax(np.abs(values)) if values.size else np.nan,
                "min": np.nanmin(values) if values.size else np.nan,
                "max": np.nanmax(values) if values.size else np.nan,
            }
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def plot_case(case: CaseSpec, out_dir: Path, args: argparse.Namespace) -> dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    gtc_input = case.gtc_input
    outputs: dict[str, object] = {"out_dir": out_dir}

    if not args.no_history:
        if case.history_path is None or not case.history_path.exists():
            raise FileNotFoundError(f"history file not found: {case.history_path}")
        history_dir = out_dir / "history"
        history_dir.mkdir(parents=True, exist_ok=True)
        df, meta = read_history(case.history_path, gtc_input, mode_order=args.history_mode_order)
        df.to_csv(history_dir / "history_parsed.csv", index=False)
        meta_path = history_dir / "history_meta.csv"
        pd.DataFrame([meta]).to_csv(meta_path, index=False)

        phys_cols = [column for column in PHYS_NAMES if column in df.columns]
        tracer_cols = [column for column in TRACER_NAMES if column in df.columns]
        plot_panel_lines(
            df,
            phys_cols,
            f"{case.label}: history physical scalars",
            history_dir / "physical_scalars.png",
            gtc_input,
            x_axis=args.x_axis,
        )
        plot_panel_lines(
            df,
            phys_cols,
            f"{case.label}: history physical scalars (absolute log scale)",
            history_dir / "physical_scalars_log_abs.png",
            gtc_input,
            log_abs=True,
            x_axis=args.x_axis,
        )
        if tracer_cols:
            plot_panel_lines(
                df,
                tracer_cols,
                f"{case.label}: history tracer scalars",
                history_dir / "tracer_scalars.png",
                gtc_input,
                x_axis=args.x_axis,
            )

        num_mode = int(meta["num_mode"])
        amp_cols = [f"phi_mode_{i}_amp" for i in range(1, num_mode + 1) if f"phi_mode_{i}_amp" in df.columns]
        amp_cols += [
            f"density_mode_{i}_amp" for i in range(1, num_mode + 1) if f"density_mode_{i}_amp" in df.columns
        ]
        plot_panel_lines(
            df,
            amp_cols,
            f"{case.label}: history mode amplitudes",
            history_dir / "mode_amplitudes.png",
            gtc_input,
            log_abs=True,
            x_axis=args.x_axis,
        )
        plot_mode_amplitudes(df, meta, gtc_input, history_dir)
        fits = density_mode_fits(df, meta, gtc_input)
        write_density_mode_fit_summary(fits, history_dir / "delta_n_growth_frequency.csv")
        plot_growth_normalized_density_modes(df, meta, fits, history_dir)
        write_history_summary(df, phys_cols, history_dir / "history_physical_summary.csv")
        outputs["history_records"] = len(df)
        outputs["history_dir"] = history_dir

    if not args.no_snap:
        snap_paths = case.snapshot_paths
        if args.max_snapshots is not None:
            snap_paths = snap_paths[: args.max_snapshots]
        snap_dir = out_dir / "snapshots"
        snap_dir.mkdir(parents=True, exist_ok=True)
        poloidal_paths: list[Path] = []
        flux_paths: list[Path] = []
        radial_paths: list[Path] = []
        eigenmode_paths: list[Path] = []
        for snap_path in snap_paths:
            snap_outputs = plot_snapshot(snap_path, snap_dir, gtc_input)
            poloidal_paths.append(snap_outputs["poloidal_phi"])
            flux_paths.append(snap_outputs["flux_phi"])
            radial_paths.append(snap_outputs["radial_1d"])
            if "eigenmode_selected_n" in snap_outputs:
                eigenmode_paths.append(snap_outputs["eigenmode_selected_n"])
        write_image_pdf(poloidal_paths, snap_dir / "poloidal_phi_snapshots.pdf")
        write_image_pdf(flux_paths, snap_dir / "flux_phi_snapshots.pdf")
        write_image_pdf(radial_paths, snap_dir / "radial_1d_snapshots.pdf")
        write_image_pdf(eigenmode_paths, snap_dir / "eigenmode_selected_n_snapshots.pdf")
        outputs["snapshots"] = len(snap_paths)
        outputs["snapshot_dir"] = snap_dir

    return outputs


def write_image_pdf(image_paths: Iterable[Path], pdf_path: Path) -> Path | None:
    image_paths = list(image_paths)
    if not image_paths:
        return None
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(pdf_path) as pdf:
        for image_path in image_paths:
            image = plt.imread(image_path)
            height, width = image.shape[:2]
            fig_width = 11.0
            fig_height = max(4.0, fig_width * height / width)
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            ax.imshow(image)
            ax.axis("off")
            fig.tight_layout(pad=0)
            pdf.savefig(fig)
            plt.close(fig)
    return pdf_path


def plot_history_compare(
    case_data: list[tuple[CaseSpec, GtcInput, pd.DataFrame]],
    columns: list[str],
    out_dir: Path,
    *,
    log_abs: bool,
    x_axis: str,
    title: str,
    filename: str,
) -> Path:
    ncols = 2
    nrows = int(np.ceil(len(columns) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 2.6 * nrows), sharex=False)
    axes = np.atleast_1d(axes).ravel()
    for ax, column in zip(axes, columns):
        for case, _gtc_input, df in case_data:
            x = time_values(df, x_axis)
            y = df[column].to_numpy(dtype=float)
            if log_abs:
                ax.semilogy(x, np.maximum(np.abs(y), np.finfo(float).tiny), linewidth=1.1, label=case.label)
                ax.set_ylabel(r"$|\mathrm{value}|$")
            else:
                ax.plot(x, y, linewidth=1.1, label=case.label)
                ax.set_ylabel("value")
        ax.set_title(display_label(column, case_data[0][1]), fontsize=10)
        ax.set_xlabel(x_axis_label(x_axis))
        ax.grid(True, alpha=0.3)
    for ax in axes[len(columns) :]:
        ax.axis("off")
    axes[0].legend(loc="best")
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_path = out_dir / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def write_compare_summary(
    case_data: list[tuple[CaseSpec, GtcInput, pd.DataFrame]], columns: list[str], path: Path
) -> Path:
    rows = []
    for case, _gtc_input, df in case_data:
        for column in columns:
            values = df[column].to_numpy(dtype=float)
            rows.append(
                {
                    "case": case.label,
                    "quantity": column,
                    "final": values[-1] if values.size else np.nan,
                    "abs_max": np.nanmax(np.abs(values)) if values.size else np.nan,
                    "min": np.nanmin(values) if values.size else np.nan,
                    "max": np.nanmax(values) if values.size else np.nan,
                }
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def compare_history(cases: list[CaseSpec], out_dir: Path, args: argparse.Namespace) -> dict[str, object]:
    if len(cases) < 2:
        raise ValueError("compare-history needs at least two cases")
    out_dir.mkdir(parents=True, exist_ok=True)
    case_data: list[tuple[CaseSpec, GtcInput, pd.DataFrame]] = []
    for case in cases:
        if case.history_path is None or not case.history_path.exists():
            raise FileNotFoundError(f"history file not found for {case.label}: {case.history_path}")
        gtc_input = case.gtc_input
        df, meta = read_history(case.history_path, gtc_input, mode_order=args.history_mode_order)
        case_dir = out_dir / "parsed" / sanitize_name(case.label)
        case_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(case_dir / "history_parsed.csv", index=False)
        pd.DataFrame([meta]).to_csv(case_dir / "history_meta.csv", index=False)
        case_data.append((case, gtc_input, df))

    if args.columns:
        requested_columns = args.columns
    else:
        requested_columns = PHYS_NAMES
    columns = [column for column in requested_columns if all(column in df.columns for _case, _input, df in case_data)]
    if not columns:
        common = sorted(set.intersection(*(set(df.columns) for _case, _input, df in case_data)))
        raise ValueError(f"no requested columns found in all cases. Common columns include: {', '.join(common[:20])}")

    png = plot_history_compare(
        case_data,
        columns,
        out_dir,
        log_abs=False,
        x_axis=args.x_axis,
        title=args.title or "History 1-D physical scalars comparison",
        filename=args.prefix + "history_1d_compare.png",
    )
    log_png = plot_history_compare(
        case_data,
        columns,
        out_dir,
        log_abs=True,
        x_axis=args.x_axis,
        title=(args.title or "History 1-D physical scalars comparison") + " (absolute log scale)",
        filename=args.prefix + "history_1d_compare_log_abs.png",
    )
    summary = write_compare_summary(case_data, columns, out_dir / (args.prefix + "history_1d_compare_summary.csv"))
    pdf = write_image_pdf([png, log_png], out_dir / (args.prefix + "history_1d_compare.pdf"))
    return {
        "out_dir": out_dir,
        "cases": len(cases),
        "columns": len(columns),
        "plot": png,
        "log_plot": log_png,
        "summary": summary,
        "pdf": pdf,
    }


def sanitize_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return value or "case"


def labels_for_count(labels: list[str], count: int) -> list[str | None]:
    if not labels:
        return [None] * count
    if len(labels) != count:
        raise ValueError(f"got {len(labels)} labels for {count} cases")
    return labels


def build_compare_cases(args: argparse.Namespace) -> list[CaseSpec]:
    cases: list[CaseSpec] = []
    if args.case:
        cases.extend(parse_case_arg(spec) for spec in args.case)
    if args.case_dir:
        labels = labels_for_count(args.label, len(args.case_dir)) if not cases else [None] * len(args.case_dir)
        for case_dir, label in zip(args.case_dir, labels):
            cases.append(case_from_dir(resolve_path(case_dir), label))
    if args.history or args.input:
        histories = [resolve_path(path) for path in args.history or []]
        inputs = [resolve_path(path) for path in args.input or []]
        if len(histories) != len(inputs):
            raise ValueError("--input and --history must have the same count when used directly")
        already_labeled = len(cases)
        labels = args.label[already_labeled:] if args.label and len(args.label) == already_labeled + len(histories) else []
        labels = labels_for_count(labels, len(histories))
        for input_path, history_path, label in zip(inputs, histories, labels):
            cases.append(
                case_from_paths(
                    label=label,
                    case_dir=history_path.parent if history_path is not None else None,
                    input_path=input_path,
                    history_path=history_path,
                    snapshot_paths=(),
                )
            )
    if args.label and cases and len(args.label) == len(cases):
        cases = [
            CaseSpec(label=label, case_dir=case.case_dir, input_path=case.input_path, history_path=case.history_path, snapshot_paths=case.snapshot_paths)
            for case, label in zip(cases, args.label)
        ]
    return cases


def build_single_case(args: argparse.Namespace) -> CaseSpec:
    case_dir = resolve_path(args.case_dir) if args.case_dir else None
    if case_dir is not None and args.input is None and args.history is None:
        snap_paths = collect_snapshot_paths(case_dir, args.snap_glob or [], args.snap or [])
        case = case_from_dir(case_dir, args.label)
        return CaseSpec(
            label=case.label,
            case_dir=case.case_dir,
            input_path=case.input_path,
            history_path=case.history_path,
            snapshot_paths=snap_paths,
        )
    input_path = resolve_path(args.input) if args.input else (case_dir / "gtc.input" if case_dir else None)
    history_path = resolve_path(args.history) if args.history else (case_dir / "history.out" if case_dir else None)
    snap_paths = collect_snapshot_paths(case_dir, args.snap_glob or [], args.snap or [])
    return case_from_paths(
        label=args.label,
        case_dir=case_dir,
        input_path=input_path,
        history_path=history_path,
        snapshot_paths=snap_paths,
    )


def add_common_plot_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--history-mode-order",
        choices=["auto", "interleaved", "field_blocks"],
        default="auto",
        help="history mode storage order. auto uses interleaved for mquantity=12 and field_blocks otherwise.",
    )
    parser.add_argument(
        "--x-axis",
        choices=["time", "step", "diagnostic_time"],
        default="time",
        help="x axis for history plots.",
    )


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot GTC history.out and snap*.out files.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    case_parser = subparsers.add_parser("case", help="Plot one case history and snapshots.")
    case_parser.add_argument("--case-dir", help="case directory containing gtc.input, history.out, and snap*.out")
    case_parser.add_argument("--input", help="path to gtc.input")
    case_parser.add_argument("--history", help="path to history.out")
    case_parser.add_argument("--snap", action="append", help="explicit snap*.out path; can be repeated")
    case_parser.add_argument("--snap-glob", action="append", help="glob pattern for snapshots; can be repeated")
    case_parser.add_argument("--out-dir", required=True, help="output plot directory")
    case_parser.add_argument("--label", help="case label used in titles")
    case_parser.add_argument("--no-history", action="store_true", help="skip history plots")
    case_parser.add_argument("--no-snap", action="store_true", help="skip snapshot plots")
    case_parser.add_argument("--max-snapshots", type=int, help="only plot the first N snapshots after sorting")
    add_common_plot_args(case_parser)

    compare_parser = subparsers.add_parser("compare-history", help="Compare history 1-D quantities from two or more cases.")
    compare_parser.add_argument(
        "--case",
        action="append",
        help="case spec, repeatable: LABEL=CASE_DIR, CASE_DIR, LABEL=GTC_INPUT:HISTORY, or GTC_INPUT:HISTORY",
    )
    compare_parser.add_argument("--case-dir", action="append", help="case directory; can be repeated")
    compare_parser.add_argument("--input", action="append", help="gtc.input path; repeat with --history")
    compare_parser.add_argument("--history", action="append", help="history.out path; repeat with --input")
    compare_parser.add_argument("--label", action="append", help="case label; repeat in case order")
    compare_parser.add_argument("--out-dir", required=True, help="output comparison directory")
    compare_parser.add_argument("--columns", nargs="+", help="history columns to compare; default is physical scalars")
    compare_parser.add_argument("--prefix", default="", help="prefix for output filenames")
    compare_parser.add_argument("--title", help="custom figure title")
    add_common_plot_args(compare_parser)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    out_dir = resolve_path(args.out_dir)
    if out_dir is None:
        raise ValueError("--out-dir is required")
    os.environ.setdefault("MPLCONFIGDIR", str(out_dir / ".matplotlib"))

    if args.command == "case":
        case = build_single_case(args)
        result = plot_case(case, out_dir, args)
        print(f"case,{case.label}")
        print(f"output_dir,{result['out_dir']}")
        if "history_records" in result:
            print(f"history_records,{result['history_records']}")
            print(f"history_dir,{result['history_dir']}")
        if "snapshots" in result:
            print(f"snapshots,{result['snapshots']}")
            print(f"snapshot_dir,{result['snapshot_dir']}")
        return 0

    if args.command == "compare-history":
        cases = build_compare_cases(args)
        result = compare_history(cases, out_dir, args)
        print(f"cases,{result['cases']}")
        print(f"columns,{result['columns']}")
        print(f"output_dir,{result['out_dir']}")
        print(f"plot,{result['plot']}")
        print(f"log_plot,{result['log_plot']}")
        print(f"summary,{result['summary']}")
        print(f"pdf,{result['pdf']}")
        return 0

    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
