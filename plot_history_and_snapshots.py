#!/usr/bin/env python3
"""Plot GTC history scalars and snapshot 2-D fields."""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

CASE_DIR = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd()
os.environ.setdefault("MPLCONFIGDIR", str(CASE_DIR / ".matplotlib"))

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

HISTORY = CASE_DIR / "history.out"
OUT_DIR = CASE_DIR / "plots"

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


def read_history(path: Path) -> tuple[pd.DataFrame, dict[str, float | int]]:
    values = [float(line.strip()) for line in path.read_text().splitlines() if line.strip()]
    irun = int(values[0])
    mquantity = int(values[1])
    mflux = int(values[2])
    num_mode = int(values[3])
    n_outputs = int(values[4])
    diagnostic_dt = values[5]
    record_len = mquantity + 4 * num_mode
    payload = np.asarray(values[6:], dtype=float)
    if payload.size % record_len != 0:
        raise ValueError(f"history payload size {payload.size} is not divisible by {record_len}")

    columns = TRACER_NAMES + PHYS_NAMES
    for field_kind in ("phi", "density"):
        for slot in range(1, num_mode + 1):
            columns.append(f"{field_kind}_mode_{slot}_real")
            columns.append(f"{field_kind}_mode_{slot}_imag")

    df = pd.DataFrame(payload.reshape((-1, record_len)), columns=columns)
    df.insert(0, "step", np.arange(1, len(df) + 1) * 5)
    df.insert(1, "diagnostic_time", np.arange(1, len(df) + 1) * diagnostic_dt)

    for field_kind in ("phi", "density"):
        for slot in range(1, num_mode + 1):
            real = df[f"{field_kind}_mode_{slot}_real"]
            imag = df[f"{field_kind}_mode_{slot}_imag"]
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
    }


def plot_lines(df: pd.DataFrame, columns: list[str], title: str, path: Path, *, logy: bool = False) -> None:
    ncols = 2
    nrows = int(np.ceil(len(columns) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 2.5 * nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, column in zip(axes, columns):
        y = df[column].to_numpy()
        if logy:
            ax.semilogy(df["step"], np.maximum(np.abs(y), np.finfo(float).tiny), linewidth=1.2)
            ax.set_ylabel("|value|")
        else:
            ax.plot(df["step"], y, linewidth=1.2)
        ax.set_title(column, fontsize=10)
        ax.grid(True, alpha=0.3)
    for ax in axes[len(columns) :]:
        ax.axis("off")
    fig.suptitle(title, fontsize=14)
    fig.supxlabel("time step")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def read_mmode() -> list[int]:
    text = (CASE_DIR / "gtc.input").read_text()
    match = re.search(r"mmode\s*=\s*([^/\n]+)", text)
    if match is None:
        raise ValueError("Could not find mmode list in gtc.input")
    return [int(value) for value in re.findall(r"[-+]?\d+", match.group(1))]


def linear_growth_rate(time_axis: np.ndarray, amplitude: np.ndarray) -> float:
    y = np.asarray(amplitude, dtype=float)
    mask = np.isfinite(time_axis) & np.isfinite(y) & (y > 0.0)
    if np.count_nonzero(mask) < 2:
        return float("nan")
    gamma, _ = np.polyfit(time_axis[mask], np.log(y[mask]), 1)
    return float(gamma)


def plot_density_mode_history_pdf(df: pd.DataFrame, meta: dict[str, float | int], hist_dir: Path) -> None:
    num_mode = int(meta["num_mode"])
    mmode = read_mmode()
    time_axis = df["step"].to_numpy(dtype=float) * 0.2

    with PdfPages(hist_dir / "history_scalars_all.pdf") as pdf:
        for slot in range(1, num_mode + 1):
            column = f"density_mode_{slot}_amp"
            amplitude = np.maximum(df[column].to_numpy(dtype=float), np.finfo(float).tiny)
            gamma = linear_growth_rate(time_axis, amplitude)
            m_value = mmode[slot - 1] if slot <= len(mmode) else slot
            if m_value > 32:
                continue

            fig, ax = plt.subplots(figsize=(10, 4))
            ax.semilogy(time_axis, amplitude, linewidth=1.2)
            ax.set_title(
                rf"$\delta n$ amp, $m={m_value}$, "
                rf"$\gamma={gamma:.4g}\,C_s/L_n$"
            )
            ax.set_xlabel(r"$t\,[(C_s/L_n)^{-1}]$")
            ax.set_ylabel(r"$|\delta n|$")
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def plot_history(df: pd.DataFrame, meta: dict[str, float | int]) -> None:
    hist_dir = OUT_DIR / "history"
    hist_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(hist_dir / "history_parsed.csv", index=False)

    plot_lines(df, PHYS_NAMES, "history physical scalars", hist_dir / "physical_scalars.png")
    plot_lines(df, PHYS_NAMES, "history physical scalars (absolute log scale)", hist_dir / "physical_scalars_log_abs.png", logy=True)
    plot_lines(df, TRACER_NAMES, "history tracer scalars", hist_dir / "tracer_scalars.png")

    num_mode = int(meta["num_mode"])
    amp_cols = [f"phi_mode_{i}_amp" for i in range(1, num_mode + 1)]
    amp_cols += [f"density_mode_{i}_amp" for i in range(1, num_mode + 1)]
    plot_lines(df, amp_cols, "history mode amplitudes", hist_dir / "mode_amplitudes.png", logy=True)

    plot_density_mode_history_pdf(df, meta, hist_dir)


def read_int(tokens: list[str], idx: int) -> tuple[int, int]:
    return int(float(tokens[idx])), idx + 1


def read_floats(tokens: list[str], idx: int, count: int) -> tuple[np.ndarray, int]:
    arr = np.asarray([float(x) for x in tokens[idx : idx + count]], dtype=float)
    return arr, idx + count


def parse_snapshot(path: Path) -> dict[str, np.ndarray | float | int]:
    tokens = path.read_text().split()
    idx = 0
    time_value = float(tokens[idx]); idx += 1
    q_mid = float(tokens[idx]); idx += 1

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
    radial_names = ["zonali", "phip00_over_rho", "marker", "fflows", "dflows", "ftem", "dtem"]
    for name in radial_names:
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
    if vmax == 0.0:
        vmax = 1.0
    return -vmax, vmax


def poloidal_delta_phi_title(step: int) -> str:
    return rf"$\delta\phi(X,Z)$, $t={0.2 * step:g}\,(C_s/L_n)^{{-1}}$"


def plot_snapshot(path: Path, out_dir: Path) -> tuple[Path, Path]:
    data = parse_snapshot(path)
    step = int(path.stem.replace("snap", ""))
    x = data["x"]
    z = data["z"]
    phi = data["phi_xz"]
    flux_phi = data["flux_phi"]

    fig, ax = plt.subplots(figsize=(7, 6))
    vmin, vmax = symmetric_limits(phi)
    triang = mtri.Triangulation(x.ravel(), z.ravel())
    mesh = ax.tricontourf(triang, phi.ravel(), levels=80, cmap="RdBu_r", vmin=vmin, vmax=vmax)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X")
    ax.set_ylabel("Z")
    ax.set_title(poloidal_delta_phi_title(step))
    fig.colorbar(mesh, ax=ax, label=r"$\delta\phi/\rho_i^2$")
    fig.tight_layout()
    xz_path = out_dir / f"poloidal_phi_{step:05d}.png"
    fig.savefig(xz_path, dpi=180)
    plt.close(fig)

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
    ax.set_xlabel("toroidal grid index")
    ax.set_ylabel("theta grid index")
    ax.set_title(f"flux-surface phi(theta,zeta), step {step}")
    fig.colorbar(mesh, ax=ax, label="phi / rho_i^2")
    fig.tight_layout()
    flux_path = out_dir / f"flux_phi_{step:05d}.png"
    fig.savefig(flux_path, dpi=180)
    plt.close(fig)

    return xz_path, flux_path


def plot_snapshots() -> None:
    snap_dir = OUT_DIR / "snapshots"
    snap_dir.mkdir(parents=True, exist_ok=True)
    snapshots = sorted(CASE_DIR.glob("snap*.out"))
    xz_paths = []
    flux_paths = []
    for snapshot in snapshots:
        xz_path, flux_path = plot_snapshot(snapshot, snap_dir)
        xz_paths.append(xz_path)
        flux_paths.append(flux_path)

    with PdfPages(snap_dir / "poloidal_phi_snapshots.pdf") as pdf:
        for snapshot in snapshots:
            data = parse_snapshot(snapshot)
            step = int(snapshot.stem.replace("snap", ""))
            fig, ax = plt.subplots(figsize=(7, 6))
            vmin, vmax = symmetric_limits(data["phi_xz"])
            triang = mtri.Triangulation(data["x"].ravel(), data["z"].ravel())
            mesh = ax.tricontourf(triang, data["phi_xz"].ravel(), levels=80, cmap="RdBu_r", vmin=vmin, vmax=vmax)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("X")
            ax.set_ylabel("Z")
            ax.set_title(poloidal_delta_phi_title(step))
            fig.colorbar(mesh, ax=ax, label=r"$\delta\phi/\rho_i^2$")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    with PdfPages(snap_dir / "flux_phi_snapshots.pdf") as pdf:
        for snapshot in snapshots:
            data = parse_snapshot(snapshot)
            step = int(snapshot.stem.replace("snap", ""))
            fig, ax = plt.subplots(figsize=(8, 5))
            vmin, vmax = symmetric_limits(data["flux_phi"])
            mesh = ax.imshow(data["flux_phi"], origin="lower", aspect="auto", cmap="RdBu_r", vmin=vmin, vmax=vmax)
            ax.set_xlabel("toroidal grid index")
            ax.set_ylabel("theta grid index")
            ax.set_title(f"flux-surface phi(theta,zeta), step {step}")
            fig.colorbar(mesh, ax=ax, label="phi / rho_i^2")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    df, meta = read_history(HISTORY)
    plot_history(df, meta)
    plot_snapshots()
    print("history outputs:", OUT_DIR / "history")
    print("snapshot outputs:", OUT_DIR / "snapshots")
    print("history records:", meta["n_outputs_parsed"])
    print("snapshots:", len(list(CASE_DIR.glob("snap*.out"))))


if __name__ == "__main__":
    main()
