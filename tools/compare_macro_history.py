#!/usr/bin/env python3
import argparse
import math
import pathlib
import re
import sys


MACRO_LABELS = [
    "ddeni",
    "eradial",
    "efield",
    "entropyi",
    "dflowi",
    "pfluxi",
    "efluxi",
    "eflux1",
    "eflux2",
    "eflux3",
    "eflux4",
    "eflux5",
]


def read_numbers(path):
    text = path.read_text(errors="ignore").replace("D", "E")
    return [float(x) for x in re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?", text)]


def read_history(path):
    values = read_numbers(path)
    if len(values) < 6:
        raise ValueError(f"{path}: history file is too short")
    header = {
        "irun": int(values[0]),
        "mquantity": int(values[1]),
        "mflux": int(values[2]),
        "num_mode": int(values[3]),
        "noutputs": int(values[4]),
        "dt": values[5],
    }
    data = values[6:]
    noutputs = header["noutputs"]
    if noutputs <= 0:
        raise ValueError(f"{path}: invalid noutputs={noutputs}")
    if len(data) % noutputs != 0:
        raise ValueError(f"{path}: data count {len(data)} is not divisible by noutputs={noutputs}")
    stride = len(data) // noutputs
    mode_count = 4 * header["num_mode"]
    macro_count = stride - mode_count
    if macro_count < 0:
        raise ValueError(f"{path}: stride={stride} is smaller than mode block={mode_count}")
    rows = [data[i * stride:(i + 1) * stride] for i in range(noutputs)]
    return header, macro_count, mode_count, rows


def macro_offset(macro_count):
    if macro_count >= 20:
        return 8
    if macro_count >= 12:
        return 0
    raise ValueError(f"macro block too short: {macro_count}")


def rel_error(a, b):
    if not (math.isfinite(a) and math.isfinite(b)):
        return math.nan
    mag = max(abs(a), abs(b))
    return abs(a - b) / mag if mag > 0.0 else 0.0


def summarize_pairs(pairs):
    max_abs = 0.0
    max_rel = 0.0
    max_rel_sig = 0.0
    max_rel_large = 0.0
    n_sig = 0
    n_large = 0
    finite = 0
    worst = None
    for idx, a, b in pairs:
        if not (math.isfinite(a) and math.isfinite(b)):
            continue
        d = abs(a - b)
        r = rel_error(a, b)
        mag = max(abs(a), abs(b))
        finite += 1
        if d > max_abs:
            max_abs = d
            worst = (idx, a, b, d, r)
        max_rel = max(max_rel, r)
        if mag > 1.0e-6:
            n_sig += 1
            max_rel_sig = max(max_rel_sig, r)
        if mag > 1.0:
            n_large += 1
            max_rel_large = max(max_rel_large, r)
    return {
        "finite": finite,
        "max_abs": max_abs,
        "max_rel": max_rel,
        "max_rel_sig": max_rel_sig,
        "n_sig": n_sig,
        "max_rel_large": max_rel_large,
        "n_large": n_large,
        "worst": worst,
    }


def print_summary(label, summary):
    worst = summary["worst"]
    if worst:
        worst_text = f" worst_step={worst[0] + 1} f={worst[1]:.8e} c={worst[2]:.8e}"
    else:
        worst_text = ""
    print(
        f"{label}: finite={summary['finite']} "
        f"max_abs={summary['max_abs']:.6e} "
        f"max_rel={summary['max_rel']:.6e} "
        f"max_rel(|x|>1e-6)={summary['max_rel_sig']:.6e} n_sig={summary['n_sig']} "
        f"max_rel(|x|>1)={summary['max_rel_large']:.6e} n_large={summary['n_large']}"
        f"{worst_text}"
    )


def compare_histories(fortran_path, c_path):
    fh, fmacro, fmode_count, frows = read_history(fortran_path)
    ch, cmacro, cmode_count, crows = read_history(c_path)
    nsteps = min(fh["noutputs"], ch["noutputs"])
    foff = macro_offset(fmacro)
    coff = macro_offset(cmacro)
    macro_len = min(len(MACRO_LABELS), fmacro - foff, cmacro - coff)
    mode_len = min(fmode_count, cmode_count)

    print("headers:")
    print(f"  fortran: mquantity={fh['mquantity']} num_mode={fh['num_mode']} noutputs={fh['noutputs']} dt={fh['dt']:.8e}")
    print(f"  c      : mquantity={ch['mquantity']} num_mode={ch['num_mode']} noutputs={ch['noutputs']} dt={ch['dt']:.8e}")
    print(f"alignment: fortran_macro_skip={foff} c_macro_skip={coff} compared_steps={nsteps}")

    print("\nmacro summaries:")
    for i in range(macro_len):
        pairs = []
        for step in range(nsteps):
            pairs.append((step, frows[step][foff + i], crows[step][coff + i]))
        print_summary(MACRO_LABELS[i], summarize_pairs(pairs))

    print("\nlast-step macro values:")
    last = nsteps - 1
    print("name                 fortran              c              rel_error")
    for i in range(macro_len):
        a = frows[last][foff + i]
        b = crows[last][coff + i]
        print(f"{MACRO_LABELS[i]:<10} {a:18.8e} {b:18.8e} {rel_error(a, b):14.6e}")

    mode_pairs = []
    fmode_start = fmacro
    cmode_start = cmacro
    for step in range(nsteps):
        for i in range(mode_len):
            mode_pairs.append((step, frows[step][fmode_start + i], crows[step][cmode_start + i]))
    print("\nmode-amplitude summary:")
    print_summary("amp_mode", summarize_pairs(mode_pairs))


def main(argv):
    parser = argparse.ArgumentParser(
        description="Compare Fortran/C GTC history.out macro quantities after removing tracking fields."
    )
    parser.add_argument("fortran", help="Fortran history.out, or a directory containing history.out")
    parser.add_argument("c", help="C history.out, or a directory containing history.out")
    args = parser.parse_args(argv)

    fpath = pathlib.Path(args.fortran)
    cpath = pathlib.Path(args.c)
    if fpath.is_dir():
        fpath = fpath / "history.out"
    if cpath.is_dir():
        cpath = cpath / "history.out"
    compare_histories(fpath, cpath)


if __name__ == "__main__":
    main(sys.argv[1:])
