#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
work="${WORK_DIR:-$root/compare_fortran_c}"
fortran_exe="${FORTRAN_EXE:-$root/gtc_fortran_ref}"
c_exe="${C_EXE:-$root/gtc_c}"
mpi_np="${MPI_NP:-1}"

mstep="${MSTEP:-2}"
msnap="${MSNAP:-1}"
ndiag="${NDIAG:-2}"
micell="${MICELL:-20}"
mpsi="${MPSI:-50}"
mthetamax="${MTHETAMAX:-300}"
mzetamax="${MZETAMAX:-10}"
rng_control="${RNG_CONTROL:-1}"
npartdom="${NPARTDOM:-1}"
nonlinear="${NONLINEAR:-1.0}"
paranl="${PARANL:-0.0}"
mode00="${MODE00:-1}"
nbound="${NBOUND:-2}"

rm -rf "$work"
mkdir -p "$work/fortran" "$work/c"

if [[ ! -x "$fortran_exe" ]]; then
  make -C "$root" TARGET="$fortran_exe"
fi
if [[ ! -x "$c_exe" ]]; then
  make -C "$root" c
fi

write_input() {
  local dir="$1"
  mkdir -p "$dir"/{restart_dir1,restart_dir2,restart_dir,OVERFLOW_restart_dir1,OVERFLOW_restart_dir2,OVERFLOW_restart_dir,phi_dir,trackp_dir}
  cat > "$dir/gtc.input" <<EOF_INPUT
&input_parameters
  irun=0,
  mstep=${mstep},
  msnap=${msnap},
  ndiag=${ndiag},
  nonlinear=${nonlinear},
  paranl=${paranl},
  mode00=${mode00},
  tstep=0.2,
  micell=${micell},
  mpsi=${mpsi},
  mthetamax=${mthetamax},
  mzetamax=${mzetamax},
  npartdom=${npartdom},
  a=0.358,
  a0=0.1,
  a1=0.9,
  q0=0.581,
  q1=1.092,
  q2=1.092,
  rc=0.5,
  rw=0.4,
  aion=1.0,
  qion=1.0,
  kappati=6.9,
  kappate=6.9,
  kappan=2.2,
  tite=1.0,
  flow0=0.0,
  flow1=0.0,
  flow2=0.0,
  r0=93.2,
  b0=19100.0,
  temperature=2500.0,
  edensity0=1.46e14,
  stdout=6,
  nbound=${nbound},
  umax=4.0,
  iload=0,
  track_particles=0,
  nptrack=0,
  rng_control=${rng_control},
  nmode=20,20,20,20,20,20,20,20,20,20,20,20,20,
  mmode=25,26,27,28,29,30,31,32,33,34,35,36,37
/
EOF_INPUT
}

write_input "$work/fortran"
write_input "$work/c"

(cd "$work/fortran" && /usr/bin/time -p mpirun -np "$mpi_np" "$fortran_exe" > run.log 2> time.log)
(cd "$work/c" && /usr/bin/time -p mpirun -np "$mpi_np" "$c_exe" > run.log 2> time.log)

python3 - "$work" <<'PY'
import math
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1])

def nums(path):
    if not path.exists():
        return []
    text = path.read_text(errors="ignore")
    return [float(x) for x in re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?", text.replace("D", "E"))]

def numeric_lines(path, stop_at_eigen=False):
    if not path.exists():
        return []
    out = []
    lines = path.read_text(errors="ignore").splitlines()
    for i, line in enumerate(lines):
        if stop_at_eigen and line.strip() == "13" and i > 100:
            tail = [x.strip() for x in lines[i:i + 15]]
            if len(tail) >= 15 and tail[1] == "9":
                break
        token = line.strip().replace("D", "E")
        if token.lower() in {"nan", "+nan", "-nan"}:
            out.append(math.nan)
            continue
        try:
            out.append(float(token))
        except ValueError:
            pass
    return out

def compare_arrays(label, a, b, tol=1e-6):
    n = min(len(a), len(b))
    max_abs = 0.0
    max_rel = 0.0
    max_rel_sig = 0.0
    max_rel_large = 0.0
    finite = 0
    skipped_nan = 0
    significant = 0
    large = 0
    for x, y in zip(a[:n], b[:n]):
        if not (math.isfinite(x) and math.isfinite(y)):
            skipped_nan += 1
            continue
        d = abs(x - y)
        mag = max(abs(x), abs(y))
        rel = d / mag if mag > 0.0 else 0.0
        max_abs = max(max_abs, d)
        max_rel = max(max_rel, rel)
        if mag > 1.0e-6:
            max_rel_sig = max(max_rel_sig, rel)
            significant += 1
        if mag > 1.0:
            max_rel_large = max(max_rel_large, rel)
            large += 1
        finite += 1
    status = "OK" if len(a) == len(b) and max_rel_sig < tol else "DIFF"
    extra = f" skipped_nan={skipped_nan}" if skipped_nan else ""
    print(
        f"{status} {label}: count_fortran={len(a)} count_c={len(b)} finite={finite} "
        f"max_abs={max_abs:.6e} max_rel={max_rel:.6e} "
        f"max_rel_mag_gt_1e-6={max_rel_sig:.6e} n_sig={significant} "
        f"max_rel_mag_gt_1={max_rel_large:.6e} n_large={large}{extra}"
    )

def compare(name):
    a = nums(root / "fortran" / name)
    b = nums(root / "c" / name)
    compare_arrays(name, a, b)
    if name.startswith("snap"):
        af = numeric_lines(root / "fortran" / name, stop_at_eigen=True)
        bf = numeric_lines(root / "c" / name, stop_at_eigen=True)
        compare_arrays(name + " pre_eigen", af, bf)

names = ["history.out", "sheareb.out", "FileExit.dat"]
snaps = sorted({p.name for p in (root / "fortran").glob("snap*.out")} | {p.name for p in (root / "c").glob("snap*.out")})
for name in names + snaps:
    compare(name)
PY
