#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")" && pwd)"
case_dir="${CASE_DIR:-$root/itg_cyclone_base_nonlinear_zonal_c}"
exe="${GTC_EXE:-$root/gtc_c}"
npartdom="${NPARTDOM:-1}"
irun="${IRUN:-0}"
mstep="${MSTEP:-400}"
msnap="${MSNAP:-4}"
ndiag="${NDIAG:-5}"
micell="${MICELL:-20}"
mpsi="${MPSI:-50}"
mthetamax="${MTHETAMAX:-300}"
mzetamax="${MZETAMAX:-10}"
rng_control="${RNG_CONTROL:-1}"
nonlinear="${NONLINEAR:-1.0}"
paranl="${PARANL:-0.0}"
mode00="${MODE00:-1}"
spectrum_mode="${SPECTRUM_MODE:-0}"
nbound="${NBOUND:-10}"
nmode_list="${NMODE_LIST:-20,20,20,20,20,20,20,20,20,20,20,20,20}"
mmode_list="${MMODE_LIST:-25,26,27,28,29,30,31,32,33,34,35,36,37}"
keep_restart="${KEEP_RESTART:-0}"
restart_source="${RESTART_SOURCE:-}"

detect_cores() {
  sysctl -n hw.ncpu 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 10
}

largest_divisor_leq() {
  local n="$1"
  local limit="$2"
  local best=1
  local d
  if (( limit < 1 )); then
    echo 1
    return
  fi
  for ((d = 1; d <= limit && d <= n; d++)); do
    if (( n % d == 0 )); then
      best="$d"
    fi
  done
  echo "$best"
}

resolve_parallel_layout() {
  local layout="$1"
  local cores="$2"
  local mzeta="$3"
  local explicit_mpi="${MPI_NP:-}"
  local explicit_omp="${OMP_NUM_THREADS:-}"

  if [[ "$layout" == "manual" || "$layout" == "env" ||
        ( -z "${GTC_LAYOUT+x}" && ( -n "$explicit_mpi" || -n "$explicit_omp" ) ) ]]; then
    mpi_np="${explicit_mpi:-1}"
    if [[ -n "$explicit_omp" ]]; then
      omp_threads="$explicit_omp"
    else
      omp_threads=$(( cores / mpi_np ))
      if (( omp_threads < 1 )); then omp_threads=1; fi
    fi
    layout_requested="manual"
    return
  fi

  case "$layout" in
    omp|openmp|single-rank|single_rank|single_rank_multi_thread)
      mpi_np=1
      omp_threads="$cores"
      ;;
    mpi|pure-mpi|multi-rank|multi_rank|multi_rank_single_thread)
      mpi_np="$(largest_divisor_leq "$mzeta" "$cores")"
      omp_threads=1
      ;;
    hybrid|auto)
      local best_mpi=1
      local best_omp="$cores"
      local best_total=0
      local best_balance=2147483647
      local d omp total balance
      for ((d = 1; d <= cores && d <= mzeta; d++)); do
        if (( mzeta % d != 0 )); then continue; fi
        omp=$(( cores / d ))
        if (( omp < 1 )); then omp=1; fi
        total=$(( d * omp ))
        balance=$(( d > omp ? d - omp : omp - d ))
        if (( total > best_total ||
              (total == best_total && balance < best_balance) ||
              (total == best_total && balance == best_balance && d < best_mpi) )); then
          best_mpi="$d"
          best_omp="$omp"
          best_total="$total"
          best_balance="$balance"
        fi
      done
      mpi_np="$best_mpi"
      omp_threads="$best_omp"
      ;;
    *)
      echo "unknown GTC_LAYOUT=$layout; use auto, omp, mpi, hybrid, or manual" >&2
      exit 2
      ;;
  esac
}

host_cores="$(detect_cores)"
parallel_cores="${GTC_CORES:-$host_cores}"
layout_requested="${GTC_LAYOUT:-auto}"
mpi_np=""
omp_threads=""
resolve_parallel_layout "$layout_requested" "$parallel_cores" "$mzetamax"
linear_case="$(awk -v nonlinear="$nonlinear" 'BEGIN { print ((nonlinear + 0.0) < 0.5) ? 1 : 0 }')"
parallel_total_workers=$(( mpi_np * omp_threads ))
if (( mpi_np == 1 && omp_threads > 1 )); then
  layout_resolved="single_rank_multi_thread"
elif (( mpi_np > 1 && omp_threads == 1 )); then
  layout_resolved="multi_rank_single_thread"
elif (( mpi_np > 1 && omp_threads > 1 )); then
  layout_resolved="multi_rank_multi_thread"
else
  layout_resolved="serial"
fi
if (( mzetamax % mpi_np == 0 )); then
  mzetamax_decomposition="exact"
else
  mzetamax_decomposition="rounded_by_solver"
  echo "warning: MZETAMAX=$mzetamax is not divisible by MPI_NP=$mpi_np; solver will round mzetamax" >&2
fi

if [[ ! -x "$exe" ]]; then
  echo "missing executable: $exe" >&2
  echo "build it with: make c" >&2
  exit 1
fi

export OMP_NUM_THREADS="$omp_threads"
export OMP_PROC_BIND="${OMP_PROC_BIND:-spread}"
export OMP_PLACES="${OMP_PLACES:-cores}"
echo "parallel layout: requested=$layout_requested resolved=$layout_resolved cores=$parallel_cores mpi_np=$mpi_np omp_threads=$OMP_NUM_THREADS total_workers=$parallel_total_workers"

if [[ -e "$case_dir" ]]; then
  echo "case directory already exists: $case_dir" >&2
  echo "choose another name with CASE_DIR=/path/to/case $0" >&2
  exit 1
fi

mkdir -p "$case_dir"
cd "$case_dir"
mkdir -p restart_dir1 restart_dir2 restart_dir OVERFLOW_restart_dir1 OVERFLOW_restart_dir2 OVERFLOW_restart_dir phi_dir trackp_dir
if [[ -n "$restart_source" ]]; then
  cp -R "$restart_source"/restart_dir1/. restart_dir1/ 2>/dev/null || true
  cp -R "$restart_source"/restart_dir2/. restart_dir2/ 2>/dev/null || true
  cp -R "$restart_source"/restart_dir/. restart_dir/ 2>/dev/null || true
  cp "$restart_source"/FileExit.dat . 2>/dev/null || true
fi

{
cat <<EOF
&input_parameters
  irun=${irun},
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
  ! spectrum_mode: 0 keeps all toroidal n in the field, 1 filters to nmode
  ! nmode selects toroidal mode numbers; linear selected-mode history also pairs with mmode
spectrum_mode=${spectrum_mode},
nmode=${nmode_list}
EOF
if [[ "$linear_case" == "1" ]]; then
  cat <<EOF
mmode=${mmode_list}
EOF
fi
cat <<EOF
/
EOF
} > gtc.input

set +e
/usr/bin/time -p mpirun -np "$mpi_np" "$exe" > run.log 2> time.log
run_status=$?
set -e
if [[ "$run_status" -ne 0 ]] && ! grep -q 'Program ends' run.log; then
  exit "$run_status"
fi

{
  echo "case,cyclone_base_nonlinear_electrostatic_zonal_c"
  echo "mpirun_exit_code,$run_status"
  if grep -q 'Program ends' run.log; then
    echo "gtc_completed,1"
  else
    echo "gtc_completed,0"
  fi
  echo "parallel_layout_requested,$layout_requested"
  echo "parallel_layout_resolved,$layout_resolved"
  echo "parallel_cores,$parallel_cores"
  echo "parallel_total_workers,$parallel_total_workers"
  echo "mzetamax_decomposition,$mzetamax_decomposition"
  echo "mpi_np,$mpi_np"
  echo "omp_num_threads,$OMP_NUM_THREADS"
  echo "omp_proc_bind,$OMP_PROC_BIND"
  echo "omp_places,$OMP_PLACES"
  echo "mstep,$mstep"
  echo "msnap,$msnap"
  echo "ndiag,$ndiag"
  echo "micell,$micell"
  echo "mpsi,$mpsi"
  echo "mthetamax,$mthetamax"
  echo "mzetamax,$mzetamax"
  echo "nonlinear,$nonlinear"
  echo "mode00,$mode00"
  echo "spectrum_mode,$spectrum_mode"
  echo "linear_case,$linear_case"
  echo "nbound,$nbound"
  echo "nmode,$nmode_list"
  if [[ "$linear_case" == "1" ]]; then
    echo "mmode,$mmode_list"
  fi
  awk '/^real /{print "time_real_sec," $2}' time.log
  awk '/^user /{print "time_user_sec," $2}' time.log
  awk '/^sys /{print "time_sys_sec," $2}' time.log
  awk -F: '/MAIN LOOP TIME/{gsub(/^[ \t]+|[ \t]+$/,"",$2); print "main_loop_wall_sec," $2}' run.log
  awk -F: '/TOTAL WALL CLOCK TIME/{gsub(/^[ \t]+|[ \t]+$/,"",$2); print "program_wall_sec," $2}' run.log
  find . -maxdepth 1 -name 'snap*.out' | wc -l | awk '{print "snapshot_count," $1}'
  find restart_dir1 restart_dir2 restart_dir -type f -name 'restart_*.bp' | wc -l | awk '{print "restart_file_count," $1}'
} > summary.csv

if [[ "$keep_restart" != "1" ]]; then
  rm -rf restart_dir1 restart_dir2 restart_dir OVERFLOW_restart_dir1 OVERFLOW_restart_dir2 OVERFLOW_restart_dir phi_dir trackp_dir
  echo "restart_files_retained,0" >> summary.csv
else
  echo "restart_files_retained,1" >> summary.csv
fi

if command -v python3 >/dev/null 2>&1; then
  python3 "$root/gtc_plot.py" case --case-dir "$case_dir" --out-dir "$case_dir/plots" > plot.log 2>&1 || true
fi

echo "case output: $case_dir"
echo "summary: $case_dir/summary.csv"
