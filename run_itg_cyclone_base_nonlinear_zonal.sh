#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")" && pwd)"
case_dir="${CASE_DIR:-$root/itg_cyclone_base_nonlinear_zonal}"
exe="${GTC_EXE:-$root/gtc_local_num_mode13}"
mpi_np="${MPI_NP:-10}"
mstep="${MSTEP:-400}"
msnap="${MSNAP:-4}"
ndiag="${NDIAG:-5}"
keep_restart="${KEEP_RESTART:-0}"

if [[ ! -x "$exe" ]]; then
  echo "missing executable: $exe" >&2
  echo "build it with: make" >&2
  exit 1
fi

if [[ -e "$case_dir" ]]; then
  echo "case directory already exists: $case_dir" >&2
  echo "choose another name with CASE_DIR=/path/to/case $0" >&2
  exit 1
fi

mkdir -p "$case_dir"
cd "$case_dir"
mkdir -p restart_dir1 restart_dir2 restart_dir OVERFLOW_restart_dir1 OVERFLOW_restart_dir2 OVERFLOW_restart_dir phi_dir trackp_dir

cat > gtc.input <<EOF
&input_parameters
  irun=0,
  mstep=${mstep},
  msnap=${msnap},
  ndiag=${ndiag},
  nonlinear=1.0,
  paranl=0.0,
  mode00=1,
  tstep=0.2,
  micell=20,
  mpsi=50,
  mthetamax=300,
  mzetamax=10,
  npartdom=1,
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
  nbound=10,
  umax=4.0,
  iload=0,
  track_particles=0,
  nptrack=0,
  rng_control=0,
  nmode=20,20,20,20,20,20,20,20,20,20,20,20,20,
  mmode=25,26,27,28,29,30,31,32,33,34,35,36,37
/
EOF

/usr/bin/time -p mpirun -np "$mpi_np" "$exe" > run.log 2> time.log

{
  echo "case,cyclone_base_nonlinear_electrostatic_zonal"
  echo "mpi_np,$mpi_np"
  echo "mstep,$mstep"
  echo "msnap,$msnap"
  echo "ndiag,$ndiag"
  awk '/^real /{print "time_real_sec," $2}' time.log
  awk '/^user /{print "time_user_sec," $2}' time.log
  awk '/^sys /{print "time_sys_sec," $2}' time.log
  awk -F: '/MAIN LOOP TIME/{gsub(/^[ \t]+|[ \t]+$/,"",$2); print "main_loop_wall_sec," $2}' run.log
  awk -F: '/TOTAL WALL CLOCK TIME/{gsub(/^[ \t]+|[ \t]+$/,"",$2); print "program_wall_sec," $2}' run.log
  find . -maxdepth 1 -name 'snap*.out' | wc -l | awk '{print "snapshot_count," $1}'
  find restart_dir1 restart_dir2 -type f -name 'restart_*.bp' | wc -l | awk '{print "restart_file_count," $1}'
} > summary.csv

if [[ "$keep_restart" != "1" ]]; then
  rm -rf restart_dir1 restart_dir2 restart_dir OVERFLOW_restart_dir1 OVERFLOW_restart_dir2 OVERFLOW_restart_dir phi_dir trackp_dir
  echo "restart_files_retained,0" >> summary.csv
else
  echo "restart_files_retained,1" >> summary.csv
fi

if command -v python3 >/dev/null 2>&1; then
  python3 "$root/plot_history_and_snapshots.py" "$case_dir" > plot.log 2>&1 || true
fi

echo "case output: $case_dir"
echo "summary: $case_dir/summary.csv"
