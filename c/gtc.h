#ifndef GTC_C_PORT_H
#define GTC_C_PORT_H

#include <mpi.h>
#ifdef GTC_USE_OPENMP
#include <omp.h>
#endif

#include <stddef.h>
#include <stdio.h>

#define GTC_MFLUX 5
#define GTC_M_POLIDAL 9
#define GTC_TINY 1.0e-20

#define GTC_PRAGMA(x) _Pragma(#x)
#ifdef GTC_USE_OPENMP
#define GTC_OMP_PARALLEL GTC_PRAGMA(omp parallel)
#define GTC_OMP_FOR_STATIC GTC_PRAGMA(omp for schedule(static))
#define GTC_OMP_FOR_DYNAMIC GTC_PRAGMA(omp for schedule(dynamic))
#define GTC_OMP_MASTER GTC_PRAGMA(omp master)
#define GTC_OMP_BARRIER GTC_PRAGMA(omp barrier)
#define GTC_OMP_PARALLEL_FOR_STATIC GTC_PRAGMA(omp parallel for schedule(static))
#define GTC_OMP_PARALLEL_FOR_DYNAMIC GTC_PRAGMA(omp parallel for schedule(dynamic))
#define GTC_OMP_PARALLEL_FOR_STATIC_REDUCTION(...) \
  GTC_PRAGMA(omp parallel for schedule(static) reduction(__VA_ARGS__))
#else
#define GTC_OMP_PARALLEL
#define GTC_OMP_FOR_STATIC
#define GTC_OMP_FOR_DYNAMIC
#define GTC_OMP_MASTER
#define GTC_OMP_BARRIER
#define GTC_OMP_PARALLEL_FOR_STATIC
#define GTC_OMP_PARALLEL_FOR_DYNAMIC
#define GTC_OMP_PARALLEL_FOR_STATIC_REDUCTION(...)
#endif

enum {
  GTC_SPECTRUM_FULL_N = 0,
  GTC_SPECTRUM_SINGLE_N = 1
};

typedef float GtcReal;
#define GTC_MPI_REAL MPI_FLOAT

enum {
  GTC_Z_PSI = 0,
  GTC_Z_THETA = 1,
  GTC_Z_ZETA = 2,
  GTC_Z_U = 3,
  GTC_Z_WEIGHT = 4,
  GTC_Z_MU = 5,
  GTC_NPARAM = 6
};

typedef struct {
  int irun, mstep, msnap, ndiag, mode00, micell, mpsi, mthetamax, mzetamax;
  int npartdom, stdout_unit, nbound, iload, rng_control, spectrum_mode;
  int num_mode;
  GtcReal nonlinear, paranl, tstep, a, a0, a1, q0, q1, q2, rc, rw;
  GtcReal aion, qion, kappati, kappate, kappan, tite, flow0, flow1, flow2;
  GtcReal r0, b0, temperature, edensity0, umax;
  int *nmode;
  int *mmode;
} GtcParameters;

typedef struct {
  unsigned long long state;
  int index;
  double array[100];
} GtcRng;

typedef struct {
  GtcParameters p;
  int rank, size, istep, irk;
  int mzeta, mgrid, mi, mimax, mstepall, isnap, irest, file_exit;
  int mtdiag, idiag1, idiag2;
  int ntoroidal;
  int particle_domain_location, toroidal_domain_location;
  int myrank_partd, myrank_toroidal, left_pe, right_pe;
  MPI_Comm partd_comm, toroidal_comm;
  GtcReal pi, deltar, deltaz, zetamin, zetamax, gyroradius;

  int *mtheta, *igrid, *itran;
  GtcReal *deltat, *qtinv, *pmarki, *phi00, *phip00, *zonali, *gradt;
  GtcReal *rden, *rtemi, *hfluxpsi, *pfluxpsi, *rdtemi;
  GtcReal *densityi, *phi, *markeri, *evector, *heatflux, *dtemper;
  GtcReal *pgyro, *tgyro;
  int *kzion, *jtion0, *jtion1, *jtp1, *jtp2;
  GtcReal *wzion, *wpion, *wtion0, *wtion1, *wpi;
  GtcReal *wtp1, *wtp2;
  GtcReal *zion, *zion0;

  GtcReal efluxi, pfluxi, ddeni, dflowi, entropyi, efield, eradial;
  GtcReal eflux[GTC_MFLUX];
  GtcReal *amp_mode;
  GtcReal *eigenmode;
  GtcReal total_field_energy[3];

  GtcRng rng;
} GtcState;

void gtc_die(GtcState *s, const char *message);
void *gtc_xcalloc(GtcState *s, size_t count, size_t size, const char *name);
void gtc_fprintf_e(FILE *fp, int precision, GtcReal value);
FILE *gtc_stdout_open(const GtcState *s, const char *mode);
void gtc_stdout_close(const GtcState *s, FILE *fp);
int gtc_mpi_init(int *argc, char ***argv);
static inline size_t gtc_grid_index(const GtcState *s, int kz, int ij) {
  return (size_t)ij * (size_t)(s->mzeta + 1) + (size_t)kz;
}

static inline size_t gtc_evector_index(const GtcState *s, int component, int kz, int ij) {
  return (((size_t)ij * (size_t)(s->mzeta + 1) + (size_t)kz) * (size_t)3) + (size_t)component;
}

static inline size_t gtc_interp_index(const GtcState *s, int side, int kz, int ij) {
  return (((size_t)ij * (size_t)(s->mzeta + 1) + (size_t)kz) * (size_t)2) + (size_t)side;
}

static inline size_t gtc_gyro_index(const GtcState *s, int larmor, int ij) {
  (void)s;
  return (size_t)ij * 4u + (size_t)larmor;
}

static inline size_t gtc_larmor_particle_index(const GtcState *s, int larmor, int m) {
  (void)s;
  return (size_t)m * 4u + (size_t)larmor;
}

static inline size_t gtc_amp_mode_index(const GtcState *s, int component, int mode, int field) {
  return (((size_t)component * (size_t)s->p.num_mode + (size_t)mode) * 2u) + (size_t)field;
}

static inline size_t gtc_eigenmode_index(const GtcState *s, int radial, int mode, int poloidal) {
  return (((size_t)(radial - 1) * (size_t)s->p.num_mode + (size_t)mode) *
          (size_t)GTC_M_POLIDAL) + (size_t)poloidal;
}

static inline int gtc_history_mode_count(const GtcState *s) {
  return (s->p.nonlinear < 0.5 && s->p.spectrum_mode == GTC_SPECTRUM_SINGLE_N) ? s->p.num_mode : 0;
}

static inline int gtc_clamp_int(int value, int lo, int hi) {
  if (value < lo) return lo;
  if (value > hi) return hi;
  return value;
}

static inline GtcReal gtc_clamp_real(GtcReal value, GtcReal lo, GtcReal hi) {
  if (value < lo) return lo;
  if (value > hi) return hi;
  return value;
}

static inline GtcReal gtc_real(GtcReal value) {
  return (GtcReal)value;
}

static inline size_t gtc_particle_index(const GtcState *s, int param, int m) {
  (void)s;
  return (size_t)m * (size_t)GTC_NPARAM + (size_t)param;
}

static inline int gtc_omp_max_threads(void) {
#ifdef GTC_USE_OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

static inline int gtc_omp_thread_num(void) {
#ifdef GTC_USE_OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

void gtc_rng_init(GtcRng *rng, unsigned long long seed);
GtcReal gtc_rng_uniform(GtcRng *rng);
GtcReal gtc_rng_normal(GtcRng *rng);

GtcReal psi2r(GtcReal pdum);
GtcReal r2psi(GtcReal rdum);
GtcReal bfield(GtcReal pdum, GtcReal tdum);
GtcReal boozer2x(GtcReal pdum, GtcReal tdum);
GtcReal boozer2z(GtcReal pdum, GtcReal tdum);

void setup(GtcState *s);
void load(GtcState *s);
void chargei(GtcState *s);
void smooth(GtcState *s, int iflag);
void poisson(GtcState *s, int iflag);
void poisson_initial(GtcState *s);
void field(GtcState *s);
void pushi(GtcState *s);
void shifti(GtcState *s);
void diagnosis(GtcState *s);
void snapshot(GtcState *s);
void restart_io(GtcState *s, const char *iop);
void set_random_zion(GtcState *s);
void rand_num_gen_init(GtcState *s);
void fftr1d(int isign, int irank, GtcReal scale, GtcReal *x, GtcReal *y, int icount);
void fftc1d(int isign, int irank, GtcReal scale, GtcReal *x);

void gtc_free(GtcState *s);

#endif
