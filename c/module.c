#include "gtc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void gtc_die(GtcState *s, const char *message) {
  int rank = 0;
  if (s) rank = s->rank;
  fprintf(stderr, "gtc_c error on rank %d: %s\n", rank, message);
  MPI_Abort(MPI_COMM_WORLD, 1);
}

int gtc_mpi_init(int *argc, char ***argv) {
#ifdef GTC_USE_OPENMP
  int provided = MPI_THREAD_SINGLE;
  MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
  if (provided < MPI_THREAD_FUNNELED) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      fprintf(stderr, "gtc_c error: MPI does not provide MPI_THREAD_FUNNELED for OpenMP build\n");
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return provided;
#else
  MPI_Init(argc, argv);
  return MPI_THREAD_SINGLE;
#endif
}

void *gtc_xcalloc(GtcState *s, size_t count, size_t size, const char *name) {
  void *ptr = calloc(count, size);
  if (!ptr) {
    char buffer[256];
    snprintf(buffer, sizeof(buffer), "cannot allocate %s", name);
    gtc_die(s, buffer);
  }
  return ptr;
}

void gtc_fprintf_e(FILE *fp, int precision, GtcReal value) {
  if (isnan(value)) {
    fputs("       NaN", fp);
    return;
  }
  if (value == 0.0) {
    fprintf(fp, "0.%0*dE+00", precision, 0);
    return;
  }
  const int negative = signbit(value);
  GtcReal ax = fabs(value);
  int exponent = (int)floor(log10(ax)) + 1;
  GtcReal scaled = ax / pow(10.0, (GtcReal)exponent);
  GtcReal unit = pow(10.0, (GtcReal)precision);
  long digits = lround(scaled * unit);
  if (digits >= (long)unit) {
    digits /= 10;
    exponent++;
  }
  if (negative) {
    fprintf(fp, "-.%0*ldE%+03d", precision, digits, exponent);
  } else {
    fprintf(fp, "0.%0*ldE%+03d", precision, digits, exponent);
  }
}

FILE *gtc_stdout_open(const GtcState *s, const char *mode) {
  if (!s || s->rank != 0) return NULL;
  if (s->p.stdout_unit != 0 && s->p.stdout_unit != 6) {
    FILE *fp = fopen("stdout.out", mode);
    return fp ? fp : stdout;
  }
  return stdout;
}

void gtc_stdout_close(const GtcState *s, FILE *fp) {
  if (!s || s->rank != 0 || !fp) return;
  if (fp != stdout && fp != stderr) fclose(fp);
}

void gtc_free(GtcState *s) {
  if (s->partd_comm != MPI_COMM_NULL && s->partd_comm != MPI_COMM_WORLD) MPI_Comm_free(&s->partd_comm);
  if (s->toroidal_comm != MPI_COMM_NULL && s->toroidal_comm != MPI_COMM_WORLD) MPI_Comm_free(&s->toroidal_comm);
  free(s->p.nmode);
  free(s->p.mmode);
  free(s->mtheta);
  free(s->igrid);
  free(s->itran);
  free(s->deltat);
  free(s->qtinv);
  free(s->pmarki);
  free(s->phi00);
  free(s->phip00);
  free(s->zonali);
  free(s->gradt);
  free(s->rden);
  free(s->rtemi);
  free(s->hfluxpsi);
  free(s->pfluxpsi);
  free(s->rdtemi);
  free(s->amp_mode);
  free(s->eigenmode);
  free(s->densityi);
  free(s->phi);
  free(s->markeri);
  free(s->evector);
  free(s->heatflux);
  free(s->dtemper);
  free(s->pgyro);
  free(s->tgyro);
  free(s->kzion);
  free(s->jtion0);
  free(s->jtion1);
  free(s->jtp1);
  free(s->jtp2);
  free(s->wzion);
  free(s->wpion);
  free(s->wtion0);
  free(s->wtion1);
  free(s->wpi);
  free(s->wtp1);
  free(s->wtp2);
  free(s->zion);
  free(s->zion0);
  memset(s, 0, sizeof(*s));
}
