#include "gtc.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static int read_int_line_diag(FILE *fp, int *value) {
  char line[256];
  if (!fgets(line, sizeof(line), fp)) return 0;
  return sscanf(line, "%d", value) == 1;
}

static int read_real_line_diag(FILE *fp, GtcReal *value) {
  char line[256];
  if (!fgets(line, sizeof(line), fp)) return 0;
  return sscanf(line, "%f", value) == 1;
}

static void write_history_restart_prefix(GtcState *s) {
  FILE *old = fopen("history.out", "r");
  if (!old) return;
  int old_irun = 0, mquantity = 0, mflux = 0, nmode = 0, old_outputs = 0;
  if (!read_int_line_diag(old, &old_irun) || !read_int_line_diag(old, &mquantity) ||
      !read_int_line_diag(old, &mflux) || !read_int_line_diag(old, &nmode) ||
      !read_int_line_diag(old, &old_outputs)) {
    fclose(old);
    return;
  }

  char backup[32];
  snprintf(backup, sizeof(backup), "histry%d.bak", old_irun % 10);
  FILE *bak = fopen(backup, "w");
  if (!bak) {
    fclose(old);
    return;
  }
  fprintf(bak, "%6d\n%6d\n%6d\n%6d\n%6d\n", old_irun, mquantity, mflux, nmode, old_outputs);
  const int old_real_count = (mquantity + mflux + 4 * nmode) * old_outputs + 1;
  for (int i = 0; i < old_real_count; i++) {
    GtcReal value = 0.0;
    if (!read_real_line_diag(old, &value)) break;
    gtc_fprintf_e(bak, 6, value);
    fputc('\n', bak);
  }
  fclose(bak);
  fclose(old);

  bak = fopen(backup, "r");
  FILE *history = fopen("history.out", "w");
  if (!bak || !history) {
    if (bak) fclose(bak);
    if (history) fclose(history);
    return;
  }
  if (!read_int_line_diag(bak, &old_irun) || !read_int_line_diag(bak, &mquantity) ||
      !read_int_line_diag(bak, &mflux) || !read_int_line_diag(bak, &nmode) ||
      !read_int_line_diag(bak, &old_outputs)) {
    fclose(history);
    fclose(bak);
    return;
  }
  s->p.irun = old_irun + 1;
  fprintf(history, "%6d\n%6d\n%6d\n%6d\n%6d\n", s->p.irun, mquantity, mflux, nmode,
          old_outputs + s->p.mstep / s->p.ndiag);
  const int copy_real_count = (mquantity + mflux + 4 * nmode) * old_outputs + 1;
  for (int i = 0; i < copy_real_count; i++) {
    GtcReal value = 0.0;
    if (!read_real_line_diag(bak, &value)) break;
    gtc_fprintf_e(history, 6, value);
    fputc('\n', history);
  }
  s->mstepall = old_outputs * s->p.ndiag;
  fclose(history);
  fclose(bak);
}

void diagnosis(GtcState *s) {
  const int first = s->istep == s->p.ndiag;
  const int history_modes = gtc_history_mode_count(s);
  const GtcReal vthi_diag = s->gyroradius * fabs(s->p.qion) / s->p.aion;
  const GtcReal tem_inv_diag = 1.0 / (s->p.aion * vthi_diag * vthi_diag);
  GtcReal local[15] = {0.0};
  local[0] = gtc_real(s->efield);
  local[1] = gtc_real(s->entropyi);
  local[2] = gtc_real(s->dflowi / vthi_diag);
  local[3] = gtc_real(s->pfluxi / vthi_diag);
  local[4] = gtc_real(s->efluxi * tem_inv_diag / vthi_diag);
  for (int i = 0; i < GTC_MFLUX; i++) local[5 + i] = gtc_real(s->eflux[i] * tem_inv_diag / vthi_diag);
  local[10] = gtc_real((GtcReal)s->mi);
  GtcReal global[15] = {0.0};
  MPI_Reduce(local, global, 15, GTC_MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  if (first) {
    if (s->rank == 0) {
      if (s->p.irun == 0) {
        s->mstepall = 0;
      } else {
        write_history_restart_prefix(s);
      }
    }
    MPI_Bcast(&s->mstepall, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  GtcReal *ddum_local = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*ddum_local), "diagnosis ddum local");
  GtcReal *ddum_global = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*ddum_global), "diagnosis ddum global");
  for (int i = 1; i <= s->p.mpsi; i++) {
    GtcReal sum = 0.0;
    for (int k = 1; k <= s->mzeta; k++) {
      for (int j = 1; j <= s->mtheta[i]; j++) {
        const GtcReal phi = s->phi[gtc_grid_index(s, k, s->igrid[i] + j)];
        sum += phi * phi;
      }
    }
    ddum_local[i] = gtc_real(sum);
  }
  MPI_Reduce(ddum_local, ddum_global, s->p.mpsi + 1, GTC_MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  if (s->rank != 0) {
    MPI_Barrier(MPI_COMM_WORLD);
    free(ddum_local);
    free(ddum_global);
    return;
  }

  FILE *history = fopen("history.out", first && s->p.irun == 0 ? "w" : "a");
  FILE *shear = fopen("sheareb.out", first && s->p.irun == 0 ? "w" : "a");
  if (!history || !shear) gtc_die(s, "cannot open diagnosis output");

  if (first && s->p.irun == 0) {
    fprintf(history, "%6d\n%6d\n%6d\n%6d\n%6d\n", s->p.irun, 12, GTC_MFLUX,
            history_modes, s->p.mstep / s->p.ndiag);
    gtc_fprintf_e(history, 6, s->p.tstep * (GtcReal)s->p.ndiag);
    fputc('\n', history);
    fprintf(shear, "%6d\n%6d\n", 4, s->p.mpsi);
  }
  if (first) {
    FILE *out = gtc_stdout_open(s, "a");
    for (int i = 1; i <= s->p.mpsi; i++) {
      const GtcReal rbin = s->p.a0 + (s->p.a1 - s->p.a0) * ((GtcReal)i - 0.5) / (GtcReal)s->p.mpsi;
      GtcReal rfac = s->p.rw * (rbin - s->p.rc);
      rfac = rfac * rfac;
      rfac = rfac * rfac * rfac;
      rfac = fmax(0.1, exp(-rfac));
      GtcReal kappa = 1.0;
      if (s->p.nbound == 0) kappa = 0.0;
      kappa = 1.0 - kappa + kappa * rfac;
      s->gradt[i] = gtc_real(1.0 / (kappa * s->p.kappati * s->gyroradius));
    }
    if (out) {
      for (int i = 0; i < GTC_MFLUX; i++) {
        const GtcReal rbin = s->p.a0 + (s->p.a1 - s->p.a0) * ((GtcReal)i + 0.5) / (GtcReal)GTC_MFLUX;
        GtcReal rfac = s->p.rw * (rbin - s->p.rc);
        rfac = rfac * rfac;
        rfac = rfac * rfac * rfac;
        rfac = fmax(0.1, exp(-rfac));
        GtcReal kappa = 1.0;
        if (s->p.nbound == 0) kappa = 0.0;
        kappa = 1.0 - kappa + kappa * rfac;
        fprintf(out, " kappa_T at radial_bin= %d % .8e\n", i + 1, kappa * s->p.kappati);
      }
      gtc_stdout_close(s, out);
    }
  }
  const GtcReal norm = global[10] > 1.0 ? global[10] : 1.0;
  GtcReal row[12] = {0.0};
  const GtcReal vthi = vthi_diag;
  const GtcReal tem_inv = tem_inv_diag;
  row[0] = s->ddeni;
  row[1] = s->eradial;
  row[2] = sqrt(global[0] / (GtcReal)(s->size > 0 ? s->size : 1));
  row[3] = global[1] / norm;
  row[4] = global[2] / norm;
  row[5] = global[3] / norm;
  row[6] = global[4] / norm;
  for (int i = 0; i < GTC_MFLUX; i++) {
    const GtcReal rbin = s->p.a0 + (s->p.a1 - s->p.a0) * ((GtcReal)i + 0.5) / (GtcReal)GTC_MFLUX;
    GtcReal rfac = s->p.rw * (rbin - s->p.rc);
    rfac = rfac * rfac;
    rfac = rfac * rfac * rfac;
    rfac = fmax(0.1, exp(-rfac));
    GtcReal kappa = 1.0;
    if (s->p.nbound == 0) kappa = 0.0;
    kappa = 1.0 - kappa + kappa * rfac;
    const GtcReal xnormal = 1.0 / (kappa * s->p.kappati * s->gyroradius);
    row[7 + i] = global[5 + i] * xnormal / (GtcReal)(s->size > 0 ? s->size : 1);
  }
  for (int i = 0; i < 12; i++) {
    gtc_fprintf_e(history, 6, row[i]);
    fputc('\n', history);
  }

  for (int mode = 0; mode < history_modes; mode++) {
    gtc_fprintf_e(history, 6, s->amp_mode[gtc_amp_mode_index(s, 0, mode, 0)]);
    fputc('\n', history);
    gtc_fprintf_e(history, 6, s->amp_mode[gtc_amp_mode_index(s, 1, mode, 0)]);
    fputc('\n', history);
    gtc_fprintf_e(history, 6, s->amp_mode[gtc_amp_mode_index(s, 0, mode, 1)]);
    fputc('\n', history);
    gtc_fprintf_e(history, 6, s->amp_mode[gtc_amp_mode_index(s, 1, mode, 1)]);
    fputc('\n', history);
  }
  fflush(history);

  for (int i = 1; i <= s->p.mpsi; i++) {
    gtc_fprintf_e(shear, 6, s->phip00[i] / vthi);
    fputc(i == s->p.mpsi ? '\n' : ' ', shear);
  }
  for (int i = 1; i <= s->p.mpsi; i++) {
    const GtcReal ddum = tem_inv * sqrt(ddum_global[i] / (GtcReal)(s->p.mzetamax * s->mtheta[i]));
    gtc_fprintf_e(shear, 6, ddum);
    fputc(i == s->p.mpsi ? '\n' : ' ', shear);
  }
  for (int i = 1; i <= s->p.mpsi; i++) {
    gtc_fprintf_e(shear, 6, s->hfluxpsi[i] * s->gradt[i] * tem_inv / vthi);
    fputc(i == s->p.mpsi ? '\n' : ' ', shear);
  }
  FILE *out = gtc_stdout_open(s, "a");
  if (out) {
    fprintf(out, " %d % .8e % .8e % .8e % .8e % .8e % .8e\n",
            s->istep + s->mstepall, row[2], s->eradial, row[9],
            s->total_field_energy[0], s->total_field_energy[1], s->total_field_energy[2]);
    gtc_stdout_close(s, out);
  }

  fclose(history);
  fclose(shear);
  MPI_Barrier(MPI_COMM_WORLD);
  free(ddum_local);
  free(ddum_global);
}
