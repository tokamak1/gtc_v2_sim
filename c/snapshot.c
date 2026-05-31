#include "gtc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void snapshot(GtcState *s) {
  restart_io(s, "write");
  const int mbin_psi = 5;
  const int mbin_u = 41;
  const int jm = s->mtheta[s->p.mpsi / 2];
  const size_t ubin_count = (size_t)mbin_u * (size_t)mbin_psi;
  GtcReal *ubin = gtc_xcalloc(s, ubin_count, sizeof(*ubin), "snapshot ubin");
  GtcReal *dubin = gtc_xcalloc(s, ubin_count, sizeof(*dubin), "snapshot dubin");
  GtcReal *pbin = gtc_xcalloc(s, ubin_count, sizeof(*pbin), "snapshot pbin");
  GtcReal *dpbin = gtc_xcalloc(s, ubin_count, sizeof(*dpbin), "snapshot dpbin");
  GtcReal *ebin = gtc_xcalloc(s, ubin_count, sizeof(*ebin), "snapshot ebin");
  GtcReal *debin = gtc_xcalloc(s, ubin_count, sizeof(*debin), "snapshot debin");
  GtcReal *marker = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*marker), "snapshot marker");
  GtcReal *fflows = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*fflows), "snapshot fflows");
  GtcReal *dflows = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*dflows), "snapshot dflows");
  GtcReal *ftem = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*ftem), "snapshot ftem");
  GtcReal *dtem = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*dtem), "snapshot dtem");

  const GtcReal uright = s->p.umax;
  const GtcReal uleft = -uright;
  const GtcReal dela = (GtcReal)mbin_psi / (s->p.a1 - s->p.a0);
  const GtcReal delu = 1.0 / (uright - uleft);
  const GtcReal eright = 1.0 / (s->p.umax * s->p.umax / 4.0);
  const GtcReal vthi_inv = s->p.aion / (s->gyroradius * fabs(s->p.qion));
  const GtcReal aion_inv = 1.0 / s->p.aion;
  const GtcReal delr = 1.0 / s->deltar;

  for (int m = 0; m < s->mi; m++) {
    const GtcReal psi = s->zion[gtc_particle_index(s, GTC_Z_PSI, m)];
    const GtcReal theta = s->zion[gtc_particle_index(s, GTC_Z_THETA, m)];
    const GtcReal u = s->zion[gtc_particle_index(s, GTC_Z_U, m)];
    const GtcReal delta_weight = s->zion[gtc_particle_index(s, GTC_Z_WEIGHT, m)];
    const GtcReal mu = s->zion[gtc_particle_index(s, GTC_Z_MU, m)];
    const GtcReal weight = s->zion0[gtc_particle_index(s, GTC_Z_MU, m)];
    const GtcReal r = gtc_real(sqrt(gtc_real(2.0 * psi)));
    const int ibin = gtc_clamp_int(1 + (int)((r - s->p.a0) * dela), 1, mbin_psi) - 1;
    const int ip = gtc_clamp_int((int)((r - s->p.a0) * delr + 0.5), 0, s->p.mpsi);
    const GtcReal b = gtc_real(1.0 / gtc_real(1.0 + gtc_real(r * cos(theta))));
    const GtcReal upara = gtc_real(gtc_real(gtc_real(u * b) * s->p.qion) * gtc_real(aion_inv * vthi_inv));
    const GtcReal energy = fmax(1.0e-20, gtc_real(gtc_real(0.5 * upara * upara) +
                               gtc_real(gtc_real(mu * mu) * gtc_real(b * aion_inv * vthi_inv * vthi_inv))));
    const GtcReal pitch = gtc_real(upara / gtc_real(sqrt(gtc_real(2.0 * energy))));
    int iu = 1 + (int)((GtcReal)(mbin_u - 1) * (upara - uleft) * delu);
    if (iu >= 1 && iu <= mbin_u) {
      ubin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] = gtc_real(ubin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] + weight);
      dubin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] = gtc_real(dubin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] + delta_weight);
    }
    iu = 1 + (int)((GtcReal)(mbin_u - 1) * energy * eright);
    if (iu >= 1 && iu <= mbin_u) {
      ebin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] = gtc_real(ebin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] + weight);
      debin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] = gtc_real(debin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] + delta_weight);
    }
    iu = 1 + (int)((GtcReal)(mbin_u - 1) * 0.5 * (pitch + 1.0));
    if (iu >= 1 && iu <= mbin_u) {
      pbin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] = gtc_real(pbin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] + weight);
      dpbin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] = gtc_real(dpbin[(size_t)ibin * mbin_u + (size_t)(iu - 1)] + delta_weight);
    }
    marker[ip] = gtc_real(marker[ip] + weight);
    fflows[ip] = gtc_real(fflows[ip] + upara * weight);
    dflows[ip] = gtc_real(dflows[ip] + delta_weight * upara);
    ftem[ip] = gtc_real(ftem[ip] + energy * weight);
    dtem[ip] = gtc_real(dtem[ip] + delta_weight * energy);
  }

  GtcReal *tmp_u = gtc_xcalloc(s, ubin_count, sizeof(*tmp_u), "snapshot reduce u");
  GtcReal *tmp_r = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*tmp_r), "snapshot reduce r");
#define REDUCE_ARRAY(local, count)                         \
  do {                                                     \
    memset(tmp_u, 0, ubin_count * sizeof(*tmp_u));         \
    MPI_Reduce((local), tmp_u, (int)(count), GTC_MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD); \
    if (s->rank == 0) memcpy((local), tmp_u, (count) * sizeof(*(local))); \
  } while (0)
  REDUCE_ARRAY(ubin, ubin_count);
  REDUCE_ARRAY(dubin, ubin_count);
  REDUCE_ARRAY(pbin, ubin_count);
  REDUCE_ARRAY(dpbin, ubin_count);
  REDUCE_ARRAY(ebin, ubin_count);
  REDUCE_ARRAY(debin, ubin_count);
#undef REDUCE_ARRAY
#define REDUCE_R(local)                                    \
  do {                                                     \
    memset(tmp_r, 0, (size_t)(s->p.mpsi + 1) * sizeof(*tmp_r)); \
    MPI_Reduce((local), tmp_r, s->p.mpsi + 1, GTC_MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD); \
    if (s->rank == 0) memcpy((local), tmp_r, (size_t)(s->p.mpsi + 1) * sizeof(*(local))); \
  } while (0)
  REDUCE_R(marker);
  REDUCE_R(fflows);
  REDUCE_R(dflows);
  REDUCE_R(ftem);
  REDUCE_R(dtem);
#undef REDUCE_R

  GtcReal *midphi_local = gtc_xcalloc(s, (size_t)jm * (size_t)s->mzeta, sizeof(*midphi_local), "snapshot midphi local");
  GtcReal *midphi_all = NULL;
  for (int j = 1; j <= jm; j++) {
    const int ij = s->igrid[s->p.mpsi / 2] + j;
    for (int k = 1; k <= s->mzeta; k++) midphi_local[(size_t)(j - 1) * s->mzeta + (size_t)(k - 1)] = s->phi[gtc_grid_index(s, k, ij)];
  }
  if (s->myrank_toroidal == 0) {
    midphi_all = gtc_xcalloc(s, (size_t)jm * (size_t)s->mzeta * (size_t)s->ntoroidal,
                             sizeof(*midphi_all), "snapshot midphi all");
  }
  MPI_Gather(midphi_local, jm * s->mzeta, GTC_MPI_REAL, midphi_all, jm * s->mzeta, GTC_MPI_REAL,
             0, s->toroidal_comm);

  if (s->rank != 0) {
    free(ubin); free(dubin); free(pbin); free(dpbin); free(ebin); free(debin);
    free(marker); free(fflows); free(dflows); free(ftem); free(dtem);
    free(tmp_u); free(tmp_r); free(midphi_local); free(midphi_all);
    return;
  }

  for (int i = 0; i <= s->p.mpsi; i++) {
    const GtcReal denom = marker[i] < 1.0 ? 1.0 : marker[i];
    fflows[i] = gtc_real(fflows[i] / denom);
    dflows[i] = gtc_real(dflows[i] / denom);
    ftem[i] = gtc_real(ftem[i] / denom);
    dtem[i] = gtc_real(dtem[i] / denom);
    marker[i] = marker[i] < 1.0 ? 1.0 : marker[i];
    marker[i] = gtc_real(marker[i] * s->pmarki[i]);
  }
  GtcReal ubin_max = 0.0, pbin_max = 0.0, ebin_max = 0.0;
  for (size_t i = 0; i < ubin_count; i++) {
    if (ubin[i] > ubin_max) ubin_max = ubin[i];
    if (pbin[i] > pbin_max) pbin_max = pbin[i];
    if (ebin[i] > ebin_max) ebin_max = ebin[i];
  }
  for (size_t i = 0; i < ubin_count; i++) {
    dubin[i] = gtc_real(dubin[i] / ubin_max);
    dpbin[i] = gtc_real(dpbin[i] / pbin_max);
    debin[i] = gtc_real(debin[i] / ebin_max);
  }

  char name[64];
  snprintf(name, sizeof(name), "snap%05d.out", s->mstepall + s->istep);
  FILE *out = fopen(name, "w");
  if (!out) gtc_die(s, "cannot open snapshot output");

#define WRITE_E10(value)       \
  do {                         \
    gtc_fprintf_e(out, 4, value); \
    fputc('\n', out);          \
  } while (0)
  WRITE_E10((GtcReal)(s->mstepall + s->istep) * s->p.tstep);
  WRITE_E10(s->p.q0 + 0.5 * s->p.q1 + 0.25 * s->p.q2);
  fprintf(out, "%6d\n%6d\n%6d\n%6d\n%6d\n%6d\n", mbin_u, mbin_psi, s->p.mpsi, jm,
          s->p.mzetamax, s->p.mpsi / 2);

  fprintf(out, "%6d\n", 8);
  for (int i = 0; i < mbin_u; i++) {
    WRITE_E10(-s->p.umax + 2.0 * s->p.umax * (GtcReal)i / (GtcReal)(mbin_u - 1));
  }
  for (int i = 0; i < mbin_u; i++) WRITE_E10(-1.0 + 2.0 * (GtcReal)i / (GtcReal)(mbin_u - 1));
  for (int i = 0; i < mbin_u; i++) WRITE_E10(0.25 * s->p.umax * s->p.umax * (GtcReal)i / (GtcReal)(mbin_u - 1));
  for (int block = 0; block < mbin_psi * 6; block++) {
    const int b = block / 6;
    const int q = block % 6;
    for (int i = 0; i < mbin_u; i++) {
      GtcReal value = 0.0;
      const size_t idx = (size_t)b * mbin_u + (size_t)i;
      if (q == 0) value = ubin[idx];
      if (q == 1) value = dubin[idx];
      if (q == 2) value = pbin[idx];
      if (q == 3) value = dpbin[idx];
      if (q == 4) value = ebin[idx];
      if (q == 5) value = debin[idx];
      WRITE_E10(value);
    }
  }

  fprintf(out, "%6d\n", 8);
  for (int i = 1; i <= s->p.mpsi; i++) WRITE_E10((s->p.a0 + s->deltar * (GtcReal)i) / s->gyroradius);
  for (int block = 0; block < 7; block++) {
    for (int i = 1; i <= s->p.mpsi; i++) {
      GtcReal value = s->phi00[i];
      if (block == 0) value = s->zonali[i];
      if (block == 1) value = s->phip00[i] / s->gyroradius;
      if (block == 2) value = marker[i];
      if (block == 3) value = fflows[i];
      if (block == 4) value = dflows[i];
      if (block == 5) value = ftem[i];
      if (block == 6) value = dtem[i];
      WRITE_E10(value);
    }
  }

  fprintf(out, "%6d\n", 3);
  for (int j = 0; j <= jm; j++) {
    int jj = j;
    if (j > jm / 2) jj = j - jm;
    for (int i = 1; i <= s->p.mpsi; i++) {
      const GtcReal pdum = r2psi(s->p.a0 + s->deltar * (GtcReal)i);
      const GtcReal tdum = s->pi + 2.0 * s->pi * (GtcReal)jj / (GtcReal)jm;
      WRITE_E10(boozer2x(pdum, tdum));
    }
  }
  for (int j = 0; j <= jm; j++) {
    int jj = j;
    if (j > jm / 2) jj = j - jm;
    for (int i = 1; i <= s->p.mpsi; i++) {
      const GtcReal pdum = r2psi(s->p.a0 + s->deltar * (GtcReal)i);
      const GtcReal tdum = s->pi + 2.0 * s->pi * (GtcReal)jj / (GtcReal)jm;
      WRITE_E10(boozer2z(pdum, tdum));
    }
  }
  for (int j = 0; j <= jm; j++) {
    int jj = j;
    if (j > jm / 2) jj = j - jm;
    for (int i = 1; i <= s->p.mpsi; i++) {
      const GtcReal tdum = s->pi + 2.0 * s->pi * (GtcReal)jj / (GtcReal)jm;
      int jt = 1 + (int)(tdum / s->deltat[i]);
      jt = gtc_clamp_int(jt, 1, s->mtheta[i]);
      const GtcReal wt = tdum / s->deltat[i] - (GtcReal)(jt - 1);
      const GtcReal value = (wt * s->phi[gtc_grid_index(s, 0, s->igrid[i] + jt)] +
                            (1.0 - wt) * s->phi[gtc_grid_index(s, 0, s->igrid[i] + jt - 1)]) /
                           (s->gyroradius * s->gyroradius);
      WRITE_E10(value);
    }
  }

  fprintf(out, "%6d\n", 1);
  for (int j = 1; j <= jm; j++) {
    for (int kz = 1; kz <= s->ntoroidal; kz++) {
      for (int k = 1; k <= s->mzeta; k++) {
        const size_t idx = (size_t)(kz - 1) * (size_t)jm * (size_t)s->mzeta +
                           (size_t)(j - 1) * (size_t)s->mzeta + (size_t)(k - 1);
        WRITE_E10(midphi_all[idx] / (s->gyroradius * s->gyroradius));
      }
    }
  }

  const int snapshot_modes = s->eigenmode ? s->p.num_mode : 0;
  fprintf(out, "%6d\n%6d\n", snapshot_modes, GTC_M_POLIDAL);
  for (int mode = 0; mode < snapshot_modes; mode++) {
    fprintf(out, "%6d\n", s->p.nmode[mode]);
  }
  for (int i = 1; i <= s->p.mpsi; i++) {
    for (int mode = 0; mode < snapshot_modes; mode++) {
      for (int m = 0; m < GTC_M_POLIDAL; m++) {
        WRITE_E10(s->eigenmode[gtc_eigenmode_index(s, i, mode, m)]);
      }
    }
  }
#undef WRITE_E10
  fclose(out);
  free(ubin); free(dubin); free(pbin); free(dpbin); free(ebin); free(debin);
  free(marker); free(fflows); free(dflows); free(ftem); free(dtem);
  free(tmp_u); free(tmp_r); free(midphi_local); free(midphi_all);

  if (s->istep == s->p.mstep) {
    time_t now = time(NULL);
    struct tm tmv;
    struct tm *tmp = localtime(&now);
    if (tmp) tmv = *tmp;
    else memset(&tmv, 0, sizeof(tmv));
    char date[16], clock_time[16];
    strftime(date, sizeof(date), "%Y%m%d", &tmv);
    strftime(clock_time, sizeof(clock_time), "%H%M%S.000", &tmv);
    FILE *fp = gtc_stdout_open(s, "a");
    if (fp) {
      fprintf(fp, " Program ends at DATE=%s TIME=%s\n", date, clock_time);
      gtc_stdout_close(s, fp);
    }
  }
}
