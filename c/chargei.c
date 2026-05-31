#include "gtc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void chargei(GtcState *s) {
  const size_t field_size = (size_t)(s->mzeta + 1) * (size_t)s->mgrid;
  memset(s->densityi, 0, field_size * sizeof(*s->densityi));
  s->ddeni = 0.0;

  const GtcReal delr = 1.0 / s->deltar;
  const GtcReal delz = 1.0 / s->deltaz;
  const GtcReal smu_inv = sqrtf(s->p.aion) / (fabsf(s->p.qion) * s->gyroradius);
  const GtcReal pi2_inv = 0.5 / s->pi;
  const int mzeta1 = s->mzeta + 1;
  const int single_zeta_cell = s->mzeta == 1;
  GtcReal *restrict densityi = s->densityi;
  GtcReal *delt = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*delt), "charge delt");
  for (int i = 0; i <= s->p.mpsi; i++) delt[i] = gtc_real(2.0 * s->pi / s->deltat[i]);

GTC_OMP_PARALLEL_FOR_STATIC
  for (int m = 0; m < s->mi; m++) {
    const GtcReal *zion = &s->zion[(size_t)m * (size_t)GTC_NPARAM];
    const GtcReal psitmp = zion[GTC_Z_PSI];
    const GtcReal thetatmp = zion[GTC_Z_THETA];
    const GtcReal zetatmp = zion[GTC_Z_ZETA];
    const GtcReal rhoi = zion[GTC_Z_MU] * smu_inv;
    const GtcReal r = sqrtf(2.0f * psitmp);
    const int ip = gtc_clamp_int((int)((r - s->p.a0) * delr + 0.5f), 0, s->p.mpsi);
    const int jt = gtc_clamp_int((int)(thetatmp * pi2_inv * delt[ip] + 0.5f), 0, s->mtheta[ip]);
    const int ipjt = s->igrid[ip] + jt;
    GtcReal wz1 = (zetatmp - s->zetamin) * delz;
    if (single_zeta_cell) {
      s->kzion[m] = 0;
      s->wzion[m] = wz1;
    } else {
      int kk = gtc_clamp_int((int)wz1, 0, s->mzeta - 1);
      s->kzion[m] = kk;
      s->wzion[m] = wz1 - (GtcReal)kk;
    }

    for (int larmor = 0; larmor < 4; larmor++) {
      const GtcReal rdum = delr * fmaxf(0.0f, fminf(s->p.a1 - s->p.a0,
                                 r + rhoi * s->pgyro[gtc_gyro_index(s, larmor, ipjt)] - s->p.a0));
      int ii = gtc_clamp_int((int)rdum, 0, s->p.mpsi - 1);
      GtcReal wp1 = rdum - (GtcReal)ii;
      const size_t o = gtc_larmor_particle_index(s, larmor, m);
      s->wpion[o] = wp1;

      const GtcReal tflr = thetatmp + rhoi * s->tgyro[gtc_gyro_index(s, larmor, ipjt)];
      int im = ii;
      GtcReal tdum = pi2_inv * (tflr - zetatmp * s->qtinv[im]) + 10.0;
      tdum = (tdum - (GtcReal)((int)tdum)) * delt[im];
      int j00 = gtc_clamp_int((int)tdum, 0, s->mtheta[im] - 1);
      s->jtion0[o] = s->igrid[im] + j00;
      s->wtion0[o] = tdum - (GtcReal)j00;

      im = ii + 1;
      tdum = pi2_inv * (tflr - zetatmp * s->qtinv[im]) + 10.0;
      tdum = (tdum - (GtcReal)((int)tdum)) * delt[im];
      int j01 = gtc_clamp_int((int)tdum, 0, s->mtheta[im] - 1);
      s->jtion1[o] = s->igrid[im] + j01;
      s->wtion1[o] = tdum - (GtcReal)j01;
    }
  }

  if (s->istep == 0) {
    free(delt);
    return;
  }

  const GtcReal *restrict zion = s->zion;
  const GtcReal *restrict wzion = s->wzion;
  const GtcReal *restrict wpion = s->wpion;
  const GtcReal *restrict wtion0 = s->wtion0;
  const GtcReal *restrict wtion1 = s->wtion1;
  const int *restrict kzion = s->kzion;
  const int *restrict jtion0 = s->jtion0;
  const int *restrict jtion1 = s->jtion1;

  GtcReal *density_private = NULL;
  const int nthreads = gtc_omp_max_threads();
  if (nthreads > 1) {
    density_private = gtc_xcalloc(s, (size_t)nthreads * field_size, sizeof(*density_private),
                                  "charge OpenMP density");
  }

  if (density_private) {
GTC_OMP_PARALLEL
    {
      GtcReal *local_density =
          density_private + (size_t)gtc_omp_thread_num() * field_size;
      if (single_zeta_cell) {
GTC_OMP_FOR_STATIC
        for (int m = 0; m < s->mi; m++) {
          const GtcReal weight = zion[(size_t)m * (size_t)GTC_NPARAM + (size_t)GTC_Z_WEIGHT];
          const GtcReal wz1 = weight * wzion[m];
          const GtcReal wz0 = weight - wz1;
          for (int larmor = 0; larmor < 4; larmor++) {
            const size_t o = gtc_larmor_particle_index(s, larmor, m);
            const GtcReal wp1 = wpion[o];
            const GtcReal wp0 = 1.0 - wp1;
            const GtcReal wt10 = wp0 * wtion0[o];
            const GtcReal wt00 = wp0 - wt10;
            const GtcReal wt11 = wp1 * wtion1[o];
            const GtcReal wt01 = wp1 - wt11;
            size_t idx = (size_t)jtion0[o] * 2u;
            local_density[idx] += wz0 * wt00;
            local_density[idx + 1] += wz1 * wt00;
            idx += 2u;
            local_density[idx] += wz0 * wt10;
            local_density[idx + 1] += wz1 * wt10;
            idx = (size_t)jtion1[o] * 2u;
            local_density[idx] += wz0 * wt01;
            local_density[idx + 1] += wz1 * wt01;
            idx += 2u;
            local_density[idx] += wz0 * wt11;
            local_density[idx + 1] += wz1 * wt11;
          }
        }
      } else {
GTC_OMP_FOR_STATIC
        for (int m = 0; m < s->mi; m++) {
          const GtcReal weight = zion[(size_t)m * (size_t)GTC_NPARAM + (size_t)GTC_Z_WEIGHT];
          const int kk = kzion[m];
          const GtcReal wz1 = weight * wzion[m];
          const GtcReal wz0 = weight - wz1;
          for (int larmor = 0; larmor < 4; larmor++) {
            const size_t o = gtc_larmor_particle_index(s, larmor, m);
            const GtcReal wp1 = wpion[o];
            const GtcReal wp0 = 1.0 - wp1;
            const GtcReal wt10 = wp0 * wtion0[o];
            const GtcReal wt00 = wp0 - wt10;
            const GtcReal wt11 = wp1 * wtion1[o];
            const GtcReal wt01 = wp1 - wt11;
            size_t idx = (size_t)jtion0[o] * (size_t)mzeta1 + (size_t)kk;
            local_density[idx] += wz0 * wt00;
            local_density[idx + 1] += wz1 * wt00;
            idx += (size_t)mzeta1;
            local_density[idx] += wz0 * wt10;
            local_density[idx + 1] += wz1 * wt10;
            idx = (size_t)jtion1[o] * (size_t)mzeta1 + (size_t)kk;
            local_density[idx] += wz0 * wt01;
            local_density[idx + 1] += wz1 * wt01;
            idx += (size_t)mzeta1;
            local_density[idx] += wz0 * wt11;
            local_density[idx + 1] += wz1 * wt11;
          }
        }
      }
    }
GTC_OMP_PARALLEL_FOR_STATIC
    for (size_t idx = 0; idx < field_size; idx++) {
      GtcReal value = 0.0;
      for (int tid = 0; tid < nthreads; tid++) {
        value += density_private[(size_t)tid * field_size + idx];
      }
      densityi[idx] = gtc_real(value);
    }
    free(density_private);
  } else if (single_zeta_cell) {
    for (int m = 0; m < s->mi; m++) {
      const GtcReal weight = zion[(size_t)m * (size_t)GTC_NPARAM + (size_t)GTC_Z_WEIGHT];
      const GtcReal wz1 = weight * wzion[m];
      const GtcReal wz0 = weight - wz1;
      for (int larmor = 0; larmor < 4; larmor++) {
        const size_t o = gtc_larmor_particle_index(s, larmor, m);
        const GtcReal wp1 = wpion[o];
        const GtcReal wp0 = 1.0 - wp1;
        const GtcReal wt10 = wp0 * wtion0[o];
        const GtcReal wt00 = wp0 - wt10;
        const GtcReal wt11 = wp1 * wtion1[o];
        const GtcReal wt01 = wp1 - wt11;
        size_t idx = (size_t)jtion0[o] * 2u;
        densityi[idx] += wz0 * wt00;
        densityi[idx + 1] += wz1 * wt00;
        idx += 2u;
        densityi[idx] += wz0 * wt10;
        densityi[idx + 1] += wz1 * wt10;
        idx = (size_t)jtion1[o] * 2u;
        densityi[idx] += wz0 * wt01;
        densityi[idx + 1] += wz1 * wt01;
        idx += 2u;
        densityi[idx] += wz0 * wt11;
        densityi[idx + 1] += wz1 * wt11;
      }
    }
  } else {
    for (int m = 0; m < s->mi; m++) {
      const GtcReal weight = zion[(size_t)m * (size_t)GTC_NPARAM + (size_t)GTC_Z_WEIGHT];
      const int kk = kzion[m];
      const GtcReal wz1 = weight * wzion[m];
      const GtcReal wz0 = weight - wz1;
      for (int larmor = 0; larmor < 4; larmor++) {
        const size_t o = gtc_larmor_particle_index(s, larmor, m);
        const GtcReal wp1 = wpion[o];
        const GtcReal wp0 = 1.0 - wp1;
        const GtcReal wt10 = wp0 * wtion0[o];
        const GtcReal wt00 = wp0 - wt10;
        const GtcReal wt11 = wp1 * wtion1[o];
        const GtcReal wt01 = wp1 - wt11;
        size_t idx = (size_t)jtion0[o] * (size_t)mzeta1 + (size_t)kk;
        densityi[idx] += wz0 * wt00;
        densityi[idx + 1] += wz1 * wt00;
        idx += (size_t)mzeta1;
        densityi[idx] += wz0 * wt10;
        densityi[idx + 1] += wz1 * wt10;
        idx = (size_t)jtion1[o] * (size_t)mzeta1 + (size_t)kk;
        densityi[idx] += wz0 * wt01;
        densityi[idx + 1] += wz1 * wt01;
        idx += (size_t)mzeta1;
        densityi[idx] += wz0 * wt11;
        densityi[idx + 1] += wz1 * wt11;
      }
    }
  }

  if (s->p.npartdom > 1) {
    GtcReal *dnitmp = gtc_xcalloc(s, field_size, sizeof(*dnitmp), "charge partd reduce");
    MPI_Allreduce(s->densityi, dnitmp, (int)field_size, GTC_MPI_REAL, MPI_SUM, s->partd_comm);
    memcpy(s->densityi, dnitmp, field_size * sizeof(*s->densityi));
    free(dnitmp);
  }

  if (single_zeta_cell) {
    for (int i = 0; i <= s->p.mpsi; i++) {
      const size_t edge = (size_t)(s->igrid[i] + s->mtheta[i]) * 2u;
      const size_t zero = (size_t)s->igrid[i] * 2u;
      densityi[edge] += densityi[zero];
      densityi[edge + 1] += densityi[zero + 1];
    }
  } else {
    for (int i = 0; i <= s->p.mpsi; i++) {
      const int edge = s->igrid[i] + s->mtheta[i];
      const int zero = s->igrid[i];
      for (int k = 0; k <= s->mzeta; k++) {
        densityi[(size_t)edge * (size_t)mzeta1 + (size_t)k] += densityi[(size_t)zero * (size_t)mzeta1 + (size_t)k];
      }
    }
  }

  GtcReal *sendl = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*sendl), "charge sendl");
  GtcReal *recvr = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*recvr), "charge recvr");
  if (single_zeta_cell) {
    for (int ij = 0; ij < s->mgrid; ij++) sendl[ij] = densityi[(size_t)ij * 2u];
  } else {
    for (int ij = 0; ij < s->mgrid; ij++) sendl[ij] = densityi[(size_t)ij * (size_t)mzeta1];
  }
  MPI_Sendrecv(sendl, s->mgrid, GTC_MPI_REAL, s->left_pe, s->myrank_toroidal,
               recvr, s->mgrid, GTC_MPI_REAL, s->right_pe, s->right_pe,
               s->toroidal_comm, MPI_STATUS_IGNORE);
  if (s->myrank_toroidal == s->ntoroidal - 1) {
    for (int i = 0; i <= s->p.mpsi; i++) {
      const int ii = s->igrid[i];
      const int jtmax = s->mtheta[i];
      for (int j = 1; j <= jtmax; j++) {
        const int shifted = ii + 1 + ((j - 1 + s->itran[i]) % jtmax);
        if (single_zeta_cell) {
          densityi[(size_t)(ii + j) * 2u + 1u] += recvr[shifted];
        } else {
          densityi[(size_t)(ii + j) * (size_t)mzeta1 + (size_t)s->mzeta] += recvr[shifted];
        }
      }
    }
  } else if (single_zeta_cell) {
    for (int ij = 0; ij < s->mgrid; ij++) densityi[(size_t)ij * 2u + 1u] += recvr[ij];
  } else {
    for (int ij = 0; ij < s->mgrid; ij++) densityi[(size_t)ij * (size_t)mzeta1 + (size_t)s->mzeta] += recvr[ij];
  }
  free(sendl);
  free(recvr);

  for (int i = 0; i < s->p.nbound; i++) {
    const GtcReal scale = s->p.nbound == 0 ? 1.0 : (GtcReal)i / (GtcReal)s->p.nbound;
    int il = i;
    int ir = s->p.mpsi - i;
    if (il < 0 || ir < 0 || il > s->p.mpsi || ir > s->p.mpsi) continue;
    if (single_zeta_cell) {
      for (int j = 0; j <= s->mtheta[il]; j++) {
        size_t idx = (size_t)(s->igrid[il] + j) * 2u;
        densityi[idx] *= scale;
        densityi[idx + 1] *= scale;
      }
      for (int j = 0; j <= s->mtheta[ir]; j++) {
        size_t idx = (size_t)(s->igrid[ir] + j) * 2u;
        densityi[idx] *= scale;
        densityi[idx + 1] *= scale;
      }
    } else {
      for (int k = 0; k <= s->mzeta; k++) {
        for (int j = 0; j <= s->mtheta[il]; j++) densityi[(size_t)(s->igrid[il] + j) * (size_t)mzeta1 + (size_t)k] *= scale;
        for (int j = 0; j <= s->mtheta[ir]; j++) densityi[(size_t)(s->igrid[ir] + j) * (size_t)mzeta1 + (size_t)k] *= scale;
      }
    }
  }

  for (int i = 0; i <= s->p.mpsi; i++) s->zonali[i] = 0.0;
  if (single_zeta_cell) {
    for (int i = 0; i <= s->p.mpsi; i++) {
      GtcReal zonal = 0.0;
      for (int j = 1; j <= s->mtheta[i]; j++) {
        size_t idx = (size_t)(s->igrid[i] + j) * 2u + 1u;
        zonal += 0.25 * densityi[idx];
        densityi[idx] = 0.25 * densityi[idx] * s->markeri[idx];
      }
      s->zonali[i] = zonal;
    }
  } else {
    for (int i = 0; i <= s->p.mpsi; i++) {
      GtcReal zonal = 0.0;
      for (int j = 1; j <= s->mtheta[i]; j++) {
        int ij = s->igrid[i] + j;
        for (int k = 1; k <= s->mzeta; k++) {
          size_t idx = (size_t)ij * (size_t)mzeta1 + (size_t)k;
          zonal += 0.25 * densityi[idx];
          densityi[idx] = 0.25 * densityi[idx] * s->markeri[idx];
        }
      }
      s->zonali[i] = zonal;
    }
  }

  GtcReal *adum = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*adum), "charge zonal reduce");
  MPI_Allreduce(s->zonali, adum, s->p.mpsi + 1, GTC_MPI_REAL, MPI_SUM, s->toroidal_comm);
  for (int i = 0; i <= s->p.mpsi; i++) s->zonali[i] = adum[i] * s->pmarki[i];
  free(adum);

  if (single_zeta_cell) {
    for (int i = 0; i <= s->p.mpsi; i++) {
      for (int j = 1; j <= s->mtheta[i]; j++) {
        densityi[(size_t)(s->igrid[i] + j) * 2u + 1u] -= s->zonali[i];
      }
      densityi[(size_t)s->igrid[i] * 2u + 1u] =
          densityi[(size_t)(s->igrid[i] + s->mtheta[i]) * 2u + 1u];
    }
  } else {
    for (int i = 0; i <= s->p.mpsi; i++) {
      for (int j = 1; j <= s->mtheta[i]; j++) {
        int ij = s->igrid[i] + j;
        for (int k = 1; k <= s->mzeta; k++) {
          densityi[(size_t)ij * (size_t)mzeta1 + (size_t)k] -= s->zonali[i];
        }
      }
      for (int k = 1; k <= s->mzeta; k++) {
        densityi[(size_t)s->igrid[i] * (size_t)mzeta1 + (size_t)k] =
            densityi[(size_t)(s->igrid[i] + s->mtheta[i]) * (size_t)mzeta1 + (size_t)k];
      }
    }
  }

  GtcReal rdum = 0.0;
  GtcReal tdum = 0.0;
  for (int i = 1; i < s->p.mpsi; i++) {
    const GtcReal r = s->p.a0 + s->deltar * (GtcReal)i;
    rdum += r;
    tdum += r * s->zonali[i];
  }
  tdum = rdum != 0.0 ? tdum / rdum : 0.0;
  s->ddeni = tdum;
  for (int i = 1; i < s->p.mpsi; i++) s->zonali[i] -= tdum;
  free(delt);
}
