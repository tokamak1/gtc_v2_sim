#include "gtc.h"

#include <math.h>
#include <stdlib.h>

void field(GtcState *s) {
  const size_t n = 3 * (size_t)(s->mzeta + 1) * (size_t)s->mgrid;
GTC_OMP_PARALLEL_FOR_STATIC
  for (size_t idx = 0; idx < n; idx++) s->evector[idx] = 0.0;

  const GtcReal diffr = 0.5 / s->deltar;
  const GtcReal diffz = 0.5 / s->deltaz;

GTC_OMP_PARALLEL_FOR_STATIC
  for (int k = 1; k <= s->mzeta; k++) {
    for (int i = 1; i < s->p.mpsi; i++) {
      const GtcReal r = s->p.a0 + s->deltar * (GtcReal)i;
      const GtcReal drdp = 1.0 / r;
      for (int j = 1; j <= s->mtheta[i]; j++) {
        const int ij = s->igrid[i] + j;
        size_t up = gtc_interp_index(s, 0, k, ij);
        size_t dn = gtc_interp_index(s, 1, k, ij);
        const int jup = s->jtp1[up];
        const int jdn = s->jtp1[dn];
        const GtcReal wup = s->wtp1[up];
        const GtcReal wdn = s->wtp1[dn];
        const GtcReal phi_up = (1.0 - wup) * s->phi[gtc_grid_index(s, k, jup)] +
                              wup * s->phi[gtc_grid_index(s, k, jup + 1)];
        const GtcReal phi_dn = (1.0 - wdn) * s->phi[gtc_grid_index(s, k, jdn)] +
                              wdn * s->phi[gtc_grid_index(s, k, jdn + 1)];
        s->evector[gtc_evector_index(s, 0, k, ij)] = gtc_real(drdp * diffr * (phi_up - phi_dn));
      }
    }
  }

GTC_OMP_PARALLEL_FOR_STATIC
  for (int i = 1; i < s->p.mpsi; i++) {
    const GtcReal difft = 0.5 / s->deltat[i];
    for (int k = 1; k <= s->mzeta; k++) {
      for (int j = 1; j <= s->mtheta[i]; j++) {
        const int ij = s->igrid[i] + j;
        const int jt = j + 1 - s->mtheta[i] * (j / s->mtheta[i]);
        s->evector[gtc_evector_index(s, 1, k, ij)] =
            gtc_real(difft * (s->phi[gtc_grid_index(s, k, s->igrid[i] + jt)] -
                              s->phi[gtc_grid_index(s, k, s->igrid[i] + j - 1)]));
      }
    }
  }

  GtcReal *sendr = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*sendr), "field sendr");
  GtcReal *sendl = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*sendl), "field sendl");
  GtcReal *recvl = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*recvl), "field recvl");
  GtcReal *recvr = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*recvr), "field recvr");
  for (int ij = 0; ij < s->mgrid; ij++) {
    sendr[ij] = s->phi[gtc_grid_index(s, s->mzeta, ij)];
    sendl[ij] = s->phi[gtc_grid_index(s, 1, ij)];
  }
  const int left = s->left_pe;
  const int right = s->right_pe;
  MPI_Sendrecv(sendr, s->mgrid, GTC_MPI_REAL, right, s->myrank_toroidal,
               recvl, s->mgrid, GTC_MPI_REAL, left, left,
               s->toroidal_comm, MPI_STATUS_IGNORE);
  MPI_Sendrecv(sendl, s->mgrid, GTC_MPI_REAL, left, s->myrank_toroidal,
               recvr, s->mgrid, GTC_MPI_REAL, right, right,
               s->toroidal_comm, MPI_STATUS_IGNORE);

GTC_OMP_PARALLEL_FOR_STATIC
  for (int i = 1; i < s->p.mpsi; i++) {
    const int ii = s->igrid[i];
    const int jtmax = s->mtheta[i];
    for (int j = 1; j <= jtmax; j++) {
      const int ij = ii + j;
      int left_j = j;
      int right_j = j;
      if (s->myrank_toroidal == 0) left_j = 1 + ((j - 1 - s->itran[i] + jtmax) % jtmax);
      if (s->myrank_toroidal == s->ntoroidal - 1) right_j = 1 + ((j - 1 + s->itran[i]) % jtmax);
      const GtcReal pleft = recvl[ii + left_j];
      const GtcReal pright = recvr[ii + right_j];
      if (s->mzeta == 1) {
        s->evector[gtc_evector_index(s, 2, 1, ij)] = gtc_real((pright - pleft) * diffz);
      } else if (s->mzeta == 2) {
        s->evector[gtc_evector_index(s, 2, 1, ij)] =
            gtc_real((s->phi[gtc_grid_index(s, 2, ij)] - pleft) * diffz);
        s->evector[gtc_evector_index(s, 2, 2, ij)] =
            gtc_real((pright - s->phi[gtc_grid_index(s, 1, ij)]) * diffz);
      } else {
        s->evector[gtc_evector_index(s, 2, 1, ij)] =
            gtc_real((s->phi[gtc_grid_index(s, 2, ij)] - pleft) * diffz);
        s->evector[gtc_evector_index(s, 2, s->mzeta, ij)] =
            gtc_real((pright - s->phi[gtc_grid_index(s, s->mzeta - 1, ij)]) * diffz);
        for (int k = 2; k < s->mzeta; k++) {
          s->evector[gtc_evector_index(s, 2, k, ij)] =
              gtc_real((s->phi[gtc_grid_index(s, k + 1, ij)] - s->phi[gtc_grid_index(s, k - 1, ij)]) * diffz);
        }
      }
    }
  }
  free(sendr);
  free(sendl);
  free(recvl);
  free(recvr);

GTC_OMP_PARALLEL_FOR_STATIC
  for (int i = 1; i < s->p.mpsi; i++) {
    const GtcReal r = s->p.a0 + s->deltar * (GtcReal)i;
    const GtcReal q = s->p.q0 + s->p.q1 * r / s->p.a + s->p.q2 * r * r / (s->p.a * s->p.a);
    const GtcReal delq = 1.0 / q - s->qtinv[i];
    for (int j = 1; j <= s->mtheta[i]; j++) {
      const int ij = s->igrid[i] + j;
      for (int k = 1; k <= s->mzeta; k++) {
        s->evector[gtc_evector_index(s, 2, k, ij)] =
            gtc_real(s->evector[gtc_evector_index(s, 2, k, ij)] + delq * s->evector[gtc_evector_index(s, 1, k, ij)]);
      }
    }
  }

  if (s->p.mode00 == 1) {
GTC_OMP_PARALLEL_FOR_STATIC
    for (int i = 1; i < s->p.mpsi; i++) {
      const GtcReal r = s->p.a0 + s->deltar * (GtcReal)i;
      for (int j = 1; j <= s->mtheta[i]; j++) {
        const int ij = s->igrid[i] + j;
        for (int k = 1; k <= s->mzeta; k++) {
          s->evector[gtc_evector_index(s, 0, k, ij)] =
              gtc_real(s->evector[gtc_evector_index(s, 0, k, ij)] + s->phip00[i] / r);
        }
      }
    }
  }

  GtcReal *sendrs = gtc_xcalloc(s, 3 * (size_t)s->mgrid, sizeof(*sendrs), "field sendrs");
  GtcReal *recvls = gtc_xcalloc(s, 3 * (size_t)s->mgrid, sizeof(*recvls), "field recvls");
  for (int ij = 0; ij < s->mgrid; ij++) {
    for (int c = 0; c < 3; c++) {
      sendrs[(size_t)c * (size_t)s->mgrid + (size_t)ij] = s->evector[gtc_evector_index(s, c, s->mzeta, ij)];
    }
  }
  MPI_Sendrecv(sendrs, 3 * s->mgrid, GTC_MPI_REAL, right, s->myrank_toroidal,
               recvls, 3 * s->mgrid, GTC_MPI_REAL, left, left,
               s->toroidal_comm, MPI_STATUS_IGNORE);

GTC_OMP_PARALLEL_FOR_STATIC
  for (int i = 1; i < s->p.mpsi; i++) {
    const int ii = s->igrid[i];
    const int jtmax = s->mtheta[i];
    for (int j = 1; j <= jtmax; j++) {
      const int ij0 = ii + j;
      const int srcj = s->myrank_toroidal == 0 ? 1 + ((j - 1 - s->itran[i] + jtmax) % jtmax) : j;
      const int ijs = ii + srcj;
      for (int c = 0; c < 3; c++) {
        s->evector[gtc_evector_index(s, c, 0, ij0)] =
            recvls[(size_t)c * (size_t)s->mgrid + (size_t)ijs];
      }
    }
  }

GTC_OMP_PARALLEL_FOR_STATIC
  for (int i = 1; i < s->p.mpsi; i++) {
    for (int k = 0; k <= s->mzeta; k++) {
      const int edge = s->igrid[i] + s->mtheta[i];
      const int zero = s->igrid[i];
      for (int c = 0; c < 3; c++) {
        s->evector[gtc_evector_index(s, c, k, zero)] = s->evector[gtc_evector_index(s, c, k, edge)];
      }
    }
  }
  free(sendrs);
  free(recvls);
}
