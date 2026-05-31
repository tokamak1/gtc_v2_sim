#include "gtc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static int initialized = 0;
static int mindex_g = 65;
static int *nindex_g = NULL;
static int *indexp_g = NULL;
static GtcReal *ring_g = NULL;

static size_t ring_base(const GtcState *s, int k, int ij) {
  return ((size_t)k * (size_t)s->mgrid + (size_t)ij) * (size_t)mindex_g;
}

static void add_ring_point(const GtcState *s, int k, int ij0, int ij, GtcReal weight) {
  size_t base = ring_base(s, k, ij0);
  if (weight < 0.001) {
    ring_g[base] = gtc_real(ring_g[base] + weight);
    return;
  }
  int n = nindex_g[(size_t)k * (size_t)s->mgrid + (size_t)ij0];
  for (int nt = 0; nt < n; nt++) {
    if (indexp_g[base + (size_t)nt] == ij) {
      ring_g[base + (size_t)nt] = gtc_real(ring_g[base + (size_t)nt] + weight);
      return;
    }
  }
  if (n >= mindex_g) {
    gtc_die((GtcState *)s, "poisson ring index overflow");
  }
  indexp_g[base + (size_t)n] = ij;
  ring_g[base + (size_t)n] = gtc_real(weight);
  nindex_g[(size_t)k * (size_t)s->mgrid + (size_t)ij0] = n + 1;
}

void poisson_initial(GtcState *s) {
  if (initialized) return;
  const int mring = 2;
  const size_t planes = (size_t)(s->mzeta + 1) * (size_t)s->mgrid;
  nindex_g = gtc_xcalloc(s, planes, sizeof(int), "poisson nindex");
  indexp_g = gtc_xcalloc(s, planes * (size_t)mindex_g, sizeof(int), "poisson indexp");
  ring_g = gtc_xcalloc(s, planes * (size_t)mindex_g, sizeof(*ring_g), "poisson ring");

  const GtcReal vring[2] = {0.9129713024553, 2.233935334042};
  const GtcReal fring[2] = {0.7193896325719, 0.2806103674281};
  const GtcReal pi2_inv = 0.5 / s->pi;
  const GtcReal delr = 1.0 / s->deltar;

GTC_OMP_PARALLEL_FOR_DYNAMIC
  for (int k = 1; k <= s->mzeta; k++) {
    const GtcReal zdum = s->zetamin + s->deltaz * (GtcReal)k;
    for (int i = 0; i <= s->p.mpsi; i++) {
      for (int j = 1; j <= s->mtheta[i]; j++) {
        const int ij0 = s->igrid[i] + j;
        size_t base = ring_base(s, k, ij0);
        nindex_g[(size_t)k * (size_t)s->mgrid + (size_t)ij0] = 1;
        indexp_g[base] = ij0;
        ring_g[base] = 0.25;

        const GtcReal rgrid = s->p.a0 + s->deltar * (GtcReal)i;
        GtcReal tgrid = s->deltat[i] * (GtcReal)j + zdum * s->qtinv[i];
        tgrid = tgrid * pi2_inv;
        tgrid = 2.0 * s->pi * (tgrid - floor(tgrid));
        int jt = gtc_clamp_int((int)(pi2_inv * (2.0 * s->pi / s->deltat[i]) * tgrid + 0.5), 0, s->mtheta[i]);
        const GtcReal b = 1.0 / (1.0 + rgrid * cos(tgrid));
        const int ipjt = s->igrid[i] + jt;

        for (int kr = 0; kr < mring; kr++) {
          for (int kp = 0; kp < 8; kp++) {
            GtcReal ddelr = 0.0;
            GtcReal ddelt = 0.0;
            GtcReal wght = 0.0625 * fring[kr];
            if (kp < 4) {
              ddelr = s->pgyro[gtc_gyro_index(s, kp, ipjt)];
              ddelt = s->tgyro[gtc_gyro_index(s, kp, ipjt)];
            } else {
              wght = 0.125 * fring[kr];
              if (kp == 4) {
                ddelr = 0.5 * (s->pgyro[gtc_gyro_index(s, 0, ipjt)] + s->pgyro[gtc_gyro_index(s, 2, ipjt)]);
                ddelt = 0.5 * (s->tgyro[gtc_gyro_index(s, 0, ipjt)] + s->tgyro[gtc_gyro_index(s, 2, ipjt)]);
              } else if (kp == 5) {
                ddelr = 0.5 * (s->pgyro[gtc_gyro_index(s, 1, ipjt)] + s->pgyro[gtc_gyro_index(s, 2, ipjt)]);
                ddelt = 0.5 * (s->tgyro[gtc_gyro_index(s, 1, ipjt)] + s->tgyro[gtc_gyro_index(s, 2, ipjt)]);
              } else if (kp == 6) {
                ddelr = 0.5 * (s->pgyro[gtc_gyro_index(s, 1, ipjt)] + s->pgyro[gtc_gyro_index(s, 3, ipjt)]);
                ddelt = 0.5 * (s->tgyro[gtc_gyro_index(s, 1, ipjt)] + s->tgyro[gtc_gyro_index(s, 3, ipjt)]);
              } else {
                ddelr = 0.5 * (s->pgyro[gtc_gyro_index(s, 0, ipjt)] + s->pgyro[gtc_gyro_index(s, 3, ipjt)]);
                ddelt = 0.5 * (s->tgyro[gtc_gyro_index(s, 0, ipjt)] + s->tgyro[gtc_gyro_index(s, 3, ipjt)]);
              }
            }
            const GtcReal factor = 2.0 * vring[kr] * sqrt(0.5 / b);
            const GtcReal r = rgrid + ddelr * factor;
            const GtcReal t = tgrid + ddelt * factor;
            const GtcReal rdum = delr * fmax(0.0, fmin(s->p.a1 - s->p.a0, r - s->p.a0));
            const int ii = gtc_clamp_int((int)rdum, 0, s->p.mpsi - 1);
            GtcReal wr = rdum - (GtcReal)ii;
            if (wr > 0.95) wr = 1.0;
            if (wr < 0.05) wr = 0.0;

            GtcReal tdum = t - zdum * s->qtinv[ii + 1];
            tdum = tdum * pi2_inv + 10.0;
            tdum = (2.0 * s->pi / s->deltat[ii + 1]) * (tdum - floor(tdum));
            int j1 = gtc_clamp_int((int)tdum, 0, s->mtheta[ii + 1] - 1);
            GtcReal wt1 = tdum - (GtcReal)j1;
            if (wt1 > 0.95) wt1 = 1.0;
            if (wt1 < 0.05) wt1 = 0.0;

            tdum = t - zdum * s->qtinv[ii];
            tdum = tdum * pi2_inv + 10.0;
            tdum = (2.0 * s->pi / s->deltat[ii]) * (tdum - floor(tdum));
            int j0 = gtc_clamp_int((int)tdum, 0, s->mtheta[ii] - 1);
            GtcReal wt0 = tdum - (GtcReal)j0;
            if (wt0 > 0.95) wt0 = 1.0;
            if (wt0 < 0.05) wt0 = 0.0;

            add_ring_point(s, k, ij0, s->igrid[ii + 1] + j1 + 1, wght * wr * wt1);
            add_ring_point(s, k, ij0, s->igrid[ii + 1] + (j1 == 0 ? s->mtheta[ii + 1] : j1), wght * wr * (1.0 - wt1));
            add_ring_point(s, k, ij0, s->igrid[ii] + j0 + 1, wght * (1.0 - wr) * wt0);
            add_ring_point(s, k, ij0, s->igrid[ii] + (j0 == 0 ? s->mtheta[ii] : j0), wght * (1.0 - wr) * (1.0 - wt0));
          }
        }
      }
    }
  }
  int maxn = 0;
  int minn = mindex_g;
  long long nsum = 0;
  GtcReal maxsum = 0.0;
  GtcReal minsum = 1.0;
  for (int k = 1; k <= s->mzeta; k++) {
    for (int ij = 0; ij < s->mgrid; ij++) {
      const int n = nindex_g[(size_t)k * (size_t)s->mgrid + (size_t)ij];
      if (n > maxn) maxn = n;
      if (n < minn) minn = n;
      nsum += n;
      const size_t base = ring_base(s, k, ij);
      GtcReal total = 0.0;
      for (int jj = 0; jj < n; jj++) total += ring_g[base + (size_t)jj];
      if (total > maxsum) maxsum = total;
      if (total < minsum) minsum = total;
    }
  }
  if (s->rank == 0) {
    FILE *out = gtc_stdout_open(s, "a");
    fprintf(out, " poisson solver= %d %d %.8e %.8e %d %lld\n",
            maxn, minn, maxsum, minsum, s->mgrid, nsum);
    gtc_stdout_close(s, out);
  }
  initialized = 1;
}

void poisson(GtcState *s, int iflag) {
  (void)iflag;
  if (!initialized) poisson_initial(s);
  const GtcReal gamma = 0.75;
  const GtcReal tmp = 1.0 / (s->p.tite + 1.0 - gamma);
  int ipartd = 0;
  int nzeta = s->mzeta;
  int izeta1 = 1;
  int izeta2 = s->mzeta;
  if (s->p.npartdom > 1 && s->mzeta % s->p.npartdom == 0) {
    ipartd = 1;
    nzeta = s->mzeta / s->p.npartdom;
    izeta1 = s->myrank_partd * nzeta + 1;
    izeta2 = (s->myrank_partd + 1) * nzeta;
  }

  GTC_OMP_PARALLEL
  {
    GtcReal *phitmp = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*phitmp), "poisson phitmp");
    GtcReal *next = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*next), "poisson next");
GTC_OMP_FOR_DYNAMIC
    for (int k = izeta1; k <= izeta2; k++) {
      for (int ij = 0; ij < s->mgrid; ij++) phitmp[ij] = s->p.qion * s->densityi[gtc_grid_index(s, k, ij)] * tmp;
      for (int it = 2; it <= 5; it++) {
        memcpy(next, phitmp, (size_t)s->mgrid * sizeof(*next));
        for (int i = 0; i <= s->p.mpsi; i++) {
          for (int j = 1; j <= s->mtheta[i]; j++) {
            const int ij = s->igrid[i] + j;
            size_t base = ring_base(s, k, ij);
            int n = nindex_g[(size_t)k * (size_t)s->mgrid + (size_t)ij];
            GtcReal ptilde = 0.0;
            for (int jj = 0; jj < n; jj++) {
              ptilde = gtc_real(ptilde + ring_g[base + (size_t)jj] * phitmp[indexp_g[base + (size_t)jj]]);
            }
            const GtcReal perr = gtc_real(ptilde - gamma * phitmp[ij]);
            next[ij] = gtc_real((s->p.qion * s->densityi[gtc_grid_index(s, k, ij)] + perr) * tmp);
          }
        }
        for (int j = 0; j <= s->mtheta[0]; j++) next[s->igrid[0] + j] = 0.0;
        for (int j = 0; j <= s->mtheta[s->p.mpsi]; j++) next[s->igrid[s->p.mpsi] + j] = 0.0;
        GtcReal *swap = phitmp; phitmp = next; next = swap;
      }
      for (int ij = 0; ij < s->mgrid; ij++) s->phi[gtc_grid_index(s, k, ij)] = phitmp[ij];
    }
    free(phitmp);
    free(next);
  }

  if (ipartd) {
    const size_t send_count = (size_t)nzeta * (size_t)s->mgrid;
    const size_t recv_count = (size_t)s->mzeta * (size_t)s->mgrid;
    GtcReal *sendbuf = gtc_xcalloc(s, send_count, sizeof(*sendbuf), "poisson sendbuf");
    GtcReal *recvbuf = gtc_xcalloc(s, recv_count, sizeof(*recvbuf), "poisson recvbuf");
    for (int k = izeta1; k <= izeta2; k++) {
      const size_t plane = (size_t)(k - izeta1) * (size_t)s->mgrid;
      for (int ij = 0; ij < s->mgrid; ij++) {
        sendbuf[plane + (size_t)ij] = s->phi[gtc_grid_index(s, k, ij)];
      }
    }
    MPI_Allgather(sendbuf, (int)send_count, GTC_MPI_REAL,
                  recvbuf, (int)send_count, GTC_MPI_REAL, s->partd_comm);
    for (int k = 1; k <= s->mzeta; k++) {
      const size_t plane = (size_t)(k - 1) * (size_t)s->mgrid;
      for (int ij = 0; ij < s->mgrid; ij++) {
        s->phi[gtc_grid_index(s, k, ij)] = recvbuf[plane + (size_t)ij];
      }
    }
    free(sendbuf);
    free(recvbuf);
  }

GTC_OMP_PARALLEL_FOR_STATIC
  for (int k = 1; k <= s->mzeta; k++) {
    for (int i = 0; i <= s->p.mpsi; i++) {
      const GtcReal unit = s->rtemi[i] * (s->p.qion * s->gyroradius) * (s->p.qion * s->gyroradius) / s->p.aion;
      for (int j = 1; j <= s->mtheta[i]; j++) {
        const int ij = s->igrid[i] + j;
        s->phi[gtc_grid_index(s, k, ij)] *= unit;
      }
      s->phi[gtc_grid_index(s, k, s->igrid[i])] = s->phi[gtc_grid_index(s, k, s->igrid[i] + s->mtheta[i])];
    }
  }
}
