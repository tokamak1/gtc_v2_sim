#include "gtc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t phiflux_index(int mtdiag, int mz, int ri, int j, int k) {
  return (((size_t)ri * (size_t)mtdiag + (size_t)j) * (size_t)mz) + (size_t)k;
}

static int spectrum_filter_enabled(const GtcState *s) {
  return s->p.spectrum_mode == GTC_SPECTRUM_SINGLE_N;
}

static int spectrum_snapshot_enabled(const GtcState *s) {
  return s->eigenmode && s->irk == 2 && s->isnap > 0 && (s->istep % s->isnap) == 0;
}

static void apply_spectrum_control(GtcState *s, int record_history_modes, int record_eigenmode) {
  const int filter_modes = spectrum_filter_enabled(s);
  record_history_modes = record_history_modes && filter_modes && s->amp_mode;
  record_eigenmode = record_eigenmode && s->eigenmode && s->p.num_mode > 0;
  const int record_modes = record_history_modes || record_eigenmode;
  if (!filter_modes && !record_modes) return;
  const int mzbig = s->mtdiag / s->p.mzetamax > 1 ? s->mtdiag / s->p.mzetamax : 1;
  const int mz = s->mzeta * mzbig;
  const int mzg = s->p.mzetamax * mzbig;
  const int diag1 = record_eigenmode ? 1 : s->idiag1;
  const int diag2 = record_eigenmode ? s->p.mpsi : s->idiag2;
  const int nr = diag2 - diag1 + 1;
  if (mz <= 0 || mzg <= 0 || nr <= 0 || mzg != s->ntoroidal * mz) return;
  if (s->mtdiag % s->ntoroidal != 0) return;
  const int meachtheta = s->mtdiag / s->ntoroidal;

  const size_t local_count = (size_t)nr * (size_t)s->mtdiag * (size_t)mz;
  const size_t each_count = (size_t)meachtheta * (size_t)nr * (size_t)mz;
  GtcReal *phiflux = gtc_xcalloc(s, local_count, sizeof(*phiflux), "smooth phiflux");
  GtcReal *eachzeta = gtc_xcalloc(s, each_count, sizeof(*eachzeta), "smooth eachzeta");
  GtcReal *allzeta = gtc_xcalloc(s, each_count * (size_t)s->ntoroidal, sizeof(*allzeta), "smooth allzeta");
  const int nc = mzg / 2 + 1;
  GtcReal *xz = gtc_xcalloc(s, (size_t)mzg, sizeof(*xz), "smooth linear xz");
  GtcReal *yz = gtc_xcalloc(s, (size_t)2 * (size_t)nc, sizeof(*yz), "smooth linear yz");
  const int mode_count = meachtheta * s->p.num_mode;
  GtcReal *y_eigen = NULL;
  GtcReal *yt = NULL;
  GtcReal *ye = NULL;
  GtcReal *mode_theta = NULL;
  GtcReal *mode_theta_all = NULL;
  if (record_history_modes) {
    y_eigen = gtc_xcalloc(s, (size_t)2 * (size_t)mode_count, sizeof(*y_eigen), "smooth y_eigen");
    yt = gtc_xcalloc(s, (size_t)2 * (size_t)mode_count * (size_t)s->ntoroidal, sizeof(*yt), "smooth yt");
  }
  if (record_modes) {
    ye = gtc_xcalloc(s, (size_t)2 * (size_t)s->mtdiag, sizeof(*ye), "smooth ye");
  }
  if (record_eigenmode) {
    mode_theta = gtc_xcalloc(s, (size_t)2 * (size_t)nr * (size_t)s->p.num_mode *
                                    (size_t)meachtheta,
                              sizeof(*mode_theta), "smooth snapshot mode theta");
    mode_theta_all = gtc_xcalloc(s, (size_t)2 * (size_t)nr * (size_t)s->p.num_mode *
                                        (size_t)s->mtdiag,
                                  sizeof(*mode_theta_all), "smooth snapshot mode theta all");
  }
  unsigned char *keep_mode = NULL;
  if (filter_modes) {
    keep_mode = gtc_xcalloc(s, (size_t)mzg, sizeof(*keep_mode), "smooth spectrum filter");
    for (int imode = 0; imode < s->p.num_mode; imode++) {
      const int mode = s->p.nmode[imode];
      if (mode < 0 || mode >= nc) continue;
      keep_mode[mode] = 1;
    }
  }

  const GtcReal dt = 2.0 * s->pi / (GtcReal)s->mtdiag;
  const GtcReal pi2_inv = 0.5 / s->pi;
GTC_OMP_PARALLEL_FOR_STATIC
  for (int k = 1; k <= s->mzeta; k++) {
    for (int kz = 1; kz <= mzbig; kz++) {
      const int lk = (k - 1) * mzbig + (kz - 1);
      const GtcReal wz = (GtcReal)kz / (GtcReal)mzbig;
      const GtcReal zdum = s->zetamin + s->deltaz * ((GtcReal)(k - 1) + wz);
      for (int i = diag1; i <= diag2; i++) {
        const int ri = i - diag1;
        const int ii = s->igrid[i];
        for (int j = 1; j <= s->mtdiag; j++) {
          GtcReal tdum = pi2_inv * (dt * (GtcReal)j - zdum * s->qtinv[i]) + 10.0;
          tdum = (tdum - floor(tdum)) * (GtcReal)s->mtheta[i];
          int jt = gtc_clamp_int((int)tdum, 0, s->mtheta[i] - 1);
          const GtcReal wt = tdum - (GtcReal)jt;
          const GtcReal phik =
              (1.0 - wt) * s->phi[gtc_grid_index(s, k, ii + jt)] +
              wt * s->phi[gtc_grid_index(s, k, ii + jt + 1)];
          const GtcReal phim =
              (1.0 - wt) * s->phi[gtc_grid_index(s, k - 1, ii + jt)] +
              wt * s->phi[gtc_grid_index(s, k - 1, ii + jt + 1)];
          phiflux[phiflux_index(s->mtdiag, mz, ri, j - 1, lk)] =
              gtc_real(wz * phik + (1.0 - wz) * phim);
        }
      }
    }
  }

  for (int jpe = 0; jpe < s->ntoroidal; jpe++) {
    memset(eachzeta, 0, each_count * sizeof(*eachzeta));
    for (int jj = 0; jj < meachtheta; jj++) {
      const int j = jpe * meachtheta + jj;
      const size_t indt = (size_t)jj * (size_t)mz;
      for (int ri = 0; ri < nr; ri++) {
        const size_t indp1 = indt + (size_t)ri * (size_t)meachtheta * (size_t)mz;
        for (int k = 0; k < mz; k++) {
          eachzeta[indp1 + (size_t)k] = phiflux[phiflux_index(s->mtdiag, mz, ri, j, k)];
        }
      }
    }
    MPI_Gather(eachzeta, (int)each_count, GTC_MPI_REAL,
               allzeta, (int)each_count, GTC_MPI_REAL, jpe, s->toroidal_comm);
  }

  for (int jj = 0; jj < meachtheta; jj++) {
    const size_t indt1 = (size_t)jj * (size_t)mz;
    for (int ri = 0; ri < nr; ri++) {
      const size_t indt = indt1 + (size_t)ri * (size_t)meachtheta * (size_t)mz;
      for (int pe = 0; pe < s->ntoroidal; pe++) {
        const size_t pe_base = (size_t)pe * each_count;
        for (int k = 0; k < mz; k++) {
          xz[pe * mz + k] = allzeta[pe_base + indt + (size_t)k];
        }
      }
      memset(yz, 0, (size_t)2 * (size_t)nc * sizeof(*yz));
      fftr1d(1, mzg, 1.0, xz, yz, 1);
      if (record_history_modes && (ri + diag1) == s->p.mpsi / 2) {
        for (int mode = 0; mode < s->p.num_mode; mode++) {
          const int n = s->p.nmode[mode];
          if (n >= 0 && n < nc) {
            const size_t idx = (size_t)mode * (size_t)meachtheta + (size_t)jj;
            y_eigen[2 * idx] = gtc_real(yz[2 * n]);
            y_eigen[2 * idx + 1] = gtc_real(yz[2 * n + 1]);
          }
        }
      }
      if (record_eigenmode) {
        for (int mode = 0; mode < s->p.num_mode; mode++) {
          const int n = s->p.nmode[mode];
          if (n >= 0 && n < nc) {
            const size_t idx = (((size_t)ri * (size_t)s->p.num_mode + (size_t)mode) *
                                (size_t)meachtheta) + (size_t)jj;
            mode_theta[2 * idx] = gtc_real(yz[2 * n]);
            mode_theta[2 * idx + 1] = gtc_real(yz[2 * n + 1]);
          }
        }
      }
      if (filter_modes) {
        for (int mode = 0; mode < nc; mode++) {
          if (!keep_mode[mode]) {
            yz[2 * mode] = 0.0;
            yz[2 * mode + 1] = 0.0;
          }
        }
        fftr1d(-1, mzg, 1.0, xz, yz, 1);
        for (int pe = 0; pe < s->ntoroidal; pe++) {
          const size_t pe_base = (size_t)pe * each_count;
          for (int k = 0; k < mz; k++) {
            allzeta[pe_base + indt + (size_t)k] = gtc_real(xz[pe * mz + k]);
          }
        }
      }
    }
  }

  if (filter_modes) {
    for (int jpe = 0; jpe < s->ntoroidal; jpe++) {
      MPI_Scatter(allzeta, (int)each_count, GTC_MPI_REAL,
                  eachzeta, (int)each_count, GTC_MPI_REAL, jpe, s->toroidal_comm);
      for (int jj = 0; jj < meachtheta; jj++) {
        const int j = jpe * meachtheta + jj;
        const size_t indt = (size_t)jj * (size_t)mz;
        for (int ri = 0; ri < nr; ri++) {
          const size_t indp1 = indt + (size_t)ri * (size_t)meachtheta * (size_t)mz;
          for (int k = 0; k < mz; k++) {
            phiflux[phiflux_index(s->mtdiag, mz, ri, j, k)] = eachzeta[indp1 + (size_t)k];
          }
        }
      }
    }

    for (int k = 1; k <= s->mzeta; k++) {
      const GtcReal zdum = s->zetamin + s->deltaz * (GtcReal)k;
      const int lk = k * mzbig - 1;
      for (int i = diag1; i <= diag2; i++) {
        const int ri = i - diag1;
        const int ii = s->igrid[i];
        for (int j = 1; j <= s->mtheta[i]; j++) {
          GtcReal tdum = pi2_inv * (s->deltat[i] * (GtcReal)j + zdum * s->qtinv[i]) + 10.0;
          tdum = (tdum - floor(tdum)) * (GtcReal)s->mtdiag;
          int jt = gtc_clamp_int((int)tdum, 0, s->mtdiag - 1);
          const GtcReal wt = tdum - (GtcReal)jt;
          const int jtp = jt + 1;
          if (jt == 0) jt = s->mtdiag;
          const GtcReal right = phiflux[phiflux_index(s->mtdiag, mz, ri, jtp - 1, lk)];
          const GtcReal leftv = phiflux[phiflux_index(s->mtdiag, mz, ri, jt - 1, lk)];
          s->phi[gtc_grid_index(s, k, ii + j)] = gtc_real(wt * right + (1.0 - wt) * leftv);
        }
      }
    }

    GtcReal *sendr = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*sendr), "smooth spectrum sendr");
    GtcReal *recvl = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*recvl), "smooth spectrum recvl");
    for (int ij = 0; ij < s->mgrid; ij++) sendr[ij] = s->phi[gtc_grid_index(s, s->mzeta, ij)];
    MPI_Sendrecv(sendr, s->mgrid, GTC_MPI_REAL, s->right_pe, s->myrank_toroidal,
                 recvl, s->mgrid, GTC_MPI_REAL, s->left_pe, s->left_pe,
                 s->toroidal_comm, MPI_STATUS_IGNORE);
    if (s->myrank_toroidal == 0) {
      for (int i = diag1; i <= diag2; i++) {
        const int ii = s->igrid[i];
        const int jtmax = s->mtheta[i];
        for (int j = 1; j <= jtmax; j++) {
          const int srcj = 1 + ((j - 1 - s->itran[i] + jtmax) % jtmax);
          s->phi[gtc_grid_index(s, 0, ii + j)] = recvl[ii + srcj];
        }
      }
    } else {
      for (int ij = 0; ij < s->mgrid; ij++) s->phi[gtc_grid_index(s, 0, ij)] = recvl[ij];
    }
    for (int i = diag1; i <= diag2; i++) {
      for (int k = 0; k <= s->mzeta; k++) {
        s->phi[gtc_grid_index(s, k, s->igrid[i])] =
            s->phi[gtc_grid_index(s, k, s->igrid[i] + s->mtheta[i])];
      }
    }
    free(sendr);
    free(recvl);
  }

  if (record_history_modes) {
    MPI_Gather(y_eigen, 2 * mode_count, GTC_MPI_REAL,
               yt, 2 * mode_count, GTC_MPI_REAL, 0, s->toroidal_comm);
    if (s->myrank_toroidal == 0) {
      const GtcReal norm = (GtcReal)(mzg * s->mtdiag) * s->gyroradius * s->gyroradius;
      for (int mode = 0; mode < s->p.num_mode; mode++) {
        for (int pe = 0; pe < s->ntoroidal; pe++) {
          const size_t pe_base = (size_t)pe * (size_t)mode_count;
          for (int jj = 0; jj < meachtheta; jj++) {
            const size_t src = pe_base + (size_t)mode * (size_t)meachtheta + (size_t)jj;
            const int j = pe * meachtheta + jj;
            ye[2 * j] = yt[2 * src];
            ye[2 * j + 1] = yt[2 * src + 1];
          }
        }
        fftc1d(1, s->mtdiag, 1.0, ye);
        if (record_history_modes) {
          const int target = s->mtdiag - s->p.mmode[mode];
          if (target >= 0 && target < s->mtdiag) {
            s->amp_mode[gtc_amp_mode_index(s, 0, mode, 1)] = gtc_real(ye[2 * target] / norm);
            s->amp_mode[gtc_amp_mode_index(s, 1, mode, 1)] = gtc_real(ye[2 * target + 1] / norm);
          } else {
            s->amp_mode[gtc_amp_mode_index(s, 0, mode, 1)] = 0.0;
            s->amp_mode[gtc_amp_mode_index(s, 1, mode, 1)] = 0.0;
          }
        }
      }
    }
  }

  if (record_eigenmode) {
    const size_t theta_count = (size_t)nr * (size_t)s->p.num_mode * (size_t)meachtheta;
    MPI_Gather(mode_theta, (int)(2 * theta_count), GTC_MPI_REAL,
               mode_theta_all, (int)(2 * theta_count), GTC_MPI_REAL, 0, s->toroidal_comm);
    if (s->myrank_toroidal == 0) {
      const GtcReal norm = (GtcReal)(mzg * s->mtdiag) * s->gyroradius * s->gyroradius;
      memset(s->eigenmode, 0, (size_t)s->p.mpsi * (size_t)s->p.num_mode *
                                  (size_t)GTC_M_POLIDAL * sizeof(*s->eigenmode));
      for (int ri = 0; ri < nr; ri++) {
        const int radial = diag1 + ri;
        for (int mode = 0; mode < s->p.num_mode; mode++) {
          memset(ye, 0, (size_t)2 * (size_t)s->mtdiag * sizeof(*ye));
          for (int pe = 0; pe < s->ntoroidal; pe++) {
            const size_t pe_base = (size_t)pe * theta_count;
            for (int jj = 0; jj < meachtheta; jj++) {
              const size_t src = pe_base +
                  (((size_t)ri * (size_t)s->p.num_mode + (size_t)mode) *
                   (size_t)meachtheta) + (size_t)jj;
              const int j = pe * meachtheta + jj;
              ye[2 * j] = mode_theta_all[2 * src];
              ye[2 * j + 1] = mode_theta_all[2 * src + 1];
            }
          }
          fftc1d(1, s->mtdiag, 1.0, ye);
          for (int m = 0; m < GTC_M_POLIDAL; m++) {
            const int mt = m == 0 ? 0 : s->mtdiag - m;
            const GtcReal re = ye[2 * mt] / norm;
            const GtcReal im = ye[2 * mt + 1] / norm;
            s->eigenmode[gtc_eigenmode_index(s, radial, mode, m)] =
                gtc_real(sqrt(re * re + im * im));
          }
        }
      }
    }
  }

  free(mode_theta_all);
  free(mode_theta);
  free(ye);
  free(yt);
  free(y_eigen);
  free(keep_mode);
  free(yz);
  free(xz);
  free(allzeta);
  free(eachzeta);
  free(phiflux);
}

void smooth(GtcState *s, int iflag) {
  const size_t field_size = (size_t)(s->mzeta + 1) * (size_t)s->mgrid;
  GtcReal *src = iflag == 0 ? s->densityi : s->phi;
  GtcReal *phitmp = gtc_xcalloc(s, field_size, sizeof(*phitmp), "smooth tmp");
  GtcReal *work = gtc_xcalloc(s, field_size, sizeof(*work), "smooth work");
  memcpy(phitmp, src, field_size * sizeof(*phitmp));

  const int ismooth = s->p.nonlinear < 0.5 ? 0 : 1;
  for (int pass = 0; pass < ismooth; pass++) {
    for (int l = 0; l < 2; l++) {
GTC_OMP_PARALLEL_FOR_STATIC
      for (int i = 1; i < s->p.mpsi; i++) {
        for (int k = 1; k <= s->mzeta; k++) {
          phitmp[gtc_grid_index(s, k, s->igrid[i])] =
              phitmp[gtc_grid_index(s, k, s->igrid[i] + s->mtheta[i])];
        }
      }

      memcpy(work, phitmp, field_size * sizeof(*work));
GTC_OMP_PARALLEL_FOR_STATIC
      for (int k = 1; k <= s->mzeta; k++) {
        for (int ij = 0; ij < s->mgrid; ij++) {
          work[gtc_grid_index(s, k, ij)] = gtc_real(0.625 * phitmp[gtc_grid_index(s, k, ij)]);
        }
      }
GTC_OMP_PARALLEL_FOR_STATIC
      for (int k = 1; k <= s->mzeta; k++) {
        for (int i = 1; i < s->p.mpsi; i++) {
          for (int j = 1; j <= s->mtheta[i]; j++) {
            const int ij = s->igrid[i] + j;
            size_t up1 = gtc_interp_index(s, 0, k, ij);
            size_t dn1 = gtc_interp_index(s, 1, k, ij);
            const GtcReal radial1 =
                (1.0 - s->wtp1[up1]) * phitmp[gtc_grid_index(s, k, s->jtp1[up1])] +
                s->wtp1[up1] * phitmp[gtc_grid_index(s, k, s->jtp1[up1] + 1)] +
                (1.0 - s->wtp1[dn1]) * phitmp[gtc_grid_index(s, k, s->jtp1[dn1])] +
                s->wtp1[dn1] * phitmp[gtc_grid_index(s, k, s->jtp1[dn1] + 1)];
            size_t up2 = up1;
            size_t dn2 = dn1;
            const GtcReal radial2 =
                (1.0 - s->wtp2[up2]) * phitmp[gtc_grid_index(s, k, s->jtp2[up2])] +
                s->wtp2[up2] * phitmp[gtc_grid_index(s, k, s->jtp2[up2] + 1)] +
                (1.0 - s->wtp2[dn2]) * phitmp[gtc_grid_index(s, k, s->jtp2[dn2])] +
                s->wtp2[dn2] * phitmp[gtc_grid_index(s, k, s->jtp2[dn2] + 1)];
            work[gtc_grid_index(s, k, ij)] =
                gtc_real(work[gtc_grid_index(s, k, ij)] + 0.25 * radial1 - 0.0625 * radial2);
          }
        }
      }
      memcpy(phitmp, work, field_size * sizeof(*phitmp));

GTC_OMP_PARALLEL_FOR_STATIC
      for (int i = 1; i < s->p.mpsi; i++) {
        const int ii = s->igrid[i];
        const int jtmax = s->mtheta[i];
        for (int k = 1; k <= s->mzeta; k++) {
          for (int j = 1; j <= jtmax; j++) {
            const int jm1 = ii + 1 + ((j - 2 + jtmax) % jtmax);
            const int jp1 = ii + 1 + (j % jtmax);
            const int jm2 = ii + 1 + ((j - 3 + jtmax) % jtmax);
            const int jp2 = ii + 1 + ((j + 1) % jtmax);
            const int ij = ii + j;
            work[gtc_grid_index(s, k, ij)] =
                0.625 * phitmp[gtc_grid_index(s, k, ij)] +
                0.25 * (phitmp[gtc_grid_index(s, k, jm1)] + phitmp[gtc_grid_index(s, k, jp1)]) -
                0.0625 * (phitmp[gtc_grid_index(s, k, jm2)] + phitmp[gtc_grid_index(s, k, jp2)]);
          }
        }
      }
      memcpy(phitmp, work, field_size * sizeof(*phitmp));
    }

    GtcReal *sendr = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*sendr), "smooth sendr");
    GtcReal *sendl = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*sendl), "smooth sendl");
    GtcReal *recvl = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*recvl), "smooth recvl");
    GtcReal *recvr = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*recvr), "smooth recvr");
    for (int ij = 0; ij < s->mgrid; ij++) {
      sendr[ij] = phitmp[gtc_grid_index(s, s->mzeta, ij)];
      sendl[ij] = phitmp[gtc_grid_index(s, 1, ij)];
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
        const int ij_left = ii + left_j;
        const int ij_right = ii + right_j;
        for (int k = 1; k <= s->mzeta; k++) {
          GtcReal left_value = k == 1 ? recvl[ij_left] : phitmp[gtc_grid_index(s, k - 1, ij)];
          GtcReal right_value = k == s->mzeta ? recvr[ij_right] : phitmp[gtc_grid_index(s, k + 1, ij)];
          work[gtc_grid_index(s, k, ij)] =
              0.5 * phitmp[gtc_grid_index(s, k, ij)] + 0.25 * (left_value + right_value);
        }
      }
    }
    free(sendr);
    free(sendl);
    free(recvl);
    free(recvr);
    memcpy(phitmp, work, field_size * sizeof(*phitmp));
  }

  GtcReal *sendr = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*sendr), "smooth bc sendr");
  GtcReal *recvl = gtc_xcalloc(s, (size_t)s->mgrid, sizeof(*recvl), "smooth bc recvl");
  for (int ij = 0; ij < s->mgrid; ij++) sendr[ij] = phitmp[gtc_grid_index(s, s->mzeta, ij)];
  const int left = s->left_pe;
  const int right = s->right_pe;
  MPI_Sendrecv(sendr, s->mgrid, GTC_MPI_REAL, right, s->myrank_toroidal,
               recvl, s->mgrid, GTC_MPI_REAL, left, left,
               s->toroidal_comm, MPI_STATUS_IGNORE);
  if (s->myrank_toroidal == 0) {
GTC_OMP_PARALLEL_FOR_STATIC
    for (int i = 1; i < s->p.mpsi; i++) {
      const int ii = s->igrid[i];
      const int jtmax = s->mtheta[i];
      for (int j = 1; j <= jtmax; j++) {
        const int srcj = 1 + ((j - 1 - s->itran[i] + jtmax) % jtmax);
        phitmp[gtc_grid_index(s, 0, ii + j)] = recvl[ii + srcj];
      }
    }
  } else {
GTC_OMP_PARALLEL_FOR_STATIC
    for (int ij = 0; ij < s->mgrid; ij++) phitmp[gtc_grid_index(s, 0, ij)] = recvl[ij];
  }
  free(sendr);
  free(recvl);

GTC_OMP_PARALLEL_FOR_STATIC
  for (int i = 1; i < s->p.mpsi; i++) {
    for (int k = 0; k <= s->mzeta; k++) {
      phitmp[gtc_grid_index(s, k, s->igrid[i])] =
          phitmp[gtc_grid_index(s, k, s->igrid[i] + s->mtheta[i])];
    }
  }
GTC_OMP_PARALLEL_FOR_STATIC
  for (int i = 0; i <= s->mtheta[0]; i++) {
    for (int k = 0; k <= s->mzeta; k++) phitmp[gtc_grid_index(s, k, s->igrid[0] + i)] = 0.0;
  }
GTC_OMP_PARALLEL_FOR_STATIC
  for (int i = 0; i <= s->mtheta[s->p.mpsi]; i++) {
    for (int k = 0; k <= s->mzeta; k++) phitmp[gtc_grid_index(s, k, s->igrid[s->p.mpsi] + i)] = 0.0;
  }

  memcpy(src, phitmp, field_size * sizeof(*src));

  if (iflag == 3) {
    GtcReal *den00 = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*den00), "den00");
    for (int i = 0; i <= s->p.mpsi; i++) s->phip00[i] = s->p.qion * s->zonali[i];
    for (int pass = 0; pass < 2; pass++) {
      den00[0] = s->phip00[0];
      den00[s->p.mpsi] = s->phip00[s->p.mpsi];
      if (s->p.mpsi > 3) {
        den00[1] = s->phip00[3];
        den00[s->p.mpsi - 1] = s->phip00[s->p.mpsi - 3];
      }
      for (int i = 2; i <= s->p.mpsi - 2; i++) den00[i] = gtc_real(s->phip00[i - 2] + s->phip00[i + 2]);
      for (int i = 1; i < s->p.mpsi; i++) {
        den00[i] = gtc_real(0.625 * s->phip00[i] + 0.25 * (s->phip00[i - 1] + s->phip00[i + 1]) -
                             0.0625 * den00[i]);
      }
      for (int i = 0; i <= s->p.mpsi; i++) s->phip00[i] = den00[i];
    }
    for (int i = 0; i <= s->p.mpsi; i++) den00[i] = s->phip00[i];
    s->phip00[0] = 0.0;
    for (int i = 1; i <= s->p.mpsi; i++) {
      const GtcReal r = s->p.a0 + s->deltar * (GtcReal)i;
      s->phip00[i] = gtc_real(s->phip00[i - 1] + 0.5 * s->deltar * ((r - s->deltar) * den00[i - 1] + r * den00[i]));
    }
    for (int i = 0; i <= s->p.mpsi; i++) {
      const GtcReal r = s->p.a0 + s->deltar * (GtcReal)i;
      s->phip00[i] = gtc_real(-s->phip00[i] / r);
      s->phi00[i] = gtc_real(den00[i] * s->gyroradius * s->gyroradius);
    }
    for (int i = 1; i < s->p.mpsi; i++) s->phip00[i] = gtc_real(s->phip00[i] + 0.5 * (s->phi00[i + 1] - s->phi00[i - 1]) / s->deltar);
    s->phi00[0] = 0.0;
    for (int i = 1; i <= s->p.mpsi; i++) {
      s->phi00[i] = gtc_real(s->phi00[i - 1] + 0.5 * s->deltar * (s->phip00[i - 1] + s->phip00[i]));
    }
    if (s->p.mode00 == 0) {
      for (int i = 0; i <= s->p.mpsi; i++) s->phip00[i] = 0.0;
    }
    free(den00);
  }

  const int idiag = ((s->irk + 1) % 2) + (s->istep % s->p.ndiag);
  if (iflag == 3 && idiag == 0) {
    const int history_modes = gtc_history_mode_count(s);
    for (int i = 0; s->amp_mode && i < history_modes; i++) {
      s->amp_mode[gtc_amp_mode_index(s, 0, i, 1)] = 0.0;
      s->amp_mode[gtc_amp_mode_index(s, 1, i, 1)] = 0.0;
    }
  }
  if (iflag > 1) {
    const int history_diag = iflag == 3 && idiag == 0;
    const int snapshot_diag = iflag == 3 && spectrum_snapshot_enabled(s);
    apply_spectrum_control(s, history_diag && spectrum_filter_enabled(s), snapshot_diag);
  }

  if (iflag == 3 && idiag == 0) {
    s->eradial = 0.0;
    for (int i = 1; i <= s->p.mpsi; i++) s->eradial += s->phip00[i] * s->phip00[i];
    s->eradial = gtc_real(sqrt(s->eradial / (GtcReal)s->p.mpsi) / s->gyroradius);

    GtcReal local_phi2 = 0.0;
    long local_count = 0;
GTC_OMP_PARALLEL_FOR_STATIC_REDUCTION(+:local_phi2,local_count)
    for (int i = 0; i <= s->p.mpsi; i++) {
      for (int j = 1; j <= s->mtheta[i]; j++) {
        const int ij = s->igrid[i] + j;
        for (int k = 1; k <= s->mzeta; k++) {
          const GtcReal v = s->phi[gtc_grid_index(s, k, ij)];
          local_phi2 += v * v;
          local_count++;
        }
      }
    }
    const GtcReal unit = s->gyroradius * s->gyroradius;
    s->efield = local_count > 0 ? gtc_real(local_phi2 / ((GtcReal)local_count * unit * unit)) : 0.0;

    const int history_modes = gtc_history_mode_count(s);
    if (s->amp_mode && history_modes > 0) {
      GtcReal *x = gtc_xcalloc(s, (size_t)s->p.mpsi, sizeof(GtcReal), "smooth fft x");
      GtcReal *y = gtc_xcalloc(s, (size_t)(s->p.mpsi / 2 + 1) * 2, sizeof(GtcReal), "smooth fft y");
      for (int i = 0; i < s->p.mpsi; i++) x[i] = s->phip00[i + 1];
      fftr1d(1, s->p.mpsi, 1.0, x, y, 1);
      for (int i = 0; i < history_modes; i++) {
        const int mode = i + 1;
        if (mode <= s->p.mpsi / 2) {
          s->amp_mode[gtc_amp_mode_index(s, 0, i, 0)] = y[2 * mode] / ((GtcReal)s->p.mpsi * s->gyroradius);
          s->amp_mode[gtc_amp_mode_index(s, 1, i, 0)] = y[2 * mode + 1] / ((GtcReal)s->p.mpsi * s->gyroradius);
        } else {
          s->amp_mode[gtc_amp_mode_index(s, 0, i, 0)] = 0.0;
          s->amp_mode[gtc_amp_mode_index(s, 1, i, 0)] = 0.0;
        }
      }
      free(x);
      free(y);
    }
  }

  free(work);
  free(phitmp);
}
