#include "gtc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  double prep;
  double interpolate;
  double advance;
  double boundary;
  double diagnostics;
  double total;
  long calls;
} PushiProfile;

static PushiProfile pushi_profile;

static int pushi_profile_enabled(void) {
  static int initialized = 0;
  static int enabled = 0;
  if (!initialized) {
    const char *value = getenv("GTC_PUSHI_PROFILE");
    enabled = value && value[0] != '\0' && strcmp(value, "0") != 0;
    initialized = 1;
  }
  return enabled;
}

static void pushi_profile_report(GtcState *s) {
  if (s->istep != s->p.mstep || s->irk != 2) return;

  const double local_times[2] = {pushi_profile.advance, pushi_profile.total};
  double sum_times[2] = {0.0, 0.0};
  double min_times[2] = {0.0, 0.0};
  double max_times[2] = {0.0, 0.0};
  MPI_Reduce(local_times, sum_times, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(local_times, min_times, 2, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(local_times, max_times, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (s->rank != 0) return;
  FILE *out = gtc_stdout_open(s, "a");
  if (!out) return;
  fprintf(out, " PUSHI PROFILE (rank0 wall sec, calls=%ld): prep %.6f interpolate %.6f advance %.6f boundary %.6f diagnostics %.6f total %.6f\n",
          pushi_profile.calls, pushi_profile.prep, pushi_profile.interpolate,
          pushi_profile.advance, pushi_profile.boundary, pushi_profile.diagnostics,
          pushi_profile.total);
  fprintf(out, " PUSHI PROFILE WORLD (ranks=%d): advance avg %.6f min %.6f max %.6f total avg %.6f min %.6f max %.6f\n",
          s->size, sum_times[0] / (double)s->size, min_times[0], max_times[0],
          sum_times[1] / (double)s->size, min_times[1], max_times[1]);
  gtc_stdout_close(s, out);
}

static GtcReal modulo_positive(GtcReal x, GtcReal period) {
  if (x >= period) {
    x -= period;
    if (x >= period) x = fmodf(x, period);
  } else if (x < 0.0) {
    x += period;
    if (x < 0.0) {
      x = fmodf(x, period);
      if (x < 0.0) x += period;
    }
  }
  return x;
}

static inline void sincos_real(GtcReal x, GtcReal *s, GtcReal *c) {
#if defined(__clang__) || defined(__GNUC__)
  __builtin_sincosf(x, s, c);
#else
  *s = sinf(x);
  *c = cosf(x);
#endif
}

void pushi(GtcState *s) {
  const int profile = pushi_profile_enabled();
  const double profile_start = profile ? MPI_Wtime() : 0.0;
  double profile_tick = profile_start;
  const GtcReal delr = 1.0 / s->deltar;
  const GtcReal pi2 = 2.0 * s->pi;
  const GtcReal psimax = 0.5 * s->p.a1 * s->p.a1;
  const GtcReal psimin = 0.5 * s->p.a0 * s->p.a0;
  const GtcReal cmratio = s->p.qion / s->p.aion;
  const GtcReal cinv = 1.0 / s->p.qion;
  const GtcReal vthi = s->gyroradius * fabs(s->p.qion) / s->p.aion;
  const GtcReal ainv = 1.0 / s->p.a;
  const GtcReal sbound = s->p.nbound == 0 ? 0.0 : 1.0;
  const GtcReal dtime = s->irk == 1 ? 0.5 * s->p.tstep : s->p.tstep;
  const int linear_orbit = s->p.nonlinear < 0.5 && s->p.paranl < 0.5 && s->p.flow0 == 0.0;
  GtcReal *temp_inv = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*temp_inv), "pushi temp inv");
  const GtcReal temp_scale = 1.0 / (s->p.aion * vthi * vthi);
GTC_OMP_PARALLEL_FOR_STATIC
  for (int i = 0; i <= s->p.mpsi; i++) temp_inv[i] = gtc_real(temp_scale / s->rtemi[i]);
  GtcReal *vdrtmp = NULL;
  if (s->p.nonlinear > 0.5) {
    vdrtmp = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*vdrtmp), "pushi vdrtmp");
  }

  if (s->irk != 1 && vdrtmp) {
GTC_OMP_PARALLEL_FOR_STATIC
    for (int i = 0; i <= s->p.mpsi; i++) vdrtmp[i] = s->pfluxpsi[i];
  }
  if (profile) {
    const double now = MPI_Wtime();
    pushi_profile.prep += now - profile_tick;
    profile_tick = now;
  }
  const int idiag = ((s->irk + 1) % 2) + (s->istep % s->p.ndiag);

  GtcReal *restrict wpi0 = s->wpi;
  GtcReal *restrict wpi1 = s->wpi + (size_t)s->mimax;
  GtcReal *restrict wpi2 = s->wpi + (size_t)2 * (size_t)s->mimax;
  const GtcReal *restrict evector = s->evector;
  const int *restrict kzion = s->kzion;
  const int *restrict jtion0 = s->jtion0;
  const int *restrict jtion1 = s->jtion1;
  const GtcReal *restrict wzion = s->wzion;
  const GtcReal *restrict wpion = s->wpion;
  const GtcReal *restrict wtion0 = s->wtion0;
  const GtcReal *restrict wtion1 = s->wtion1;

  if (linear_orbit && idiag != 0) {
#ifdef GTC_USE_METAL
    if (gtc_gpu_pushi_linear_orbit(s, temp_inv, delr, pi2, psimax, cmratio,
                                   cinv, vthi, ainv, sbound, dtime)) {
      if (profile) {
        const double now = MPI_Wtime();
        pushi_profile.advance += now - profile_tick;
        profile_tick = now;
      }
      if (s->irk == 2) {
GTC_OMP_PARALLEL_FOR_STATIC
        for (int m = 0; m < s->mi; m++) {
          GtcReal *zion = &s->zion[(size_t)m * (size_t)GTC_NPARAM];
          const GtcReal *zion0 = &s->zion0[(size_t)m * (size_t)GTC_NPARAM];
          if (zion[GTC_Z_PSI] > psimax || zion[GTC_Z_PSI] < psimin) {
            zion[GTC_Z_PSI] = zion0[GTC_Z_PSI];
            zion[GTC_Z_THETA] = gtc_real(2.0 * s->pi - zion0[GTC_Z_THETA]);
            zion[GTC_Z_ZETA] = zion0[GTC_Z_ZETA];
            zion[GTC_Z_U] = zion0[GTC_Z_U];
            zion[GTC_Z_WEIGHT] = zion0[GTC_Z_WEIGHT];
          }
        }
      }
      if (profile) {
        const double now = MPI_Wtime();
        pushi_profile.boundary += now - profile_tick;
        pushi_profile.total += now - profile_start;
        pushi_profile.calls++;
        pushi_profile_report(s);
      }
      free(vdrtmp);
      free(temp_inv);
      return;
    }
#endif
GTC_OMP_PARALLEL_FOR_STATIC
    for (int m = 0; m < s->mi; m++) {
      GtcReal *zion = &s->zion[(size_t)m * (size_t)GTC_NPARAM];
      GtcReal *zion0 = &s->zion0[(size_t)m * (size_t)GTC_NPARAM];
      if (s->irk == 1) {
        zion0[GTC_Z_PSI] = gtc_real(zion[GTC_Z_PSI]);
        zion0[GTC_Z_THETA] = gtc_real(zion[GTC_Z_THETA]);
        zion0[GTC_Z_ZETA] = gtc_real(zion[GTC_Z_ZETA]);
        zion0[GTC_Z_U] = gtc_real(zion[GTC_Z_U]);
        zion0[GTC_Z_WEIGHT] = gtc_real(zion[GTC_Z_WEIGHT]);
      }

      GtcReal epsi = 0.0, etheta = 0.0, ezeta = 0.0;
      const int kk = kzion[m];
      const GtcReal wz1 = wzion[m];
      const GtcReal wz0 = 1.0 - wz1;
      for (int larmor = 0; larmor < 4; larmor++) {
        const size_t o = gtc_larmor_particle_index(s, larmor, m);
        int ij = jtion0[o];
        const GtcReal wp0 = 1.0 - wpion[o];
        const GtcReal wt00 = 1.0 - wtion0[o];
        const GtcReal wt10 = 1.0 - wt00;
        GtcReal coeff = wp0 * wt00;
        GtcReal c0 = coeff * wz0;
        GtcReal c1 = coeff * wz1;
        size_t base0 = gtc_evector_index(s, 0, kk, ij);
        size_t base1 = base0 + 3u;
        epsi += c0 * evector[base0] + c1 * evector[base1];
        etheta += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];
        ezeta += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];
        ij++;
        coeff = wp0 * wt10;
        c0 = coeff * wz0;
        c1 = coeff * wz1;
        base0 = gtc_evector_index(s, 0, kk, ij);
        base1 = base0 + 3u;
        epsi += c0 * evector[base0] + c1 * evector[base1];
        etheta += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];
        ezeta += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];

        ij = jtion1[o];
        const GtcReal wp1 = 1.0 - wp0;
        const GtcReal wt01 = 1.0 - wtion1[o];
        const GtcReal wt11 = 1.0 - wt01;
        coeff = wp1 * wt01;
        c0 = coeff * wz0;
        c1 = coeff * wz1;
        base0 = gtc_evector_index(s, 0, kk, ij);
        base1 = base0 + 3u;
        epsi += c0 * evector[base0] + c1 * evector[base1];
        etheta += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];
        ezeta += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];
        ij++;
        coeff = wp1 * wt11;
        c0 = coeff * wz0;
        c1 = coeff * wz1;
        base0 = gtc_evector_index(s, 0, kk, ij);
        base1 = base0 + 3u;
        epsi += c0 * evector[base0] + c1 * evector[base1];
        etheta += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];
        ezeta += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];
      }
      epsi *= 0.25;
      etheta *= 0.25;
      ezeta *= 0.25;

      const GtcReal psi = zion[GTC_Z_PSI];
      const GtcReal theta = zion[GTC_Z_THETA];
      const GtcReal u = zion[GTC_Z_U];
      const GtcReal mu = zion[GTC_Z_MU];
      const GtcReal r = sqrtf(2.0f * psi);
      const GtcReal rinv = 1.0 / r;
      const int ii = gtc_clamp_int((int)floorf((r - s->p.a0) * delr), 0, s->p.mpsi - 1);
      const GtcReal wp0 = (GtcReal)(ii + 1) - (r - s->p.a0) * delr;
      const GtcReal wp1 = 1.0 - wp0;
      const GtcReal tem = wp0 * temp_inv[ii] + wp1 * temp_inv[ii + 1];
      const GtcReal q = s->p.q0 + s->p.q1 * r * ainv + s->p.q2 * r * r * ainv * ainv;
      const GtcReal qinv = 1.0 / q;
      GtcReal sint = 0.0;
      GtcReal cost = 1.0;
      sincos_real(theta, &sint, &cost);
      const GtcReal b = 1.0 / (1.0 + r * cost);
      const GtcReal dbdp = -b * b * cost * rinv;
      const GtcReal dbdt = b * b * r * sint;
      const GtcReal dedb = cinv * (u * u * s->p.qion * b * cmratio + mu * mu);
      const GtcReal upara = u * b * cmratio;
      const GtcReal energy = 0.5 * s->p.aion * upara * upara + mu * mu * b;
      GtcReal rfac = s->p.rw * (r - s->p.rc);
      rfac = rfac * rfac;
      rfac = rfac * rfac * rfac;
      rfac = expf(-rfac);
      const GtcReal kappa = ((energy * tem - 1.5) * s->p.kappati + s->p.kappan) *
                           (1.0 - sbound + sbound * rfac) * rinv;

      const GtcReal vdr = -etheta;
      const GtcReal epara = -ezeta * b;
      const GtcReal wdrive = vdr * kappa;
      const GtcReal wpara = epara * upara * s->p.qion * tem;
      const GtcReal wdrift = (dbdt * epsi - dbdp * etheta) * dedb * s->p.qion * tem;
      const GtcReal wdot = zion0[GTC_Z_MU] * (wdrive + wpara + wdrift);
      const GtcReal pdot = -dedb * dbdt;
      const GtcReal tdot = upara * b * qinv + dedb * dbdp;
      const GtcReal zdot = upara * b;
      const GtcReal rdot = -dedb * dbdt * qinv;

      zion[GTC_Z_PSI] = gtc_real(fmax(1.0e-8 * psimax, zion0[GTC_Z_PSI] + dtime * pdot));
      zion[GTC_Z_THETA] = gtc_real(modulo_positive(zion0[GTC_Z_THETA] + dtime * tdot, pi2));
      zion[GTC_Z_ZETA] = gtc_real(modulo_positive(zion0[GTC_Z_ZETA] + dtime * zdot, pi2));
      zion[GTC_Z_U] = gtc_real(zion0[GTC_Z_U] + dtime * rdot);
      zion[GTC_Z_WEIGHT] = gtc_real(zion0[GTC_Z_WEIGHT] + dtime * wdot);
    }
    if (profile) {
      const double now = MPI_Wtime();
      pushi_profile.advance += now - profile_tick;
      profile_tick = now;
    }
    if (s->irk == 2) {
GTC_OMP_PARALLEL_FOR_STATIC
      for (int m = 0; m < s->mi; m++) {
        GtcReal *zion = &s->zion[(size_t)m * (size_t)GTC_NPARAM];
        const GtcReal *zion0 = &s->zion0[(size_t)m * (size_t)GTC_NPARAM];
        if (zion[GTC_Z_PSI] > psimax || zion[GTC_Z_PSI] < psimin) {
          zion[GTC_Z_PSI] = zion0[GTC_Z_PSI];
          zion[GTC_Z_THETA] = gtc_real(2.0 * s->pi - zion0[GTC_Z_THETA]);
          zion[GTC_Z_ZETA] = zion0[GTC_Z_ZETA];
          zion[GTC_Z_U] = zion0[GTC_Z_U];
          zion[GTC_Z_WEIGHT] = zion0[GTC_Z_WEIGHT];
        }
      }
    }
    if (profile) {
      const double now = MPI_Wtime();
      pushi_profile.boundary += now - profile_tick;
      pushi_profile.total += now - profile_start;
      pushi_profile.calls++;
      pushi_profile_report(s);
    }
    free(vdrtmp);
    free(temp_inv);
    return;
  }

#ifdef GTC_USE_METAL
  if (!gtc_gpu_pushi_general(s, temp_inv, vdrtmp, vdrtmp != NULL, linear_orbit,
                             delr, pi2, psimax, cmratio, cinv, vthi, ainv,
                             sbound, dtime)) {
#endif
GTC_OMP_PARALLEL
  {
GTC_OMP_FOR_STATIC
  for (int m = 0; m < s->mi; m++) {
    GtcReal e1 = 0.0, e2 = 0.0, e3 = 0.0;
    const int kk = kzion[m];
    const GtcReal wz1 = wzion[m];
    const GtcReal wz0 = 1.0 - wz1;
    for (int larmor = 0; larmor < 4; larmor++) {
      const size_t o = gtc_larmor_particle_index(s, larmor, m);
      int ij = jtion0[o];
      const GtcReal wp0 = 1.0 - wpion[o];
      const GtcReal wt00 = 1.0 - wtion0[o];
      const GtcReal wt10 = 1.0 - wt00;
      GtcReal coeff = wp0 * wt00;
      GtcReal c0 = coeff * wz0;
      GtcReal c1 = coeff * wz1;
      size_t base0 = gtc_evector_index(s, 0, kk, ij);
      size_t base1 = base0 + 3u;
      e1 += c0 * evector[base0] + c1 * evector[base1];
      e2 += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];
      e3 += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];
      ij++;
      coeff = wp0 * wt10;
      c0 = coeff * wz0;
      c1 = coeff * wz1;
      base0 = gtc_evector_index(s, 0, kk, ij);
      base1 = base0 + 3u;
      e1 += c0 * evector[base0] + c1 * evector[base1];
      e2 += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];
      e3 += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];

      ij = jtion1[o];
      const GtcReal wp1 = 1.0 - wp0;
      const GtcReal wt01 = 1.0 - wtion1[o];
      const GtcReal wt11 = 1.0 - wt01;
      coeff = wp1 * wt01;
      c0 = coeff * wz0;
      c1 = coeff * wz1;
      base0 = gtc_evector_index(s, 0, kk, ij);
      base1 = base0 + 3u;
      e1 += c0 * evector[base0] + c1 * evector[base1];
      e2 += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];
      e3 += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];
      ij++;
      coeff = wp1 * wt11;
      c0 = coeff * wz0;
      c1 = coeff * wz1;
      base0 = gtc_evector_index(s, 0, kk, ij);
      base1 = base0 + 3u;
      e1 += c0 * evector[base0] + c1 * evector[base1];
      e2 += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];
      e3 += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];
    }
    wpi0[m] = gtc_real(0.25 * e1);
    wpi1[m] = gtc_real(0.25 * e2);
    wpi2[m] = gtc_real(0.25 * e3);
  }
    if (profile) {
GTC_OMP_MASTER
      {
        const double now = MPI_Wtime();
        pushi_profile.interpolate += now - profile_tick;
        profile_tick = now;
      }
GTC_OMP_BARRIER
    }

GTC_OMP_FOR_STATIC
  for (int m = 0; m < s->mi; m++) {
    GtcReal *zion = &s->zion[(size_t)m * (size_t)GTC_NPARAM];
    GtcReal *zion0 = &s->zion0[(size_t)m * (size_t)GTC_NPARAM];
    if (s->irk == 1) {
      zion0[GTC_Z_PSI] = gtc_real(zion[GTC_Z_PSI]);
      zion0[GTC_Z_THETA] = gtc_real(zion[GTC_Z_THETA]);
      zion0[GTC_Z_ZETA] = gtc_real(zion[GTC_Z_ZETA]);
      zion0[GTC_Z_U] = gtc_real(zion[GTC_Z_U]);
      zion0[GTC_Z_WEIGHT] = gtc_real(zion[GTC_Z_WEIGHT]);
    }
    const GtcReal psi = zion[GTC_Z_PSI];
    const GtcReal theta = zion[GTC_Z_THETA];
    const GtcReal u = zion[GTC_Z_U];
    const GtcReal weight = zion[GTC_Z_WEIGHT];
    const GtcReal mu = zion[GTC_Z_MU];
    const GtcReal r = sqrtf(2.0f * psi);
    const GtcReal rinv = 1.0 / r;
    const int ii = gtc_clamp_int((int)floorf((r - s->p.a0) * delr), 0, s->p.mpsi - 1);
    const GtcReal wp0 = (GtcReal)(ii + 1) - (r - s->p.a0) * delr;
    const GtcReal wp1 = 1.0 - wp0;
    const GtcReal tem = wp0 * temp_inv[ii] + wp1 * temp_inv[ii + 1];
    const GtcReal q = s->p.q0 + s->p.q1 * r * ainv + s->p.q2 * r * r * ainv * ainv;
    const GtcReal qinv = 1.0 / q;
    GtcReal sint = 0.0;
    GtcReal cost = 1.0;
    sincos_real(theta, &sint, &cost);
    const GtcReal b = 1.0 / (1.0 + r * cost);
    const GtcReal g = 1.0;
    const GtcReal gp = 0.0;
    const GtcReal ri = 0.0;
    const GtcReal rip = 0.0;
    const GtcReal dbdp = -b * b * cost * rinv;
    const GtcReal dbdt = b * b * r * sint;
    const GtcReal dedb = cinv * (u * u * s->p.qion * b * cmratio + mu * mu);
    const GtcReal deni = 1.0 / (g * q + ri + u * (g * rip - ri * gp));
    const GtcReal upara = u * b * cmratio;
    const GtcReal energy = 0.5 * s->p.aion * upara * upara + mu * mu * b;
    GtcReal rfac = s->p.rw * (r - s->p.rc);
    rfac = rfac * rfac;
    rfac = rfac * rfac * rfac;
    rfac = expf(-rfac);
    const GtcReal kappa = ((energy * tem - 1.5) * s->p.kappati + s->p.kappan) *
                         (1.0 - sbound + sbound * rfac) * rinv;

    const GtcReal epsi = wpi0[m];
    const GtcReal etheta = wpi1[m];
    const GtcReal ezeta = wpi2[m];
    GtcReal dptdp = epsi;
    GtcReal dptdt = etheta;
    GtcReal dptdz = ezeta - etheta * qinv;
    const GtcReal epara = -ezeta * b * q * deni;
    const GtcReal vdr = q * (ri * dptdz - g * dptdt) * deni;
    const GtcReal wdrive = vdr * kappa;
    const GtcReal wpara = epara * upara * s->p.qion * tem;
    const GtcReal wdrift = q * (g * dbdt * dptdp - g * dbdp * dptdt + ri * dbdp * dptdz) *
                          deni * dedb * s->p.qion * tem;
    const GtcReal marker_weight = zion0[GTC_Z_MU];
    const GtcReal wdot = (marker_weight - s->p.paranl * weight) * (wdrive + wpara + wdrift);

    const GtcReal profile_vdr = vdrtmp ? vdrtmp[ii] : 0.0;
    GtcReal pdot, tdot, zdot, rdot;
    if (linear_orbit) {
      pdot = q * (-dedb * dbdt) * deni;
      tdot = (upara * b + q * (dedb * dbdp)) * deni;
      zdot = upara * b * q * deni;
      rdot = -dedb * dbdt * deni;
    } else {
      dptdp = dptdp * s->p.nonlinear + s->gyroradius * s->p.flow0 * (1.0 - 0.5 * s->p.a * rinv);
      dptdt *= s->p.nonlinear;
      dptdz *= s->p.nonlinear;
      pdot = q * (-g * dedb * dbdt - g * dptdt + ri * dptdz) * deni - profile_vdr;
      tdot = (upara * b * (1.0 - q * gp * u) +
              q * g * (dedb * dbdp + dptdp)) *
             deni;
      zdot = (upara * b * q * (1.0 + rip * u) -
              q * ri * (dedb * dbdp + dptdp)) *
             deni;
      rdot = ((gp * u - 1.0) * (dedb * dbdt + s->p.paranl * dptdt) -
              s->p.paranl * q * (1.0 + rip * u) * dptdz) *
             deni;
    }

    zion[GTC_Z_PSI] = gtc_real(fmax(1.0e-8 * psimax, zion0[GTC_Z_PSI] + dtime * pdot));
    zion[GTC_Z_THETA] = gtc_real(modulo_positive(zion0[GTC_Z_THETA] + dtime * tdot, pi2));
    zion[GTC_Z_ZETA] = gtc_real(modulo_positive(zion0[GTC_Z_ZETA] + dtime * zdot, pi2));
    zion[GTC_Z_U] = gtc_real(zion0[GTC_Z_U] + dtime * rdot);
    zion[GTC_Z_WEIGHT] = gtc_real(zion0[GTC_Z_WEIGHT] + dtime * wdot);

    wpi0[m] = gtc_real(vdr);
    wpi1[m] = gtc_real(energy);
    wpi2[m] = gtc_real(b);
  }
  }
#ifdef GTC_USE_METAL
  }
#endif
  if (profile) {
    const double now = MPI_Wtime();
    pushi_profile.advance += now - profile_tick;
    profile_tick = now;
  }

  if (s->irk == 2) {
GTC_OMP_PARALLEL_FOR_STATIC
    for (int m = 0; m < s->mi; m++) {
      GtcReal *zion = &s->zion[(size_t)m * (size_t)GTC_NPARAM];
      const GtcReal *zion0 = &s->zion0[(size_t)m * (size_t)GTC_NPARAM];
      if (zion[GTC_Z_PSI] > psimax || zion[GTC_Z_PSI] < psimin) {
        zion[GTC_Z_PSI] = zion0[GTC_Z_PSI];
        zion[GTC_Z_THETA] = gtc_real(2.0 * s->pi - zion0[GTC_Z_THETA]);
        zion[GTC_Z_ZETA] = zion0[GTC_Z_ZETA];
        zion[GTC_Z_U] = zion0[GTC_Z_U];
        zion[GTC_Z_WEIGHT] = zion0[GTC_Z_WEIGHT];
      }
    }
    if (profile) {
      const double now = MPI_Wtime();
      pushi_profile.boundary += now - profile_tick;
      profile_tick = now;
    }

    if (s->p.nonlinear > 0.5 && s->p.paranl < 0.5) {
      GtcReal *dtem = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*dtem), "pushi dtem");
      GtcReal *dmark = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*dmark), "pushi dmark");
      GtcReal *dden = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*dden), "pushi dden");
      GtcReal *dtemtmp = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*dtemtmp), "pushi dtemtmp");
      GtcReal *dmarktmp = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*dmarktmp), "pushi dmarktmp");
      GtcReal *ddentmp = gtc_xcalloc(s, (size_t)(s->p.mpsi + 1), sizeof(*ddentmp), "pushi ddentmp");
      int *iir = gtc_xcalloc(s, (size_t)(s->mi > 0 ? s->mi : 1), sizeof(int), "pushi iir");
      const GtcReal tem_inv = 1.0 / (s->p.aion * vthi * vthi);
      for (int m = 0; m < s->mi; m++) {
        const GtcReal r = sqrtf(2.0f * s->zion[gtc_particle_index(s, GTC_Z_PSI, m)]);
        const int ii = gtc_clamp_int((int)floorf((r - s->p.a0) * delr), 0, s->p.mpsi);
        iir[m] = ii;
        dtem[ii] = gtc_real(dtem[ii] + s->wpi[1 * s->mimax + m] * s->zion[gtc_particle_index(s, GTC_Z_WEIGHT, m)]);
        dmark[ii] = gtc_real(dmark[ii] + s->wpi[0 * s->mimax + m]);
        dden[ii] = gtc_real(dden[ii] + 1.0);
      }
      MPI_Allreduce(dtem, dtemtmp, s->p.mpsi + 1, GTC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(dmark, dmarktmp, s->p.mpsi + 1, GTC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(dden, ddentmp, s->p.mpsi + 1, GTC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
      for (int i = 0; i <= s->p.mpsi; i++) dmark[i] = gtc_real(dmarktmp[i] / fmax(1.0, ddentmp[i]));
      int irsmooth = (int)sqrtf((GtcReal)s->p.mpsi);
      for (int iter = 0; iter < irsmooth; iter++) {
        const GtcReal rdum = dmark[1];
        const GtcReal tdum = dmark[s->p.mpsi - 1];
        for (int i = 1; i < s->p.mpsi; i++) dmarktmp[i] = gtc_real(0.5 * dmark[i] + 0.25 * (dmark[i - 1] + dmark[i + 1]));
        for (int i = 1; i < s->p.mpsi; i++) dmark[i] = dmarktmp[i];
        dmark[0] = gtc_real(0.5 * (dmark[0] + rdum));
        dmark[s->p.mpsi] = gtc_real(0.5 * (dmark[s->p.mpsi] + tdum));
      }
      for (int i = 0; i <= s->p.mpsi; i++) s->pfluxpsi[i] = 0.9 * s->pfluxpsi[i] + 0.1 * dmark[i];

      for (int i = 0; i <= s->p.mpsi; i++) dtem[i] = gtc_real(dtemtmp[i] * tem_inv / fmax(1.0, ddentmp[i]));
      irsmooth = s->p.mpsi;
      for (int iter = 0; iter < irsmooth; iter++) {
        const GtcReal rdum = dtem[1];
        const GtcReal tdum = dtem[s->p.mpsi - 1];
        for (int i = 1; i < s->p.mpsi; i++) dtemtmp[i] = gtc_real(0.5 * dtem[i] + 0.25 * (dtem[i - 1] + dtem[i + 1]));
        for (int i = 1; i < s->p.mpsi; i++) dtem[i] = dtemtmp[i];
        dtem[0] = gtc_real(0.5 * (dtem[0] + rdum));
        dtem[s->p.mpsi] = gtc_real(0.5 * (dtem[s->p.mpsi] + tdum));
      }
      for (int i = 0; i <= s->p.mpsi; i++) s->rdtemi[i] = 0.99 * s->rdtemi[i] + 0.01 * dtem[i];
GTC_OMP_PARALLEL_FOR_STATIC
      for (int m = 0; m < s->mi; m++) {
        const int ii = iir[m];
        s->zion[gtc_particle_index(s, GTC_Z_WEIGHT, m)] = gtc_real(s->zion[gtc_particle_index(s, GTC_Z_WEIGHT, m)] - (s->wpi[1 * s->mimax + m] * tem_inv - 1.5) * s->rdtemi[ii]);
      }
      free(dtem);
      free(dmark);
      free(dden);
      free(dtemtmp);
      free(dmarktmp);
      free(ddentmp);
      free(iir);
    }
  }
  if (profile && s->irk != 2) {
    const double now = MPI_Wtime();
    pushi_profile.boundary += now - profile_tick;
    profile_tick = now;
  }

  if (idiag == 0) {
    const int radial_count = s->p.mpsi + 1;
    const int nthreads = gtc_omp_max_threads();
    GtcReal *dden = gtc_xcalloc(s, (size_t)radial_count, sizeof(*dden), "pushi diag dden");
    GtcReal *hflux_private = gtc_xcalloc(s, (size_t)nthreads * (size_t)radial_count,
                                         sizeof(*hflux_private), "pushi diag hflux private");
    GtcReal *dden_private = gtc_xcalloc(s, (size_t)nthreads * (size_t)radial_count,
                                        sizeof(*dden_private), "pushi diag dden private");
    GtcReal *diag_scalars = gtc_xcalloc(s, 4u * (size_t)nthreads, sizeof(*diag_scalars),
                                        "pushi diag scalar private");
    GtcReal *hflux_tmp = gtc_xcalloc(s, (size_t)radial_count, sizeof(*hflux_tmp), "pushi diag hflux tmp");
    GtcReal *dden_tmp = gtc_xcalloc(s, (size_t)radial_count, sizeof(*dden_tmp), "pushi diag dden tmp");
    s->efluxi = 0.0;
    s->pfluxi = 0.0;
    s->dflowi = 0.0;
    s->entropyi = 0.0;
GTC_OMP_PARALLEL
    {
      const int tid = gtc_omp_thread_num();
      GtcReal *hflux_local = hflux_private + (size_t)tid * (size_t)radial_count;
      GtcReal *dden_local = dden_private + (size_t)tid * (size_t)radial_count;
      GtcReal eflux_local = 0.0;
      GtcReal pflux_local = 0.0;
      GtcReal dflow_local = 0.0;
      GtcReal entropy_local = 0.0;
GTC_OMP_FOR_STATIC
      for (int m = 0; m < s->mi; m++) {
        const GtcReal r = sqrtf(2.0f * s->zion0[gtc_particle_index(s, GTC_Z_PSI, m)]);
        const GtcReal rinv = 1.0 / r;
        const int ii = gtc_clamp_int((int)floorf((r - s->p.a0) * delr + 0.5f), 0, s->p.mpsi);
        const GtcReal weight0 = s->zion0[gtc_particle_index(s, GTC_Z_WEIGHT, m)];
        const GtcReal vdrenergy =
            rinv * s->wpi[0 * s->mimax + m] *
            (s->wpi[1 * s->mimax + m] - 1.5 * s->p.aion * vthi * vthi * s->rtemi[ii]) *
            weight0;
        hflux_local[ii] += vdrenergy;
        dden_local[ii] = gtc_real(dden_local[ii] + 1.0);
        eflux_local += vdrenergy;
        pflux_local += rinv * s->wpi[0 * s->mimax + m] * weight0;
        dflow_local += s->wpi[2 * s->mimax + m] *
                       s->zion0[gtc_particle_index(s, GTC_Z_U, m)] * weight0;
        entropy_local += weight0 * weight0;
      }
      diag_scalars[4u * (size_t)tid + 0u] = eflux_local;
      diag_scalars[4u * (size_t)tid + 1u] = pflux_local;
      diag_scalars[4u * (size_t)tid + 2u] = dflow_local;
      diag_scalars[4u * (size_t)tid + 3u] = entropy_local;
    }
    for (int i = 0; i < radial_count; i++) {
      GtcReal hsum = 0.0;
      GtcReal dsum = 0.0;
      for (int tid = 0; tid < nthreads; tid++) {
        const size_t idx = (size_t)tid * (size_t)radial_count + (size_t)i;
        hsum += hflux_private[idx];
        dsum += dden_private[idx];
      }
      s->hfluxpsi[i] = hsum;
      dden[i] = dsum;
    }
    for (int tid = 0; tid < nthreads; tid++) {
      s->efluxi += diag_scalars[4u * (size_t)tid + 0u];
      s->pfluxi += diag_scalars[4u * (size_t)tid + 1u];
      s->dflowi += diag_scalars[4u * (size_t)tid + 2u];
      s->entropyi += diag_scalars[4u * (size_t)tid + 3u];
    }
    MPI_Allreduce(s->hfluxpsi, hflux_tmp, s->p.mpsi + 1, GTC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(dden, dden_tmp, s->p.mpsi + 1, GTC_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
    for (int i = 0; i <= s->p.mpsi; i++) s->hfluxpsi[i] = hflux_tmp[i];
    for (int i = 1; i <= GTC_MFLUX; i++) {
      const int width = s->p.mpsi / GTC_MFLUX;
      const int lo = (i - 1) * width + 1;
      const int hi = i * width;
      GtcReal flux_sum = 0.0;
      GtcReal marker_sum = 0.0;
      for (int j = lo; j <= hi && j <= s->p.mpsi; j++) {
        flux_sum += s->hfluxpsi[j];
        marker_sum += dden_tmp[j];
      }
      s->eflux[i - 1] = flux_sum / marker_sum;
    }
    for (int i = 0; i <= s->p.mpsi; i++) s->hfluxpsi[i] = hflux_tmp[i] / dden_tmp[i];
    if (s->rank == 0) {
      FILE *fp = fopen("ddeni.dat", "w");
      if (fp) {
        for (int i = 1; i <= s->p.mpsi; i++) {
          gtc_fprintf_e(fp, 6, dden_tmp[i]);
          fputc(i == s->p.mpsi ? '\n' : ' ', fp);
        }
        fclose(fp);
      }
    }
    free(dden);
    free(hflux_private);
    free(dden_private);
    free(diag_scalars);
    free(hflux_tmp);
    free(dden_tmp);

    GtcReal en0 = 0.0;
    GtcReal en1 = 0.0;
    GtcReal en2 = 0.0;
GTC_OMP_PARALLEL_FOR_STATIC_REDUCTION(+:en0,en1,en2)
    for (int i = 1; i < s->p.mpsi; i++) {
      const GtcReal r = s->p.a0 + s->deltar * (GtcReal)i;
      const GtcReal drdp = 1.0 / r;
      for (int k = 1; k <= s->mzeta; k++) {
        const GtcReal zdum = s->zetamin + (GtcReal)k * s->deltaz;
        for (int j = 1; j <= s->mtheta[i]; j++) {
          const int ij = s->igrid[i] + j;
          const GtcReal tdum = (GtcReal)j * s->deltat[i] + zdum * s->qtinv[i];
          const GtcReal grad_zeta = 1.0f + r * cosf(tdum);
          const GtcReal jacobian = grad_zeta * grad_zeta;
          const GtcReal epsi = s->evector[gtc_evector_index(s, 0, k, ij)] - s->phip00[i] * drdp;
          const GtcReal etheta = s->evector[gtc_evector_index(s, 1, k, ij)];
          const GtcReal ezeta = s->evector[gtc_evector_index(s, 2, k, ij)] -
                               s->evector[gtc_evector_index(s, 1, k, ij)] * s->qtinv[i];
          en0 += (r * r * epsi * epsi + drdp * drdp * etheta * etheta +
                    ezeta * ezeta / (grad_zeta * grad_zeta)) * jacobian * r * s->deltat[i];
          en1 += (r * r * s->evector[gtc_evector_index(s, 0, k, ij)] *
                        s->evector[gtc_evector_index(s, 0, k, ij)] +
                    drdp * drdp * etheta * etheta + ezeta * ezeta / (grad_zeta * grad_zeta)) *
                   jacobian * r * s->deltat[i];
          en2 += s->phip00[i] * s->phip00[i] * jacobian * r * s->deltat[i];
        }
      }
    }
    const GtcReal scale = s->deltar * s->deltaz / (4.0 * s->pi);
    GtcReal en_field[3] = {gtc_real(en0 * scale), gtc_real(en1 * scale), gtc_real(en2 * scale)};
    GtcReal total_field_energy[3] = {0.0, 0.0, 0.0};
    MPI_Reduce(en_field, total_field_energy, 3, GTC_MPI_REAL, MPI_SUM, 0, s->toroidal_comm);
    for (int i = 0; i < 3; i++) s->total_field_energy[i] = total_field_energy[i];
  }
  if (profile) {
    const double now = MPI_Wtime();
    pushi_profile.diagnostics += now - profile_tick;
    pushi_profile.total += now - profile_start;
    pushi_profile.calls++;
    pushi_profile_report(s);
  }
  free(vdrtmp);
  free(temp_inv);
}
