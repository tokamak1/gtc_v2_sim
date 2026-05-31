#include "gtc.h"

#include <math.h>

static GtcReal sf(GtcReal x) {
  return gtc_real(x);
}

void load(GtcState *s) {
  if (s->p.irun != 0) {
    if (s->rank == 0) {
      FILE *fp = fopen("FileExit.dat", "r");
      if (fp) {
        char line[128];
        int file_exit = 0, irest = 1;
        if (fgets(line, sizeof(line), fp) && sscanf(line, "FileExit=%d", &file_exit) == 1) {
          s->file_exit = file_exit;
        }
        if (fgets(line, sizeof(line), fp) && sscanf(line, "irest   =%d", &irest) == 1) {
          s->irest = irest - 1;
        }
        (void)fgets(line, sizeof(line), fp);
        fclose(fp);
      }
    }
    MPI_Bcast(&s->irest, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    restart_io(s, "read");
    return;
  }

  rand_num_gen_init(s);

  const GtcReal c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
  const GtcReal d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
  const GtcReal rmi = 1.0 / (GtcReal)(s->mi * (s->p.npartdom > 0 ? s->p.npartdom : 1));
  const GtcReal pi2_inv = 0.5 / s->pi;
  const GtcReal delr = 1.0 / s->deltar;
  GtcReal w_initial = s->p.nonlinear < 0.5 ? 1.0e-12 : 1.0e-3;
  if (s->rank == 0) {
    FILE *out = gtc_stdout_open(s, "a");
    if (out) {
      fprintf(out, " w_initial = % .8e\n", w_initial);
      gtc_stdout_close(s, out);
    }
  }

  for (int m = 0; m < s->mi; m++) {
    const GtcReal ordinal = (GtcReal)(m * s->p.npartdom + s->myrank_partd + 1) - 0.5;
    s->zion[gtc_particle_index(s, GTC_Z_PSI, m)] = sf(sqrt(sf(s->p.a0 * s->p.a0 + ordinal * (s->p.a1 * s->p.a1 - s->p.a0 * s->p.a0) * rmi)));
  }

  set_random_zion(s);

  for (int m = 0; m < s->mi; m++) {
    s->zion[gtc_particle_index(s, GTC_Z_THETA, m)] = sf(2.0 * s->pi * (s->zion[gtc_particle_index(s, GTC_Z_THETA, m)] - 0.5));
    s->zion0[gtc_particle_index(s, GTC_Z_THETA, m)] = s->zion[gtc_particle_index(s, GTC_Z_THETA, m)];
  }
  for (int iter = 0; iter < 10; iter++) {
    for (int m = 0; m < s->mi; m++) {
      s->zion[gtc_particle_index(s, GTC_Z_THETA, m)] = sf(s->zion0[gtc_particle_index(s, GTC_Z_THETA, m)] - 2.0 * s->zion[gtc_particle_index(s, GTC_Z_PSI, m)] * sin(s->zion[gtc_particle_index(s, GTC_Z_THETA, m)]));
    }
  }
  for (int m = 0; m < s->mi; m++) {
    s->zion[gtc_particle_index(s, GTC_Z_THETA, m)] = sf(s->zion[gtc_particle_index(s, GTC_Z_THETA, m)] * pi2_inv + 10.0);
    s->zion[gtc_particle_index(s, GTC_Z_THETA, m)] = sf(2.0 * s->pi * (s->zion[gtc_particle_index(s, GTC_Z_THETA, m)] - floor(s->zion[gtc_particle_index(s, GTC_Z_THETA, m)])));
  }

  for (int m = 0; m < s->mi; m++) {
    const GtcReal z4tmp = s->zion[gtc_particle_index(s, GTC_Z_U, m)];
    GtcReal u = s->zion[gtc_particle_index(s, GTC_Z_U, m)] - 0.5;
    const GtcReal sign = u >= 0.0 ? 1.0 : -1.0;
    u = sqrt(fmax(1.0e-20, log(1.0 / fmax(1.0e-20, u * u))));
    u = u - (c0 + c1 * u + c2 * u * u) / (1.0 + d1 * u + d2 * u * u + d3 * u * u * u);
    if (u > s->p.umax) u = z4tmp;
    s->zion0[gtc_particle_index(s, GTC_Z_U, m)] = sign;
    s->zion[gtc_particle_index(s, GTC_Z_U, m)] = sf(u);
  }

  const GtcReal vthi = s->gyroradius * fabs(s->p.qion) / s->p.aion;
  for (int m = 0; m < s->mi; m++) {
    s->zion[gtc_particle_index(s, GTC_Z_ZETA, m)] = sf(s->zetamin + (s->zetamax - s->zetamin) * s->zion[gtc_particle_index(s, GTC_Z_ZETA, m)]);
    s->zion[gtc_particle_index(s, GTC_Z_U, m)] = sf(s->zion0[gtc_particle_index(s, GTC_Z_U, m)] * fmin(s->p.umax, s->zion[gtc_particle_index(s, GTC_Z_U, m)]));
    s->zion[gtc_particle_index(s, GTC_Z_WEIGHT, m)] = sf(2.0 * w_initial * (s->zion[gtc_particle_index(s, GTC_Z_WEIGHT, m)] - 0.5) * (1.0 + cos(s->zion[gtc_particle_index(s, GTC_Z_THETA, m)])));
    s->zion[gtc_particle_index(s, GTC_Z_MU, m)] = sf(fmax(1.0e-20, fmin(s->p.umax * s->p.umax, -log(fmax(1.0e-20, s->zion[gtc_particle_index(s, GTC_Z_MU, m)])))));

    const GtcReal b = 1.0 / (1.0 + s->zion[gtc_particle_index(s, GTC_Z_PSI, m)] * cos(s->zion[gtc_particle_index(s, GTC_Z_THETA, m)]));
    const GtcReal r = s->zion[gtc_particle_index(s, GTC_Z_PSI, m)];
    s->zion[gtc_particle_index(s, GTC_Z_PSI, m)] = sf(0.5 * r * r);
    s->zion[gtc_particle_index(s, GTC_Z_U, m)] = sf(vthi * s->zion[gtc_particle_index(s, GTC_Z_U, m)] * s->p.aion / (s->p.qion * b));
    s->zion[gtc_particle_index(s, GTC_Z_MU, m)] = sf(sqrt(s->p.aion * vthi * vthi * s->zion[gtc_particle_index(s, GTC_Z_MU, m)] / b));
    s->zion0[gtc_particle_index(s, GTC_Z_MU, m)] = 1.0;
  }

  if (s->p.iload != 0) {
    for (int m = 0; m < s->mi; m++) {
      const GtcReal r = sqrt(2.0 * s->zion[gtc_particle_index(s, GTC_Z_PSI, m)]);
      const int i = gtc_clamp_int((int)floor((r - s->p.a0) * delr + 0.5), 0, s->p.mpsi);
      s->zion[gtc_particle_index(s, GTC_Z_U, m)] *= sqrt(s->rtemi[i]);
      s->zion[gtc_particle_index(s, GTC_Z_MU, m)] *= sqrt(s->rtemi[i]);
      s->zion0[gtc_particle_index(s, GTC_Z_MU, m)] = fmax(0.1, fmin(10.0, s->rden[i]));
    }
  }
}
