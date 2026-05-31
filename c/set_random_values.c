#include "gtc.h"

static GtcReal sf(GtcReal x) {
  return gtc_real(x);
}

void set_random_zion(GtcState *s) {
  if (s->p.rng_control > 0) {
    GtcRng rng = s->rng;
    for (int rep = 0; rep <= s->rank; rep++) {
      for (int m = 0; m < s->mi; m++) {
        s->zion[gtc_particle_index(s, GTC_Z_THETA, m)] = sf(gtc_rng_uniform(&rng));
        s->zion[gtc_particle_index(s, GTC_Z_ZETA, m)] = sf(gtc_rng_uniform(&rng));
        s->zion[gtc_particle_index(s, GTC_Z_U, m)] = sf(gtc_rng_uniform(&rng));
        s->zion[gtc_particle_index(s, GTC_Z_WEIGHT, m)] = sf(gtc_rng_uniform(&rng));
        s->zion[gtc_particle_index(s, GTC_Z_MU, m)] = sf(gtc_rng_uniform(&rng));
      }
    }
    s->rng = rng;
    MPI_Barrier(MPI_COMM_WORLD);
  } else {
    GtcRng rng = s->rng;
    for (int param = GTC_Z_THETA; param <= GTC_Z_MU; param++) {
      for (int m = 0; m < s->mi; m++) {
        s->zion[gtc_particle_index(s, param, m)] = sf(gtc_rng_uniform(&rng));
      }
    }
    s->rng = rng;
  }
}
