#include "gtc.h"

#include <math.h>

GtcReal psi2r(GtcReal pdum) {
  return sqrt(fmax(GTC_TINY, 2.0 * pdum));
}

GtcReal r2psi(GtcReal rdum) {
  return fmax(GTC_TINY, 0.5 * rdum * rdum);
}

GtcReal bfield(GtcReal pdum, GtcReal tdum) {
  const GtcReal r = psi2r(pdum);
  return 1.0 / (1.0 + r * cos(tdum));
}

GtcReal boozer2x(GtcReal pdum, GtcReal tdum) {
  const GtcReal r = psi2r(pdum);
  return 1.0 + r * cos(tdum);
}

GtcReal boozer2z(GtcReal pdum, GtcReal tdum) {
  const GtcReal r = psi2r(pdum);
  return r * sin(tdum);
}
