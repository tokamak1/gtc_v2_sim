#include "gtc.h"
#include "pocketfft_hdronly.h"

#include <complex>
#include <exception>

namespace {
using complexd = std::complex<GtcReal>;

void run_c2c(int isign, int n, GtcReal scale, complexd *data) {
  const pocketfft::shape_t shape{static_cast<size_t>(n)};
  const pocketfft::stride_t stride{static_cast<ptrdiff_t>(sizeof(complexd))};
  const pocketfft::shape_t axes{0};
  const bool forward = isign < 0;
  pocketfft::c2c(shape, stride, stride, axes, forward, data, data, scale, 1);
}
}

extern "C" int gtc_pocketfft_c2c(int isign, int n, GtcReal scale, GtcReal *x) {
  if (n <= 0 || x == nullptr) return 0;
  try {
    const GtcReal normalized_scale = isign == -1 ? scale / static_cast<GtcReal>(n) : scale;
    run_c2c(isign, n, normalized_scale, reinterpret_cast<complexd *>(x));
    return 1;
  } catch (const std::exception &) {
    return 0;
  }
}
