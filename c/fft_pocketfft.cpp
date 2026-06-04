#include "gtc.h"

#ifndef POCKETFFT_CACHE_SIZE
#define POCKETFFT_CACHE_SIZE 32
#endif

#include "pocketfft_hdronly.h"

#include <algorithm>
#include <complex>
#include <exception>

namespace {
using complexr = std::complex<GtcReal>;

size_t thread_count(int requested, int count) {
  if (requested <= 1 || count <= 1) return 1;
  return static_cast<size_t>(std::min(requested, count));
}

void run_c2c(int isign, int n, GtcReal scale, complexr *data) {
  const pocketfft::shape_t shape{static_cast<size_t>(n)};
  const pocketfft::stride_t stride{static_cast<ptrdiff_t>(sizeof(complexr))};
  const pocketfft::shape_t axes{0};
  const bool forward = isign < 0;
  pocketfft::c2c(shape, stride, stride, axes, forward, data, data, scale, 1);
}

void run_r2c(int n, int count, GtcReal scale, const GtcReal *in,
             complexr *out, int requested_threads) {
  const int nc = n / 2 + 1;
  const pocketfft::shape_t shape{static_cast<size_t>(count), static_cast<size_t>(n)};
  const pocketfft::stride_t stride_in{
      static_cast<ptrdiff_t>((size_t)n * sizeof(GtcReal)),
      static_cast<ptrdiff_t>(sizeof(GtcReal))};
  const pocketfft::stride_t stride_out{
      static_cast<ptrdiff_t>((size_t)nc * sizeof(complexr)),
      static_cast<ptrdiff_t>(sizeof(complexr))};
  const pocketfft::shape_t axes{1};
  pocketfft::r2c(shape, stride_in, stride_out, axes, false, in, out, scale,
                 thread_count(requested_threads, count));
}

void run_c2r(int n, int count, GtcReal scale, const complexr *in,
             GtcReal *out, int requested_threads) {
  const int nc = n / 2 + 1;
  const pocketfft::shape_t shape{static_cast<size_t>(count), static_cast<size_t>(n)};
  const pocketfft::stride_t stride_in{
      static_cast<ptrdiff_t>((size_t)nc * sizeof(complexr)),
      static_cast<ptrdiff_t>(sizeof(complexr))};
  const pocketfft::stride_t stride_out{
      static_cast<ptrdiff_t>((size_t)n * sizeof(GtcReal)),
      static_cast<ptrdiff_t>(sizeof(GtcReal))};
  const pocketfft::shape_t axes{1};
  pocketfft::c2r(shape, stride_in, stride_out, axes, true, in, out,
                 scale / static_cast<GtcReal>(n), thread_count(requested_threads, count));
}
}

extern "C" int gtc_pocketfft_c2c(int isign, int n, GtcReal scale, GtcReal *x) {
  if (n <= 0 || x == nullptr) return 0;
  try {
    const GtcReal normalized_scale = isign == -1 ? scale / static_cast<GtcReal>(n) : scale;
    run_c2c(isign, n, normalized_scale, reinterpret_cast<complexr *>(x));
    return 1;
  } catch (const std::exception &) {
    return 0;
  }
}

extern "C" int gtc_pocketfft_r2c(int n, int count, GtcReal scale,
                                  const GtcReal *x, GtcReal *y,
                                  int requested_threads) {
  if (n <= 0 || count <= 0 || x == nullptr || y == nullptr) return 0;
  try {
    run_r2c(n, count, scale, x, reinterpret_cast<complexr *>(y), requested_threads);
    return 1;
  } catch (const std::exception &) {
    return 0;
  }
}

extern "C" int gtc_pocketfft_c2r(int n, int count, GtcReal scale,
                                  const GtcReal *y, GtcReal *x,
                                  int requested_threads) {
  if (n <= 0 || count <= 0 || x == nullptr || y == nullptr) return 0;
  try {
    run_c2r(n, count, scale, reinterpret_cast<const complexr *>(y), x, requested_threads);
    return 1;
  } catch (const std::exception &) {
    return 0;
  }
}
