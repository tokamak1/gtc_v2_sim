#include "gtc.h"

#ifdef GTC_USE_ACCELERATE
#include <Accelerate/Accelerate.h>
#endif

#include <math.h>
#include <stdlib.h>

#ifdef GTC_USE_POCKETFFT
extern int gtc_pocketfft_c2c(int isign, int n, GtcReal scale, GtcReal *x);
extern int gtc_pocketfft_r2c(int n, int count, GtcReal scale, const GtcReal *x,
                             GtcReal *y, int requested_threads);
extern int gtc_pocketfft_c2r(int n, int count, GtcReal scale, const GtcReal *y,
                             GtcReal *x, int requested_threads);
#endif

static size_t cidx(int index) {
  return (size_t)2 * (size_t)index;
}

static void spcpft_c(int a, int b, int c, const GtcReal *uin, GtcReal *uout, int isign) {
  const GtcReal angle = 8.0 * atan(1.0) / (GtcReal)(a * c);
  GtcReal omega_re = 1.0;
  GtcReal omega_im = 0.0;
  const GtcReal delta_re = cos(angle);
  const GtcReal delta_im = isign == 1 ? sin(angle) : -sin(angle);

  for (int ic = 0; ic < c; ic++) {
    for (int ia = 0; ia < a; ia++) {
      for (int ib = 0; ib < b; ib++) {
        int src = ib + b * ((c - 1) + c * ia);
        GtcReal sum_re = uin[cidx(src)];
        GtcReal sum_im = uin[cidx(src) + 1];
        for (int jcr = 2; jcr <= c; jcr++) {
          const int jc = c + 1 - jcr;
          src = ib + b * ((jc - 1) + c * ia);
          const GtcReal xr = sum_re;
          const GtcReal xi = sum_im;
          sum_re = uin[cidx(src)] + omega_re * xr - omega_im * xi;
          sum_im = uin[cidx(src) + 1] + omega_re * xi + omega_im * xr;
        }
        const int dst = ib + b * (ia + a * ic);
        uout[cidx(dst)] = sum_re;
        uout[cidx(dst) + 1] = sum_im;
      }
      const GtcReal next_re = delta_re * omega_re - delta_im * omega_im;
      const GtcReal next_im = delta_re * omega_im + delta_im * omega_re;
      omega_re = next_re;
      omega_im = next_im;
    }
  }
}

static int spcfft_c(GtcReal *u, int n, int isign, GtcReal *work, GtcReal interp) {
  int a = 1;
  int b = n;
  int c = 1;
  int inu = 1;

  while (b > 1) {
    a = c * a;
    c = 2;
    while (b % c != 0) c++;
    b /= c;
    if (inu) {
      spcpft_c(a, b, c, u, work, isign);
    } else {
      spcpft_c(a, b, c, work, u, isign);
    }
    inu = !inu;
  }

  if (!inu) {
    for (int i = 0; i < 2 * n; i++) u[i] = work[i];
  }
  if (isign == 1) {
    const GtcReal scale = interp / (GtcReal)n;
    for (int i = 0; i < 2 * n; i++) u[i] *= scale;
  }
  return 1;
}

static void fft_complex_transform(int isign, int n, GtcReal scale, const GtcReal *in, GtcReal *out) {
  for (int i = 0; i < 2 * n; i++) out[i] = in[i];
#ifdef GTC_USE_POCKETFFT
  if (n >= 128 && gtc_pocketfft_c2c(isign, n, scale, out)) return;
#endif
  GtcReal *work = calloc((size_t)2 * (size_t)n, sizeof(*work));
  if (!work) return;
  spcfft_c(out, n, -isign, work, scale);
  free(work);
}

static int fft_batch_threads(int count) {
#ifdef GTC_USE_OPENMP
  if (omp_in_parallel()) return 1;
  const int threads = omp_get_max_threads();
  return count < threads ? count : threads;
#else
  (void)count;
  return 1;
#endif
}

#ifdef GTC_USE_ACCELERATE
static int fft_accelerate(int isign, int n, const GtcReal *in, GtcReal *out) {
  static int setup_n = 0;
  static vDSP_DFT_SetupD setup_forward = NULL;
  static vDSP_DFT_SetupD setup_inverse = NULL;
  static GtcReal *ir = NULL;
  static GtcReal *ii = NULL;
  static GtcReal *orv = NULL;
  static GtcReal *oiv = NULL;
  static int buffer_n = 0;

  if (n <= 0) return 0;
  if (setup_n != n) {
    vDSP_DFT_DestroySetupD(setup_forward);
    vDSP_DFT_DestroySetupD(setup_inverse);
    setup_forward = vDSP_DFT_zop_CreateSetupD(NULL, (vDSP_Length)n, vDSP_DFT_INVERSE);
    setup_inverse = vDSP_DFT_zop_CreateSetupD(NULL, (vDSP_Length)n, vDSP_DFT_FORWARD);
    setup_n = n;
  }
  vDSP_DFT_SetupD setup = isign == 1 ? setup_forward : setup_inverse;
  if (!setup) return 0;

  if (buffer_n < n) {
    GtcReal *new_ir = malloc((size_t)n * sizeof(*new_ir));
    GtcReal *new_ii = malloc((size_t)n * sizeof(*new_ii));
    GtcReal *new_orv = malloc((size_t)n * sizeof(*new_orv));
    GtcReal *new_oiv = malloc((size_t)n * sizeof(*new_oiv));
    if (!new_ir || !new_ii || !new_orv || !new_oiv) {
      free(new_ir);
      free(new_ii);
      free(new_orv);
      free(new_oiv);
      return 0;
    }
    free(ir);
    free(ii);
    free(orv);
    free(oiv);
    ir = new_ir;
    ii = new_ii;
    orv = new_orv;
    oiv = new_oiv;
    buffer_n = n;
  }

  for (int i = 0; i < n; i++) {
    ir[i] = in[2 * i];
    ii[i] = in[2 * i + 1];
  }
  vDSP_DFT_ExecuteD(setup, ir, ii, orv, oiv);
  for (int i = 0; i < n; i++) {
    out[2 * i] = orv[i];
    out[2 * i + 1] = oiv[i];
  }
  return 1;
}
#endif

void fftc1d(int isign, int irank, GtcReal scale, GtcReal *x) {
  const int n = irank;
#ifdef GTC_USE_POCKETFFT
  if (n >= 128 && gtc_pocketfft_c2c(isign, n, scale, x)) return;
#endif
  GtcReal *tmp = calloc((size_t)2 * (size_t)n, sizeof(GtcReal));
  if (!tmp) return;
#ifdef GTC_USE_ACCELERATE
  if (!fft_accelerate(isign, n, x, tmp))
    fft_complex_transform(isign, n, scale, x, tmp);
#else
  fft_complex_transform(isign, n, scale, x, tmp);
#endif
  for (int i = 0; i < 2 * n; i++) x[i] = tmp[i];
  free(tmp);
}

void fftr1d(int isign, int irank, GtcReal scale, GtcReal *x, GtcReal *y, int icount) {
  const int n = irank;
  const int nc = n / 2 + 1;
  if (n <= 0 || icount <= 0) return;
#ifdef GTC_USE_POCKETFFT
  if (n >= 128) {
    const int threads = fft_batch_threads(icount);
    if (isign == 1 && gtc_pocketfft_r2c(n, icount, scale, x, y, threads)) return;
    if (isign == -1 && gtc_pocketfft_c2r(n, icount, scale, y, x, threads)) return;
  }
#endif
  GtcReal *tmp = calloc((size_t)2 * (size_t)n, sizeof(GtcReal));
  GtcReal *out = calloc((size_t)2 * (size_t)n, sizeof(GtcReal));
  if (!tmp || !out) {
    free(out);
    free(tmp);
    return;
  }
  for (int c = 0; c < icount; c++) {
    if (isign == 1) {
      for (int j = 0; j < n; j++) {
        tmp[2 * j] = x[c * n + j];
        tmp[2 * j + 1] = 0.0;
      }
      fft_complex_transform(1, n, scale, tmp, out);
      for (int k = 0; k < nc; k++) {
        y[2 * (c * nc + k)] = out[2 * k];
        y[2 * (c * nc + k) + 1] = out[2 * k + 1];
      }
    } else {
      for (int j = 0; j < 2 * n; j++) tmp[j] = 0.0;
      for (int k = 0; k < nc; k++) {
        const GtcReal re = y[2 * (c * nc + k)];
        const GtcReal im = y[2 * (c * nc + k) + 1];
        tmp[2 * k] = re;
        tmp[2 * k + 1] = im;
        if (k > 0 && k < n - k) {
          tmp[2 * (n - k)] = re;
          tmp[2 * (n - k) + 1] = -im;
        }
      }
      fft_complex_transform(-1, n, scale, tmp, out);
      for (int j = 0; j < n; j++) {
        x[c * n + j] = out[2 * j];
      }
    }
  }
  free(out);
  free(tmp);
}
