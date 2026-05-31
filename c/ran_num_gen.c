#include "gtc.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

static void portable_seed_int(int seed[8], int n) {
  const int b = 1 << 14;
  int m, sign;
  for (int i = 0; i < 8; i++) seed[i] = 0;
  if (n < 0) {
    for (int i = 0; i < 8; i++) seed[i] = b - 1;
    m = -n - 1;
    sign = -1;
  } else {
    m = n;
    sign = 1;
  }
  int i = 0;
  while (m != 0 && i < 8) {
    seed[i] += sign * (m % b);
    m /= b;
    i++;
  }
}

static void rand_axc(const int a[8], int x[8], const int c[8]) {
  const int b = 1 << 14;
  int z[8];
  for (int i = 0; i < 8; i++) z[i] = c[i];
  for (int j = 0; j < 8; j++) {
    for (int i = j; i < 8; i++) z[i] += a[j] * x[i - j];
  }
  int t = 0;
  for (int i = 0; i < 8; i++) {
    t = t / b + z[i];
    x[i] = t % b;
  }
}

static void decimal_to_seed(const char *decimal, int seed[8]) {
  const int ten[8] = {10, 0, 0, 0, 0, 0, 0, 0};
  int c[8] = {0};
  for (int i = 0; i < 8; i++) seed[i] = 0;
  for (const char *p = decimal; *p; p++) {
    if (*p >= '0' && *p <= '9') {
      c[0] = *p - '0';
      for (int i = 1; i < 8; i++) c[i] = 0;
      rand_axc(ten, seed, c);
    }
  }
}

static void seed_to_decimal(const int seed[8], char out[35]) {
  const int pow10 = 10000;
  const int b = 1 << 14;
  int z[8], k = -1;
  for (int i = 0; i < 8; i++) {
    z[i] = seed[i];
    if (z[i] > 0) k = i;
  }
  char str[37];
  memset(str, ' ', sizeof(str));
  str[36] = '\0';
  int block = 9;
  do {
    block--;
    int t = 0;
    for (int j = k; j >= 0; j--) {
      z[j] += t * b;
      t = z[j] % pow10;
      z[j] /= pow10;
    }
    if (k >= 0 && z[k] == 0) k--;
    int pos = 4 * (block + 1) - 1;
    for (int q = 0; q < 4; q++) {
      str[pos - q] = (char)('0' + (t % 10));
      t /= 10;
    }
    if (k < 0) {
      while (pos - 3 < 35 && str[pos - 3] == '0' && pos - 3 < 35) {
        str[pos - 3] = ' ';
        break;
      }
    }
  } while (k >= 0 && block > 0);
  memcpy(out, str + 2, 34);
  out[34] = '\0';
}

static void rand_next_seed(int n, const int ax[8], const int cx[8], int y[8]) {
  if (n == 0) return;
  int a[8], c[8], z[8] = {0}, t[8];
  for (int i = 0; i < 8; i++) {
    a[i] = ax[i];
    c[i] = cx[i];
  }
  int m = n;
  while (m != 0) {
    if (m % 2 > 0) rand_axc(a, y, c);
    m /= 2;
    if (m == 0) return;
    for (int i = 0; i < 8; i++) t[i] = c[i];
    rand_axc(a, c, t);
    for (int i = 0; i < 8; i++) t[i] = a[i];
    rand_axc(t, a, z);
  }
}

static void portable_step_seed(int seed[8], int n0) {
  static const int af0[8] = {15741, 8689, 9280, 4732, 12011, 7130, 6824, 12302};
  static const int cf0[8] = {16317, 10266, 1198, 331, 10769, 8310, 2779, 13880};
  static const int ab0[8] = {9173, 9894, 15203, 15379, 7981, 2280, 8071, 429};
  static const int cb0[8] = {8383, 3616, 597, 12724, 15663, 9639, 187, 4866};
  if (n0 > 0) rand_next_seed(n0, af0, cf0, seed);
  else if (n0 < 0) rand_next_seed(-n0, ab0, cb0, seed);
}

static void portable_seed_time(int seed[8]) {
  time_t now = time(NULL);
  struct tm tmv;
  struct tm *tmp = localtime(&now);
  if (tmp) tmv = *tmp;
  else memset(&tmv, 0, sizeof(tmv));
  char decimal[32];
  snprintf(decimal, sizeof(decimal), "%04d%02d%02d0%03d%02d%02d%02d%03d",
           tmv.tm_year + 1900, tmv.tm_mon + 1, tmv.tm_mday, 0,
           tmv.tm_hour, tmv.tm_min, tmv.tm_sec, 0);
  decimal_to_seed(decimal, seed);
}

static void portable_rand_batch(GtcRng *rng) {
  double w[1009 - 100];
  for (int i = 0; i < 63; i++) {
    double tmp = rng->array[i] + rng->array[i + 100 - 63];
    w[i] = tmp - floor(tmp);
  }
  for (int i = 63; i < 100; i++) {
    double tmp = rng->array[i] + w[i - 63];
    w[i] = tmp - floor(tmp);
  }
  for (int i = 100; i < 1009 - 100; i++) {
    double tmp = w[i - 100] + w[i - 63];
    w[i] = tmp - floor(tmp);
  }
  for (int i = 1009 - 100; i < 1009 - 100 + 63; i++) {
    double tmp = w[i - 100] + w[i - 63];
    rng->array[i - 1009 + 100] = tmp - floor(tmp);
  }
  for (int i = 1009 - 100 + 63; i < 1009; i++) {
    double tmp = w[i - 100] + rng->array[i - 1009 + 100 - 63];
    rng->array[i - 1009 + 100] = tmp - floor(tmp);
  }
  rng->index = 0;
}

static void portable_init(GtcRng *rng, const int seed_in[8]) {
  const int b = 1 << 14;
  const double del = pow(2.0, -14);
  const double ulp = pow(2.0, -47);
  const int a0 = 15661, a1 = 678, a2 = 724, a3 = 5245;
  const int a4 = 13656, a5 = 11852, a6 = 29, c0 = 1;
  int s[8], z[8];
  for (int i = 0; i < 8; i++) s[i] = seed_in[i];
  int odd = (s[7] % 2) != 0;
  rng->array[0] = (((s[7] * del + s[6]) * del + s[5]) * del + (s[4] / 512)) * 512 * del;
  for (int j = 1; j < 100; j++) {
    z[0] = c0 + a0 * s[0];
    z[1] = a0 * s[1] + a1 * s[0];
    z[2] = a0 * s[2] + a1 * s[1] + a2 * s[0];
    z[3] = a0 * s[3] + a1 * s[2] + a2 * s[1] + a3 * s[0];
    z[4] = a0 * s[4] + a1 * s[3] + a2 * s[2] + a3 * s[1] + a4 * s[0];
    z[5] = a0 * s[5] + a1 * s[4] + a2 * s[3] + a3 * s[2] + a4 * s[1] + a5 * s[0];
    z[6] = a0 * s[6] + a1 * s[5] + a2 * s[4] + a3 * s[3] + a4 * s[2] + a5 * s[1] + a6 * s[0];
    z[7] = a0 * s[7] + a1 * s[6] + a2 * s[5] + a3 * s[4] + a4 * s[3] + a5 * s[2] + a6 * s[1];
    int t = 0;
    for (int i = 0; i < 8; i++) {
      t = t / b + z[i];
      s[i] = t % b;
    }
    odd = odd || ((s[7] % 2) != 0);
    rng->array[j] = (((s[7] * del + s[6]) * del + s[5]) * del + (s[4] / 512)) * 512 * del;
  }
  rng->index = 100;
  if (!odd) {
    z[0] = c0 + a0 * s[0];
    z[1] = a0 * s[1] + a1 * s[0];
    z[2] = a0 * s[2] + a1 * s[1] + a2 * s[0];
    z[3] = a0 * s[3] + a1 * s[2] + a2 * s[1] + a3 * s[0];
    z[4] = a0 * s[4] + a1 * s[3] + a2 * s[2] + a3 * s[1] + a4 * s[0];
    z[5] = a0 * s[5] + a1 * s[4] + a2 * s[3] + a3 * s[2] + a4 * s[1] + a5 * s[0];
    z[6] = a0 * s[6] + a1 * s[5] + a2 * s[4] + a3 * s[3] + a4 * s[2] + a5 * s[1] + a6 * s[0];
    z[7] = a0 * s[7] + a1 * s[6] + a2 * s[5] + a3 * s[4] + a4 * s[3] + a5 * s[2] + a6 * s[1];
    int t = 0;
    for (int i = 0; i < 8; i++) {
      t = t / b + z[i];
      s[i] = t % b;
    }
    int j = (s[7] * 100) / b;
    rng->array[j] += ulp;
  }
}

void gtc_rng_init(GtcRng *rng, unsigned long long seed) {
  rng->state = seed ? seed : 0x9e3779b97f4a7c15ULL;
  rng->index = 100;
  for (int i = 0; i < 100; i++) rng->array[i] = 0.0;
  int portable_seed[8];
  portable_seed_int(portable_seed, (int)seed);
  portable_init(rng, portable_seed);
}

static unsigned long long gtc_rng_next(GtcRng *rng) {
  unsigned long long x = rng->state;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  rng->state = x;
  return x * 2685821657736338717ULL;
}

GtcReal gtc_rng_uniform(GtcRng *rng) {
  if (rng->index >= 0) {
    const double ulps = pow(2.0, -23);
    const double mult = pow(2.0, 23);
    if (rng->index >= 100) portable_rand_batch(rng);
    GtcReal value = (GtcReal)(((int)(mult * rng->array[rng->index]) + 0.5) * ulps);
    rng->index++;
    return value;
  }
  const unsigned long long x = gtc_rng_next(rng);
  return (GtcReal)(x >> 11) * (1.0 / 9007199254740992.0);
}

GtcReal gtc_rng_normal(GtcRng *rng) {
  const GtcReal u1 = fmax(gtc_rng_uniform(rng), 1.0e-12);
  const GtcReal u2 = gtc_rng_uniform(rng);
  return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
}

void rand_num_gen_init(GtcState *s) {
  if (s->p.rng_control > 0) {
    int seed[8];
    portable_seed_int(seed, s->p.rng_control - 1);
    if (s->rank == 0) {
      char printed[35];
      seed_to_decimal(seed, printed);
      fprintf(stderr, "Seed is set to %s\n", printed);
    }
    portable_init(&s->rng, seed);
  } else if (s->p.rng_control < 0) {
    int seed[8];
    portable_seed_time(seed);
    portable_step_seed(seed, s->rank);
    portable_init(&s->rng, seed);
  } else {
    unsigned long long base = 0x475443ULL + (unsigned long long)(s->rank + 1) * 111111ULL;
    if (s->p.irun == 0) base += 1ULL;
    else base += (unsigned long long)time(NULL);
    gtc_rng_init(&s->rng, base);
    s->rng.index = -1;
    if (s->rank == 0) {
      FILE *out = gtc_stdout_open(s, "a");
      if (out) {
        fprintf(out, " random_seed= C intrinsic-compatible seed path %llu\n", base);
        gtc_stdout_close(s, out);
      }
    }
  }
}
