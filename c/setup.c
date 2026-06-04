#include "gtc.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum {
  INPUT_IRUN = 0,
  INPUT_MSTEP,
  INPUT_MSNAP,
  INPUT_NDIAG,
  INPUT_NONLINEAR,
  INPUT_PARANL,
  INPUT_MODE00,
  INPUT_TSTEP,
  INPUT_MICELL,
  INPUT_MPSI,
  INPUT_MTHETAMAX,
  INPUT_MZETAMAX,
  INPUT_NPARTDOM,
  INPUT_A,
  INPUT_A0,
  INPUT_A1,
  INPUT_Q0,
  INPUT_Q1,
  INPUT_Q2,
  INPUT_RC,
  INPUT_RW,
  INPUT_AION,
  INPUT_QION,
  INPUT_KAPPATI,
  INPUT_KAPPATE,
  INPUT_KAPPAN,
  INPUT_TITE,
  INPUT_FLOW0,
  INPUT_FLOW1,
  INPUT_FLOW2,
  INPUT_R0,
  INPUT_B0,
  INPUT_TEMPERATURE,
  INPUT_EDENSITY0,
  INPUT_STDOUT,
  INPUT_NBOUND,
  INPUT_UMAX,
  INPUT_ILOAD,
  INPUT_RNG_CONTROL,
  INPUT_SPECTRUM_MODE,
  INPUT_NMODE,
  INPUT_MMODE,
  INPUT_COUNT
};

static const char *const input_names[INPUT_COUNT] = {
    "irun",       "mstep",    "msnap",       "ndiag",     "nonlinear", "paranl",
    "mode00",     "tstep",    "micell",      "mpsi",      "mthetamax", "mzetamax",
    "npartdom",   "a",        "a0",          "a1",        "q0",        "q1",
    "q2",         "rc",       "rw",          "aion",      "qion",      "kappati",
    "kappate",    "kappan",   "tite",        "flow0",     "flow1",     "flow2",
    "r0",         "b0",       "temperature", "edensity0", "stdout",    "nbound",
    "umax",       "iload",    "rng_control", "spectrum_mode", "nmode", "mmode"};

#define INPUT_BIT(name) (1ull << (INPUT_##name))

static const unsigned long long required_input_mask =
    INPUT_BIT(IRUN) | INPUT_BIT(MSTEP) | INPUT_BIT(MSNAP) | INPUT_BIT(NDIAG) |
    INPUT_BIT(NONLINEAR) | INPUT_BIT(PARANL) | INPUT_BIT(MODE00) | INPUT_BIT(TSTEP) |
    INPUT_BIT(MICELL) | INPUT_BIT(MPSI) | INPUT_BIT(MTHETAMAX) | INPUT_BIT(MZETAMAX) |
    INPUT_BIT(NPARTDOM) | INPUT_BIT(A) | INPUT_BIT(A0) | INPUT_BIT(A1) |
    INPUT_BIT(Q0) | INPUT_BIT(Q1) | INPUT_BIT(Q2) | INPUT_BIT(RC) | INPUT_BIT(RW) |
    INPUT_BIT(AION) | INPUT_BIT(QION) | INPUT_BIT(KAPPATI) | INPUT_BIT(KAPPATE) |
    INPUT_BIT(KAPPAN) | INPUT_BIT(TITE) | INPUT_BIT(FLOW0) | INPUT_BIT(FLOW1) |
    INPUT_BIT(FLOW2) | INPUT_BIT(R0) | INPUT_BIT(B0) | INPUT_BIT(TEMPERATURE) |
    INPUT_BIT(EDENSITY0) | INPUT_BIT(STDOUT) | INPUT_BIT(NBOUND) | INPUT_BIT(UMAX) |
    INPUT_BIT(ILOAD) | INPUT_BIT(RNG_CONTROL) | INPUT_BIT(SPECTRUM_MODE);

static void initialize_parameters(GtcParameters *p) {
  memset(p, 0, sizeof(*p));
}

static void mark_input(unsigned long long *seen, int input_id) {
  if (seen && input_id >= 0 && input_id < INPUT_COUNT) *seen |= 1ull << input_id;
}

static char *skip_space(char *p) {
  while (*p && (isspace((unsigned char)*p) || *p == ',')) p++;
  return p;
}

static void normalize_key(char *key) {
  char *src = key;
  while (*src && (isspace((unsigned char)*src) || *src == '&')) src++;
  char *dst = key;
  while (*src) {
    if (!isspace((unsigned char)*src)) *dst++ = (char)tolower((unsigned char)*src);
    src++;
  }
  *dst = '\0';
}

static void normalize_number_token(char *token) {
  for (char *p = token; *p; p++) {
    if (*p == 'd' || *p == 'D') *p = 'e';
  }
}

static int spectrum_mode_from_token(const char *token, int *value) {
  if (strcmp(token, "full") == 0 || strcmp(token, "full_n") == 0 || strcmp(token, "all") == 0 ||
      strcmp(token, "all_n") == 0) {
    *value = GTC_SPECTRUM_FULL_N;
    return 1;
  }
  if (strcmp(token, "single") == 0 || strcmp(token, "single_n") == 0 || strcmp(token, "selected") == 0 ||
      strcmp(token, "selected_n") == 0 || strcmp(token, "list") == 0) {
    *value = GTC_SPECTRUM_SINGLE_N;
    return 1;
  }
  return 0;
}

static int parse_int_list(char *value, int **dest) {
  free(*dest);
  *dest = NULL;
  int count = 0;
  char *cursor = value;
  while ((cursor = skip_space(cursor)) && *cursor && *cursor != '/') {
    char *end = cursor;
    long parsed = strtol(cursor, &end, 10);
    if (end == cursor) break;
    int *next = realloc(*dest, (size_t)(count + 1) * sizeof(**dest));
    if (!next) {
      free(*dest);
      *dest = NULL;
      return 0;
    }
    *dest = next;
    (*dest)[count++] = (int)parsed;
    cursor = end;
  }
  return count;
}

static void assign_values(GtcParameters *p, const char *key, char *value, int *num_mmode,
                          unsigned long long *seen) {
  if (strcmp(key, "nmode") == 0) {
    p->num_mode = parse_int_list(value, &p->nmode);
    if (p->num_mode > 0) mark_input(seen, INPUT_NMODE);
    return;
  }
  if (strcmp(key, "mmode") == 0) {
    *num_mmode = parse_int_list(value, &p->mmode);
    if (*num_mmode > 0) mark_input(seen, INPUT_MMODE);
    return;
  }

  GtcReal reals[64];
  int ints[64];
  int nr = 0, ni = 0;
  char *cursor = value;
  while ((cursor = skip_space(cursor)) && *cursor && *cursor != '/') {
    char token[128];
    int n = 0;
    while (cursor[n] && cursor[n] != ',' && cursor[n] != '/' && cursor[n] != '=' &&
           !isspace((unsigned char)cursor[n]) && n < (int)sizeof(token) - 1) {
      token[n] = cursor[n];
      n++;
    }
    token[n] = '\0';
    if (n == 0) break;
    normalize_number_token(token);
    if (strcmp(key, "spectrum_mode") == 0 && ni == 0 && nr == 0) {
      int mode = 0;
      if (spectrum_mode_from_token(token, &mode)) {
        p->spectrum_mode = mode;
        mark_input(seen, INPUT_SPECTRUM_MODE);
        return;
      }
    }
    char *end = token;
    GtcReal v = strtod(token, &end);
    if (end == token) break;
    if (nr < 64) reals[nr++] = v;
    if (ni < 64) ints[ni++] = (int)v;
    cursor += n;
  }

#define SETI(name, flag)   \
  if (strcmp(key, #name) == 0 && ni > 0) { \
    p->name = ints[0];     \
    mark_input(seen, flag); \
    return;                \
  }
#define SETR(name, flag)   \
  if (strcmp(key, #name) == 0 && nr > 0) { \
    p->name = reals[0];    \
    mark_input(seen, flag); \
    return;                \
  }
  SETI(irun, INPUT_IRUN) SETI(mstep, INPUT_MSTEP) SETI(msnap, INPUT_MSNAP)
  SETI(ndiag, INPUT_NDIAG) SETI(mode00, INPUT_MODE00) SETI(micell, INPUT_MICELL)
  SETI(mpsi, INPUT_MPSI) SETI(mthetamax, INPUT_MTHETAMAX)
  SETI(mzetamax, INPUT_MZETAMAX) SETI(npartdom, INPUT_NPARTDOM)
  SETI(nbound, INPUT_NBOUND) SETI(iload, INPUT_ILOAD)
  SETI(rng_control, INPUT_RNG_CONTROL) SETI(spectrum_mode, INPUT_SPECTRUM_MODE)
  if (strcmp(key, "stdout") == 0 && ni > 0) {
    p->stdout_unit = ints[0];
    mark_input(seen, INPUT_STDOUT);
    return;
  }
  SETR(nonlinear, INPUT_NONLINEAR) SETR(paranl, INPUT_PARANL)
  SETR(tstep, INPUT_TSTEP) SETR(a, INPUT_A) SETR(a0, INPUT_A0)
  SETR(a1, INPUT_A1) SETR(q0, INPUT_Q0) SETR(q1, INPUT_Q1)
  SETR(q2, INPUT_Q2) SETR(rc, INPUT_RC) SETR(rw, INPUT_RW)
  SETR(aion, INPUT_AION) SETR(qion, INPUT_QION)
  SETR(kappati, INPUT_KAPPATI) SETR(kappate, INPUT_KAPPATE)
  SETR(kappan, INPUT_KAPPAN) SETR(tite, INPUT_TITE)
  SETR(flow0, INPUT_FLOW0) SETR(flow1, INPUT_FLOW1) SETR(flow2, INPUT_FLOW2)
  SETR(r0, INPUT_R0) SETR(b0, INPUT_B0) SETR(temperature, INPUT_TEMPERATURE)
  SETR(edensity0, INPUT_EDENSITY0) SETR(umax, INPUT_UMAX)
#undef SETI
#undef SETR
}

static void append_missing_input(char *buffer, size_t buffer_size, const char *name) {
  if (buffer_size == 0) return;
  const size_t used = strlen(buffer);
  if (used >= buffer_size - 1) return;
  if (used > 0) strncat(buffer, ", ", buffer_size - strlen(buffer) - 1);
  strncat(buffer, name, buffer_size - strlen(buffer) - 1);
}

static void validate_required_inputs(GtcState *s, unsigned long long seen) {
  const unsigned long long missing_mask = required_input_mask & ~seen;
  if (missing_mask) {
    char missing[1024] = "";
    for (int i = 0; i < INPUT_COUNT; i++) {
      if (missing_mask & (1ull << i)) append_missing_input(missing, sizeof(missing), input_names[i]);
    }
    char message[1200];
    snprintf(message, sizeof(message), "missing required gtc.input parameter(s): %s", missing);
    gtc_die(s, message);
  }
  if (s->p.spectrum_mode < GTC_SPECTRUM_FULL_N || s->p.spectrum_mode > GTC_SPECTRUM_SINGLE_N) {
    gtc_die(s, "spectrum_mode must be 0/full_n or 1/selected_n");
  }
  if (s->p.spectrum_mode == GTC_SPECTRUM_SINGLE_N && !(seen & INPUT_BIT(NMODE))) {
    gtc_die(s, "missing required gtc.input parameter: nmode when spectrum_mode=1");
  }
  if (s->p.nonlinear < 0.5 && s->p.spectrum_mode == GTC_SPECTRUM_SINGLE_N &&
      !(seen & INPUT_BIT(MMODE))) {
    gtc_die(s, "missing required gtc.input parameter: mmode for linear selected-mode history");
  }
}

static void read_input_params(GtcState *s, unsigned long long *seen) {
  GtcParameters *p = &s->p;
  FILE *fp = fopen("gtc.input", "r");
  if (!fp) gtc_die(s, "cannot open required gtc.input");
  int num_mmode = 0;

  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    gtc_die(s, "cannot seek required gtc.input");
  }
  long len = ftell(fp);
  if (len < 0) {
    fclose(fp);
    gtc_die(s, "cannot determine required gtc.input size");
  }
  rewind(fp);
  char *raw = calloc((size_t)len + 2u, 1);
  char *buf = calloc((size_t)len + 2u, 1);
  if (!raw || !buf) {
    free(buf);
    free(raw);
    fclose(fp);
    gtc_die(s, "cannot allocate gtc.input parser buffer");
  }
  size_t nread = fread(raw, 1, (size_t)len, fp);
  raw[nread] = '\0';
  fclose(fp);

  int in_comment = 0;
  size_t out = 0;
  for (size_t i = 0; i < nread; i++) {
    char ch = raw[i];
    if (in_comment) {
      if (ch == '\n') {
        in_comment = 0;
        buf[out++] = ' ';
      }
      continue;
    }
    if (ch == '!') {
      in_comment = 1;
      continue;
    }
    if (ch == '\n' || ch == '\r' || ch == '\t') ch = ' ';
    buf[out++] = ch;
  }
  buf[out] = '\0';

  char *cursor = buf;
  while (*cursor) {
    while (*cursor && (*cursor == '/' || *cursor == ',' || isspace((unsigned char)*cursor))) cursor++;
    if (*cursor == '&') {
      cursor++;
      while (*cursor && (isalnum((unsigned char)*cursor) || *cursor == '_')) cursor++;
      continue;
    }
    if (!*cursor) break;
    char *eq = strchr(cursor, '=');
    if (!eq) break;
    char key[128];
    size_t key_len = (size_t)(eq - cursor);
    if (key_len >= sizeof(key)) key_len = sizeof(key) - 1;
    memcpy(key, cursor, key_len);
    key[key_len] = '\0';
    normalize_key(key);

    char *value = eq + 1;
    char *next = value;
    while (*next) {
      if (*next == '/') break;
      if (*next == '=') {
        char *probe = next;
        while (probe > value && (isspace((unsigned char)probe[-1]) || probe[-1] == ',')) probe--;
        while (probe > value && (isalnum((unsigned char)probe[-1]) || probe[-1] == '_')) probe--;
        next = probe;
        break;
      }
      next++;
    }
    char saved = *next;
    *next = '\0';
    assign_values(p, key, value, &num_mmode, seen);
    *next = saved;
    cursor = next;
    if (saved == '=') {
      while (cursor > buf && cursor[-1] != ',' && cursor[-1] != '/' && !isspace((unsigned char)cursor[-1])) cursor--;
    }
  }

  free(buf);
  free(raw);

  if (p->mmode) {
    int *source = p->mmode;
    p->mmode = NULL;
    if (p->num_mode > 0 && num_mmode >= p->num_mode) {
      int *mmode = calloc((size_t)p->num_mode, sizeof(*mmode));
      if (mmode) {
        for (int i = 0; i < p->num_mode; i++) mmode[i] = source[i];
        p->mmode = mmode;
      } else {
        free(source);
        gtc_die(s, "cannot allocate mmode input list");
      }
    }
    free(source);
  }
}

static void write_input_parameters(FILE *out, const GtcParameters *p) {
  fprintf(out, " &input_parameters\n");
  fprintf(out, " irun=%d, mstep=%d, msnap=%d, ndiag=%d,\n", p->irun, p->mstep, p->msnap, p->ndiag);
  fprintf(out, " nonlinear=% .8e, paranl=% .8e, mode00=%d, tstep=% .8e,\n",
          p->nonlinear, p->paranl, p->mode00, p->tstep);
  fprintf(out, " micell=%d, mpsi=%d, mthetamax=%d, mzetamax=%d, npartdom=%d,\n",
          p->micell, p->mpsi, p->mthetamax, p->mzetamax, p->npartdom);
  fprintf(out, " a=% .8e, a0=% .8e, a1=% .8e, q0=% .8e, q1=% .8e, q2=% .8e,\n",
          p->a, p->a0, p->a1, p->q0, p->q1, p->q2);
  fprintf(out, " rc=% .8e, rw=% .8e, aion=% .8e, qion=% .8e,\n",
          p->rc, p->rw, p->aion, p->qion);
  fprintf(out, " kappati=% .8e, kappate=% .8e, kappan=% .8e, tite=% .8e,\n",
          p->kappati, p->kappate, p->kappan, p->tite);
  fprintf(out, " flow0=% .8e, flow1=% .8e, flow2=% .8e,\n",
          p->flow0, p->flow1, p->flow2);
  fprintf(out, " r0=% .8e, b0=% .8e, temperature=% .8e, edensity0=% .8e,\n",
          p->r0, p->b0, p->temperature, p->edensity0);
  fprintf(out, " stdout=%d, nbound=%d, umax=% .8e, iload=%d,\n",
          p->stdout_unit, p->nbound, p->umax, p->iload);
  fprintf(out, " ! spectrum_mode: 0 keeps all toroidal n in the field, 1 filters evolution to nmode\n");
  fprintf(out, " ! nmode selects toroidal mode numbers for history and snapshot diagnostics\n");
  fprintf(out, " rng_control=%d, spectrum_mode=%d,\n", p->rng_control, p->spectrum_mode);
  fprintf(out, " nmode=");
  for (int i = 0; i < p->num_mode; i++) fprintf(out, "%s%d", i == 0 ? "" : ",", p->nmode[i]);
  if (p->mmode) {
    fprintf(out, "\n mmode=");
    for (int i = 0; i < p->num_mode; i++) fprintf(out, "%s%d", i == 0 ? "" : ",", p->mmode[i]);
  }
  fprintf(out, "\n /\n");
}

void setup(GtcState *s) {
  initialize_parameters(&s->p);
  s->partd_comm = MPI_COMM_NULL;
  s->toroidal_comm = MPI_COMM_NULL;
  MPI_Comm_rank(MPI_COMM_WORLD, &s->rank);
  MPI_Comm_size(MPI_COMM_WORLD, &s->size);
  unsigned long long seen = 0;
  if (s->rank == 0) {
    read_input_params(s, &seen);
    validate_required_inputs(s, seen);
  }

  int *root_nmode = s->p.nmode;
  int *root_mmode = s->p.mmode;
  MPI_Bcast(&s->p, sizeof(s->p), MPI_BYTE, 0, MPI_COMM_WORLD);
  if (s->rank == 0) {
    s->p.nmode = root_nmode;
    s->p.mmode = root_mmode;
  } else {
    s->p.nmode = NULL;
    s->p.mmode = NULL;
  }
  if (s->p.num_mode < 0) s->p.num_mode = 0;
  const int filter_spectrum = s->p.spectrum_mode == GTC_SPECTRUM_SINGLE_N;
  const int need_history_modes = s->p.nonlinear < 0.5 && filter_spectrum;
  if (!filter_spectrum && s->p.num_mode < 1) {
    free(s->p.nmode);
    free(s->p.mmode);
    s->p.nmode = NULL;
    s->p.mmode = NULL;
    s->p.num_mode = 0;
  }
  if (filter_spectrum && s->p.num_mode < 1) {
    gtc_die(s, "nmode list is required when spectrum_mode filters the field evolution");
  }
  if (need_history_modes && s->rank == 0 && !s->p.mmode) {
    gtc_die(s, "mmode list is required for linear selected-mode history diagnostics");
  }
  if (s->p.num_mode > 0) {
    if (s->rank != 0) s->p.nmode = gtc_xcalloc(s, (size_t)s->p.num_mode, sizeof(*s->p.nmode), "nmode");
    MPI_Bcast(s->p.nmode, s->p.num_mode, MPI_INT, 0, MPI_COMM_WORLD);
  }
  if (need_history_modes) {
    if (s->rank != 0) {
      s->p.mmode = gtc_xcalloc(s, (size_t)s->p.num_mode, sizeof(*s->p.mmode), "mmode");
    }
    MPI_Bcast(s->p.mmode, s->p.num_mode, MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    free(s->p.mmode);
    s->p.mmode = NULL;
  }
  s->pi = acos(-1.0);
  s->p.a0 *= s->p.a;
  s->p.a1 *= s->p.a;
  s->p.mstep = s->p.mstep < 2 ? 2 : s->p.mstep;
  s->p.ndiag = s->p.ndiag < 1 ? 1 : s->p.ndiag;
  s->idiag1 = s->p.mpsi / 2;
  s->idiag2 = s->p.mpsi / 2;
  if (s->p.nonlinear < 0.5) {
    s->p.paranl = 0.0;
    s->p.mode00 = 0;
  }
  if (filter_spectrum) {
    s->idiag1 = 1;
    s->idiag2 = s->p.mpsi;
  }
  if (s->p.npartdom < 1) s->p.npartdom = 1;
  while (s->size % s->p.npartdom != 0) {
    s->p.npartdom--;
    if (s->p.npartdom == 1) break;
  }
  s->ntoroidal = s->size / s->p.npartdom;
  if (s->ntoroidal < 1) s->ntoroidal = 1;
  if (s->rank == 0) {
    FILE *out = gtc_stdout_open(s, "w");
    if (out) {
      time_t now = time(NULL);
      struct tm tmv;
      localtime_r(&now, &tmv);
      char stamp[64];
      strftime(stamp, sizeof(stamp), "%Y%m%d TIME=%H%M%S", &tmv);
      fprintf(out, " Program starts at DATE=%s\n", stamp);
      write_input_parameters(out, &s->p);
      fprintf(out, "\n===================================\n");
      fprintf(out, " Run with MPI ranks=%d, toroidal ranks=%d, OpenMP threads/rank=%d\n",
              s->size, s->ntoroidal, gtc_omp_max_threads());
      fprintf(out, "===================================\n\n");
      fprintf(out, "*******************************************************\n");
      fprintf(out, "  Using npartdom= %d and ntoroidal= %d\n", s->p.npartdom, s->ntoroidal);
      fprintf(out, "  Requested mzetamax= %d\n", s->p.mzetamax);
      fprintf(out, "*******************************************************\n\n");
      gtc_stdout_close(s, out);
    }
  }
  s->p.mzetamax = s->ntoroidal * (int)fmax(1.0, floor((GtcReal)s->p.mzetamax / (GtcReal)s->ntoroidal + 0.5));
  s->p.mpsi = 2 * (s->p.mpsi / 2);
  s->particle_domain_location = s->rank % s->p.npartdom;
  s->toroidal_domain_location = s->rank / s->p.npartdom;
  MPI_Comm_split(MPI_COMM_WORLD, s->toroidal_domain_location, s->particle_domain_location, &s->partd_comm);
  MPI_Comm_split(MPI_COMM_WORLD, s->particle_domain_location, s->toroidal_domain_location, &s->toroidal_comm);
  MPI_Comm_rank(s->partd_comm, &s->myrank_partd);
  MPI_Comm_rank(s->toroidal_comm, &s->myrank_toroidal);
  s->left_pe = (s->myrank_toroidal + s->ntoroidal - 1) % s->ntoroidal;
  s->right_pe = (s->myrank_toroidal + 1) % s->ntoroidal;
  if (s->p.msnap < 1) s->p.msnap = 1;
  if (s->p.msnap > s->p.mstep / s->p.ndiag) s->p.msnap = s->p.mstep / s->p.ndiag;
  if (s->p.msnap < 1) s->p.msnap = 1;

  s->mzeta = s->p.mzetamax / s->ntoroidal;
  if (s->mzeta < 1) s->mzeta = 1;
  s->zetamin = 2.0 * s->pi * (GtcReal)s->toroidal_domain_location / (GtcReal)s->ntoroidal;
  s->zetamax = 2.0 * s->pi * (GtcReal)(s->toroidal_domain_location + 1) / (GtcReal)s->ntoroidal;
  s->deltaz = (s->zetamax - s->zetamin) / (GtcReal)s->mzeta;
  const GtcReal ulength = s->p.r0;
  s->gyroradius = 102.0 * sqrt(s->p.aion * s->p.temperature) /
                  (fabs(s->p.qion) * s->p.b0) / ulength;
  s->p.tstep = s->p.tstep * s->p.aion / (fabs(s->p.qion) * s->gyroradius * s->p.kappati);
  s->p.rc *= s->p.a0 + s->p.a1;
  s->p.rw = 1.0 / (s->p.rw * (s->p.a1 - s->p.a0));
  s->deltar = (s->p.a1 - s->p.a0) / (GtcReal)s->p.mpsi;

  const int npsi = s->p.mpsi + 1;
  s->mtheta = gtc_xcalloc(s, npsi, sizeof(int), "mtheta");
  s->igrid = gtc_xcalloc(s, npsi, sizeof(int), "igrid");
  s->itran = gtc_xcalloc(s, npsi, sizeof(int), "itran");
  s->deltat = gtc_xcalloc(s, npsi, sizeof(*s->deltat), "deltat");
  s->qtinv = gtc_xcalloc(s, npsi, sizeof(*s->qtinv), "qtinv");

  const GtcReal theta_unit = s->pi * s->p.a / (GtcReal)s->p.mthetamax;
  for (int i = 0; i <= s->p.mpsi; i++) {
    const GtcReal r = s->p.a0 + s->deltar * (GtcReal)i;
    s->mtheta[i] = 2 * (int)fmax(1.0, floor(s->pi * r / theta_unit + 0.5));
    s->deltat[i] = gtc_real(2.0 * s->pi / (GtcReal)s->mtheta[i]);
    const GtcReal q = s->p.q0 + s->p.q1 * r / s->p.a + s->p.q2 * r * r / (s->p.a * s->p.a);
    int itran = (int)floor((GtcReal)s->mtheta[i] / q + 0.5);
    if (itran < 1) itran = 1;
    s->qtinv[i] = gtc_real((GtcReal)itran / (GtcReal)s->mtheta[i]);
    s->itran[i] = itran - s->mtheta[i] * (itran / s->mtheta[i]);
  }
  s->mtdiag = (s->p.mthetamax / s->p.mzetamax) * s->p.mzetamax;
  if (s->mtdiag < s->p.mzetamax) s->mtdiag = s->p.mzetamax;
  s->p.mthetamax = s->mtheta[s->p.mpsi];
  s->igrid[0] = 0;
  for (int i = 1; i <= s->p.mpsi; i++) s->igrid[i] = s->igrid[i - 1] + s->mtheta[i - 1] + 1;
  s->mgrid = s->igrid[s->p.mpsi] + s->mtheta[s->p.mpsi] + 1;

  const int mi_local = s->p.micell * (s->mgrid - s->p.mpsi) * s->mzeta;
  s->mi = mi_local / s->p.npartdom;
  if (s->mi < mi_local % s->p.npartdom) s->mi++;
  if (s->mi < 1) s->mi = 1;
  s->mimax = s->mi + 100 * (int)ceil(sqrt((GtcReal)s->mi));
  s->isnap = s->p.mstep / s->p.msnap;
  s->irest = 1;
  s->file_exit = 0;

  s->pmarki = gtc_xcalloc(s, npsi, sizeof(*s->pmarki), "pmarki");
  s->phi00 = gtc_xcalloc(s, npsi, sizeof(*s->phi00), "phi00");
  s->phip00 = gtc_xcalloc(s, npsi, sizeof(*s->phip00), "phip00");
  s->zonali = gtc_xcalloc(s, npsi, sizeof(*s->zonali), "zonali");
  s->gradt = gtc_xcalloc(s, npsi, sizeof(*s->gradt), "gradt");
  s->rden = gtc_xcalloc(s, npsi, sizeof(*s->rden), "rden");
  s->rtemi = gtc_xcalloc(s, npsi, sizeof(*s->rtemi), "rtemi");
  s->hfluxpsi = gtc_xcalloc(s, npsi, sizeof(*s->hfluxpsi), "hfluxpsi");
  s->pfluxpsi = gtc_xcalloc(s, npsi, sizeof(*s->pfluxpsi), "pfluxpsi");
  s->rdtemi = gtc_xcalloc(s, npsi, sizeof(*s->rdtemi), "rdtemi");
  if (gtc_history_mode_count(s) > 0) {
    s->amp_mode = gtc_xcalloc(s, (size_t)2 * (size_t)gtc_history_mode_count(s) * (size_t)2,
                              sizeof(*s->amp_mode), "amp_mode");
  }
  if (s->p.num_mode > 0) {
    s->eigenmode = gtc_xcalloc(s, (size_t)s->p.mpsi * (size_t)s->p.num_mode *
                                      (size_t)GTC_M_POLIDAL,
                                sizeof(*s->eigenmode), "eigenmode");
  }
  for (int i = 0; i <= s->p.mpsi; i++) {
    s->pmarki[i] = 1.0;
    s->rden[i] = 1.0;
    s->rtemi[i] = 1.0;
    const GtcReal r = s->p.a0 + (s->p.a1 - s->p.a0) * ((GtcReal)i - 0.5) / (GtcReal)s->p.mpsi;
    GtcReal rfac = s->p.rw * (r - s->p.rc);
    rfac = rfac * rfac;
    rfac = rfac * rfac * rfac;
    rfac = fmax(0.1, exp(-rfac));
    GtcReal kappa = 1.0;
    if (s->p.nbound == 0) kappa = 0.0;
    kappa = 1.0 - kappa + kappa * rfac;
    s->gradt[i] = gtc_real(1.0 / (kappa * s->p.kappati * s->gyroradius));
  }

  const size_t field_size = (size_t)(s->mzeta + 1) * (size_t)s->mgrid;
  s->densityi = gtc_xcalloc(s, field_size, sizeof(*s->densityi), "densityi");
  s->phi = gtc_xcalloc(s, field_size, sizeof(*s->phi), "phi");
  s->markeri = gtc_xcalloc(s, field_size, sizeof(*s->markeri), "markeri");
  s->evector = gtc_xcalloc(s, 3 * field_size, sizeof(*s->evector), "evector");
  s->heatflux = gtc_xcalloc(s, field_size, sizeof(*s->heatflux), "heatflux");
  s->dtemper = gtc_xcalloc(s, field_size, sizeof(*s->dtemper), "dtemper");
  s->pgyro = gtc_xcalloc(s, 4 * (size_t)s->mgrid, sizeof(*s->pgyro), "pgyro");
  s->tgyro = gtc_xcalloc(s, 4 * (size_t)s->mgrid, sizeof(*s->tgyro), "tgyro");
  s->zion = gtc_xcalloc(s, (size_t)GTC_NPARAM * (size_t)s->mimax, sizeof(*s->zion), "zion");
  s->zion0 = gtc_xcalloc(s, (size_t)GTC_NPARAM * (size_t)s->mimax, sizeof(*s->zion0), "zion0");
  s->kzion = gtc_xcalloc(s, (size_t)s->mimax, sizeof(int), "kzion");
  s->jtion0 = gtc_xcalloc(s, 4 * (size_t)s->mimax, sizeof(int), "jtion0");
  s->jtion1 = gtc_xcalloc(s, 4 * (size_t)s->mimax, sizeof(int), "jtion1");
  s->jtp1 = gtc_xcalloc(s, 2 * (size_t)s->mgrid * (size_t)(s->mzeta + 1), sizeof(int), "jtp1");
  s->jtp2 = gtc_xcalloc(s, 2 * (size_t)s->mgrid * (size_t)(s->mzeta + 1), sizeof(int), "jtp2");
  s->wzion = gtc_xcalloc(s, (size_t)s->mimax, sizeof(*s->wzion), "wzion");
  s->wpion = gtc_xcalloc(s, 4 * (size_t)s->mimax, sizeof(*s->wpion), "wpion");
  s->wtion0 = gtc_xcalloc(s, 4 * (size_t)s->mimax, sizeof(*s->wtion0), "wtion0");
  s->wtion1 = gtc_xcalloc(s, 4 * (size_t)s->mimax, sizeof(*s->wtion1), "wtion1");
  s->wpi = gtc_xcalloc(s, 3 * (size_t)s->mimax, sizeof(*s->wpi), "wpi");
  s->wtp1 = gtc_xcalloc(s, 2 * (size_t)s->mgrid * (size_t)(s->mzeta + 1), sizeof(*s->wtp1), "wtp1");
  s->wtp2 = gtc_xcalloc(s, 2 * (size_t)s->mgrid * (size_t)(s->mzeta + 1), sizeof(*s->wtp2), "wtp2");

  for (int i = 0; i <= s->p.mpsi; i++) {
    const GtcReal r = s->p.a0 + s->deltar * (GtcReal)i;
    GtcReal shell = 0.0;
    for (int j = 0; j <= s->mtheta[i]; j++) {
      const int ij = s->igrid[i] + j;
      const GtcReal theta = s->deltat[i] * (GtcReal)j;
      const GtcReal b = 1.0 / (1.0 + r * cos(theta));
      const GtcReal rhoi = sqrt(2.0 / b) * s->gyroradius;
      s->pgyro[gtc_gyro_index(s, 0, ij)] = gtc_real(-rhoi);
      s->pgyro[gtc_gyro_index(s, 1, ij)] = gtc_real(rhoi);
      s->pgyro[gtc_gyro_index(s, 2, ij)] = gtc_real(0.5 * rhoi * rhoi / r);
      s->pgyro[gtc_gyro_index(s, 3, ij)] = gtc_real(0.5 * rhoi * rhoi / r);
      s->tgyro[gtc_gyro_index(s, 0, ij)] = 0.0;
      s->tgyro[gtc_gyro_index(s, 1, ij)] = 0.0;
      s->tgyro[gtc_gyro_index(s, 2, ij)] = gtc_real(-rhoi / r);
      s->tgyro[gtc_gyro_index(s, 3, ij)] = gtc_real(rhoi / r);
    }
    for (int j = 1; j <= s->mtheta[i]; j++) {
      const int ij = s->igrid[i] + j;
      for (int k = 1; k <= s->mzeta; k++) {
        const GtcReal z = s->zetamin + s->deltaz * (GtcReal)k;
        const GtcReal td = s->deltat[i] * (GtcReal)j + z * s->qtinv[i];
        const GtcReal marker = (1.0 + r * cos(td)) * (1.0 + r * cos(td));
        s->markeri[gtc_grid_index(s, k, ij)] = gtc_real(marker);
        shell += marker;
      }
    }
    const GtcReal rmax = fmin(s->p.a1, r + 0.5 * s->deltar);
    const GtcReal rmin = fmax(s->p.a0, r - 0.5 * s->deltar);
    const GtcReal tdum = (GtcReal)(s->mi * (s->p.npartdom > 0 ? s->p.npartdom : 1)) *
                        (rmax * rmax - rmin * rmin) /
                        (s->p.a1 * s->p.a1 - s->p.a0 * s->p.a0);
    for (int j = 1; j <= s->mtheta[i]; j++) {
      const int ij = s->igrid[i] + j;
      for (int k = 1; k <= s->mzeta; k++) {
        size_t idx = gtc_grid_index(s, k, ij);
        GtcReal value = shell > 0.0 ? tdum * s->markeri[idx] / shell : 1.0;
        s->markeri[idx] = value != 0.0 ? gtc_real(1.0 / value) : 1.0;
      }
    }
    for (int k = 1; k <= s->mzeta; k++) {
      s->markeri[gtc_grid_index(s, k, s->igrid[i])] =
          s->markeri[gtc_grid_index(s, k, s->igrid[i] + s->mtheta[i])];
    }
    s->pmarki[i] = tdum != 0.0 ? gtc_real(1.0 / ((GtcReal)s->ntoroidal * tdum)) : 1.0;
  }

  for (int k = 1; k <= s->mzeta; k++) {
    const GtcReal zdum = s->zetamin + s->deltaz * (GtcReal)k;
    for (int i = 1; i < s->p.mpsi; i++) {
      for (int ip = 1; ip <= 2; ip++) {
        const int indp = i + ip < s->p.mpsi ? i + ip : s->p.mpsi;
        const int indt = i - ip > 0 ? i - ip : 0;
        for (int j = 1; j <= s->mtheta[i]; j++) {
          const int ij = s->igrid[i] + j;
          GtcReal tdum = ((GtcReal)j * s->deltat[i] + zdum * (s->qtinv[i] - s->qtinv[indp])) / s->deltat[indp];
          int jt = (int)floor(tdum);
          GtcReal wt = tdum - (GtcReal)jt;
          jt = ((jt % s->mtheta[indp]) + s->mtheta[indp]) % s->mtheta[indp];
          size_t base = gtc_interp_index(s, 0, k, ij);
          if (ip == 1) {
            s->wtp1[base] = gtc_real(wt);
            s->jtp1[base] = s->igrid[indp] + jt;
          } else {
            s->wtp2[base] = gtc_real(wt);
            s->jtp2[base] = s->igrid[indp] + jt;
          }

          tdum = ((GtcReal)j * s->deltat[i] + zdum * (s->qtinv[i] - s->qtinv[indt])) / s->deltat[indt];
          jt = (int)floor(tdum);
          wt = tdum - (GtcReal)jt;
          jt = ((jt % s->mtheta[indt]) + s->mtheta[indt]) % s->mtheta[indt];
          base = gtc_interp_index(s, 1, k, ij);
          if (ip == 1) {
            s->wtp1[base] = gtc_real(wt);
            s->jtp1[base] = s->igrid[indt] + jt;
          } else {
            s->wtp2[base] = gtc_real(wt);
            s->jtp2[base] = s->igrid[indt] + jt;
          }
        }
      }
    }
  }

  if (s->rank == 0) {
    FILE *out = gtc_stdout_open(s, "a");
    if (out) {
      fprintf(out, " C GTC port starts: mstep=%d mpsi=%d mgrid=%d local_ions=%d local_mzeta=%d global_mzetamax=%d\n",
              s->p.mstep, s->p.mpsi, s->mgrid, s->mi, s->mzeta, s->p.mzetamax);
      gtc_stdout_close(s, out);
    }
  }
}
