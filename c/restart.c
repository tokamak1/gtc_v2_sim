#include "gtc.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>

static void restart_name(const GtcState *s, char *name, size_t n) {
  const char *dir = (s->irest % 2) == 0 ? "restart_dir1" : "restart_dir2";
  snprintf(name, n, "%s/restart_%05d.bp", dir, s->myrank_toroidal);
}

static void write_file_exit_restart(const GtcState *s) {
  if (s->rank != 0) return;
  FILE *fp = fopen("FileExit.dat", "w");
  if (!fp) return;
  fprintf(fp, "FileExit=%1d\n", s->file_exit);
  fprintf(fp, "irest   =%5d\n", s->irest + 1);
  fprintf(fp, "%s\n", (s->irest % 2) == 0 ? "restart_dir1" : "restart_dir2");
  fclose(fp);
}

static void copy_file(const char *src, const char *dst) {
  FILE *in = fopen(src, "r");
  if (!in) return;
  FILE *out = fopen(dst, "w");
  if (!out) {
    fclose(in);
    return;
  }
  char buf[8192];
  size_t nread;
  while ((nread = fread(buf, 1, sizeof(buf), in)) > 0) fwrite(buf, 1, nread, out);
  fclose(out);
  fclose(in);
}

static int read_int_line(FILE *fp, int *value) {
  char line[256];
  if (!fgets(line, sizeof(line), fp)) return 0;
  return sscanf(line, "%d", value) == 1;
}

static int read_real_line(FILE *fp, GtcReal *value) {
  char line[256];
  if (!fgets(line, sizeof(line), fp)) return 0;
  return sscanf(line, "%f", value) == 1;
}

static void write_restart_history(const GtcState *s, const char *dst) {
  FILE *in = fopen("history.out", "r");
  if (!in) return;
  FILE *out = fopen(dst, "w");
  if (!out) {
    fclose(in);
    return;
  }
  int irun = 0, mquantity = 0, mflux = 0, nmode = 0, mstepfinal = 0;
  if (!read_int_line(in, &irun) || !read_int_line(in, &mquantity) ||
      !read_int_line(in, &mflux) || !read_int_line(in, &nmode) ||
      !read_int_line(in, &mstepfinal)) {
    fclose(out);
    fclose(in);
    return;
  }
  const int noutputs = mstepfinal - s->p.mstep / s->p.ndiag + s->istep / s->p.ndiag;
  fprintf(out, "%6d\n%6d\n%6d\n%6d\n%6d\n", irun, mquantity, mflux, nmode, noutputs);
  const int nreal = (mquantity + 4 * nmode) * noutputs + 1;
  for (int i = 0; i < nreal; i++) {
    GtcReal value = 0.0;
    if (!read_real_line(in, &value)) break;
    gtc_fprintf_e(out, 6, value);
    fputc('\n', out);
  }
  fclose(out);
  fclose(in);
}

static void checked_write(GtcState *s, const void *ptr, size_t size, size_t count, FILE *fp, const char *what) {
  if (fwrite(ptr, size, count, fp) != count) {
    char msg[256];
    snprintf(msg, sizeof(msg), "restart write failed for %s: %s", what, strerror(errno));
    gtc_die(s, msg);
  }
}

static void checked_read(GtcState *s, void *ptr, size_t size, size_t count, FILE *fp, const char *what) {
  if (fread(ptr, size, count, fp) != count) {
    char msg[256];
    snprintf(msg, sizeof(msg), "restart read failed for %s", what);
    gtc_die(s, msg);
  }
}

void restart_io(GtcState *s, const char *iop) {
  char name[128];
  restart_name(s, name, sizeof(name));
  const size_t npsi = (size_t)s->p.mpsi + 1;
  const size_t field_size = (size_t)(s->mzeta + 1) * (size_t)s->mgrid;

  if (strcmp(iop, "write") == 0) {
    FILE *fp = fopen(name, "wb");
    if (!fp) {
      char msg[256];
      snprintf(msg, sizeof(msg), "cannot open restart output %s", name);
      gtc_die(s, msg);
    }
    int header[7] = {GTC_MFLUX, s->p.mpsi, s->mzeta, s->mi, s->mimax, s->mgrid, GTC_NPARAM};
    checked_write(s, header, sizeof(header[0]), 7, fp, "header");
    checked_write(s, s->rdtemi, sizeof(*s->rdtemi), npsi, fp, "rdtemi");
    checked_write(s, s->pfluxpsi, sizeof(*s->pfluxpsi), npsi, fp, "pfluxpsi");
    checked_write(s, s->phi00, sizeof(*s->phi00), npsi, fp, "phi00");
    checked_write(s, s->phip00, sizeof(*s->phip00), npsi, fp, "phip00");
    checked_write(s, s->zonali, sizeof(*s->zonali), npsi, fp, "zonali");
    for (int m = 0; m < s->mimax; m++) {
      const GtcReal marker = m < s->mi ? s->zion0[gtc_particle_index(s, GTC_Z_MU, m)] : 0.0;
      checked_write(s, &marker, sizeof(marker), 1, fp, "zion0 marker");
    }
    checked_write(s, s->phi, sizeof(*s->phi), field_size, fp, "phi");
    checked_write(s, s->zion, sizeof(*s->zion), (size_t)GTC_NPARAM * (size_t)s->mimax, fp, "zion");
    fclose(fp);

    if (s->rank == 0 && s->istep <= s->p.mstep) {
      const char *dir = (s->irest % 2) == 0 ? "restart_dir1" : "restart_dir2";
      char dst[128];
      snprintf(dst, sizeof(dst), "%s/history_restart.out", dir);
      write_restart_history(s, dst);
      snprintf(dst, sizeof(dst), "%s/sheareb_restart.out", dir);
      copy_file("sheareb.out", dst);
    }
    write_file_exit_restart(s);
    s->irest++;
  } else if (strcmp(iop, "read") == 0) {
    FILE *fp = fopen(name, "rb");
    if (!fp) {
      char msg[256];
      snprintf(msg, sizeof(msg), "cannot open restart input %s", name);
      gtc_die(s, msg);
    }
    int header[7] = {0};
    checked_read(s, header, sizeof(header[0]), 7, fp, "header");
    if (header[0] != GTC_MFLUX || header[1] != s->p.mpsi || header[2] != s->mzeta ||
        header[3] != s->mi || header[4] != s->mimax || header[5] != s->mgrid ||
        header[6] != GTC_NPARAM) {
      gtc_die(s, "restart file dimensions do not match current run");
    }
    checked_read(s, s->rdtemi, sizeof(*s->rdtemi), npsi, fp, "rdtemi");
    checked_read(s, s->pfluxpsi, sizeof(*s->pfluxpsi), npsi, fp, "pfluxpsi");
    checked_read(s, s->phi00, sizeof(*s->phi00), npsi, fp, "phi00");
    checked_read(s, s->phip00, sizeof(*s->phip00), npsi, fp, "phip00");
    checked_read(s, s->zonali, sizeof(*s->zonali), npsi, fp, "zonali");
    for (int m = 0; m < s->mimax; m++) {
      GtcReal marker = 0.0;
      checked_read(s, &marker, sizeof(marker), 1, fp, "zion0 marker");
      if (m < s->mi) s->zion0[gtc_particle_index(s, GTC_Z_MU, m)] = marker;
    }
    checked_read(s, s->phi, sizeof(*s->phi), field_size, fp, "phi");
    checked_read(s, s->zion, sizeof(*s->zion), (size_t)GTC_NPARAM * (size_t)s->mimax, fp, "zion");
    fclose(fp);
    s->irest++;
  } else {
    gtc_die(s, "restart_io iop must be read or write");
  }
}
