#include "gtc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { PACK_FIELDS = 2 * GTC_NPARAM };

static void pack_particle(const GtcState *s, int m, GtcReal *buf) {
  for (int p = 0; p < GTC_NPARAM; p++) {
    buf[p] = s->zion[gtc_particle_index(s, p, m)];
    buf[GTC_NPARAM + p] = s->zion0[gtc_particle_index(s, p, m)];
  }
}

static void unpack_particle(GtcState *s, int m, const GtcReal *buf) {
  for (int p = 0; p < GTC_NPARAM; p++) {
    s->zion[gtc_particle_index(s, p, m)] = buf[p];
    s->zion0[gtc_particle_index(s, p, m)] = buf[GTC_NPARAM + p];
  }
}

static void copy_particle(GtcState *s, int dst, int src) {
  for (int p = 0; p < GTC_NPARAM; p++) {
    s->zion[gtc_particle_index(s, p, dst)] = s->zion[gtc_particle_index(s, p, src)];
    s->zion0[gtc_particle_index(s, p, dst)] = s->zion0[gtc_particle_index(s, p, src)];
  }
}

void shifti(GtcState *s) {
  const int left = s->left_pe;
  const int right = s->right_pe;
  const GtcReal two_pi = 2.0 * s->pi;
  const GtcReal pi_inv = 1.0 / s->pi;
  int m0 = 0;

  if (s->ntoroidal == 1) return;

  for (int iteration = 1; iteration <= s->ntoroidal; iteration++) {
    int *move = gtc_xcalloc(s, (size_t)(s->mi > 0 ? s->mi : 1), sizeof(int), "shifti move");
    int *ileft = gtc_xcalloc(s, (size_t)(s->mi > 0 ? s->mi : 1), sizeof(int), "shifti ileft");
    int *iright = gtc_xcalloc(s, (size_t)(s->mi > 0 ? s->mi : 1), sizeof(int), "shifti iright");
    int nmove = 0, nleft = 0, nright = 0;

    for (int m = m0; m < s->mi; m++) {
      const GtcReal zeta = s->zion[gtc_particle_index(s, GTC_Z_ZETA, m)];
      GtcReal zetaright = fmin(two_pi, zeta) - s->zetamax;
      GtcReal zetaleft = zeta - s->zetamin;
      if (zetaright * zetaleft > 0.0) {
        GtcReal zr = zetaright * 0.5 * pi_inv;
        zr = zr - floor(zr);
        move[nmove++] = m;
        if (zr < 0.5) {
          iright[nright++] = m;
        } else {
          ileft[nleft++] = m;
        }
      }
    }

    if (iteration > 1) {
      int send_total = nleft + nright;
      int global_total = 0;
      MPI_Allreduce(&send_total, &global_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (global_total == 0) {
        free(move);
        free(ileft);
        free(iright);
        return;
      }
    }

    GtcReal *send_right = gtc_xcalloc(s, (size_t)(nright > 0 ? nright : 1) * PACK_FIELDS, sizeof(*send_right), "shifti send right");
    GtcReal *send_left = gtc_xcalloc(s, (size_t)(nleft > 0 ? nleft : 1) * PACK_FIELDS, sizeof(*send_left), "shifti send left");
    for (int i = 0; i < nright; i++) pack_particle(s, iright[i], &send_right[(size_t)i * PACK_FIELDS]);
    for (int i = 0; i < nleft; i++) pack_particle(s, ileft[i], &send_left[(size_t)i * PACK_FIELDS]);

    int mtop = s->mi - 1;
    const int new_mi = s->mi - nmove;
    int lasth = nmove - 1;
    for (int i = 0; i < nmove; i++) {
      const int m = move[i];
      if (m >= new_mi) break;
      while (lasth >= 0 && mtop == move[lasth]) {
        mtop--;
        lasth--;
      }
      if (mtop < new_mi) break;
      copy_particle(s, m, mtop);
      mtop--;
      if (mtop < new_mi) break;
    }
    s->mi = new_mi;

    int recv_left_count = 0;
    int recv_right_count = 0;
    MPI_Sendrecv(&nright, 1, MPI_INT, right, s->myrank_toroidal,
                 &recv_left_count, 1, MPI_INT, left, left,
                 s->toroidal_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&nleft, 1, MPI_INT, left, s->myrank_toroidal,
                 &recv_right_count, 1, MPI_INT, right, right,
                 s->toroidal_comm, MPI_STATUS_IGNORE);

    GtcReal *recv_left = gtc_xcalloc(s, (size_t)(recv_left_count > 0 ? recv_left_count : 1) * PACK_FIELDS, sizeof(*recv_left), "shifti recv left");
    GtcReal *recv_right = gtc_xcalloc(s, (size_t)(recv_right_count > 0 ? recv_right_count : 1) * PACK_FIELDS, sizeof(*recv_right), "shifti recv right");

    MPI_Sendrecv(send_right, nright * PACK_FIELDS, GTC_MPI_REAL, right, s->myrank_toroidal,
                 recv_left, recv_left_count * PACK_FIELDS, GTC_MPI_REAL, left, left,
                 s->toroidal_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(send_left, nleft * PACK_FIELDS, GTC_MPI_REAL, left, s->myrank_toroidal,
                 recv_right, recv_right_count * PACK_FIELDS, GTC_MPI_REAL, right, right,
                 s->toroidal_comm, MPI_STATUS_IGNORE);

    if (s->mi + recv_left_count + recv_right_count > s->mimax) {
      gtc_die(s, "need bigger particle array in shifti");
    }
    const int first_new = s->mi;
    for (int i = 0; i < recv_left_count; i++) {
      unpack_particle(s, s->mi, &recv_left[(size_t)i * PACK_FIELDS]);
      s->mi++;
    }
    for (int i = 0; i < recv_right_count; i++) {
      unpack_particle(s, s->mi, &recv_right[(size_t)i * PACK_FIELDS]);
      s->mi++;
    }
    m0 = first_new;

    free(move);
    free(ileft);
    free(iright);
    free(send_right);
    free(send_left);
    free(recv_left);
    free(recv_right);
  }
  gtc_die(s, "endless particle sorting loop");
}
