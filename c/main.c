#include "gtc.h"

#include <math.h>
#include <string.h>
#include <time.h>

static GtcReal wall_now(void) {
  return MPI_Wtime();
}

static GtcReal cpu_now(void) {
  return (GtcReal)clock() / (GtcReal)CLOCKS_PER_SEC;
}

static void write_file_exit(const GtcState *state) {
  FILE *fp = fopen("FileExit.dat", "w");
  if (!fp) return;
  const int irest = state->irest;
  fprintf(fp, "FileExit=%1d\n", 1);
  fprintf(fp, "irest   =%5d\n", irest);
  fprintf(fp, "%s\n", ((irest + 1) % 2) == 0 ? "restart_dir1" : "restart_dir2");
  fclose(fp);
}

int main(int argc, char **argv) {
  (void)gtc_mpi_init(&argc, &argv);
  GtcState state;
  memset(&state, 0, sizeof(state));
  GtcReal cpu_time[8] = {0.0};
  GtcReal wall_time[8] = {0.0};
  const GtcReal cpu_start = cpu_now();
  const GtcReal wall_start = wall_now();
  GtcReal cpu0 = cpu_start;
  GtcReal wall0 = wall_start;

  setup(&state);
  load(&state);
  chargei(&state);
  cpu_time[6] += cpu_now() - cpu0;
  wall_time[6] += wall_now() - wall0;
  cpu0 = cpu_now();
  wall0 = wall_now();
  const GtcReal loop_start = wall0;

  for (state.istep = 1; state.istep <= state.p.mstep; state.istep++) {
    for (state.irk = 1; state.irk <= 2; state.irk++) {
      const int idiag = ((state.irk + 1) % 2) + (state.istep % state.p.ndiag);

      smooth(&state, 3);
      cpu_time[4] += cpu_now() - cpu0;
      wall_time[4] += wall_now() - wall0;
      cpu0 = cpu_now();
      wall0 = wall_now();

      field(&state);
      cpu_time[5] += cpu_now() - cpu0;
      wall_time[5] += wall_now() - wall0;
      cpu0 = cpu_now();
      wall0 = wall_now();

      pushi(&state);
      cpu_time[0] += cpu_now() - cpu0;
      wall_time[0] += wall_now() - wall0;
      cpu0 = cpu_now();
      wall0 = wall_now();

      shifti(&state);
      cpu_time[1] += cpu_now() - cpu0;
      wall_time[1] += wall_now() - wall0;
      cpu0 = cpu_now();
      wall0 = wall_now();

      chargei(&state);
      cpu_time[2] += cpu_now() - cpu0;
      wall_time[2] += wall_now() - wall0;
      cpu0 = cpu_now();
      wall0 = wall_now();

      smooth(&state, 0);
      cpu_time[4] += cpu_now() - cpu0;
      wall_time[4] += wall_now() - wall0;
      cpu0 = cpu_now();
      wall0 = wall_now();

      poisson(&state, 0);
      cpu_time[3] += cpu_now() - cpu0;
      wall_time[3] += wall_now() - wall0;
      cpu0 = cpu_now();
      wall0 = wall_now();

      if (idiag == 0) {
        diagnosis(&state);
      }
    }

    if (state.isnap > 0 && state.istep % state.isnap == 0) snapshot(&state);
  }

  cpu_time[7] = cpu_now() - cpu_start;
  wall_time[7] = wall_now() - wall_start;
  const GtcReal loop_time = wall_now() - loop_start;
  if (state.rank == 0) {
    FILE *out = gtc_stdout_open(&state, "a");
    if (out) {
      fprintf(out, " CPU TIME USAGE (in SEC):\n");
      fprintf(out, " pusher     shift      charge     poisson    smooth     field      load       total\n");
      for (int i = 0; i < 8; i++) fprintf(out, "%10.3e%c", cpu_time[i], i == 7 ? '\n' : ' ');
      fprintf(out, " WALL CLOCK TIMES (in SEC):\n");
      fprintf(out, " pusher     shift      charge     poisson    smooth     field      load       total\n");
      for (int i = 0; i < 8; i++) fprintf(out, "%10.3e%c", wall_time[i], i == 7 ? '\n' : ' ');
      fprintf(out, "MAIN LOOP TIME(SEC):%12.3f\n", loop_time);
      fprintf(out, "TOTAL CPU TIME USAGE (SEC):%12.3f\n", cpu_time[7]);
      fprintf(out, "TOTAL WALL CLOCK TIME(SEC):%12.3f\n", wall_time[7]);
      gtc_stdout_close(&state, out);
    }
    state.file_exit = 1;
    write_file_exit(&state);
  }

  gtc_free(&state);
  MPI_Finalize();
  return 0;
}
