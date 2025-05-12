//
#include <math.h>
#include <stdio.h>
#include <time.h>
//
#include "my_log.h"

void progress_print(double time) {
  int hours_ = (int)time;
  double minutesRemainder = (time - hours_) * 60.0;
  int minutes_ = (int)minutesRemainder;
  if (hours_ == 0) {
    fprintf(stdout, "%d minutes\n", minutes_);
  } else {
    fprintf(stdout, "%d hours & %d minutes\n", hours_, minutes_);
  }
  fflush(stdout);
}

clock_t start_time;
static int progress_last = 0;
double time_last = 0.;
int WARN = 0;
double time_per_percent_last = 1.0e6;
double t_prev_progress = 0.;

static void progress_check(int i, double t, double tend) {
  if (i == 0) {
    progress_last = 0;
    time_last = 0.;
    WARN = 0;
    time_per_percent_last = 1.0e6;
    t_prev_progress = 0.;
#if _OPENMP
    start_time = omp_get_wtime();
#else
    start_time = clock();
#endif
  }
  double dt1 = t - t_prev_progress;
#if _OPENMP
  double time = (double)(omp_get_wtime() - start_time) / 3600.0;
#else
  double time = (double)(clock() - start_time) / CLOCKS_PER_SEC / 3600.0;
#endif
  int triggered = 0;
  int triggered_global = 0;

  /// time might not be exactly the same for all processes, 
  // hence the use of the MPI reduced 
  double time_per_percent = time - time_last;
  if (time_per_percent > 2.0 * time_per_percent_last && WARN == 1) {
    triggered = 1;
  }

#if _MPI
  MPI_Allreduce (&triggered, &triggered_global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

  triggered_global = max(triggered_global, triggered);

  if (triggered_global==1){
    fprintf(stderr, "triggered %d triggered_global %d pid %d\n",triggered, triggered_global, pid());
#ifdef SLOW_DUMP
    dump(file = "dump_slow");
#endif
    WARN = 0;
    fprintf(stderr, "pid %d\n",pid());
    if (pid()==0){
      fprintf(stdout, "Warning: solver has slowed down a lot\n");
      fflush(stdout);
    }
  }

  int progress = (int)(floor(100.0 * t / (tend - dt1)));
  if (t > (tend - 1.5 * dt1)) {
    progress = 100;
  }
  if (progress != progress_last) {
#if _OPENMP
    double time = (double)(omp_get_wtime() - start_time) / 3600.0;
#else
    double time = (double)(clock() - start_time) / CLOCKS_PER_SEC / 3600.0;
#endif
    double time_est = ((time - time_last) / 1.0 * (100.0 - progress));
    fprintf(stdout, "%d %%, ETR ", progress);
    fflush(stdout);
    progress_print(time_est);
    fflush(stdout);

    progress_last = progress;
    time_last = time;
    if (progress == 100) {
      fprintf(stdout, "Total duration: ");
      fflush(stdout);
      progress_print(time);
    }
    time_per_percent_last = time_per_percent;
    WARN = 1;
  }
  fflush(stdout);
  t_prev_progress = t;
}

