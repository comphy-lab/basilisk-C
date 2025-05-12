/**
This sample code shows a bug with the time marching event. Time
marching events sometimes incomplete at the end of time when events
with the step-marching and time-marching are defined and set to finish
at the same time. For instance below, this test case finishes without
running a sequence at t_end (t=1e-3), but a sequence at i_end (i=10)
is done. However, when i_end and t_end are set to 5 and 5e-4,
respectively, a sequence at t_end (t=5e-4) is correctly done. */

#include "grid/cartesian.h"
#include "run.h"

int i_end = 10;
double t_end = 10e-4;

int main(){
  init_grid(16);
  L0 = 1.0;
  dt = 1e-4;
  run(); 
}

/** Step-marching event */
event ievent (i = 0; i++; i <= i_end){
  dt = (double)dtnext(dt);
  printf("ievent i = %d t = %5.3e\n", i, t);
}

/**
Time-marching event

SP: this is not really a bug. Due to round-off errors, the test *t
*<= t_end* dos not produce the intended result. The solution is to
*remove the test or add an epsilon i.e. *t <= t_end + 1e-5*. */

event tevnet (t = 0; t += 1e-4 /*t <= t_end*/ ){
  printf("tevent i = %d t = %5.3e\n", i, t);
}
