// Generic predictor/corrector time-integration

#include "utils.h"

// Required from solver
// fields updated by time-integration
extern scalar * evolving;
// how to compute updates
double (* update) (scalar * evolving, scalar * updates, double dtmax) = NULL;

// User-provided functions
// gradient
double (* gradient)  (double, double, double) = minmod2;

// the timestep
double dt = 0.;

trace
static void advance_generic (scalar * output, scalar * input, scalar * updates,
			     double dt)
{
  if (input != output)
    trash (output);
  foreach() {
    scalar o, i, u;
    for (o,i,u in output,input,updates)
      o[] = i[] + dt*u[];
  }
}

static void (* advance) (scalar * output, scalar * input, scalar * updates,
			 double dt) = advance_generic;

event defaults (i = 0)
{
  // limiting
  for (scalar s in all)
    s.gradient = gradient;
  
  display ("box();");
}

trace
void run()
{
  t = 0., iter = 0; dimensional (dt == DT);
  init_grid (N);

  // main loop
  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {
    // list of updates
    scalar * updates = list_clone (evolving);
    dt = dtnext (update (evolving, updates, DT));
    if (gradient != zero) {
      /* 2nd-order time-integration */
      scalar * predictor = list_clone (evolving);
      /* predictor */
      advance (predictor, evolving, updates, dt/2.);
      /* corrector */
      update (predictor, updates, dt);
      delete (predictor);
      free (predictor);
    }
    advance (evolving, evolving, updates, dt);
    delete (updates);
    free (updates);
    update_perf();
    iter = inext, t = tnext;
  }
  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
}
