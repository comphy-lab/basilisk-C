/**
# Generic predictor/corrector time-integration
*/
#include "utils.h"

// Maximum number of source terms
#define NSource 10

// Required from solver
// fields updated by time-integration
extern scalar * evolving;

// how to compute updates
double (* update) (scalar * evolving, scalar * updates, double dtmax) = NULL;

// User-provided functions
// gradient
double (* gradient)  (double, double, double) = minmod2;

double t = 0., dt = 0.;

static void advance_generic (scalar * output, scalar * input, scalar * updates,
			     double dt){
  if (input != output)
    trash (output);
  foreach() {
    scalar o, i, u;
    for (o,i,u in output,input,updates)
      o[] = i[] + dt*u[];
  }
  boundary (output);
}


static void (* advance) (scalar * output, scalar * input, scalar * updates,
			 double dt) = advance_generic;

// how to compute source terms
void (* updatesource[NSource] ) (scalar * evolving,
				scalar * sources,
				double dtmax,
				int numbersource);

// Fixing minor bug when a pointer is calling NULL
static void fnull(){
}

// Because of the additivity of source terms, we need to reset them at each step
static void resetsource(scalar * evolving,
				scalar * sources,
				double dtmax,
				int numbersource){
  scalar dsh = sources[0];
  vector dshu = { sources[1], sources[2] }; 
  foreach(){
      dsh[] = 0;
      foreach_dimension() dshu.x[] = 0;
    }
  numbersource++;
  updatesource[numbersource](evolving, sources,
				dtmax, numbersource);
}

event defaults (i = 0)
{
  // limiting
  for (scalar s in all)
    s.gradient = gradient;

  // default values
  foreach()
    for (scalar s in all)
      s[] = 0.;
  boundary (all);
}


void run()
{
  t = 0.;
  init_grid (N);
  
  // main loop
  timer start = timer_start();
  int i = 0; long tnc = 0;
 
  while (events (i, t, true)) {
    // list of updates
    scalar * updates = list_clone (evolving);
    scalar * sources = list_clone (evolving);
    dt = dtnext (t, update (evolving, updates, DT));
    if (gradient != zero) {
      /* 2nd-order time-integration */
      scalar * predictor = list_clone (evolving);
      /* predictor */
      advance (predictor, evolving, updates, dt/2.);

      updatesource[0](predictor, sources, dt/2.,0);
      advance (predictor, predictor, sources, dt/2.);
      
      /* corrector */
      update (predictor, updates, dt);
      
      delete (predictor);
      free (predictor);
    }

    advance (evolving, evolving, updates, dt);

    updatesource[0] (evolving, sources, dt,0); 
    advance (evolving, evolving, sources, dt);
    
    delete (updates);
    free (updates);
    delete(sources);
    free(sources);
    
   long nc = 0;
    foreach (reduction(+:nc)) 
      nc++;
    tnc += nc;
    i++; t = tnext;
  }
  timer_print (start, i, tnc);

  free_grid();
}

// Overloading the first updatesource function with resetsource
int numbersource;
event initres(i = 0 ){
  numbersource = 0;
  updatesource[numbersource]=resetsource;
  numbersource++;
  updatesource[numbersource] = fnull;
}
