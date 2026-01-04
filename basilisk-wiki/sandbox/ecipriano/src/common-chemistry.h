/**
# Chemistry structure

Empty abstract class for chemistry interfaces.
*/

#define CHEMISTRY 1

typedef struct {
  int NS;
  double f;
  double rho;
  double cp;
  double P;
  double T;
  double * sources;
} UserDataODE;

typedef void(*ode_function)(const double * y, // initial conditions
                            const double t,   // integration time step
                            double * dy,      // dydt equation
                            void * args);     // additional arguments


void stiff_ode_solver_dummy (ode_function batch,
    unsigned int neq, double dt, double * y, void * args)
{
  return;
}

void (* stiff_ode_solver)(ode_function batch,
    unsigned int, double, double *, void *) = stiff_ode_solver_dummy;

