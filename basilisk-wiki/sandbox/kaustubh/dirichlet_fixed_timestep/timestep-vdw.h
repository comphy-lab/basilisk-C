/**
# Timestep for the vdw mac-cormack 
*/

#ifndef growthTimestep // the user can define its own growth factor for the timestep relative to the previous one
#define growthTimestep 0.1
#endif

double timestep (vector u, vector q, double dtmax){
  static double previous = 0.;
	dtmax /= CFL;
	double dtmax2 =HUGE;
	foreach(reduction(min:dtmax) reduction(min:dtmax2))
	foreach_dimension(){
		if (u.x[] != 0.) {
			double dt = Delta/fabs(u.x[]);
			if (dt < dtmax) dtmax = dt;
		}
		if (q.x[] != 0.) {
			double dt = Delta/fabs(q.x[]);
			if (dt < dtmax2) dtmax2 = dt;
		}
	}
	dtmax = min(dtmax, dtmax2);
	dtmax *=CFL;
	if (dtmax > previous)
		dtmax = (previous + growthTimestep*dtmax)/(1.+growthTimestep);
	previous = dtmax;
	return dtmax;
}