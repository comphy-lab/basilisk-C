/**
  What is the smallest cells size (CS) in the simulation? 
  This will be the characteristic size with which I model
  the presence of my microswimmer.*/
#define ZITA (L0/(N))

#define radius sqrt(sq(PS(x,p().x+Shift.x))+sq(PS(y,p().y+Shift.y)))
#define SmoothedProfile(RADIUS) 0.5*(tanh((RADIUS-radius)/ZITA)+1.)

/**
	Define the scalar _continous_ fields that will account
	for the presence of body, flagella ... to which
	variable viscosity or variable forces will be applied.*/
/**
	*SP*, smoothed profile, will tend to be 1 inside the particle
				and zero elsewhere. The integral on the volume should tend
				to be equal to the area covered by the particle.
	*SPN*, smoothed profile NORMALIZED, so that the integral is equal
				to one, this is required in case I want to apply a given
				force to the body.*/
scalar SP0[];  
scalar SPN0[]; 

void compute_SP(Particles p){
	foreach()
		SPN0[] = 0.;
	foreach()
		SP0[] = 0.;

	foreach_particle_in(p){
		coord Shift = {0.,0.};
/**
	### Compute SPN */
/**
	Determine the _area_ covered by each blob,
	so to normalize all quantities wrt the present configuration.*/
		double integral = 0.;
		foreach(reduction(+:integral)) // should work in parallel as well...
			integral += SmoothedProfile(p().r)*sq(Delta);
/**
	Now I can normalize, so that the integral of blob for the
	present particle is equal to 1 .*/
		foreach()
			SPN0[] += SmoothedProfile(p().r)*sq(Delta)/integral;
/**
	### Compute SP, must be normalized so that the max is equal to 1.*/
  	double maxi = - 1e100;
  	foreach (reduction(max:maxi))
  	  if (fabs(SmoothedProfile(p().r)) > maxi)
  	    maxi = fabs(SmoothedProfile(p().r));
		foreach()
			SP0[] += SmoothedProfile(p().r)/maxi;
	}
}


/**
	# CORE of Fluid Particle Dynamics.
	i.e. the computation of variable viscosity. */
extern face vector muv;

void compute_variable_viscosity()
{
	foreach_face()
		muv.x[] = fm.x[] + ETA*face_value(SP0,0);
}
