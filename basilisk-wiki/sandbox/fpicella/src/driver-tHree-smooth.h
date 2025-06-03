/**
	# 3-forces microswimmer model
	Build up a forcing field...*/

face vector av[];
/*
  What is the smallest cells size (CS) in the simulation? 
  This will be the characteristic size with which I model
  the presence of my microswimmer.*/
#define CS (L0/(N))
// This is right if I have a fixed mesh, but should check when using
// variable refinement (adaptive...)

/**
	Blob on which I compute the force on the center of the microswimmer. 
	It is a spherical blob of diameter one mesh cell.*/
#define blob(xp,yp,zp) (sq (PS(x,xp)) + sq (PS(y,yp)) + sq(PS(z,zp)) - sq (2.*CS))

#define radius sqrt(sq(PS(x,p().x+Shift.x))+sq(PS(y,p().y+Shift.y)))
/**
	Smoothed Profile, SP .*/
#define SP(RADIUS) 0.5*(tanh((RADIUS-radius)/CS)+1.)


#define sinTheta     sin(p().theta)
#define cosTheta     cos(p().theta)

#define ALPHA alpha*p().r
#define BETA   beta*p().r

#define ShiftxCenter   -0.*BETA*cosTheta+0.*ALPHA*sinTheta
#define ShiftyCenter   -0.*BETA*sinTheta-0.*ALPHA*cosTheta

#define ShiftxTrans    -1.*BETA*cosTheta+1.*ALPHA*sinTheta
#define ShiftyTrans    -1.*BETA*sinTheta-1.*ALPHA*cosTheta

#define ShiftxCis      -1.*BETA*cosTheta-1.*ALPHA*sinTheta
#define ShiftyCis      -1.*BETA*sinTheta+1.*ALPHA*cosTheta

vector bodyForce[];
void compute_microswimmer_forcing_smooth(Particles p)
{
	foreach()
		foreach_dimension()
			bodyForce.x[] = 0.;
	foreach_particle_in(p){
		coord Shift = {0.,0.};
	// Center
		Shift.x = ShiftxCenter; Shift.y = ShiftyCenter;
		foreach(){
				bodyForce.x[] += SP(p().r)*(p().B.x+cos(p().theta)*p().Thrust)/(p().r*p().r*M_PI);
				bodyForce.y[] += SP(p().r)*(p().B.y+sin(p().theta)*p().Thrust)/(p().r*p().r*M_PI);
		}
	// Trans
		Shift.x = ShiftxTrans; Shift.y = ShiftyTrans;
		foreach(){
				bodyForce.x[] += SP(p().r)*(-cos(p().theta)*p().Thrust/2.)/(p().r*p().r*M_PI);
				bodyForce.y[] += SP(p().r)*(-sin(p().theta)*p().Thrust/2.)/(p().r*p().r*M_PI);
		}
	// Cis
		Shift.x = ShiftxCis; Shift.y = ShiftyCis;
		foreach(){
				bodyForce.x[] += SP(p().r)*(-cos(p().theta)*p().Thrust/2.)/(p().r*p().r*M_PI);
				bodyForce.y[] += SP(p().r)*(-sin(p().theta)*p().Thrust/2.)/(p().r*p().r*M_PI);
		}
	}
	foreach_face()
		av.x[] = face_value(bodyForce.x,0);
}

scalar sp[];
void compute_sp(Particles p){
  foreach()
    sp[] = 0.;
	foreach_particle_in(p){
		coord Shift = {0.,0.};
		foreach()
				sp[] += SP(p().r);
	}
}

/** Include _steric_ interactions.
To do so, use some potential to tune body force
on the particle when close to boundaries or other particles.*/
/** For the moment, I've got a purely repulsive potential, identical
to the 12-term in Lennard-Jones.*/
/** Forces are applied directly at the center of the body, as if 
they where some sort of _body-force_.*/

void compute_repulsion(Particles p, double repulsion_distance, bool top, bool bottom, bool right, bool left, double repulsion_strength)
{
	foreach()
		foreach_dimension()
			bodyForce.x[] = 0.;
	foreach_particle_in(p){
		
		coord Shift = {0.,0.}; // legacy from previous functions...keep it to zero
													 // if you want the force to be applied around the 
													 // microswimmer body.
		coord Intensity = {0.,0.}; // intensity of the force to apply.
		// Intensity will be built upon body-wall and body-body interactions.
		/** What is the distance from which I want the repulsion force to kick in?*/
		double sigma = repulsion_distance*p().r;
		/** Set up a variable for the particle-surface distance, for convenience...*/
		double distance = 0.;
	/** Compute the intensity of repulsion, between the particle and the bottom wall.*/
		if(bottom==true){
			distance = L0/2.+y;
			Intensity.y += pow(sigma/distance,repulsion_strength);//-pow(sigma/distance,6); 
		}
	/** Compute the intensity of repulsion, between the particle and the top wall.*/
		if(top==true){
			distance = L0/2.-y;
			Intensity.y -= pow(sigma/distance,repulsion_strength);//-pow(sigma/distance,6); 
		}
	/** Compute the intensity of repulsion, between the particle and the right wall.*/
		if(right==true){
			distance = L0/2.-x;
			Intensity.x -= pow(sigma/distance,repulsion_strength);//-pow(sigma/distance,6); 
		}
	/** Compute the intensity of repulsion, between the particle and the bottom wall.*/
		if(left==true){
			distance = L0/2.+x;
			Intensity.x += pow(sigma/distance,repulsion_strength);//-pow(sigma/distance,6); 
		}
	/** Compute particle-particle repelling force. */
	/** This clearly is not the most efficient way to do it:
			I'm computing twice for the same k-j particles.
			Still, it is super simple and it works. For the moment
			I keep it like this.
			Honestly, I don't think that for my purposes the computational
			overhead will be even noticeable.
			*/
  	for (int _k_particle = 0; _k_particle < pn[_l_particle]; _k_particle++) {
			if(_k_particle != _j_particle){
				coord DISTANCE;
				foreach_dimension()
					DISTANCE.x = pl[_l_particle][_k_particle].x - pl[_l_particle][_j_particle].x;
				distance = sqrt(sq(DISTANCE.x)+sq(DISTANCE.y))-sigma;
				double INTENSITY = pow(sigma/distance,12);//-pow(sigma/distance,6);
				double angle = atan2(DISTANCE.y/distance,DISTANCE.x/distance);
				//fprintf(stderr,"particle-particle repulsion %d %d %+6.5e %+6.5e %+6.5e \n",_j_particle,_k_particle,angle,distance,INTENSITY);
				Intensity.x -= INTENSITY*cos(angle);
				Intensity.y -= INTENSITY*sin(angle);
			
			}
		}
	// Center
		foreach(){
				bodyForce.x[] += Intensity.x*SP(p().r)/(p().r*p().r*M_PI);
				bodyForce.y[] += Intensity.y*SP(p().r)/(p().r*p().r*M_PI);
		}
	foreach_face()
		av.x[] += face_value(bodyForce.x,0);
	}
}

/** For plotting purposes only.*/
scalar bodyPlot[];
void compute_bodyPlot(Particles p){
	foreach()
		bodyPlot[] = 0.;
	foreach_particle_in(p){
		coord Shift = {0.,0.}; // legacy from previous functions...keep it to zero
													 // if you want the force to be applied around the 
													 // microswimmer body.
	// Center
		Shift.x = ShiftxCenter; Shift.y = ShiftyCenter;
		foreach()
				bodyPlot[] += SP(p().r);
	// Trans
		Shift.x = ShiftxTrans; Shift.y = ShiftyTrans;
		foreach()
				bodyPlot[] += SP(p().r);
	// Cis
		Shift.x = ShiftxCis; Shift.y = ShiftyCis;
		foreach()
				bodyPlot[] += SP(p().r);
	}
}
