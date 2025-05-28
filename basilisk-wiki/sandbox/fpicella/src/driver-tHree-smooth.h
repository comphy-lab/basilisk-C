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
