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

#define ShiftxCenter      p().x-0.*BETA*cosTheta+0.*ALPHA*sinTheta
#define ShiftyCenter      p().y-0.*BETA*sinTheta-0.*ALPHA*cosTheta

#define ShiftxTrans       p().x-1.*BETA*cosTheta+1.*ALPHA*sinTheta
#define ShiftyTrans       p().y-1.*BETA*sinTheta-1.*ALPHA*cosTheta

#define ShiftxCis         p().x-1.*BETA*cosTheta-1.*ALPHA*sinTheta
#define ShiftyCis         p().y-1.*BETA*sinTheta+1.*ALPHA*cosTheta


/**
	Next step, improve the _blob_ model.*/




/**
	A field for visualisation purposes only. */
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
				bodyForce.x[] += SP(p().r)*(p().B.x+cos(p().theta)*p().Thrust);
				bodyForce.y[] += SP(p().r)*(p().B.y+sin(p().theta)*p().Thrust);
		}
	// Trans
		Shift.x = ShiftxTrans; Shift.y = ShiftyTrans;
		foreach(){
				bodyForce.x[] += SP(p().r)*(-cos(p().theta)*p().Thrust/2.);
				bodyForce.y[] += SP(p().r)*(-sin(p().theta)*p().Thrust/2.);
		}
	// Cis
		Shift.x = ShiftxCis; Shift.y = ShiftyCis;
		foreach(){
				bodyForce.x[] += SP(p().r)*(-cos(p().theta)*p().Thrust/2.);
				bodyForce.y[] += SP(p().r)*(-sin(p().theta)*p().Thrust/2.);
		}
	}
	foreach_face()
		av.x[] = face_value(bodyForce.x,0);
}
