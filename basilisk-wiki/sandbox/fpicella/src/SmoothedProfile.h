/**
# Implementation of the Smoothed Profile method
 [Nakayana Yamamoto 2005](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.71.036707)
*/
scalar csFP[];
face vector fsFP[];
face vector ufFP[]; // target fluid velocity under the particle
face vector ufDIFF[];// difference between the target and the actual value
scalar divDIFF[]; // divergence of ufDIFF
scalar pFP[]; // pressure correction owing to the presence of the particle
face vector fFP[]; // force required to ensure RBM

face vector av[]; // acceleration face vector, will contain
									// - fFP, the forcing of the fluid particle
									// - whatever I want

/**
### Define the shape of a cylinder */
#define circleShape(x,y,r) -(sq (x) + sq (y) - sq (r))
#define r sqrt(sq(xPWP)+sq(yPWP))
#define SmoothedProfile(r,RADIUS,ZITA) 0.5*(tanh((RADIUS-r)/ZITA)+1.)


void solidParticles()
{
	scalar csLOCAL[];
	face vector fsLOCAL[];

	foreach()
		csFP[] = 0.;
	foreach_face()
		fsFP.x[] = 0.;
	foreach_face()
		ufFP.x[] = 0.;

	foreach_particle(){
/**
	use of PWP's macros from periodic-shift-treatment
	to aboid clumsy periodicity definitions...
*/
		//fraction(csLOCAL,circleShape(xPWP,yPWP,RADIUS)); 
		//solid(csLOCAL,fsLOCAL,circleShape(xPWP,yPWP,RADIUS)); 
		foreach()
			csLOCAL[] = SmoothedProfile(r,RADIUS,ZITA);
		foreach_face()
			fsLOCAL.x[] = face_value(csLOCAL,0);
/**
Compute csFP[]*/
		foreach()
			csFP[]   += csLOCAL[];
		foreach_face()
			fsFP.x[] += fsLOCAL.x[];
/**
Compute the fluid's velocity, covered by a RBM particle...*/
		coord sgn ={-1.,1.};
		//foreach(){
		//	coord r = {x,y,z}; // coordinate x,y,z are not permuted with foreach_dimension()
		//	foreach_dimension()
		//		uFP.x[] = ((p().u.x) + sgn.x*p().w.x*PS(r.y,p().y))*csLOCAL[];
		//}
		foreach_face()
			ufFP.x[] = ((p().u.x) + sgn.x*p().w.x*PS(y,p().y))*fsLOCAL.x[];
	}
}
void compute_divDIFF()
{
	foreach_face()
		ufDIFF.x[] = (ufFP.x[]-uf.x[])*fsFP.x[];
// straight from src/poisson.h
  foreach() {
    divDIFF[] = 0.;
    foreach_dimension()
      divDIFF[] += ufDIFF.x[1] - ufDIFF.x[];
    divDIFF[] /= dt*Delta;
	}
}



event acceleration(i++){
	solidParticles();
	compute_divDIFF();
	poisson(pFP,divDIFF,fs);
	/**
	compute the force (acceleration) required to enforce RBM */
	foreach_face()
		fFP.x[] = (ufDIFF.x[]/dt - face_gradient_x (pFP, 0))*fs.x[];
	
  foreach_face()
    av.x[]=fFP.x[]*force_tune;
}

