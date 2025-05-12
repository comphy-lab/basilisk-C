//#define PERIODIC 1. // so to make it faster on some loops?
//#define RADIUS 1.0
//#define ZITA   0.1 // width of diffused interface
//#define NPARTICLES 1
//#define ALPHA 1.0 // how much off-axis sideways are the forces?
//#define BETA  1.0 // how much puller/pusher is the particle?
//#define STRENGTH 10. // strenght of the swimmer
//#define B 1000.      // reorientation owing to gravity. Large value = set it off.

scalar omega[];   // field that will contain vorticity...

/*
Fields that will contain the phase fractions of all microswimmers */
scalar csFP[];
face vector fsFP[];
scalar bodyforce[]; // a nice scalar field, to "visualize" how forces are applied on the field.
face vector forceFP[];

/*
Some practical shortcuts
*/
#define circle(px,py,pr) (sq (PS(x,px)) + sq (PS(y,py)) - sq (pr))

#define sinTheta     sin(p().theta)
#define cosTheta     cos(p().theta)

#define xCenter      p().x+0.*ALPHA*cosTheta+0.*BETA*sinTheta
#define yCenter      p().y+0.*ALPHA*sinTheta-0.*BETA*cosTheta

#define xTrans       p().x+1.*ALPHA*cosTheta+1.*BETA*sinTheta
#define yTrans       p().y+1.*ALPHA*sinTheta-1.*BETA*cosTheta

#define xCis         p().x-1.*ALPHA*cosTheta+1.*BETA*sinTheta
#define yCis         p().y-1.*ALPHA*sinTheta-1.*BETA*cosTheta

#define ED(px,py) sqrt(sq(x-px) + sq(y-py)) // Euclidean Distance

#define EPSILON 0.5

#define RS(px,py) (15.*pow(EPSILON,4))/pow(8.*M_PI*(sq(ED(px,py))+sq(EPSILON)),7./2.)

/*
### Microswimmer reorientation equation
For the moment, simple gyrotaxis owing to gravity.
But could be anything!
*/
// BEGIN -- DEFINE RE-ORIENTATION EQUATION (i.e. gyrotaxis owing to rheology, gravity, light...whatever!) // FP 20250308 
#define OMEGA_EQUATION() -sin(p().theta)/(2.*B) + interpolate_linear(locate (x, y, z), omega, x, y, z)/2.  // FP 20250308 



static FILE *singleParticleFile[NPARTICLES] = {NULL}; // for the moment only one set of particles...?

#include "fpicella/src/tracer-particles-FP.h"

int alternating_series(int n) {
    if (n == 0) return 0; // Base case
    int value = (n + 1) / 2;
    return (n % 2 == 1) ? value : -value;
}


void locate_particle_in_zero(){
  int _l_particle = 0;
    while (pn[_l_particle] != terminate_int) {
      for (int _j_particle = 0; _j_particle < pn[_l_particle]; _j_particle++) {
        p().x = 3.*RADIUS*alternating_series(_j_particle);
        p().y = L0/2.*RADIUS*alternating_series(_j_particle);
				p().z = 0.; 
				p().theta=+_j_particle*M_PI*1.;
/*
First, barebones version, all particles have the same radius...
...but provided the numerical method improved, I could play around with 
different shapes, size...
For the moment, I stick with cylinders and spheres
*/
				p().r = RADIUS;
      }
      _l_particle++;
    }
}

// // // // // // // // //
Particles ParticleList; // // Call the particles!
// // // // // // // // //


// INITIALIZE PARTICLE LOCATION
event init (t = 0){
  ParticleList = init_tp_circle(NPARTICLES);
  locate_particle_in_zero();
  //fprintf(stderr,"Number of particles in list ParticleList: %lu\n", pn[ParticleList]); // ok, this one tells me how many particle's I've got here.
}

// Particle output, super-compact (and working in serial!) way :)
// FP, 20250212 11h37
event output_particle_initialize(i=0){ // first iteration, define name and open files...
  foreach_particle(){
    char filename[100];
    sprintf(filename, "particle_%03d.dat", _j_particle);
    //fprintf(stderr,"%d %s \n",_j_particle,filename); 
    singleParticleFile[_j_particle] = fopen(filename,"w");
  }
}
#define FORMAT_SPEC_7 ("%+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n") // to avoid writing every time it...
event output_particle(i++){
  foreach_particle(){
    //fprintf(stderr,"%d \n",_j_particle);
    fprintf(singleParticleFile[_j_particle],FORMAT_SPEC_7,t,p().x,p().y,p().theta,p().u.x,p().u.y,p().omega);
    //fprintf(singleParticleFile[_j_particle],FORMAT_SPEC_7,t,p().x,p().y,0.000000000,p().u.x,p().u.y,0.0000000);

    fflush(singleParticleFile[_j_particle]);
  }
}

/*
### Compute the _variable_ viscosity field associated to the presence of each microswimmer
*/
#define circle(px,py,pr) (sq (PS(x,px)) + sq (PS(y,py)) - sq (pr))
event compute_variable_viscosity_field(i=-1){ // I'll explicitly call this from the main
	foreach()
		csFP[] = 1.; // reset field to one;
	scalar csLOCAL[];
	foreach_particle(){
		fraction(csLOCAL,circle(xCenter,yCenter,p().r));
		foreach()
			csFP[] *= csLOCAL[]; // intersection...
	}
	foreach_face()
		fsFP.x[] = (face_value(csFP,0) + face_value(csFP,1))*0.5; // "naive" face interpolation

  foreach_face()
#ifndef EMBED
    muv.x[] = (nu)*1.0000+(1.-fsFP.x[])*nuRatio*nu - 1.0000*(1.-fsFP.x[])*nu;
#else
    muv.x[] = (nu)*fs.x[]+(1.-fsFP.x[])*nuRatio*nu - fs.x[]*(1.-fsFP.x[])*nu;
#endif
}

/*
### Compute the body-forcing fields.
From a theoretical point of view, I should have a point force.
Which can not be possible, provided I'm using a finite volume method.
The approach I'll apply instead, will be applying "regularized forces"
on a non-zero volume.
For the moment, those forces are regularized on a volume that is identical
to the microswimmer's body.
*/
event compute_propulsion_forces(i=-1){
	scalar csLOCAL[];
	face vector fsLOCAL[];

	foreach(){
		bodyforce[] = 0.;
	}
	foreach_face()
		forceFP.x[] = 0.;
	foreach_particle(){
		p().areaCenter = 0.;
		p().areaTrans  = 0.;
		p().areaCis    = 0.;
/*
### Regularized Stokeslet
Idea, replace point force with a blob...
*/
	
///*
//	Blob associated to the presence of the microswimmers's body.
//*/
//		foreach(){
//			bodyforce[] += RS(xCenter,yCenter);
//			forceFP.y[] += RS(xCenter,yCenter)*SEDIMENTATION;
//		}
		solid(csLOCAL,fsLOCAL,circle(xCenter,yCenter,p().r));
		foreach(){
			bodyforce[] += -(1.-csLOCAL[]);
			p().areaCenter += (1.-csLOCAL[])*sq(Delta); // MUST BE MODIFIED IN 3D!
		}
		foreach_face(y)
			forceFP.y[] -= (1.-fsLOCAL.y[])*(+SEDIMENTATION)*sq(Delta)/p().areaCenter; 

		double TOTAL_FORCE_CHECK = 0.;
		foreach_face(y)
			TOTAL_FORCE_CHECK += forceFP.y[];
		fprintf(stderr,"TOTAL_FORCE_CHECK %+6.5e\n",TOTAL_FORCE_CHECK);
			
	}

///*
//	Blob associated to the presence of the microswimmers's body.
//*/
//		fraction(bodyforceLOCAL,circle(xCenter,yCenter,p().r));
//		foreach(){
//			bodyforce[] += -(1.-bodyforceLOCAL[]);
//			p().areaCenter += (1.-bodyforceLOCAL[])*sq(Delta); // MUST BE MODIFIED IN 3D!
//		}
//		// Since I'm still here, I can apply the forces associated to the body of the microswimmer
//		foreach(){
//			forceFP.x[] -= (1.-bodyforceLOCAL[])*sq(Delta)*(+sinTheta*2.*PROPULSION)                 / p().areaCenter; 
//			forceFP.y[] -= (1.-bodyforceLOCAL[])*sq(Delta)*(-cosTheta*2.*PROPULSION + SEDIMENTATION) / p().areaCenter;
//		}
//
///*
//	Blob associated to the presence of the microswimmers's trans flagella
//*/
//		fraction(bodyforceLOCAL,circle(xTrans,yTrans,p().r));
//		foreach(){
//			bodyforce[] += +(1.-bodyforceLOCAL[]);
//			p().areaTrans += (1.-bodyforceLOCAL[])*sq(Delta); // MUST BE MODIFIED IN 3D!
//		}
//		// Since I'm still here, I can apply the forces associated to the Trans flagella of the microswimmer
//		foreach(){
//			forceFP.x[] += (1.-bodyforceLOCAL[])*sq(Delta)*(+sinTheta*1.*PROPULSION) / p().areaTrans; 
//			forceFP.y[] += (1.-bodyforceLOCAL[])*sq(Delta)*(-cosTheta*1.*PROPULSION) / p().areaTrans;
//		}
///*
//	Blob associated to the presence of the microswimmers's trans flagella
//*/
//		fraction(bodyforceLOCAL,circle(xCis,yCis,p().r));
//		foreach(){
//			bodyforce[] += +(1.-bodyforceLOCAL[]);
//			p().areaCis += (1.-bodyforceLOCAL[])*sq(Delta); // MUST BE MODIFIED IN 3D!
//		}
//		// Since I'm still here, I can apply the forces associated to the Cis flagella of the microswimmer
//		foreach(){
//			forceFP.x[] += (1.-bodyforceLOCAL[])*sq(Delta)*(+sinTheta*1.*PROPULSION) / p().areaCis; 
//			forceFP.y[] += (1.-bodyforceLOCAL[])*sq(Delta)*(-cosTheta*1.*PROPULSION) / p().areaCis;
//		}
//	}

}
/*
#Apply acceleration */
event acceleration (i++) {
	foreach_face()
		av.x[] = forceFP.x[];
	a = av;
}

// // // // On top of that now, compute the FORCING FIELD, associated to each particle.
// // // face vector forcesFP[]; // forces on fluid particle...
// // // scalar forceMAG[];
// // // scalar totalForce[]; // for plotting purposes only
// // // 
// // // void compute_forcesFP(face vector forcesFP){
// // // //#define RCENTER double rCenter = sqrt(sq(x-(plx+xp+0.000))+sq(y-(ply+yp+0.000)))
// // // #define sinTheta     sin(p().theta)
// // // #define cosTheta     cos(p().theta)
// // // #define xCenter      plx+xp+0.*ALPHA*cosTheta+0.*BETA*sinTheta
// // // #define yCenter      ply+yp+0.*ALPHA*sinTheta-0.*BETA*cosTheta
// // // #define rCenter      sqrt(sq(x-xCenter)+sq(y-yCenter))
// // // 
// // // #define xTrans       plx+xp+1.*ALPHA*cosTheta+1.*BETA*sinTheta
// // // #define yTrans       ply+yp+1.*ALPHA*sinTheta-1.*BETA*cosTheta
// // // #define rTrans       sqrt(sq(x-xTrans )+sq(y-yTrans))
// // // 
// // // #define xCis         plx+xp-1.*ALPHA*cosTheta+1.*BETA*sinTheta
// // // #define yCis         ply+yp-1.*ALPHA*sinTheta-1.*BETA*cosTheta
// // // #define rCis         sqrt(sq(x-xCis )+sq(y-yCis))
// // // // // // Interator...over all particles...barebones.
// // //   foreach_face(){
// // //     forcesFP.x[] = 0.;
// // //   }
// // //   foreach()
// // //     totalForce[] = 0.;
// // //   int _l_particle = 0;
// // //     while (pn[_l_particle] != terminate_int) {
// // //       for (int _j_particle = 0; _j_particle < pn[_l_particle]; _j_particle++) {
// // //         double plx = p().x; double ply = p().y; //double plz = p().z;
// // //         double forceCellCounter = 0.;
// // // 				double forceCellTotalArea = 0.;
// // //         foreach()
// // //           forceMAG[] = 0.;
// // // 
// // //         for (double xp =  -L0; xp <= L0; xp += L0)
// // //           for (double yp = -L0; yp <= L0; yp += L0)
// // // //            for (double zp = -L0; zp <= L0; zp += L0) // for the moment, 2D ONLY 202503
// // //               foreach(){ // aligned in y direction...
// // // 								forceMAG[] += 1.0*(0.5*(1.0+tanh((RADIUS-rCenter)/ZITA)));
// // // 				  forceCellCounter += 1.0*(0.5*(1.0+tanh((RADIUS-rCenter)/ZITA)));
// // // 								forceMAG[] -= 0.5*(0.5*(1.0+tanh((RADIUS-rTrans )/ZITA)));
// // // 								forceMAG[] -= 0.5*(0.5*(1.0+tanh((RADIUS-rCis   )/ZITA)));
// // // 				  forceCellTotalArea += 1.0*(0.5*(1.0+tanh((RADIUS-rCenter)/ZITA)))*Delta*Delta; // DOVREI APPROSSIMARE L'AREA DEL CERCHIO DI RAGGIO RADIUS...per RADIUS = 1, dovrei avere /pi... e ce l'ho!
// // //                 //forceMAG[] += rCenter < RADIUS ? +2.0 : 0.;// force applied in the particle's center 
// // //                 //forceMAG[] += rTrans  < RADIUS ? -1.0 : 0.;// force applied at flagella...Trans one? (left) 
// // //                 //forceMAG[] += rCis    < RADIUS ? -1.0 : 0.;// force applied at flagella...Cis ones? (right)
// // //         				//forceCellCounter += rCenter < RADIUS ? 1.0 : 0.;// so to count all the cells UNDER the particle's body
// // //               }
// // // 							//fprintf(stderr,"forceCellCounter %f, forceCellTotalArea %4.5e \n", forceCellCounter, forceCellTotalArea);
// // //         			foreach()
// // //         			totalForce[] += forceMAG[]; // for plotting purposes only...set aside so not to influence the eventual MPI loop
// // // 
// // // 
// // //         //value is divided by forceCellCounter, the total number of cells where force is defined.
// // //         //so to ensure that force integral is constant, regardless of mesh refinement.
// // //         foreach_face(x)
// // //           forcesFP.x[] +=STRENGTH*-sinTheta* forceMAG[]/forceCellCounter/Delta/Delta; // divided by sq(Delta) to be
// // //                                                                                       // independent on mesh size
// // //                                                                                       // FP 20250308 15h51
// // //         foreach_face(y)
// // //           forcesFP.y[] +=STRENGTH*+cosTheta* forceMAG[]/forceCellCounter/Delta/Delta;
// // // //                }
// // // // BEGIN // // ADD GRAVITY FORCES
// // //           for (double xp =  -L0; xp <= L0; xp += L0)
// // //             for (double yp = -L0; yp <= L0; yp += L0)
// // // //              for (double zp = -L0; zp <= L0; zp += L0)
// // //                 foreach_face(y) // aligned in y direction...
// // // 								  forcesFP.y[] -= GRAVITY*1.0*(0.5*(1.0+tanh((RADIUS-rCenter)/ZITA)))/forceCellCounter/Delta/Delta; // 20250308 17h43, validated as well
// // // // END
// // // // BEGIN // Add Lennard-Jones like repulsive potential
// // // #if LJ 1
// // // 					// Loop on all particles...but the one I'm working on now (_j_particle)
// // // 					#define x_j_particle pl[_l_particle][_j_particle].x+xp
// // // 					#define y_j_particle pl[_l_particle][_j_particle].y+yp
// // // 					#define x_J_particle pl[_l_particle][_J_particle].x
// // // 					#define y_J_particle pl[_l_particle][_J_particle].y
// // // 					coord potentialLJ; // the quantity I'm looking for... 
// // //           for (double xp =  -L0; xp <= L0; xp += L0){
// // //             for (double yp = -L0; yp <= L0; yp += L0){
// // // // PERIODICITY
// // //       				for (int _J_particle = 0; _J_particle < pn[_l_particle]; _J_particle++) { // now the index is _J_particle...)
// // // 								//fprintf(stderr,"_j_particle %d, _J_particle %d \n", _j_particle, _J_particle);
// // // 								if (_J_particle != _j_particle){
// // // 									double jJdistance = sqrt(sq(x_j_particle-x_J_particle)+sq(y_j_particle-y_J_particle));
// // // 									double jJdistanceX= sqrt(sq(x_j_particle-x_J_particle));
// // // 									double jJdistanceY= sqrt(sq(y_j_particle-y_J_particle));
// // // 								//	fprintf(stderr,"DIFFERENT PARTICLES, jJdistance = %4.5e \n",jJdistance);	
// // // 									double x_sign = sign(x_J_particle-x_j_particle);
// // // 									double y_sign = sign(y_J_particle-y_j_particle);
// // // 									//fprintf(stderr,"_j_particle %d, _J_particle %d, jJdistance %+4.5e, x_sign %+4.5e, y_sign %+4.5e  \n", _j_particle, _J_particle, jJdistance, x_sign, y_sign);
// // // 									double xProjection = -jJdistanceX*x_sign; // minus is important, because it is opposite wrt the direction vector...
// // // 									double yProjection = -jJdistanceY*y_sign;
// // // 									//fprintf(stderr,"_j_particle %d, _J_particle %d, jJdistance %+4.5e, x_sign %+4.5e, y_sign %+4.5e, xProjection %+4.5e, yProjection %+4.5e \n", _j_particle, _J_particle, jJdistance, x_sign, y_sign, xProjection, yProjection);
// // // 									double potentialLJ_modulus = 4.0*LJ_epsilon*(pow(LJ_sigma/jJdistance,12.0)-pow(LJ_sigma/jJdistance,6.0)); 
// // // 									potentialLJ.x = potentialLJ_modulus * xProjection;
// // // 									potentialLJ.y = potentialLJ_modulus * yProjection;
// // // 									//fprintf(stderr,"_j_particle %d, _J_particle %d, jJdistance %+4.5e, potentialLJ.x %+4.5e, potentialLJ.y %+4.5e \n", _j_particle, _J_particle, jJdistance, potentialLJ.x, potentialLJ.y);
// // // 									// Un po lungo, ma almeno capisco pezzo a pezzo cosa faccia...ed ha senso per me adesso. FP 20250310 11h38
// // // 									// Now I apply the potential as an additional force to the system...
// // // 									foreach_face()	// add potentialLJ to the forces acting on the BODY of the microswimmer...i.e. on rCenter.	
// // // 																	// identical as for applying a body-force, as I did for GRAVITY
// // // 								  	forcesFP.x[] -= potentialLJ.x*1.0*(0.5*(1.0+tanh((RADIUS-rCenter)/ZITA)))/forceCellCounter/Delta/Delta;
// // // 																	// MINUS SIGN, because it is a REPULSIVE POTENTIAL!!!
// // // 								}
// // // 								//else{
// // // 								//	fprintf(stderr,"SAME PARTICLE \n");
// // // 								//}
// // // 							}
// // // // PERIODICITY
// // // 						}
// // // 					}
// // // // END   // Lennar - Jones
// // // #endif
// // //       }
// // //       _l_particle++;
// // //     }
// // // }
// // // 
// // // 
// // // event properties(i++){
// // // 	compute_forcesFP(forcesFP);
// // //   fprintf(stderr,"COMPUTE FORCING FIELD \n");
// // // //Consider for re-orientation too!
// // // 	vorticity(u,omega);
// // // 	
// // // }
// // // 
// // // // Apply force owing to the presence of the microswimmer...
// // // event acceleration (i++) {
// // // 	foreach_face()
// // // 		//av.x[] = forcesFP.x[];
// // // 		av.x[] = forcesFP.x[]*fs.x[]; // times fs so to turn it off when it gets into an embbeded boundary!
// // // }
