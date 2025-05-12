//#include "LS_reinit.h"
#include "redistance2.h"
#include "diffusion.h"
#include "tracer.h"


bool advect_diff_phase_field = 0;
int reinit_skip_steps = 100;
bool surfactant = 1;


scalar c1[], gamma2[],d2[], pfield[]; //volumetric surfactant conc, area surfactant conc, signed distance, phase-field respectively

scalar * tracers = {c1,d2,pfield};

double zeta = 1. [*]; // Velocity scale parameter >= max|u| for CFL/boundedness criteria
//Epsilon is the interface thickness parameter > 0.5 delta x. Cannot avoid macro here because I need the smallest delta x at all times
#define EPSILON ((L0/(1 << grid->maxdepth))*0.75) 
double D_s = 0.01 [*];  //Surface Diffusivity of surfactant
double varepsilon = 1.e-6; // \varepsilon in paper 
double sharpening_coefficient = 2. ; // a in Jain 2023 


struct HDiffusion {
  face vector D;
  face vector beta;
};

static void h_relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct HDiffusion * p = (struct HDiffusion *) data;
  face vector D = p->D, beta = p->beta;

  scalar c = a;
  foreach_level_or_leaf (l) {
    double n = - sq(Delta)*b[], d = cm[]/dt*sq(Delta);
    foreach_dimension() {  
      n += D.x[1]*a[1] + D.x[]*a[-1] +
	Delta*(beta.x[1]*a[1] - beta.x[]*a[-1])/2.;
      d += D.x[1] + D.x[] - Delta*(beta.x[1] - beta.x[])/2.;
    }
    c[] = n/d;
  }
}

static double h_residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct HDiffusion * p = (struct HDiffusion *) data;
  face vector D = p->D, beta = p->beta;
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = D.x[]*face_gradient_x (a, 0) + beta.x[]*face_value (a, 0);
  foreach (reduction(max:maxres)) {
    res[] = b[] + cm[]/dt*a[];
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta;
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] + cm[]/dt*a[];
    foreach_dimension()
      res[] -= (D.x[1]*face_gradient_x (a, 1) -
		D.x[0]*face_gradient_x (a, 0) +
		beta.x[1]*face_value (a, 1) -
		beta.x[0]*face_value (a, 0))/Delta;  	  
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif // !TREE    
  return maxres;
}

event stability(i++){
	
 double maxvel = 1.e-6 [*];

  if (surfactant ){
    foreach_face(reduction(max:maxvel))
    if (u.x[] != 0.) {
      // double positivity_criteria = 2.*D_s/(fabs(u.x[]) + D_s*sharpening_coefficient/2./epsilon); Delta/fabs(u.x[]);
      // assert(Delta < positivity_criteria);
          double deltas = (pfield[]*(1. - pfield[]))/EPSILON;

      if (deltas > 0.01){
        if (fabs(u.x[]) > maxvel) maxvel = fabs(u.x[]);
      }
    }
    //smallest grid size
    
    double deltaxmin = 1. [*] ;
    deltaxmin = L0/(1 << grid->maxdepth) ;
    
  //   double positivity_criteria = 2.*D_s/(maxvel + D_s*sharpening_coefficient/2./EPSILON);
//    if (deltaxmin > positivity_criteria){
//	    printf("Warning: positivity criteria not fullfilled, increase resolution, dxmin = %e pc = %e level = %d",deltaxmin,positivity_criteria,grid->maxdepth);
//    }
     // assert(deltaxmin < positivity_criteria);
    zeta = 1.1*maxvel; //zeta must be more than maxvel
    // assert(zeta > maxvel);
    double dt = dtmax;
  // time restriction for phase-field transport
  if (advect_diff_phase_field){
     dt = 0.75*sq(deltaxmin)/zeta/EPSILON;
      if (dt < dtmax)
        dtmax = dt;
  }
    // time restriction for surfactant-transport  
    dt = 0.75*sq(deltaxmin)/4./D_s;
    #if dimension == 3
      dt =  0.75*sq(deltaxmin)/6./D_s;
    #endif
      // if (dt < dtmax)
      //   dtmax = dt;
  }

}
int counter = 0;

//small values pfield anywhere else in domain cause
//c1 to leak hence clamp2 keeps it in check
//mass conservation of c1 isn't affected by this
//vof takes care of fluid mass
//conservation of pfield 
double clamp2 (double a){
  if (a < 1.e-6)
    return 0.;
  else if (a > 1. - 1.e-6)
    return 1.;
  else
    return a;
}
  
 
//"properties" event is called by other routines too, we named it differently so that it is not called
// more than once in a time step because redistancing is an expensive step
event properties2 (i++)
{ 
  // If not redistancing at each step phase field must be transported
  if (reinit_skip_steps > 1)
   advect_diff_phase_field = 1;
   else
   advect_diff_phase_field = 0;
  // Avoid redistancing at the 0th step because d2 is initialized accurately at cheap cost
  // Skip redistancing and let the phase field advect and diffuse to save the cost
  if (counter >=0  && counter%reinit_skip_steps == 0){
	scalar d2temp[];
    double deltamin = L0/(1 << grid->maxdepth);
    foreach()
      d2temp[] = (2.*f[] - 1.)*deltamin*0.75;
    #if TREE
      restriction({d2temp});
    #endif
    // LS_reinit (d2, dt = 0.5*L0/(1 << grid->maxdepth),
    //     it_max = 0.5*(1 << grid->maxdepth));
      redistance (d2temp, imax = 0.5*(1 << grid->maxdepth));
     printf("\n %d redistancing",i); fflush(stdout);
	double d2max = statsf(d2temp).max;
	double d2min = statsf(d2temp).min;
	bool signed_distance_faulty = 0;	
	if (d2max > 6. || d2min < -6.){

		signed_distance_faulty = 1; //something went wrong in computing distance
		counter = counter - 2; //try again after two computational steps
	}	 	
	if (!signed_distance_faulty){
		foreach(){
			d2[] = d2temp[];
			pfield[] = 0.5*(1. - tanh((d2[])/2./EPSILON));
			pfield[] = clamp2(pfield[]);
		}
		boundary({pfield});
	} 
  }

  counter++;

}
event tracer_diffusion (i++)
{
if (surfactant){

  scalar psi[];
  face vector cflux[],geta[];
       
  face vector  beta[];
  scalar r[];
  //~ scalar pfield[];
  face vector D[];
if (advect_diff_phase_field){
  foreach(){
    //~ pfield[] = 0.5*(1. - tanh((d[])/2./epsilon));
        pfield[] = clamp2(pfield[]);

    psi[] = EPSILON*log((pfield[] + varepsilon)/(1. - pfield[] + varepsilon));
  }
  boundary({psi});

    foreach_face(){      
      cflux.x[] = 0.;
      double psif = (psi[] + psi[-1])/2.;
      
      double gradpsi = (psi[] - psi[-1])/Delta;
      if (fabs(gradpsi) > varepsilon){
        #if dimension == 2

         double psiup = 0.25*(psi[] + psi[-1] + psi[0,1] + psi[-1,1]);
         double psidown = 0.25*(psi[] + psi[-1] + psi[0,-1] + psi[-1,-1]);
          double mag_grad_psi = sqrt(sq(psi[] - psi[-1]) + sq(psiup - psidown))/Delta;
        #endif

        #if dimension == 3
         
         double psiup = 0.25*(psi[0,0,0] + psi[-1,0,0] + psi[0,1,0] + psi[-1,1,0]);
         double psidown = 0.25*(psi[0,0,0] + psi[-1,0,0] + psi[0,-1,0] + psi[-1,-1,0]);
         
          double psifront = 0.25*(psi[0,0,0] + psi[-1,0,0] + psi[0,0,1] + psi[-1,0,1]);
         double psiback = 0.25*(psi[0,0,0] + psi[-1,0,0] + psi[0,0,-1] + psi[-1,0,-1]);

          double mag_grad_psi = sqrt(sq(psi[0,0,0] - psi[-1,0,0]) + sq(psiup - psidown) + sq(psifront - psiback))/Delta;

        #endif

        
        cflux.x[] = fm.x[]*(1. - sq(tanh(psif/2./EPSILON)))*gradpsi/mag_grad_psi;
      }

      D.x[] = fm.x[]*zeta*EPSILON;
      beta.x[] = 0.;
      
    }

    foreach(){
      #if dimension == 2
        r[] = - cm[]*pfield[]/dt + 0.25*zeta*(cflux.x[1] - cflux.x[] + (cflux.y[0,1] - cflux.y[]))/Delta;
      #endif

      #if dimension == 3
        r[] = - cm[]*pfield[]/dt + 0.25*zeta*(cflux.x[1,0,0] - cflux.x[0,0,0] + cflux.y[0,1,0] - cflux.y[0,0,0] + cflux.z[0,0,1] - cflux.z[0,0,0]  )/Delta;

      #endif
    }


   
   
    restriction ({D,beta, cm});

    struct HDiffusion q1;
    q1.D = D;
    q1.beta = beta;
    
    mg_solve ({pfield}, {r}, h_residual, h_relax, &q1);
   
}


  
    foreach() {
      r[] = - cm[]*c1[]/dt;
    }

     foreach(){
	    pfield[] = clamp2(pfield[]);

      psi[] = EPSILON*log((clamp(pfield[],0.,1.) + varepsilon)
              /(1. - clamp(pfield[],0.,1.) + varepsilon));
    }

 
  foreach_face() {

    double phif = (pfield[] + pfield[-1])/2.;
    D.x[] = fm.x[]*D_s;
    double gradpsi = (psi[] - psi[-1])/Delta;
    beta.x[] = 0.;
    
        if (fabs(gradpsi) > varepsilon){

        #if dimension == 2  
        // The derivative of Tanh is too steep and linear interpolation messes it up
          double psiup = 0.5*(psi[0,1] + psi[-1,1]);
          double psidown = 0.5*( psi[0,-1] + psi[-1,-1]);
          double mag_grad_psi = sqrt(sq(psi[] - psi[-1]) + sq(psiup - psidown)/4.)/Delta;

        #endif

        #if dimension == 3

         double psiup = 0.5*(psi[0,1,0] + psi[-1,1,0]);
         double psidown = 0.5*( psi[0,-1,0] + psi[-1,-1,0]);

         double psifront = 0.5*(psi[0,0,1] + psi[-1,0,1]);
         double psiback = 0.5*( psi[0,0,-1] + psi[-1,0,-1]);

         double mag_grad_psi = sqrt(sq(psi[] - psi[-1]) + sq(psiup - psidown)/4.  + sq(psifront - psiback)/4.)/Delta;

        #endif
        beta.x[] = - D.x[]* sharpening_coefficient*(0.5 - phif)/EPSILON * gradpsi/mag_grad_psi ;
      }
    }

    restriction ({D, beta, cm});    
    struct HDiffusion q;
    q.D = D;
    q.beta = beta;
    
    
    mg_solve ({c1}, {r}, h_residual, h_relax, &q);

    
}
}

/**
## References

~~~bib
@article{jain2023model,
  title={A model for transport of interface-confined scalars and insoluble surfactants in two-phase flows},
  author={Jain, Suhas S},
  journal={arXiv preprint arXiv:2311.11076},
  year={2023}
}
@article{farsoiya2024vof,
  title={A coupled Volume of Fluid - Phase Field method for direct numerical simulation of insoluble surfactant-laden interfacial flows and application to rising bubbles},
  author={Farsoiya, Palas Kumar and Popinet, Stephane and Stone, Howard A. and Deike, Luc},
  journal={prep},
  year={2024}
}
~~~
*/
