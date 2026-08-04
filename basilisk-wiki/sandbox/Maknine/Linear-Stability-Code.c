#include "stabilityRT.h"
#include "amr.h"

/*Initialization of the fields to be solved:*/
scalar REPsi[];
scalar IMPsi[];
scalar REOmega[];
scalar IMOmega[];


scalar f[];

double omtol = 0.01;

double lambda_max;

/*Initialization of the coefficients of the differential equation:*/


/*Boundary Conditions:*/

REPsi[top]=dirichlet(0.);
REPsi[bottom]=dirichlet(0.);

IMPsi[top]=dirichlet(0.);
IMPsi[bottom]=dirichlet(0.);

REOmega[top]=dirichlet(0.);
REOmega[bottom]=dirichlet(0.);

IMOmega[top]=dirichlet(0);
IMOmega[bottom]=dirichlet(0);

/*REPsi[top]=neumann(0.);
REPsi[bottom]=neumann(0.);

IMPsi[top]=neumann(0.);
IMPsi[bottom]=neumann(0.);

REOmega[top]=neumann(0.);
REOmega[bottom]=neumann(0.);

IMOmega[top]=neumann(0);
IMOmega[bottom]=neumann(0);*/

double fraction_profile(double x, double y, double DeltaI) {

    if (y>DeltaI/2) return 1.;              
    
    if (y<-DeltaI/2) return 0.;

    return 0.5 + y/DeltaI;


}
int main() {

    /*Initialisation of the domain*/
    L0 =10;
    
    Y0 = -L0/2.;
    
    X0=0;

    periodic(right);


    /*Input Parametres*/

    double DeltaI=0.1; //Interface Thickness

    rho1 = 1.; rho2 = 0.7; //Density Ratio (Light/Heavy)
    mu1 = 1e-5; mu2 = 1e-5;

    TOLERANCE = 1e-7;//Tolerance of the Solver

    double k0=0.1; //Initial Wave Numer
    double kmax=10; //Max Wave Numer
    double k;

    double sigma=0.; //Surface Tension
    g = 1.;
    h = 1.;


    lambda_max = 0.23776144401031843;

    

    int mylev=7; 
    init_grid (1 << (mylev-2));;
    int Nobj = pow(2,mylev*dimension);  // initial objective number of element


    //initial guess
    foreach(){
	    REPsi[] = 1. ;
	    IMPsi[] = 1. ;
	    REOmega[] = 1. ;
	    IMOmega[] = 1. ;
    }
    om= get_omega (interpolate (REPsi, 0, 0), interpolate (IMPsi, 0, 0), k0);
		    

    double complex omold = om;
    double complex omerr;

    /* Loop over the wave number to get a complete dispersion relation */

    //for (k=k0 ;k<=kmax; k+=0.1) {
    k=(2*pi/lambda_max)*1;


    for (mylev=7 ;mylev<=12; mylev++) {
        Nobj = pow(2,mylev*dimension);

    

	    for (int ki=0; ( fabs((double)(grid->tn - Nobj)/Nobj)) > 0.1 ; ki++){
		    /* Definition of the viscosity profile */

		    foreach ()
			    f[] = fraction_profile(x, y, DeltaI);

		    double complex omnew = solve_stability_RT (REPsi, IMPsi, REOmega, IMOmega, f, rho1, rho2, sigma, g, k);
			    
		    double complex omerr = (omnew-omold)/omnew;
           
//		    fprintf(stdout, " sol %g %g N/Nobj %g\n",creal(om),cimag(om), (double) grid->tn/Nobj);

		    if (ki < 3 ) {  // estimate AMReps the first run with uniform refinement 
			    struct PreFactorData cd = compute_prefactors (2, {REOmega} );
			    AMReps =  cd.cuniform*pow(Nobj,-1)/pow(Nobj,1./2.); 
		    } else {
			    AMReps = update_epsilon_control(Nobj);
		            if ( cabs(omerr) < omtol ) break; //tolerance frequency calculation
		    }

		    adapt_metric({REOmega});

		    omold=omnew;

		    if (ki > 20) break;
	    }


	    fprintf(stdout, "  %g %g %g %lg %lg \n", k, creal(om),cimag(om),creal(omerr),cimag(omerr));
    }

 

}
