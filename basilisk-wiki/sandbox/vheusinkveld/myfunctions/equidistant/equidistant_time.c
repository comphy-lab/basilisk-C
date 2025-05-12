/** 
# Equidistant averaging in time  

Testing the equidistant averaging functions for convergence in time.
First we define a function, cubed in space and squared in time. Moving at a speed of 1/tau.
De time integration of this function is define afterwards. The choice for a more general function would be better, but this is just to see if the functions actually do what they say. Convergence in space and time will be very specific to the problem in question. 
*/ 

#define func(s,d,f) (cube(s-sq(t/tau)) + cube(d) + cube(f))
#define Ifunc(s,d,f) (cube(d) + cube(f) + cube(s) - sq(s)*sq(t)/sq(tau) + 3*s*pow(t, 4)/5./pow(tau,4) - pow(t,6)/(7*pow(tau,6)))

double tau;
double ave_err; 		// Diagnosing total error 

#include "grid/octree.h"	// 3D
#include "run.h"

int diagi = 0;
int diaglvl = 4;		// Level at which we want equidistant diagnostisation

int minlevel = 4; 	 
int maxlevel = 4;		// set equal to diag level such that spacial error = 0
scalar b[];

double *equifield = NULL; 	// The diagnostics field

/** Used functions which are define at the end of the document */
int equi_diag (scalar s, int lvl, int diagii);
void equi_output (scalar s, FILE * fp, int lvl, int diagii);

int main(){    
    FILE * fp3 = fopen("output_time", "w"); // Specify te error output file
    L0=1.; 				
    int pp = 0;
/** Loop over different DT to see if average converges */
    for(DT = 0.1; DT > 0.001; DT = DT/2.){
        init_grid(1<<maxlevel);
        
        tau = 0.01; 		// 'fast' moving function 

        run();
	
	if(pid() == 0) {	
	    if(pp==0){
  	        fprintf(fp3, "err\tdt\tL0\n");
  	    }
            fprintf(fp3, "%g\t%g\t%g\n", ave_err, DT, L0); // Write away error
	}
	
 	free(equifield);  	// we dont want memory leaks
	equifield = NULL; 	// reset the field pointer 
        diagi=0;
	pp++;
    }	
}

/** Ability to test for non uniform grid */ 
event init(i=0){
	//unrefine(y <= L0/4. && level > minlevel-1); 
}

/** Set values of scalar b based on defined function */
event update(i++) {
    dt = dtnext(DT);
    foreach() {
	b[] = func(x, y, z);
    }
    boundary({b});
}

/** Keeping track of the field b (averaging), here we cheat a bit since the error determination is defined in the function itself. */    
event diag(i++){
    fprintf(stderr, "time=%g\n",t);
    diagi = equi_diag(b, diaglvl, diagi); 
}

/** Write away output of the averaging function */
event end(t=1){
    FILE * fp2 = fopen("output", "w");
    equi_output(b, fp2, diaglvl, diagi);
    fclose(fp2);
}
/** ## Results 
The average error seems to nicely converge in a linear fasion.

  ~~~gnuplot
    set xr [0.:0.11]
    set yr [0.:0.3]
    set xlabel 'DT'
    set ylabel 'Error
    plot 'output_time' using 2:1 title 'DT vs err' 
  ~~~

*/


/** ## Used functions 

NOTE: Very much in progress.

These functions are geared towards equidistant diagnostics. The level for which diagnostics are done can be specified. Values will be linearly interpolated if necesarry. 

So far this function works with MPI and only works in 3D. 


Short description of what they do:

//void equi_diag(scalar s){
    adds values to the equifield array
}

//void equi_output(scalar s, FILE * fp){
    outputs equifield to file fp in the format: 
    place = x*n*n + y*n + z*n
}


TODO implement 2D
TODO walking average since values in equifield get big, or divide before summing?
*/

int equi_diag (scalar s, int lvl, int diagii){  
    int len = 1;		// TODO implement multiple scalar field averaging
    int n = 1<<lvl;
    if(!equifield) {
        equifield = calloc(len*n*n*n, sizeof(double)); // Init with zeros
    } 	

    double dDelta = 0.9999999*L0/n;  

    restriction({s}); 		// Make sure that coarse level values exist
    boundary({s}); 		// Update boundaries
    
    // Iterate over diag level of leaf cell
    foreach_level_or_leaf(lvl) {
  	if(level < lvl){
	    int pn = 1<<(lvl-level);
	    for(int i = 0; i < pn; i++){
		for(int j = 0; j < pn; j++){
		    for(int k = 0; k < pn; k++){
			
			double ctrans = dDelta/2. - Delta/2.;	// Translation to corner of cell
			double xx = x + ctrans + i*dDelta;   	// Coordinates in equidistant grid
			double yy = y + ctrans + j*dDelta;
 			double zz = z + ctrans + k*dDelta;

		        int ii = round((xx - dDelta/2.)/dDelta); // Iteration in equidistant grid
	         	int jj = round((yy - dDelta/2.)/dDelta);
            	  	int kk = round((zz - dDelta/2.)/dDelta);

	    		int place = ii*n*n + jj*n + kk*len;	// Location in equidistant grid 
			double temp = interpolate(s, xx, yy, zz); // Linear interpolation 

 			equifield[place] += temp;
		    }
		}
	    }

        } else {
 
            int ii = round((x - dDelta/2.)/dDelta); // x adaptive is x equidistant 
       	    int jj = round((y - dDelta/2.)/dDelta); // Derive place in equidistant grid
     	    int kk = round((z - dDelta/2.)/dDelta);

	    int place = ii*n*n + jj*n + kk*len;
	    equifield[place] += s[];
	}
    }
    return diagii + 1;
}

void equi_output (scalar s, FILE * fp, int lvl, int diagii){
    int len = 1;
    int n = 1<<lvl;
    double dDelta = 0.9999999*L0/n;  

    if(pid() == 0){
// Make sure all thread values are summed up using MPI
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, equifield, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
@endif
    double temperr = 0;
    int len = 1;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
		
		double x = (i + 0.5)*dDelta;   	// Coordinates in equidistant grid
		double y = (j + 0.5)*dDelta;
 		double z = (k + 0.5)*dDelta;


		int place = i*n*n + j*n + k*len;
		double temp = (equifield[place]/((double)diagii)); // Divide by number of additions		
		double val = Ifunc(x,y,z);

	        fprintf(fp, "%g\t", temp);
		temperr += fabs((temp - val)/val);

	    }
	}    
    }
    ave_err = temperr/(cube(n));
    }
@if _MPI
    else // slave
    MPI_Reduce (equifield, NULL, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
@endif
}

