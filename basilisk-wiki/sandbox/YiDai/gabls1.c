/**
# The First GABLS Intercomparison
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "diffusion.h"

/**
## The declaration of the parameters 
*/
#define GV 9.81            // Gravitational constant m s-2
#define IN_VER 0.01        // Inversion 0.01 K m-1
#define REF_TEMP 263.5     // Reference surface potential temperature K
#define INI_TEMP 265.      // Initial potential temperature K
#define U_G 8.             // geostraphic wind m/s
#define TEMP_init(z) (z<100.) ? (GV*INI_TEMP/REF_TEMP) : (GV/REF_TEMP*(INI_TEMP+IN_VER*(z-100.)))

#define COOLING_RATE 0.25  // Surface cooling rate   
#define SUR_FORCE(t) GV/REF_TEMP*(INI_TEMP-COOLING_RATE*(t/3600.))     // surface cooling 0.25K/h  -->  0.25/3600 * DT = surface cooling rate at specific time


face vector av[];
#define damp(s) -s*(exp(y-300)-1)*(y>300)/2.

/* randomnizer function adapt from DALES */
int randomnizer(scalar field, int level, double amplitude, double ir);


/* damping function for u and b */
/* up[] - = (u[]-u_average[])*tsc[]
   tsc[]  = rnuo * sin(0.5* M_PI * (y-300)/100)**2
   rnuo   = 2.75e-3
*/

int minlevel = 3, maxlevel = 4;
double ue=0.01;
double be=0.01;

/* double T_total = 32400.; */
double T_total = 600.;  // for testing 

scalar b[];
scalar * tracers = {b};


int main(){
        L0 = 400.;  // domain size

        DT = 1;


        init_grid(1<<3);

        run();
}


/**
## Initial condition */

event init(t = 0){

        /* Boundary condition */

        u.n[bottom] = dirichlet(0.);
        u.t[bottom] = dirichlet(0.);

        periodic(left);

        b[bottom] = dirichlet(GV/REF_TEMP*INI_TEMP);
        b[top] = neumann(GV/REF_TEMP*IN_VER);



        #if dimension ==3
            /* boundary condition*/
            u.r[bottom] = dirichlet(0.);
            periodic(left);
        #endif

        /* wind profile and temperature profile */
        do {
        foreach(){
                u.x[] = U_G;
                b[] = TEMP_init(y);
                /* randomnizer(b[], 1, 0.1, 0.); */
        }
        boundary ({b, u});
        } while(adapt_wavelet({u,b}, (double []){ue,ue,ue,be}, maxlevel, minlevel).nf);
        boundary ({b, u});
}


/* the boundary condition of b changing with time; MOST theory? */
event Sur_forcing(i++){
        b[bottom] = dirichlet(SUR_FORCE(t));
        periodic(left);
}


/**
## Sponge layer higher than 300m */



event acceleration (i++) { // how to add damping to tendency of u and b (inside of center solver) 
  foreach_face(y)
    av.y[] = (b[] + b[0,-1])/2.;
}

event tracer_diffusion (i++) {
  diffusion (b, dt);
}

/**
## Coriolis force */



/**
## adaptive grid */
event adapt(i++){
    adapt_wavelet((scalar *){u,b}, (double []){ue,ue,ue,be},maxlevel,minlevel);
}


/**
## output data */
event dumpfield(t = 300; t+=300){
        char datadump[180];
        sprintf(datadump, "dump_t=%dT",(int)((t/300)));
        dump(file = datadump, list = all);
}

event end(t = T_total){
}

/*
int randomnizer(scalar field, int level, double amplitude, double ir){
        int imm = 134456;
        int ia = 8121; 
        int ic = 28411;
        scalar ran[];
        foreach(){
        ir = (ir*ia+ic)%(imm);
        ran[] = ir/imm;
        field[] = field[] + (ran[]-0.5)*2.0*amplitude;
        }
        return 0;
}

*/


