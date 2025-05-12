/**
This is a standard file to reproduce the simulations presented in:

[Rivière, Mostert, Perrard and Deike, Sub-Hinze scale bubble production in turbulent bubble break up. Journal of Fluid Mechanics, 917:A40, 2021](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/subhinze-scale-bubble-production-in-turbulent-bubble-breakup/7F8DADA87ABE117F5FE39183417D650B).

[Perrard, Rivière, Mostert and Deike, Bubble deformation by a turbulent flow. Journal of Fluid Mechanics, 920:A15, 2021](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/bubble-deformation-by-a-turbulent-flow/D6E47109EEE9595E52FDD2E1B9E83A1D).

# Deformation and breakup of a bubble in a turbulent environment.

A fully developed turbulent flow is generated as a precursor simulation, based on the example case
isotropic.c, itself based on the case investigated by
[Rosales & Meneveau, 2005](/src/references.bib#rosales2005); once fully developed, a gas bubble is inserted.

This simulation is in two parts. First, an equivalent version of 
isotropic.c is run. Second, the simulation is dumped and restarted, and the bubble is initialized in the flow. The turbulent flow may be forced or unforced in the liquid phase after bubble insertion (usually forced), but in any case there is no forcing inside the bubble.

Both parts are incorporated into this one problem file using preprocessor commands.

Include the Output_facets to have a good description of the interface

Include a way to restart the files. Should have a file dump file called dump.

Level in the bulk fixed at 8.

Write several dump files to avoid problems when you want to restart your simulations.

Initially flow inside the bubble put to 0.

 */

#include "grid/octree.h"
#include "navier-stokes/centered.h"

/**
We use the $\lambda_2$ criterion and Basilisk View for visualisation
of vortices. */

#include "lambda2.h"
#include "view.h"

/**
Include surface tension and two-phase schemes.*/
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"

/**
We monitor performance statistics and control the maximum runtime. */

#include "navier-stokes/perfs.h"
#include "maxruntime.h"

/**
   Include tag functionalities for droplet/bubbles calculation. */
#include "tag.h"

/**
   Include density and viscosity ratios.*/
#define RHOR 850.0
#define MUR 100.0

/**
   Define domain size. */
#define WIDTH 120.0

/**
We need to store the variable forcing term. */
face vector av[];

/**
The code takes the level of refinement as optional command-line
argument (as well as an optional maximum runtime). */

int maxlevel = 8;	//only for the compilation
double MAXTIME = 20;	//idem
int FORCED = 0;		//equals 0 if it is not forced, 1 otherwise
double amp_force = 0.1; //amplitude of the forcing
double SIGMA = 18.;
double R0 = 8.0;


int main (int argc, char * argv[]) {
  if (argc > 1)
  {
    maxlevel = atoi(argv[1]); //maximum refinement level
  	MAXTIME = atoi(argv[2]); //maximum runtime
  	SIGMA = atof(argv[3]); //surface tension strength, which sets the Weber number
	R0 = atof(argv[4]); //bubble diameter, if inserted
        FORCED = atoi(argv[5]); //equals 1 if forced, 0 otherwise
        amp_force = atof(argv[6]);  //amplitude of the forcing
  }
  size (WIDTH);
  foreach_dimension()
    periodic (right);

  /**
     The acceleration is defined below.
     The level of refinement is *maxlevel*. */
  mu1 = 0.01 * sq(WIDTH/(2.0*pi))/(2.0);
  mu2 = mu1/MUR;
  /**
     Use reference density for fluid, and define gas by density ratio.*/
  rho1 = 1.0;
  rho2 = rho1/RHOR;
  f.sigma = SIGMA;
  a = av;
#if TREE
  N = 1 << (maxlevel-3);
#else
  N = 1 << maxlevel;
#endif

  /**
     Set tolerance on flow divergence. */
  TOLERANCE = 1e-4;
  run();
}


/**
## Initial conditions

The initial condition is "ABC" flow. This is a laminar base flow that 
is easy to implement in both Basilisk and a spectral code. 
The init event is structured so that if there are no dump files from which to restart, a precursor simulation (i.e. without a bubble) is begun. If a dump file is present, its name will determine what kind of simulation is to be done. If it is called "restart", then a bubble is injected into the already turbulent flow. If it is called "dump", then it is merely a checkpoint from which to restart, and can be a precursor or not.
*/

event init (i = 0) {
  /**
     If there are no existing dump files from which to restart, initialize the flow (without a bubble). */
  if (!restore (file = "restart") && !restore(file="dump"))
    {
      double waveno = WIDTH/(2.0*pi);
      double amp = WIDTH/(2.0*pi);
      foreach() {
        f[] = 1.0;
        double xw = x/waveno;
        double yw = y/waveno;
        double zw = z/waveno;
        u.x[] = amp*( cos(yw) + sin(zw) );
        u.y[] = amp*( sin(xw) + cos(zw) );
        u.z[] = amp*( cos(xw) + sin(yw) );
    }
    }
  /**
     If there is a dump file named "dump", do not do anything.*/
  else if(restore (file = "dump")){
  }
  /**
     If there is a dump file named "restart", then insert the bubble. */
	else
    {

      fraction (f, sq(x-WIDTH/2.0) + sq(y-WIDTH/2.0) + sq(z-WIDTH/2.0) - sq(R0));
	foreach(){      
	foreach_dimension(){
          	u.x[] = f[]*u.x[];
      }
}
    }
}

/**
   Include accelerations:*/
event acceleration (i++) {
/**
## Linear forcing term
We compute the average velocity and add the corresponding linear
forcing term, if desired at runtime by FORCED variable. */
  if (FORCED)
{
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
  foreach_face()
      av.x[] += f[]*amp_force*((u.x[] + u.x[-1])/2. - ubar.x);
}
}
/**
## Outputs

We log the evolution of viscous dissipation, kinetic energy, and
microscale Reynolds number. This is as in the isotropic.c example case.*/

event logfile (i++) {
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }

  double ke = 0., vd = 0., vol = 0.;
  foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
    vol += dv();
    foreach_dimension() {
      // mean fluctuating kinetic energy
      ke += dv()*sq(u.x[] - ubar.x);
      // viscous dissipation
      vd += dv()*(sq(u.x[1] - u.x[-1]) +
                  sq(u.x[0,1] - u.x[0,-1]) +
                  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
  }
  ke /= 2.*vol;
  vd *= mu1/vol;
  static FILE * fd = fopen("stats.dat","a");//"w" before
  if (i == 0) {
    fprintf (ferr, "t dissipation energy Reynolds\n");          //screen
    fprintf (fd, "t dissipation energy Reynolds\n");
  }
  fprintf (ferr, "%g %g %g %g\n",
           t, vd, ke, 2./3.*ke/mu1*sqrt(15.*mu1/vd));
  fprintf (fd, "%g %g %g %g\n",
           t, vd, ke, 2./3.*ke/mu1*sqrt(15.*mu1/vd));
}

/**
   We also want to count the drops and bubbles in the flow. */
event countDropsBubble(i++)
{
  scalar m1[]; //droplets
  scalar m2[]; //bubbles
  foreach(){
    m1[] = f[] > 0.5; //i.e. set m true if f[] is close to unity (droplets)
    m2[] = f[] < 0.5; //m true if f[] close to zero (bubbles)
  }
  int n1 = tag(m1);
  int n2 = tag(m2);
  /**
     Having counted the bubbles, now we find their size. This example
     is similar to the jet atomization problem. We are interested in
     the volumes and positions of each droplet/bubble.*/
  double v1[n1]; //droplet
  coord b1[n1];  //droplet
  double v2[n2]; //bubble
 coord b2[n2];  //bubble
  /**
     We initialize: */
  for (int j=0; j<n1; j++)
    {
      v1[j] = b1[j].x = b1[j].y = b1[j].z = 0.0;
    }
  for (int j=0; j<n2; j++)
    {
      v2[j] = b2[j].x = b2[j].y = b2[j].z = 0.0;
    }
  /**
     We proceed with calculation. */
  foreach_leaf() //droplets
    {
      if (m1[] > 0) {
      int j = m1[] - 1;
      v1[j] += dv()*f[]; //increment the volume of the droplet
      coord p = {x,y,z};
      foreach_dimension()
        b1[j].x += dv()*f[]*p.x;
      }
    }
  foreach_leaf() //bubbles
    {
      if (m2[] > 0) {
        int j = m2[] - 1;
        v2[j] += dv()*(1.0-f[]);
        coord p = {x,y,z};
        foreach_dimension()
          b2[j].x += dv()*(1.0-f[])*p.x;
      }
    }
  /**
     Reduce for MPI. */
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v1, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b1, 3*n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, v2, n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b2, 3*n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  /**
     Output the volume and position of each droplet to file. */
  static FILE * fdrop = fopen("droplets.dat","w");
  static FILE * fbubb = fopen("bubbles.dat","w");
  for (int j=0; j<n1; j++)
    {
      fprintf (fdrop, "%d %g %d %g %g %g %g\n", i, t,
               j, v1[j], b1[j].x/v1[j], b1[j].y/v1[j], b1[j].z/v1[j]);
    }
  for (int j=0; j<n2; j++)
    {
      fprintf (fbubb, "%d %g %d %g %g %g %g\n", i, t,
               j, v2[j], b2[j].x/v2[j], b2[j].y/v2[j], b2[j].z/v2[j]);
    }
}

/**
We generate a movie of the vortices. Commented lines are optional alternative outputs.*/
#define POPEN(name, mode) fopen (name ".ppm", mode)

event movie (t += 0.25)
{
  view (fov = 44, camera = "iso", ty = .2,
        width = 600, height = 600, bg = {1,1,1}, samples = 4);
  clear();
  squares ("u.y", linear = true, n = {0,1,0});
  squares ("u.x", linear = true, n = {1,0,0});
  squares ("u.z", linear = true, n = {0,0,1});
  //scalar omega[];
  //vorticity (u, omega);
  //squares ("omega", linear = true, n = {0,1,0});
  //scalar l2[];
  //lambda2 (u, l2);
  //isosurface ("l2", -1);
  draw_vof ("f");
 {
    static FILE * fp = POPEN ("movie", "w");
    save (fp = fp);
  }
}

/**
   Gas-liquid interface data is also outputted. */

event getInterface(i++)
{
  char fixname[100];
  char fiyname[100];
  char fizname[100];
  char fiintname[100];
  sprintf(fixname, "fix_%d.dat", pid() );
  sprintf(fiyname, "fiy_%d.dat", pid() );
  sprintf(fizname, "fiz_%d.dat", pid() );
  sprintf(fiintname, "fint_%d.dat", pid() );
  FILE * fix = fopen(fixname, "a");
  FILE * fiy = fopen(fiyname, "a");
  FILE * fiz = fopen(fizname, "a");
  FILE * fint = fopen(fiintname, "a");
  fprintf (fix, "%g ", t);
  fprintf (fiy, "%g ", t);
  fprintf (fiz, "%g ", t);
  fprintf (fint, "%g\n", t);
  output_facets(f,fint);

  foreach()
    {
      if (f[] > 0.4 && f[] < 0.6)
        {
          fprintf (fix, "%g ", x);
          fprintf (fiy, "%g ", y);
          fprintf (fiz, "%g ", z);
        }
    }
  fprintf (fix, "\n");
  fprintf (fiy, "\n");
  fprintf (fiz, "\n");
  fprintf (fint, "\n");
  fclose(fix);
  fclose(fiy);
  fclose(fiz);
  fclose(fint);

  return 0;
}

/**
We refine on the velocity field in the precursor, and on the VOF field in the full simulation.
 */
#if TREE
event adapt (i++) {
  double uemax = 0.2*normf(u.x).avg;
  double femax = 0.01;
  adapt_wavelet ((scalar *){u}, (double[]){uemax,uemax,uemax}, maxlevel);
  //adapt_wavelet ((scalar *){f, u}, (double[]){femax, uemax, uemax, uemax},maxlevel);
}
#endif


/**
   We output a full snapshot every time unit. */
int j = 0;

event snapshot (t=1; t <=MAXTIME; i+=10)
{
  //dump(NAME);
  //
  char namedump[80];
sprintf(namedump,"dump_%d",j);
  dump(file = namedump);
j++;
if (j==2)
j=0;		  
}

event allsnapshots (t=60; t<=MAXTIME; t +=1)
{
char step[100];
p.nodump = false;
sprintf(step, "time=%g",t);
dump(step); 
}

/**
   End the simulation after 300 time units. */
event end (t=MAXTIME) {
  dump ("end");
}

