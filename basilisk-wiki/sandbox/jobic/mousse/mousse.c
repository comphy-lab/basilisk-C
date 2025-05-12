
/**
# Howto run

Sequential run           : CFLAGS='-DCLUSTER=1' make clean mousse.tst
                           This creats a dump file with the calculated levelset from the distance function
			   Should be created for each maxlevel.

For parallel             : The file Dump should be created with a sequential run, then 
                           CC='mpicc -D_MPI=8' CFLAGS='-DCLUSTER=1' make clean mousse.tst
                           cd mousse; ./mousse -restoreDump
			   
Run on workstation       : make clean mousse.tst
                           then use the navigator to open the html file to run it
*/

/**
# ToDo

  *) Restart not working from 1s physical time of the simulation

  *) Parallel not always working on 16 processors (mainly smaller discretisations).
    
*/

/**
# Example file
It's inspired from the example [Here](http://basilisk.fr/src/examples/porous3D.c)
*/

#include <libgen.h> /* for basename func */

#include "grid/octree.h"
#include "embed.h"

#include "navier-stokes/centered.h"
//#include "navier-stokes/double-projection.h"
//#include "navier-stokes/perfs.h"

//Borrowed from acastillo's sandbox
#include "output_vtu_foreach.h"
#include "curvature.h"
#include "distance.h"

#if !defined(CLUSTER)
#include "run.h"
#include "display.h"
#endif

#include "utils.h"
#include "view.h"

#if defined(_MPI)
  int isparallelrun=1;
#else
  int isparallelrun=0;
#endif

int Usage(char *); 
int CommandLineArgs(char *,int, char *argv[]);
int porous (scalar, face vector, int);

/**
  Global variables */

face vector muc[];
int maxlevel  = 6;
int minlevel = 6;
long int nbCells=0;
int refine=0;
int restoreDump=0;
int restartfrom1=0;
const face vector grav[] = {1, 0, 0};
scalar un[];

/**
  main function
*/
int main(int argc, char *argv[]){
  char *progname=  basename(argv[0]);
  if (CommandLineArgs(progname,argc,argv)) {
    return 0;
  }

  a = grav;
  mu = muc;

  N=1<<maxlevel; /* default discretization of the domain */
  NITERMIN = 2;
  NITERMAX = 300;
  //TOLERANCE=1e-4;
  DT=0.2;
  periodic(right);
  run();
}

int porous (scalar cs, face vector fs, int isadapted) {
  fprintf(stdout,"We consider %d maxlevel\n",maxlevel);
  vertex scalar phi[];

  if (!restoreDump) {
    if (isparallelrun) {
      fprintf(stdout,"We cannot use distance function in parallel\n");
      Usage("");
      return 1;
    }
    scalar d[];
    fprintf(stdout,"We load the stl file\n");
    coord * p = input_stl (fopen ("../outside_rot.stl", "r"));
    coord min, max;
    bounding_box (p, &min, &max);  
    double maxl = -HUGE;
    foreach_dimension()
      if (max.x - min.x > maxl)
    	maxl = max.x - min.x;
    size (maxl);

    X0 = (max.x + min.x)/2. - L0/2;
    Z0 = (max.z + min.z)/2. - L0/2;
    Y0 = (max.y + min.y)/2. - L0/2;

    distance (d, p);
  
    boundary ({d}); /* should not be necessary, as it's sequential */
    
    foreach_vertex()
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1] + d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
    
    boundary ({phi});
    char fdump[10];
    sprintf (fdump, "dump-%1d", maxlevel);
    dump(file = fdump, list={phi}); /* we save only phi */
    
  } else {
  
    char fdump[10];
    sprintf (fdump, "dump-%1d", maxlevel);
    fprintf(stdout,"In restore of %s dump file\n",fdump);
    if (!restore (file = fdump)){
      fprintf(stderr,"dump file for phi not present\n");
      fprintf(stderr,"Simulation aborded\n");
      fprintf(stderr,"\n");
      fflush(stderr);fflush(stdout);
      return 1;
    }  
    
    boundary ({phi}); /* initialise the ghost-cell values */
    
    /** 
    After, we load the rest for a restart at 1s 
    */
    if (restartfrom1) {
      fprintf(stdout,"In restore from 1s\n");
      if (!restore (file = "dump_at_1s", list={p,u.x,u.y,u.z,g.x,g.y,g.z})){
        fprintf(stderr,"dump file for simulation at 1s not present\n");
        fprintf(stderr,"Simulation aborded\n");
        fprintf(stderr,"\n");
        fflush(stderr);fflush(stdout);
        return 1;
      }
    }  
  }
  
  stats s = statsf (phi);
  fprintf (stdout, "phi max=%g min=%g\n", s.max, s.min);
  fflush(stdout);


  if (isadapted) {
    unrefine(level>=minlevel && fabs(phi[])>0.15);
    unrefine(level>=1 && phi[]<-0.1);
    //adapt_wavelet ({phi}, (double[]){1e-2}, maxlevel);
  }

  view (fov = 19.1765, camera="left", bg = {1,1,1},
        width = 600, height = 600);
  squares ("phi", linear = true, n={1,0,0}, alpha=0.);
  //cells(n={1,0,0}, alpha=0.);
  char name[80];
  
  if (isparallelrun) 
    sprintf (name, "phi-MPI-%d-%d.png", maxlevel,minlevel);
  else
    sprintf (name, "phi-SEQ-RESTART-%d-%d.png", maxlevel,minlevel);
  save (name);    


  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
  boundary ({cs});   

  
  update_perf();
  fprintf(stdout,"Total number of leaf cells : %ld\n",perf.tnc);
  nbCells=perf.tnc;
  fflush(stdout);

  scalar kappa[];
  curvature (cs, kappa);
  view (fov=40, theta = 0.5, phi = 0.5);
  draw_vof("cs", "fs", color = "kappa");
  box();
  sprintf (name, "kelvin-%d-%d.png", maxlevel,minlevel);
  save (name);

  clear();
  view (fov = 22, camera="left", bg = {1,1,1},
        width = 600, height = 600);
  //X=0
  squares ("cs", linear = true, n={1,0,0}, alpha=0.);
  cells(n={1,0,0}, alpha=0.);
  sprintf (name, "cs-plan-%d-%d.png", maxlevel,minlevel);
  save (name);
  
  return 0;
}

event init (i = 0) {

  if (porous(cs,fs,1)) return 1;

  u.n[embed] = dirichlet (0.);
  u.t[embed] = dirichlet (0.);
  u.r[embed] = dirichlet (0.);
  
}

event defaults (i++) {
  foreach_face()
    muc.x[] = 1.*fm.x[];
  boundary((scalar*){muc});
}

/**
  Get statistic on fluid field 
*/
event stats(i++){
  stats s = statsf (u.x);

  fprintf(stdout, "Stat 1 : %g %g %g %g %g %g %g\n", t, s.volume,
    normf(u.z).max, normf(u.z).avg, normf(p).max, normf(p).avg, s.sum/s.volume);
  fflush(stdout);
}

/*
event adapt (i++)
  adapt_wavelet ({cs,u}, (double[]){0.001,1e-7,1e-7,1e-7}, maxlevel);
*/

event init_un (i = 0) {
  foreach()
    un[] = u.x[];
}

event logfile (i++; i <= 500) {
  double avg = normf(u.x).avg, du = change (u.x, un)/(avg + SEPS);
  scalar div[];
  foreach() {
    div[] = 0.;
    /**
    if we are in the solid, we do not compute the divergence */
    if (cs[]>0.5 && cs[1]>0.5) {
      foreach_dimension() 
        div[] += u.x[1] - u.x[];
      div[] /= Delta;
    }
  }
  stats s0 = statsf (div);
  fprintf (stdout, "Stat 2 : %g %.9g %.9g %.9g \n", t, du, s0.sum/s0.volume, s0.max);
  fflush(stdout);
  
  /**
  If the relative change of the velocity is small enough we stop this
  simulation. */
  
  if (i>3 && du < 1e-3) {
  
    fprintf (stdout, "i=%d, t=%g, i=%g\n",i,t,du);
    view (fov = 32.2073, quat = {-0.309062,0.243301,0.0992085,0.914026},
	  tx = 0.0122768, ty = 0.0604286, bg = {1,1,1},
	  width = 600, height = 600);
    box();
    draw_vof("cs", "fs", fc = {0.5,0.5,0.5});
    char name[80];
    sprintf (name, "cs-%d-%d.png", maxlevel,minlevel);
    save (name);

    scalar nu[];
    foreach()
      nu[] = norm(u);
    stats s = statsf (u.x);
     
    view (fov = 25, camera="left", bg = {1,1,1},
	  width = 600, height = 600);
    squares ("nu", linear = true, n={1,0,0}, alpha=0.);
    cells(n={1,0,0}, alpha=0.);
    sprintf (name, "cross-%d-%d.png", maxlevel,minlevel);
    char str[99];
    sprintf (str, "U=%1.7f  NbCells=%ld", s.sum/s.volume, nbCells);
    draw_string (str, 1, lw = 3, lc = {1, 0, 1});
    save (name);    
        
    /*
    sprintf(name, "./field_time%g_at_%dmax-%dmin",t,maxlevel,minlevel);
    output_vtu((scalar *) {p,cs}, (vector *) {u}, name);
    */
    
    if (refine) {
      maxlevel++;
      refine=0;
      porous(cs,fs,1);
    } else {
      return 1; /* stop */
    }
      
  }
}

/*
event stop (t = 1.0){
  char name[80];
  sprintf(name, "./paraview/field_time%g_cpu",t);
  output_vtu((scalar *) {p,cs}, (vector *) {u}, name);
  p.nodump = false;
  dump(file = "dump_at_1s"); 
  if (refine) {
    maxlevel++;
    porous(cs,fs,1);
  } else {
    return 1; 
  }
}
*/

int Usage(char * prog) {
    fprintf(stdout,"Options suivantes :\n");
    fprintf(stdout,"   -maxlevel level : nb max level, must be >= to 6\n");
    fprintf(stdout,"   -minlevel level : nb min level\n");
    fprintf(stdout,"   -restoreDump    : load the levelset dump file for parallel runs\n");
    fprintf(stdout,"   -init level     : run maxlevel-1 for initialization of run at maxlevel level\n");
  return 1;
}

int CommandLineArgs(char * progname, int argc, char *argv[]) {
  int valRetour=0;
  /* En entree au clavier, pratique ! */
  while (--argc>0 ) {
    if ( (*++argv) [0] == '-' ){
      if ( !strcmp (*argv, "-h")){
        valRetour=Usage(progname);
      }else if ( !strcmp (*argv, "-maxlevel")){   
        maxlevel = atoi(*++argv);
        fprintf (stdout,"maxlevel set to : %d\n",maxlevel);
	fflush(stdout);
        argc--;
      }else if ( !strcmp (*argv, "-minlevel")){   
        minlevel = atoi(*++argv);
        fprintf (stdout,"minlevel set to : %d\n",minlevel);
	fflush(stdout);
        argc--;
      }else if ( !strcmp (*argv, "-init")){   
        maxlevel = atoi(*++argv);
	maxlevel--;
	refine=1;
        fprintf (stdout,"We compute first with level%d\n",maxlevel);
	fflush(stdout);
        argc--;
      }else if ( !strcmp (*argv, "-restoreDump")){   
        restoreDump = 1;
        fprintf (stdout,"For parallel, we load the precalculated levelset function\n");
	fflush(stdout);
      }else {
       printf("Option:%s non reconnu\n\n",*argv);
       valRetour=Usage(progname);
     }
    }  /* fin option */
  }    /* fin while */
  return valRetour;
}
