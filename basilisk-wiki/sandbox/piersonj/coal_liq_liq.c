/**
Coalescence of a drop at a liquid-liquid interface : comparison with Mohamed-Kassim & LongMire (2003) */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
#include "tag.h"
#include "output_mpi.h"



#define LEVEL 12
#define INLEVEL 7
#define MINLEVEL 7

/**Fluids Parameters */

/* Initial Diameter */
#define D0 1.
#define g 10.
#define rhoInt 1

/* Dimensionless parameters as in Mohamed-Kassim & LongMire (2003)*/
#define Re 71.022
#define muRatio 0.3315
#define rhoRatio 0.841
#define We 6.4018

/*Characteristic parameters*/
#define uC sqrt((1.-rhoRatio)*g*D0)
#define tC D0/uC
#define tEnd 22*tC


/**
Size of the domain .*/
#define L0 24*D0
/* Initial drop distance from the west boundary*/
#define L 4*D0
/* Initial interface distance from the west boundary */
#define Li 20*D0


/**
The boundary conditions are slip walls*/


int main() {

  /**
  The domain */
  size (L0);
  init_grid (1 << INLEVEL);
  
  /* Dimensionless criterion for mass conservation */
  TOLERANCE = 1e-5;

  /* We use the reduced gravity approach (advantage -> exact hydrostatic balance) */
  G.x = -g;


 /*
  We define the parameters of the fluids. The internal fluid (the oil) is 
  the fluid 1. The external fluid (the water) is the fluid 2.
  All fluids parameters are define at the begining of this source code */

  rho1 = rhoInt;
  rho2 = rhoInt/rhoRatio;
  mu2 = D0*uC*rho2/(Re); 
  mu1 = mu2*muRatio;
  f.sigma = D0*uC*uC*rho2/We;
    
  run();
}


event init (t = 0) {
    if (!restore (file = "dump-")){

       /**
       we refine the mesh locally to the maximum level, in a sphere of diameter 2D around the bubble */
       refine (sq(x-L) + sq(y) - sq(D0) < 0 && level < LEVEL);

       /**
       We initialise the geometry of the drop...*/
       fraction(f, -sq(x-L)-sq(y)+sq(D0/2));
       
       /**
       ... and the pool.*/      
       scalar p[];
       fraction(p, x-Li);
       foreach()
         f[] += p[];
    }
}


/**
Adaptative mesh refinement.
We adapte the mesh by controling the error of volume fraction and velocity field*/
event adapt (i++) {
  double uemax = 1e-4*uC;
  adapt_wavelet ({f,u}, (double[]){0.001,uemax,uemax,uemax}, LEVEL, MINLEVEL); 

}


/*
We output restart files, interface and fields . */

event saveDatas (t += 0.2*tC; t <= tEnd) {

  //We output the interface
  char name[80];
  sprintf (name, "interface/interface-%g-%d.txt", t, pid());
  FILE* fp = fopen (name, "w");
  output_facets (f, fp);
  fclose (fp);
 
  // We output the Vorticity field in binary files
  scalar vort[];
  vorticity (u, vort);
  char namev[80];
  sprintf (namev, "vort/bin-vort-%g", t);
  FILE* fvor = fopen(namev,"w");
  output_matrix_part_mpi(vort,fvor,n=2400,linear =true,box = {{0.,0.},{L0,2*D0}});
  fclose(fvor);

  // binary file , dump
  char named[80];
  sprintf (named, "dump/dump-%d", i);
  dump (file = named);

}


event logfile(i++) {
  printf ("i = %d t = %g\n", i,t);
  fflush(stdout);
  
}


/*
Output position, volume, velocity ... of the drop.
We do not save it during each iterations because the foreach loops can be expensive*/
event posVit (i+=10) {

  //We output the drop and pool interfaces
  //char nameInt[80];
  //sprintf (nameInt, "./interfacedp/drop-interface-%d-%d.txt", i, pid());
  //FILE* fpInt = fopen (nameInt, "w");
  //output_facets (f, fpInt);
  //fclose (fpInt);

  /*
  Tag each drop which are separated*/
  scalar m[];
  foreach()
    m[] = f[] > 1e-6;
  int n = tag (m);

  vector h[];
  heights (f, h);

  /*
  Velocity, position ... of the drop */
  double ud = 0., xd = 0., vd = 0., ecd = 0.;
  double xdMax = -HUGE;;
  double ydMax = -HUGE;;
  double xdMin = +HUGE;;
  double xpMin = -HUGE;;


  // we assume that the tag of the drop is larger than 0
  bool testDrop=true;
  int j=1;

  /* while loop to save only information related to the drop
     since the reduction operators is limited to 10
     could be improved since it can do two foreach loop */
  while (testDrop && j<n+1)
  {

    if (j > 1)
    {
      double ud = 0., xd = 0., vd = 0., ecd = 0.;
      double xdMax = -HUGE;;
      double ydMax = -HUGE;;
      double xdMin = +HUGE;;
      double xpMin = -HUGE;;
    }

    // Reduction operations does not exist yet for arrays, but here we use reduction for scalars 
    foreach(reduction(+:ud) reduction(+:xd) reduction(+:vd) reduction(+:ecd)
            reduction(max:xdMax) reduction(max:ydMax) reduction(min:xdMin) reduction(max:xpMin)) 
    {

      // we assume that the drop is the fluid tag j
      if (m[] ==j) 
      {

        double dv = f[]*dv();
        vd += dv;
        xd += x*dv;
        ud += u.x[]*dv;
        ecd += 0.5*rho1*u.x[]*u.x[]*dv;
        
        // drop x max position
        if (h.x[] != nodata) 
        {
          double xj = x + height(h.x[])*Delta; 
          if (xj > xdMax)
            xdMax = xj;
        }
        // drop x min position
        if (h.x[] != nodata) 
        {
          double xjm = x + height(h.x[])*Delta; 
          if (xjm < xdMin)
            xdMin = xjm;
        }
        // drop y max
        if (h.y[] != nodata) 
        {
          double yj = y + height(h.y[])*Delta;
          //double xyj = x + height(h.x[])*Delta;  
          if (yj > ydMax)
            ydMax = yj;
            //xydMax = xyj;
        }
         
      }

      // pool x min position (max position for the classical void fraction)
      if (h.x[] != nodata) 
      {
        double xm = x + height(h.x[])*Delta; 
        if (xm > xpMin)
          xpMin = xm;
      }

    }


    if(vd < 0.1)
    {
      testDrop=false;
    } 
    else
    {
      j+=1;
    }

  }


  //fprintf(stderr, "%d %g %g %g %g %g %g %g %g %g\n", i, t, vd, xd/vd, ud/vd, ecd/vd/vd, xdMin, xdMax, ydMax, ydMax/(xdMax-xdMin));
  //fflush (stderr);

  char namePosVitd[50];
  sprintf (namePosVitd, "pos_velocity_drop.txt");
  static FILE * fnamePosVitd = fopen (namePosVitd, "w");
  fprintf (fnamePosVitd, "%g %g %g %g %g  %g %g %g %g\n", t, vd, xd/vd, ud/vd, ecd/vd/vd, xdMin, xdMax, ydMax, ydMax/(xdMax-xdMin));
  //fclose (fnamePosVitd);
  fflush (fnamePosVitd);

  char dist[50];
  sprintf (dist, "dist.txt");
  static FILE * fdist = fopen (dist, "w");
  fprintf (fdist, "%g %g\n", t, fabs(xdMax - xpMin));
  //fclose (fdist);
  fflush (fdist);

}

event movies (t += 0.2*tC; t <= tEnd){

  scalar vort[];
  vorticity (u, vort);
  scalar m[];
  foreach() {
    	if (f[]<0.90 && f[] >0.10)
     	 m[] = -1;
    	else
    	  m[] = 0;}

  char name_vort[100];
  sprintf (name_vort, "vort/vort-%g.png", t);// save of the picture 
  FILE* fvort = fopen (name_vort, "w");
  double om = 1;
  output_ppm (vort, fvort, min=-om,max=om,mask=m, map=cool_warm, n = 1000,box = {{0.,0.},{L0,2*D0}});
  fclose (fvort);

  scalar le[];
  foreach(){
    le[] = level;
  }
  static FILE * lev = fopen ("grid.gif", "w");
  output_ppm (le, lev,min=MINLEVEL,max=LEVEL,n = 1000,box = {{0.,0.},{L0,2*D0}});


}


///////////////////////////////////////////////// End fonction ////////////////////////////
event end (t = tEnd) {}

