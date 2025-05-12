/**
# Bubble rise and coalescence (film drainage only) in a Newtonian surrounding

This code includes a film identification method.
*/
  
#include "axi.h"         // Axisymmetric simulation
#include "navier-stokes/centered.h"
#include "vof.h"         // The header two-phase could be used instead
#include "tension.h"
#include "view.h"

/**Here, we define the dimensionless numbers governing the flow:*/

#define Ar 25.0         // Archimedes number
#define Bo 5.0           // Bond Number
#define MUR 100.          // Viscosity Ratio
#define RHOR 100.         // Density Ratio

/** Some flow geometries*/
/** $Ar = \sqrt{\rho_1|\Delta \rho|g D^3} / \mu_1$; $Bo = |\Delta \rho|g D^2 / \sigma$; $\mu_r = \mu_1 / \mu_2$; $\rho_r = \rho_1 / \rho_2$*/
/** The indexes 1 and 2 are for the surrounding and bubble phases, respectively.*/

#define R0 0.5           // Drop Radius
#define HEIGHT 50*R0     // Domain Size
#define Xo 20.*R0        // Drop Initial Position
#define Xi 30.*R0         // Interface Position

/** The variables of the flow */
scalar rhov[];
face vector alphav[], muv[], muv2[];

/** Two volume fraction fields are used, one for the bubble and one for the top layer interface. This is done to facilitate the film identification later*/
scalar f1[], f2[], * interfaces = {f1,f2};

double xdd = 0.;         // Drop center of mass position
double vdd = 0.;         // Drop center o mass velocity
double uddx = 0.;        // Drop volume

double rho1, rho2, mu1, mu2, U;

/** By balacing the bouyant ($\tau_b = |\Delta\rho| g D$) and inertial ($\tau_i = \rho_1 U^2$) stress, the characteristic velocity $U = \sqrt{|\Delta\rho| g D / \rho_1}$ is obtained.
*/

/* Simulations time */
double tEnd = 6.81;

/* Mesh refinement levels */
int LEVEL = 10;
int MINLEVEL = 6;


int main ()
{
  size (HEIGHT);
  init_grid (1 << MINLEVEL);

  rho1 = 1.;
  rho2 = rho1/RHOR;
  U = sqrt(rho1 - rho2);
  mu1 = U/Ar;
  mu2 = mu1/MUR;  
  f1.sigma = f2.sigma = (rho1 - rho2)/Bo;

  alpha = alphav;
  rho = rhov;
  mu = muv;

  TOLERANCE = 1e-4;  // Tolerance
  NITERMAX = 100;    // Maximum number of iterations
  run();
}


/** Gravity acceleration */
event acceleration (i++)
{
  face vector av = a;
  foreach_face(x)
    av.x[] = -1.;
}


/** Initialization */
event init (t = 0)
{
  if (!restore (file = "restart"))
  {
    refine (sq(x - Xo) + sq(y) - sq(1.2*R0) < 0 && level < LEVEL);
    fraction (f1, -sq(x - Xo) - sq(y) + sq(R0));

    fraction (f2, x - Xi);
  }
}


/** Calculation of the vicosity and density in each phase */
event properties (i++) 
{
  scalar eta[];    // Viscosity function
  scalar fa1[], fa2[];     // For a smoothing of the viscosity jump
  foreach()
  {
    fa1[] =  (4.*f1[] + 2.*(f1[-1,0] + f1[1,0] + f1[0,-1] + f1[0,1]) + f1[1,1] + f1[-1,1] + f1[1,-1] + f1[-1,-1])/16.;
    fa2[] =  (4.*f2[] + 2.*(f2[-1,0] + f2[1,0] + f2[0,-1] + f2[0,1]) + f2[1,1] + f2[-1,1] + f2[1,-1] + f2[-1,-1])/16.;
   
    eta[] = mu1;
    
    rhov[] = cm[]*((1 - f1[] - f2[])*rho1 + f1[]*rho2 + f2[]*rho2);
  }
//  boundary ({fa, eta, rhov});

// Here we pass the viscosity, specific volume and density.
  foreach_face()
  {
    double fm1 = (fa1[] + fa1[-1])/2.;
    double fm2 = (fa2[] + fa2[-1])/2.;
    double etam = (eta[] + eta[-1])/2.;   // Viscosity is passed at the faces of the cells
    
    muv.x[] = fm.x[] / ((1. - fm1 - fm2)/etam + fm1/mu2 + fm2/mu2);
    alphav.x[] = fm.x[]/ ((1 - fm1 - fm2)*rho1 + fm1*rho2 + fm2*rho2);
  }
//  boundary ({muv,alphav,rhov});
}




// Print
/*
event logfile(i+=1)
{
  scalar le[];
  int j = 0;
  foreach(reduction(+:j))
  {
    le[] = level;
    j++;
  }
  stats s = statsf(le);

  printf ("i = %d t = %g, # = %d, l.mi = %g, l.ma = %g\n", i, t, j, s.min, s.max);
  fflush(stdout);
}*/



/** This part is used to identify the film region and to apply a mesh refinement there*/

/** My restriction and prolongation function, which forces the variable to always be equal to zero in the mesh next and prior refinement levels*/
static inline void myrestrict(Point point, scalar s){
  s[] = 0;
}

static inline void myprolongation(Point point, scalar s){
  foreach_child()
    s[] = 0;
}



double thick_ms = 0.3;
scalar sign6[];
int q = 0;
event thin_film (i++)
{
  scalar drox[], droy[], intx[], inty[];
  position (f1, drox, {1,0});
  position (f1, droy, {0,1});
  position (f2, intx, {1,0});
  position (f2, inty, {0,1});
  double posit1ma, posit2ma, posityma, posit1mi;
  posit1ma = statsf(drox).max;
  posit2ma = statsf(intx).max;
  posityma = statsf(droy).max;
  posit1mi = statsf(drox).min;

  if (posit1ma > Xi || q == 1)
  {
    q = 1;
    int i, j;
    i = j = 0;

  //reduction(+:i) reduction(+:j)
    foreach(serial)
    {
      if (drox[] < 90. && drox[] > xdd*vdd-1.)
      {
        i++;
      }

      if (inty[] < 2. && inty[] > 0.)
      {
        j++;
      }
      sign6[] = 0.;
    }


    double dropx[i], interx[j];
    double dropy[i], intery[j];
    for (int k = 0; k < i; k++)
    {
      dropx[k] = 0.;
      dropy[k] = 0.;
    }
    for (int k = 0; k < j; k++)
    {
      interx[k] = 0.;
      intery[k] = 0.;
    }
    i = j = 0;
  
    foreach(serial)
    {
      if (drox[] < 90. && drox[] > xdd*vdd-1.)
      {
        dropx[i] = drox[];
        dropy[i] = droy[];
        i++;
      }

      if (inty[] < 2. && inty[] > 0.)
      {
        interx[j] = intx[];
        intery[j] = inty[];
        j++;
      }
    }


    double max_film = 0.3;
    double thick, thick_m;
    thick_ms = max_film;

    for (int k = 0; k < i; k++)
    {
      thick_m = max_film;
      for (int o = 0; o < j; o++)
      {
        if (dropx[k] > 0. && interx[o] > 0.)
        {
          thick = sqrt( sq(dropx[k] - interx[o]) + sq(dropy[k] - intery[o]) );
          if(thick < max_film && thick < thick_m)
          {
            thick_m = thick;
          }
          if(thick < max_film && thick < thick_ms)
          {
            thick_ms = thick;
          }
        }
      }
    } 

    double thin_max = max_film;
    thin_max = thick_ms*5. < max_film ? thick_ms*3 : max_film;


    if (thick_ms < thin_max)
    {  
      foreach(serial)
      {
        double s4 = 0.;
        double s5 = 0.;

        for (int k = 0; k < i; k++){
            {
              double thin1 = sqrt( sq(dropx[k] - x) + sq(dropy[k] - y));
              if(thin1 < thin_max && (f1[]+f2[]) < 0.99 && x > Xi)
              {
                s4 = 1.;
                break;
              }
            }
        }
      

        for (int k = 0; k < j; k++){
              double thin2 = sqrt( sq(interx[k] - x) + sq(intery[k] - y));
              if(thin2 < thin_max && (f1[]+f2[]) < 0.99 && x > Xi)
              {
                s5 = 1.;
                break;
              }
        }
 
        if ((s4 + s5) == 2.)
          sign6[] = 1.;
        else
          sign6[] = 0.; 
      }
    }
  //boundary((scalar *){sign6});


    double num_cell = 3;
    if (thick_ms < max_film)
      LEVEL = (int)(log(num_cell*HEIGHT/thick_ms)/log(2));
    if (LEVEL <= 10)
      LEVEL = 10;
    if (LEVEL > 13)
      LEVEL = 13;

    sign6.prolongation = myprolongation;
    sign6.restriction = myrestrict;
  }


  /** Calculating drop position and velocity */

  double xd1 = 0.;
  double vd1 = 0.;
  double udx1 = 0.;
  foreach(reduction(+:vd1) reduction(+:xd1) reduction(+:udx1))
  {
     vd1 += dv()*f1[];
     xd1 += dv()*f1[]*x;
     udx1 += u.x[]*dv()*f1[];  
  } 

  xdd = xd1;
  vdd = vd1;
  uddx = udx1;

  char nameu[50];
  sprintf (nameu, "pos.txt");
  static FILE * fname = fopen (nameu, "w");
  fprintf (fname, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n", t, vd1, xd1/vd1, udx1/vd1, posit1mi, posit1ma, posit2ma, posityma, thick_ms);
  fflush (fname);
}


/** Mesh adaptation  */
event adapt (i++)
{
    adapt_wavelet ({f1, f2, u, sign6}, (double[]){1e-3, 1e-3, 1e-3, 1e-3, 1e-1}, LEVEL, MINLEVEL);
}


/** Here, we identify the film region again, but this time to write a file with the variables in the film. We also write for the fim shape*/
/*
event outputs_fields (t = 0; t <= tEnd; t += 0.05)
{
  scalar drox[], droy[], intx[], inty[];

///We identify the interface positon on the interfacil cells and it pass to a array.

  position (f1, drox, {1,0});
  position (f1, droy, {0,1});
  position (f2, intx, {1,0});
  position (f2, inty, {0,1});
  double posit1ma, posit2ma;
  posit1ma = statsf(drox).max;
  posit2ma = statsf(intx).max;


if(thick_ms < 0.3)
{

  int i, j;
  i = j = 0;

  foreach(serial)
  {
    if (drox[] < 90. && drox[] > 0. && x > xdd/vdd-0.1)
    //if (drox[] < 10. && drox[] > 0.)
    {
      i++;
    }

    if (inty[] < 2. && inty[] > 0.)
    {
      j++;
    }
  }

    double dropx[i+1], interx[j+1];
    double dropy[i+1], intery[j+1];
    for (int k = 0; k <= i; k++)
    {
      dropx[k] = 0.;
      dropy[k] = 0.;
    }
    for (int k = 0; k <= j; k++)
    {
      interx[k] = 0.;
      intery[k] = 0.;
    }

    i = j = 0;

    foreach(serial)
    {
      if (drox[] < 90. && drox[] > 0. && x > xdd/vdd-0.1)
      {
        dropx[i] = drox[];
        dropy[i] = droy[];
        i++;
      }

      if (inty[] < 2. && inty[] > 0.)
      {
        interx[j] = intx[];
        intery[j] = inty[];
        j++;
      }
    }


    double xx, yy;
    interx[0] = posit2ma;
    intery[0] = 0.0;
    for(int k = 0.; k  <= j; k++)
    {
      for(int o = k+1; o <= j; o++)
      {
        if(interx[k] < interx[o])
        {
          xx = interx[k];
          yy = intery[k];
          interx[k] = interx[o];
          intery[k] = intery[o];
          interx[o] = xx;
          intery[o] = yy;
        }
      }
    }



    dropx[0] = posit1ma;
    dropy[0] = 0.0;
    for(int k = 0.; k  <= i; k++)
    {
      for(int o = k+1; o <= i; o++)
      {
        if(dropx[k] < dropx[o])
        {
          xx = dropx[k];
          yy = dropy[k];
          dropx[k] = dropx[o];
          dropy[k] = dropy[o];
          dropx[o] = xx;
          dropy[o] = yy;
        }
      }
    }




    double max_film = 0.3;
    double thick, thick_m, thick_max;
    double thick_p[i+1];
    thick_ms = max_film;
    thick_max = 0.;

    char film_t[50];
    sprintf (film_t, "film-thickness-%.2f.txt", t);
    FILE *filmt = fopen(film_t, "w");  

    double ss[i];
    ss[0] = 0.;
    thick_p[0] = posit2ma - posit1ma;
    fprintf(filmt, "%g %g %g %g\n", ss[0], thick_p[0], dropx[0], dropy[0]);
    for(int k=1; k<i; k++) 
    {
      thick_m = max_film;
      thick_p[k] = max_film;
      for (int o = 0; o < j; o++)
      {
        if (dropx[k] > 0. && interx[o] > 0.)
        {
          thick = sqrt( sq(dropx[k] - interx[o]) + sq(dropy[k] - intery[o]) );
          if(thick < max_film && thick < thick_m)
          {
            thick_m = thick;
            thick_p[k] = thick;
          }
          if(thick < max_film && thick < thick_ms)
          {
            thick_ms = thick;
          }
          if(thick < max_film && thick > thick_max)
          {
            thick_max = thick;
          }
        }
      }
      ss[k] = sqrt( sq(dropx[k] - dropx[k-1]) + sq(dropy[k] - dropy[k-1]) ) + ss[k-1];
      fprintf(filmt, "%g %g %g %g\n", ss[k], thick_p[k], dropx[k], dropy[k]);     // Film shape
    }
    fclose(filmt);


  //  Field variables in the film
    char field[80];
    sprintf (field, "fields-film-%.2f.txt", t);
    FILE* fld = fopen (field, "w");
  
    fprintf(fld, "%g 0.0 %g %g %g\n", posit1ma, interpolate(u.x, posit1ma, 0), interpolate(u.y, posit1ma, 0), interpolate(p, posit1ma, 0));
    fprintf(fld, "%g 0.0 %g %g %g\n", posit2ma, interpolate(u.x, posit2ma, 0), interpolate(u.y, posit2ma, 0), interpolate(p, posit2ma, 0)); 
    foreach(serial)
    {
      if (sign6[] == 1.)
      {
        fprintf(fld, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]); 
     }
    }
    fclose (fld);
  }
}*/





/** Generating the bubble shape and variable field */
/*
event outputs (t += 1.0)
{

  char namei1[80];
  sprintf (namei1, "interface-bubble-%.2f.txt", t);
  FILE* fp1 = fopen (namei1, "w");
  output_facets (f, fp1);
  fclose (fp1);

  char field2[80];
  sprintf (field2, "fields-output-%.2f.txt", t);
  FILE* fld2 = fopen (field2, "w");
  output_field ((scalar *){u, ,f}, fld2, n = 500, linear = true,  box = {{xdd/vdd - 5*R0, 0.},{xdd/vdd + 5*R0, 5*R0}});
  fclose (fld2);
}
*/


/** Backup */
/*
event snapshot (t = 0; t <= tEnd; t += 1.0)
{
  char named[80];
  sprintf (named, "dump-%.2f", t);
  dump (file = named);
}
*/


event movie_mesh (t+=0.02)
{
  clear(); 
  view(fov = 2, tx = 0., ty = -0.4, psi = -pi/2);
  
  translate (x =  Xo - xdd/vdd)
  {
    draw_vof ("f1");
    draw_vof ("f2");
    cells();
    mirror (n = {0,-1})
            {
            draw_vof ("f1");
            draw_vof ("f2");
            cells();
          }
  }
  save("movie_mesh.mp4");

  translate (x =  Xo - xdd/vdd)
  {
    draw_vof ("f1");
    draw_vof ("f2");
    squares("sign6", linear = true);
    mirror (n = {0,-1})
          {
            draw_vof ("f1");
            draw_vof ("f2");
            squares("sign6", linear = true);
          }
  }
  save("movie_film.mp4");
}

event movie_vof (t+=0.02)
{
  scalar ff[];
  foreach()
  {
    ff[] = f1[] + f2[];
  }

  clear(); 
  view(fov = 6.0, tx = 0., ty = -0.5, psi = -pi/2);
   
  draw_vof ("f1");
  draw_vof ("f2");
  squares("ff", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f1");
          draw_vof ("f2");
          squares("ff", linear = true);
        }
  save("movie_vf.mp4");
}



event end (t = tEnd) {}


/**
## Results
### Volume fraction field
![Volume fraction field](Bubble-rise-coalescence-Newtonian/movie_vf.mp4)

### Mesh
![Mesh](Bubble-rise-coalescence-Newtonian/movie_mesh.mp4)

### Film region
![Film region](Bubble-rise-coalescence-Newtonian/movie_film.mp4)
*/

/**
~~~gnuplot Velocity profile: comparison between numerical and analytical solutions
reset
set title "Center of mass velocity"
set xlabel 'Time'
set ylabel 'Velocity'
set xrange [0:12]
set yrange [0:1]
plot 'pos.txt' u 1:4 w l lw 3 lc 7 notitle
~~~
*/