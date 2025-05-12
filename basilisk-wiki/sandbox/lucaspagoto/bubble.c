/**
# Bubble Coalescence

My thin structure function using the restriction and prolongation (it works) and new interface profile method (it also workrs), with drop[i][2] and inter[i][2]

*/


#include "axi.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"
#include "view.h"

//#include "tracer.h"
#include "diffusion.h"
scalar cf[];
mgstats mgd1, mgd2;
face vector Diff[];
scalar r1[], beta1[], theta1[];
//scalar * tracers = {cf};

#define Ar 5.0
#define Bo 5.0

#define MUR 10.
#define RHOR 10.

#define R0 0.5
#define HEIGHT 20*R0
#define X0 5*R0
#define MPLD 7.5*R0

scalar rhov[], rhov2[];
face vector alphav[], muv[], muv2[];
scalar f1[], f2[], * interfaces = {f1, f2};
scalar Lambda1[], sour[];

double tEnd = 10.01;

double xdd = 0.;
double vdd = 0.;
double uddx = 0.;

double rho1, rho2, mu1, mu2, U;

int LEVEL = 9;
int MIDLEVEL = 6;
int MINLEVEL = 5;



static inline void myrestrict(Point point, scalar s){
  s[] = 0;
}

static inline void myprolongation(Point point, scalar s){
  foreach_child()
    s[] = 0;
}




int main ()
{
  size (HEIGHT);
  init_grid (1 << MIDLEVEL);

  rho1 = 1.;
  mu1 = sqrt(rho1 - rho2)/(Ar);
  rho2 = rho1/RHOR;
  mu2 = mu1/MUR;
  U = sqrt(rho1 - rho2);
  f1.sigma = f2.sigma = (rho1 - rho2)/Bo;

  alpha = alphav;
  rho = rhov;
  mu = muv;

//  f1.tracers = {cf};
//  cf.inverse = true;

  TOLERANCE = 1e-4;
  NITERMAX = 100;
  run();
}



event acceleration (i++)
{
  face vector av = a;
  foreach_face(x)
    av.x[] = -1.;
}




event init (t = 0)
{
  if (!restore (file = "restart"))
  {
    refine (sq(x - X0) + sq(y) - sq(1.2*R0) < 0 && level < LEVEL);
    fraction (f1, -sq(x - X0) - sq(y) + sq(R0));
 
    fraction (cf, -sq(x - X0) - sq(y) + sq(R0));

//    foreach()
//      cf[] = f1[]*0.5;
    boundary ((scalar *){cf});

    scalar cfconcentration[];
    foreach()
      cfconcentration[] = (f1[] > 1e-10 ? cf[]/f1[] : 0.);

    fraction (f2, x - MPLD);
  }
}




static scalar * interfaces1 = NULL;
event vof (i++)
{
  cf.gradient = zero;  // where is this used?
  cf.refine = refine_linear;
  cf.restriction = restriction_volume_average;

  f1.tracers = (scalar *){cf};
  vof_advection ({f1}, i);

  interfaces1 = interfaces, interfaces = NULL;
}


event tracer_advection (i++) {
  interfaces = interfaces1;
}


event tracer_diffusion (i++)
{
  foreach_face()
    Diff.x[] = 0.*fm.x[];
  boundary ((scalar *){Diff});

  double t_eq = 1e-1;
  double cf_eq = 0.5;
  foreach()
  {
    beta1[] = (-1/(cf_eq*t_eq))*cm[];
    r1[] = (1/t_eq)*cm[];
    theta1[] = rho1*cm[];
  }
  boundary ((scalar *){beta1, r1, theta1});
 
  foreach()
    cf[] = (f1[] > 1e-10 ? cf[]/f1[] : 0.);
  boundary({cf});
  
//  int j = 10;
//  double dt_d = dt/j;
//  printf("dt_d = %g dt = %g\n", dt_d, dt);

//  for(int k = 0; k < j; k++)
  mgd1 = diffusion (cf, dt, Diff, r = r1, beta = beta1, theta = theta1);

  foreach()
    cf[] *= f1[];
  boundary({cf});
}








event properties (i++) 
{
// Here we pass the viscosity, specific volume and density. I could also include smearing of the jump.
  foreach_face()
  {
    double fm1 = (f1[] + f1[-1])/2.;
    double fm2 = (f2[] + f2[-1])/2.;
    
    muv.x[] = fm.x[] / ( (1. - fm1 - fm2)/mu1 + fm1/mu2 + fm2/mu2  );
    
    alphav.x[] = fm.x[]/((1 - fm1 - fm2)*rho1 + fm1*rho2 + fm2*rho2);
  }

  foreach()
  {
    rhov[] = cm[]*((1 - f1[] - f2[])*rho1 + f1[]*rho2 + f2[]*rho2);
  }
  boundary ({muv,alphav,rhov});
}





double thick_ms = 0.;
scalar le[];
event logfile(i+=1)
{
  int j = 0;
  foreach()
  {
    le[] = level;
    j++;
  }
  stats s = statsf(le);

  printf ("i = %d t = %g, # = %d, l.mi = %g, l.ma = %g fm = %g\n", i, t, j, s.min, s.max, thick_ms);
  fflush(stdout);
}





scalar drox[], droy[], intx[], inty[];
scalar sign4[], sign5[], sign6[];
event thin_film (i++)
{

/**
We identify the interface positon on the interfacil cells and it pass to a array.
*/
  position (f1, drox, {1,0});
  position (f1, droy, {0,1});
  position (f2, intx, {1,0});
  position (f2, inty, {0,1});

  int i, j;
  i = j = 0;

  foreach(reduction(+:i) reduction(+:j))
  {
    if (drox[] < 10. && drox[] > 0.)
    {
      i++;
    }

    if (inty[] < 5. && inty[] > 0.)
    {
      j++;
    }
    sign6[] = 0.;
  }



  double drop[i][2], inter[j][2];
  for (int k = 0; k < i; k++)
  {
    drop[k][0] = 0.;
    drop[k][1] = 0.;
  }
  for (int k = 0; k < j; k++)
  {
    inter[k][0] = 0.;
    inter[k][1] = 0.;
  }


  i = j = 0;

  foreach(reduction(+:i) reduction(+:j) serial)
  {
    if (drox[] < 10. && drox[] > 0.)
    {
      drop[i][0] = drox[];
      drop[i][1] = droy[];
      i++;
    }

    if (inty[] < 5. && inty[] > 0.)
    {
      inter[j][0] = intx[];
      inter[j][1] = inty[];
      j++;
    }
  }

/*
  char drops[50];
  sprintf (drops, "drops.txt");
  FILE *dropp = fopen(drops, "w");

  char inters[50];
  sprintf (inters, "inter.txt");
  FILE *interr = fopen(inters, "w");

  for (int k = 0; k <= i; k++)
  {
    fprintf(dropp, "%g %g\n", drop[k][0], drop[k][1]);
  }
  for (int k = 0; k <= j; k++)
  {
    fprintf(interr, "%g %g\n", inter[k][0], inter[k][1]);
  }
  fclose(dropp);
  fclose(interr);
*/



/**
Here we write the film profile in a file
*/
  char film[50];
  sprintf (film, "film.txt");
  FILE *filmm = fopen(film, "w");  

  double max_film = 0.3;
  double thick, thick_m;
  thick_ms = max_film;
  int q = 1;

  for (int k = 0; k < i; k++)
  {
    thick_m = max_film;
    for (int o = 0; o < j; o++)
    {
      if (drop[k][0] > 0. && inter[o][0] > 0.)
      {
        thick = sqrt( sq(drop[k][0] - inter[o][0]) + sq(drop[k][1] - inter[o][1]) );
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
    if (thick_m < max_film && drop[k][0] > 0.)
    {
      fprintf(filmm, "%d %g %g %g\n",q ,thick_m, drop[k][0], drop[k][1]);
      q++;
    }
  }
  fclose(filmm);



  /**
  Here we identify the thin structure region
  */

  double ymax;
  ymax = statsf(droy).max;
  double thin_max = max_film;
  thin_max = thick_ms*2 < max_film ? thick_ms*2 : max_film;
  double thin1, thin2;
  
  if (thick_ms < thin_max)
  {  
    foreach(serial)
    {
      sign4[] = 0.;
      sign5[] = 0.;

      for (int k = 0; k < i; k++)
      {
          {
            thin1 = sqrt( sq(drop[k][0] - x) + sq(drop[k][1] - y));
            if(thin1 < thin_max && (f1[]+f2[]) < 0.99 && y < 1.1*ymax)
              sign4[] = 1.;
          }
      }

      for (int k = 0; k < j; k++)
      {
          {
            thin2 = sqrt( sq(inter[k][0] - x) + sq(inter[k][1] - y));
            if(thin2 < thin_max && (f1[]+f2[]) < 0.99 && y < 1.1*ymax)
              sign5[] = 1.;
          }
      }
 
      if ((sign4[] + sign5[]) == 2.)
        sign6[] = 1.;
      else
        sign6[] = 0.; 
    }
  }
  boundary((scalar *){sign6});

   /**
  May I do this?
   */
  double num_cell = 10;
  if (thick_ms < max_film)
    LEVEL = (int)(log(num_cell*HEIGHT/thick_ms)/log(2));
  if (LEVEL > 12)
    LEVEL = 12;

  sign6.prolongation = myprolongation;
  sign6.restriction = myrestrict;
}






//vector h1[], h2[];
event dropData (i++)
{
  double xd1 = 0.;
  double vd1 = 0.;
  double udx1 = 0.;
  double vdc = 0.;

  foreach(reduction(+:vd1) reduction(+:xd1) reduction(+:udx1) reduction(+:vdc))
  {
     vd1 += dv()*f1[];
     vdc += dv()*cf[];
     xd1 += dv()*f1[]*x;
     udx1 += u.x[]*dv()*f1[];
  } 

//  heights (f1, h1);
//  heights (f2, h2);

  xdd = xd1;
  vdd = vd1;
  uddx = udx1;


  char nameu2[50];
  sprintf (nameu2, "pos.txt");
  static FILE * fname2 = fopen (nameu2, "w");
  fprintf (fname2, "%.8f %.8f %.8f %.8f %.8f %.8f\n", t, vdc, vd1, xd1/vd1, udx1/vd1, thick_ms);
  fflush (fname2);
}





#if 0
vector uu[];
event outputs (i++)
{
  foreach()
  {
    uu.x[] = (sign6[]  == 1. ? u.x[] : nodata);
    uu.y[] = (sign6[]  == 1. ? u.y[] : nodata);
  }
}





event outputs2 (t = 0; t <= tEnd; t += 0.5)
{
  char field2[80];
  sprintf (field2, "fields-output-%.2f.txt", t);
  FILE* fld2 = fopen (field2, "w");
  
  output_field ((scalar *){uu, p, sign6}, fld2, n = 500, linear = true,  box = {{3,0},{4,1}});

//  foreach()
//    fprintf(fld2, "%g %g %g %g %g\n", x, y, uu.x[], uu.y[], p[]);
  
  fclose (fld2);
}
#endif





#if 1
event snapshot (t = 0; t <= tEnd; t += 0.5)
{
  char named[80];
  sprintf (named, "dump-%.2f", t);
  dump (file = named);
}
#endif





#if 0
event interface (t = 0; t <= tEnd; t += 0.5)
{
  char namei1[80];
  sprintf (namei1, "interface-drop-%.2f.txt", t);
  FILE* fp1 = fopen (namei1, "w");
  output_facets (f1, fp1);
  fclose (fp1);

  char namei2[80];
  sprintf (namei2, "interface-mpl-%.2f.txt", t);
  FILE* fp2 = fopen (namei2, "w");
  output_facets (f2, fp2);
  fclose (fp2);
}
#endif





#if 0
static inline void myprolong(Point point, scalar s){
  if(is_refined(cell))
  {
    foreach_child()
    {
      s[] = 1.;
      if(is_refined(cell))
        myprolong(point, s);
    }
  }
}

static inline void myrestrict1 (Point point, scalar s)
{
  double sum = 0.;
  scalar ss[];
  foreach_child()
  {
    if (f1[] > 0.0 && f1[] < 1.0)
    {
      s[] = 1.;
      sum++;
    }
  }
  ss[] = sum;

  if (ss[] > 0.)
  {
    s[] = 1.;
    foreach_child()
    {
      s[] = 1.;
      myprolong(point, s);
    }
  }
}



static inline void myrestrict2 (Point point, scalar s)
{
  double sum = 0.;
  scalar ss[];
  foreach_child()
  {
    if (f2[] > 0.0 && f2[] < 1.0)
    {
      s[] = 1.;
      sum++;
    }
  }
  ss[] = sum;

  if (ss[] > 0.)
  {
    s[] = 1.;
    foreach_child()
    {
      s[] = 1.;
      myprolong(point, s);
    }
  }
}



scalar sign1[], sign2[], sign3[];

event compute_thin_structure(i++)
{
  sign3.prolongation = myprolongation;
  sign3.restriction = myrestrict;
  sign3.refine = refine_injection;

  foreach()
  {
    sign1[] = 0.;
    sign2[] = 0.;
    sign3[] = 0.;
  }

   int l_sign = 4;
   for (int ilev = depth() - 1; ilev >= l_sign; ilev--)
    foreach_level(ilev){
      if(is_refined(cell))
      {
        myrestrict1(point, sign1);
        myrestrict2(point, sign2);
      }
    }

  foreach()
  {
    if ((sign1[] + sign2[])==2 && f1[] < 1. && f2[] < 1.)
      sign3[] = 1.;
  }
}
#endif



event adapt (i++)
{
  adapt_wavelet ({f1, f2, u, sign6, cf}, (double[]){1e-3, 1e-3, 1e-3, 1e-3, 1e-1, 1e-3}, LEVEL, MINLEVEL);
}



#if 1
event movies (t += 0.1)
{
  clear(); 
  view(fov = 5., tx = -0.35, ty = 0.001);

  draw_vof ("f1", lc = {1,0,0});
  draw_vof ("f2", lc = {1,0,0});
  cells();
  squares("le");
  mirror (n = {0,-1})
        {
          draw_vof ("f1", lc = {1,0,0});
          draw_vof ("f2", lc = {1,0,0});
          cells();
          squares("le");
        }
  save("movie_mesh.mp4");

  scalar vel[];
  foreach()
    vel[] = pow(sq(u.x[])+sq(u.y[]), 0.5);
  boundary({vel});
  draw_vof ("f1", lc = {1,0,0});
  draw_vof ("f2", lc = {1,0,0});
  squares("vel", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f1", lc = {1,0,0});
          draw_vof ("f2", lc = {1,0,0});
          squares("vel", linear = true);
        }
  save("movie_vel.mp4");

  scalar ff[];
  foreach()
  {
     ff[] = f1[] + f2[];
  }
  boundary({ff}); 
  draw_vof ("f1", lc = {1,0,0});
  draw_vof ("f2", lc = {1,0,0});
  squares("ff", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f1", lc = {1,0,0});
          draw_vof ("f2", lc = {1,0,0});
          squares("ff", linear = true);
        }
  save("movie_vf.mp4");

  draw_vof ("f1", lc = {1,0,0});
  draw_vof ("f2", lc = {1,0,0});
  squares("sign6", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f1", lc = {1,0,0});
          draw_vof ("f2", lc = {1,0,0});
          squares("sign6", linear = true);
        }
  save("movie_film.mp4");


  draw_vof ("f1", lc = {1,0,0});
  draw_vof ("f2", lc = {1,0,0});
  squares("cf", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f1", lc = {1,0,0});
          draw_vof ("f2", lc = {1,0,0});
          squares("cf", linear = true);
        }
  save("movie_tracer.mp4");
}
#endif


event end (t = tEnd) {}