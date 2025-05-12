/**
# Bubble Rise and Coalescence in Oldroyd-B material
*/

#define ELASTIC 1

#include "axi.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"
#include "view.h"
#include "log-conform.h"

#define Fr 200.0         // Archimedes number
#define Bo 2.0          // Bond Number
#define Wi 4.0           // Weisenberg number

#define BETA 0.0         // Viscosity ratio of the solvent to the solute
#define MUR 10.          // Viscosity Ratio
#define RHOR 10.         // Density Ratio

#define R0 0.5           // Drop Radius
#define HEIGHT 50*R0     // Domain Size
#define Xo 20.*R0          // Drop Initial Position
#define Xi 40.*R0         // Interface Position

scalar rhov[];
face vector alphav[], muv[];
scalar f1[], f2[], * interfaces = {f1, f2};
scalar lambdav[], mupv[];


double tEnd = 20.01;

double xdd = 0.;         // Drop center of mass position
double vdd = 0.;         // Drop center o mass velocity
double uddx = 0.;        // Drop volume

double rho1, rho2, mu1, mu2, U, tau_y, mu_p;

int LEVEL = 12;
int MIDLEVEL = 6;
int MINLEVEL = 6;


int main ()
{
  size (HEIGHT);
  init_grid (1 << MIDLEVEL);

  rho1 = 1.;
  rho2 = rho1/RHOR;
  U = sqrt(Fr*(rho1 - rho2));
  mu1 = BETA*U/Fr;
  mu2 = U/(Fr*MUR);
  
  f1.sigma = f2.sigma = (rho1 - rho2)/Bo;

  alpha = alphav;
  rho = rhov;
  mu = muv;
  mup = mupv;
  lambda = lambdav;

  TOLERANCE = 1e-4;
  NITERMAX = 100;
  run();
}



/** Functions to refine the mesh in the film*/
static inline void myrestrict(Point point, scalar s){
  s[] = 0;
}

static inline void myprolongation(Point point, scalar s){
  foreach_child()
    s[] = 0;
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


event properties (i++)
{ 
  scalar fa1[], fa2[];
  foreach()
  {
    fa1[] =  (4.*f1[] + 2.*(f1[-1,0] + f1[1,0] + f1[0,-1] + f1[0,1]) + f1[1,1] + f1[-1,1] + f1[1,-1] + f1[-1,-1])/16.;
    fa2[] =  (4.*f2[] + 2.*(f2[-1,0] + f2[1,0] + f2[0,-1] + f2[0,1]) + f2[1,1] + f2[-1,1] + f2[1,-1] + f2[-1,-1])/16.;
  }
  boundary ({fa1, fa2});

  // Here we pass the viscosity, specific volume and density.
  foreach_face()
  {
    double fm1 = (fa1[] + fa1[-1])/2.;
    double fm2 = (fa2[] + fa2[-1])/2.;
    
//    muv.x[] = fm.x[] / ((1. - fm1 - fm2)/mu1 + fm1/mu2 + fm2/mu2);
    muv.x[] = fm.x[] * ((1. - fm1 - fm2)*mu1 + fm1*mu2 + fm2*mu2);
    alphav.x[] = fm.x[]/ ((1 - fm1 - fm2)*rho1 + fm1*rho2 + fm2*rho2);
  }

#if ELASTIC
  foreach()
  {
    rhov[] = cm[]*((1 - f1[] - f2[])*rho1 + f1[]*rho2 + f2[]*rho2);
    mupv[] = (1. - BETA)*(1 - fa1[] - fa2[])*U/Fr;
    lambdav[] = Wi*(1 - fa1[] - fa2[])/U;
  }
  boundary ({muv,alphav,rhov,mupv,lambdav});
#else
  foreach()
  {
    rhov[] = cm[]*((1 - f1[] - f2[])*rho1 + f1[]*rho2 + f2[]*rho2);
  }
  boundary ({muv,alphav,rhov});
#endif
}




double thick_ms = 0.3;
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

  printf ("i = %d t = %g, # = %d, l.mi = %g, l.ma = %g fm = %g\n", i, t, j, s.min, s.max, thick_ms);
  printf("U = %g Wi = %g Fr  = %g Bo = %g Mur = %g\n", U, Wi, Fr, Bo, MUR);
  fflush(stdout);
}






scalar sign6[];
int q = 0;
#if 1
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

if (posit1ma > (Xi - 0.3) || q == 1)
{
  q = 1;
  int i, j;
  i = j = 0;

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

/**
Here we write the film profile in a file
*/
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
  thin_max = thick_ms*3. < max_film ? thick_ms*3 : max_film;



  if (thick_ms < thin_max)
  {  
    foreach(serial)
    {
      double s4 = 0.;
      double s5 = 0.;

      for (int k = 0; k < i; k++){
          {
            double thin1 = sqrt( sq(dropx[k] - x) + sq(dropy[k] - y));
            if(thin1 < thin_max && (f1[]+f2[]) < 0.99)
            {
              s4 = 1.;
              break;
            }
          }
      }
      

      for (int k = 0; k < j; k++){
            double thin2 = sqrt( sq(interx[k] - x) + sq(intery[k] - y));
            if(thin2 < thin_max && (f1[]+f2[]) < 0.99)
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
  boundary((scalar *){sign6});


  double num_cell = 10;
  if (thick_ms < max_film)
    LEVEL = (int)(log(num_cell*HEIGHT/thick_ms)/log(2));
  if (LEVEL <= 12)
    LEVEL = 12;
  if (LEVEL > 17)
    LEVEL = 17;

  sign6.prolongation = myprolongation;
  sign6.restriction = myrestrict;
}

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
#endif




// Mesh adapt
event adapt (i++)
{
  #if ELASTIC
    adapt_wavelet ({f1, f2, u, sign6}, (double[]){1e-3, 1e-3, 1e-3, 1e-3, 1e-1}, LEVEL, MINLEVEL);
  #else
    adapt_wavelet ({f1, f2, u, sign6}, (double[]){1e-3, 1e-3, 1e-3, 1e-3, 1e-1}, LEVEL, MINLEVEL);
  #endif
}






#if 1
event outputs_fields (t = 0; t <= tEnd; t += 0.05)
{
  scalar drox[], droy[], intx[], inty[];
/**
We identify the interface positon on the interfacil cells and it pass to a array.
*/
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
    fprintf(filmt, "%g %g %g %g\n", ss[k], thick_p[k], dropx[k], dropy[k]);
  }
  fclose(filmt);



  char field[80];
  sprintf (field, "fields-film-%.2f.txt", t);
  FILE* fld = fopen (field, "w");
  
  fprintf(fld, "%g 0.0 %g %g %g %g %g  %g %g %g %g\\n", posit1ma, interpolate(u.x, posit1ma, 0), interpolate(u.y, posit1ma, 0), interpolate(p, posit1ma, 0), interpolate(f1, posit1ma, 0), interpolate(f2, posit1ma, 0), interpolate(tau_p.x.x, posit2ma, 0), interpolate(tau_p.y.y, posit2ma, 0), interpolate(tau_p.x.y, posit2ma, 0), interpolate(tau_qq, posit2ma, 0));
  fprintf(fld, "%g 0.0 %g %g %g %g %g  %g %g %g %g\\n", posit2ma, interpolate(u.x, posit2ma, 0), interpolate(u.y, posit2ma, 0), interpolate(p, posit2ma, 0), interpolate(f1, posit2ma, 0), interpolate(f2, posit2ma, 0), interpolate(tau_p.x.x, posit2ma, 0), interpolate(tau_p.y.y, posit2ma, 0), interpolate(tau_p.x.y, posit2ma, 0), interpolate(tau_qq, posit2ma, 0)); 
  foreach(serial)
  {
    if (sign6[] == 1.)
    {
      fprintf(fld, "%g %g %g %g %g %g %g %g %g %g %g\\n", x, y, u.x[], u.y[], p[], f1[], f2[], tau_p.x.x[], tau_p.y.y[], tau_p.x.y[], tau_qq[]); 
    }
  }
  fclose (fld);
 }
}
#endif





#if 1
event outputs_shape (t += 1.00)
{
  char namei1[80];
  sprintf (namei1, "interface-drop-%.2f.txt", t);
  FILE* fp1 = fopen (namei1, "w");
  output_facets (f1, fp1);
  fclose (fp1);

  char namei2[80];
  sprintf (namei2, "interface-inter-%.2f.txt", t);
  FILE* fp2 = fopen (namei2, "w");
  output_facets (f2, fp2);
  fclose (fp2);


  char field2[80];
  sprintf (field2, "fields-output-%.2f.txt", t);
  FILE* fld2 = fopen (field2, "w");

  scalar ff[];
  foreach()
  {
    ff[] = f1[] + f2[];
  }

  output_field ((scalar *){u, p, ff, sign6, tau_p}, fld2, n = 500, linear = true,  box = {{xdd/vdd - 5*R0, 0.},{xdd/vdd + 5*R0, 5*R0}});
  fclose (fld2);
}
#endif




// Backup
#if 1
event snapshot (t = 0; t <= tEnd; t += 1.00)
{
  char named[80];
  sprintf (named, "dump-%.2f", t);
  dump (file = named);
}
#endif





#if 1
event movies (t+=0.01)
{
  clear(); 
  view(fov = 2.0, tx = 0., ty = -0.4, psi = -pi/2);
  
  translate (x =  Xo - xdd/vdd)
  {
  draw_vof ("f1", lw = 4, lc = {1,0,0});
  draw_vof ("f2", lw = 4, lc = {1,0,0});
  cells();
  mirror (n = {0,-1})
        {
          draw_vof ("f1", lw = 4, lc = {1,0,0});
          draw_vof ("f2", lw = 4, lc = {1,0,0});
          cells();
        }
  }
  save("movie_movie.mp4");

  scalar ff[];
  foreach()
  {
    ff[] = f1[] + f2[];
  }

  boundary({ff});
  translate (x =  Xo - xdd/vdd)
  { 
  draw_vof ("f1");
  draw_vof ("f2");
  squares("ff", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f1");
          draw_vof ("f2");
          squares("ff", linear = true);
        }
  }
  save("movie_vf.mp4");

}
#endif



event end (t = tEnd) {}
