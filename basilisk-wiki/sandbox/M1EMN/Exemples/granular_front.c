/**
# Granular front

We propose an implementation of the Jop Pouliquen Forterre µ(I) rheology for an avalanche along a slope.
We focus on the shape of the front. We use a moving frame as the process requieres a long domain


# Code 
Includes and definitions
*/
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "heights.h"
// Domain extent
#define LDOMAIN 10.
// heap definition
double  H0,Lb,R0,D;
double dy;
//film mpg size
double  Lz;
double theta,u0;

/** the $\mu(I)$ parameters
*/
#define mu_s 0.32  // 17.7446716250569
#define Delta_mu 0.28
#define I0 0.4

/** the analytical Bagnold solution
$$u= \sqrt{cos(\theta)} (2/3) Ia/D (1-(1.-y/H)^{3/2})$$ 
*/
#define Ia fmax((-mu_s*I0 + I0*tan(theta))/(mu_s+Delta_mu - tan(theta)),0.01) 
double  U0(double y, double H){
//  return (sqrt(cos(theta))*Ia/D)*((y<=H)*( (1-pow(fmax(1.-y/H,0),1.5) ) - 3./5*0.9) + 0.*(y>H))*(2./3.) ; }
    return .0;}
/**
Initial heap
*/
double  Hi(double x){
 // return H0*sqrt(1-(x>Lb)*(x-Lb)*(x-Lb)/R0/R0);
   return H0*(1-(x>Lb)*(x-Lb)/R0);
  }
  
  
double xdeh(double h){
 double d = (tan(theta)-mu_s)/Delta_mu;
 double x = ((-1+d)*h-(2*atan((1+2*sqrt(h))/sqrt(3)))/sqrt(3)-(2*log(1 - sqrt(h)))/3.+log(1+sqrt(h)+h)/3.)/(-1+d) ;
return x/(tan(theta)-mu_s);
}

//passive fluid
#define RHOF 1e-3
#define mug  1e-4
// Maximum refinement
#define LEVEL 7 // 8 is OK
#define LEVELmax 7
#define LEVELmin 3
char s[80];
FILE * fpf;
scalar f[];
scalar lev[];
scalar nub[];
scalar det[];
scalar * interfaces = {f}; 
face vector alphav[];
face vector muv[];
scalar rhov[];
scalar eta[];
// note the v
/**
Boundary conditions for granular flow, pressure must be zero at the surface
*/
p[top] = dirichlet(-RHOF*LDOMAIN);
p[right] = neumann(0);
p[left] =  neumann(0);

u.n[top] = neumann(0);
u.t[bottom] =  dirichlet(0);

u.t[bottom] = dirichlet(U0(0,H0)-u0);
u.n[bottom] = dirichlet(0);
//u.x[left] = (y < 1 ? neumann(0) : dirichlet(0));
//u.x[left] = (y < .5 ? dirichlet(1.5*y) : dirichlet(0));
//u.n[left] = (y <= H0 ? dirichlet( U0(y,H0)) : dirichlet(0));
u.n[left] = neumann(0);
u.t[left] = neumann(0);


u.n[right] = neumann(0);
u.t[right] = neumann(0);

u.t[top] = dirichlet(0);

f[left]= (y <= H0 ? 1: 0 );
//f[left]= neumann(0);
f[right]= neumann(0);
/**
The main
*/
int main() {
  L0 = LDOMAIN;
  // number of grid points
  N  = 1 << LEVEL ;
  // maximum timestep
  DT = 0.5e-2;
  TOLERANCE = 1e-3;
  
// Initial condition
  H0=1.000001;
  Lb=1;
  R0=1.0001;
  // Grain size
  D=1./30;
  theta=0.33;
  const face vector g[] = {sin(theta),-cos(theta)};
  a = g;
  alpha =  alphav;
  mu = muv;
  rho = rhov;
// slip volcity
  u0 = 0.;
  
  FILE *fe;
  fe = fopen ("out_ex", "w");
  for(double h=.01;h<.95;h+=.01)
     fprintf (fe, "%g %g \n",xdeh(h)-xdeh(0),h);
  fclose (fe);   
  
  fprintf (stderr, "# U0(1,1) = %g \n",U0(1,1));
  fpf = fopen ("out_f", "w");
  Lz = LDOMAIN;
  run();
  fclose (fpf);
}
/**
initial heap and velocity, 
*/
event init (t = 0) {
 
 foreach()
    u.x[] = 0;
 
  foreach()
   lev[]= (y<2 ? LEVEL : 3); 
 
  scalar phi[];
  foreach_vertex()
    phi[] =  ( - y + Hi(x));
//    phi[] = min( - y + H0*sqrt(1-x*x/R0/R0), R0-x);
  fractions (phi, f);
/**
initialisation of hydrostatic pressure for granular phase  
*/
  foreach()
    p[] = f[]*(H0-y);
}
/**
total density 
*/
#define rho(f) ((f) + RHOF*(1. - (f)))
/**
Viscosity computing $D_2=D_{ij}D_{ji}$; the inertial number $I$ is 
and $\mu = \mu_s+ \frac{\Delta \mu}{1+I/I_0}$ 
the viscosity is $\eta = \mu(I)p/D_2$:

attention, we increase $\mu_g$ 

attention etamin = 10*sqrt(D*D*D);
	eta[] = max((muI*p[])/D2, etamin);
	eta[] = min(eta[],100); 
*/
event properties (i++) {
  trash ({alpha});
  foreach() {
    eta[] = mug;
    nub[]=0;
    if (p[] > 0.) {
      double D2 = 0.;
      foreach_dimension() {
	double dxx = u.x[1,0] - u.x[-1,0];
	double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
	D2 += sq(dxx) + sq(dxy);
      }
      if (D2 > 0.) {
	D2 = sqrt(2.*D2)/(2.*Delta);
	double In = D2*D/sqrt(p[]);
	double muI = mu_s + Delta_mu*In/(I0 + In);
	double etamin = 1.*sqrt(D*D*D);
	eta[] = max((muI*p[])/D2, etamin);
	eta[] = min(eta[],100);   
    eta[]+=0.5;
    nub[]=Delta_mu*(I0/In)/(1+I0/In)/(mu_s*(I0/In)+mu_s+Delta_mu);
    det[]=4*nub[]*(nub[]-1)+muI*muI*sq(1-nub[]/2)*f[];
  //  if(det[]>0){ eta[]=1;}
      }
    }
  }
  boundary ({eta});
  
  
  scalar fa[];
  foreach()
    fa[] = (4.*f[] + 
	    2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
	    f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
  boundary ({fa});
  
  foreach_face() {
    double fm = (fa[] + fa[-1,0])/2.;
    muv.x[] = (fm*(eta[] + eta[-1,0])/2. + (1. - fm)*(mug*(y<1.5*H0) + 10*(y>=1.5*H0) ));
    // mu.x[] = 1./(2.*fm/(eta[] + eta[-1,0]) + (1. - fm)/mug);
    alphav.x[] = 1./rho(fm);
  }
    
    foreach()
    rhov[] = rho(fa[]);
 // boundary_normal ({mu,alpha});
 boundary ({muv,alphav,rhov});
}
/**
convergence outputs
*/
void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,  
	         mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0);
}
/**
convergence outputs
*/
event logfile (i++) {
  stats s = statsf (f);
  fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  fflush (stderr);
}

/**
interface outputs
*/
//event interface (t = {0,1.,2.,3.,04.,05., 06.,07.,08.,09.,10.}) {
event interface (t = {0,1.,2.,3.,04.,05., 06.,07.,08.,09.,10., 15.,20.,25,30,35}) {
  output_facets (f, fpf); 
  char s[80];
  sprintf (s, "field-%4.2f.txt", t);
  FILE * fp = fopen (s, "w");
  output_field ({f,p,u,f}, fp, linear = true);
  fclose (fp);
  sprintf (s, "field.txt");
  fp = fopen (s, "w");
  output_field ({f,p,u,f}, fp, linear = true);
  fclose (fp);
  sprintf (s, "u1.txt");
  fp = fopen (s, "w");
  for (double yy = 0 ; yy < 4; yy += LDOMAIN/(pow(2,LEVEL)))
     fprintf (fp,"%g %g %g %g\n",yy,interpolate(u.x,1,yy),interpolate(p,1,yy),interpolate(f,1,yy));
  fclose (fp);
  sprintf (s, "u5.txt");
  fp = fopen (s, "w");
  for (double yy = 0 ; yy < 4; yy += LDOMAIN/(pow(2,LEVEL)))
     fprintf (fp,"%g %g %g %g\n",yy,interpolate(u.x,5,yy),interpolate(p,5,yy),interpolate(f,5,yy));
  fclose (fp);
  sprintf (s, "u10.txt");
  fp = fopen (s, "w");
  for (double yy = 0 ; yy < 4; yy += LDOMAIN/(pow(2,LEVEL)))
     fprintf (fp,"%g %g %g %g\n",yy,interpolate(u.x,10,yy),interpolate(p,10,yy),interpolate(f,10,yy));
  fclose (fp);
  
}
/**
 Saving the top position and runout positionfor slump measurements
 */
vector h[];
event timeseries (t += 0.1 ) {
    heights (f, h);
    double maxy = - HUGE,maxx = - HUGE;;
    foreach()
    if ((h.y[] != nodata) && (h.x[] != nodata)) {
        double yi = y + height(h.y[])*Delta;
        double xi = x + height(h.x[])*Delta;
        if (yi > maxy)
            maxy = yi;
        if (xi > maxx)
            maxx = xi;
    }
    char s[80];
    sprintf (s, "hmax");
    static FILE * fp0 = fopen (s, "w");
    fprintf (fp0, "%g %g %g\n", t, maxx, maxy);
    fflush (fp0);
}
/**
film output
*/
#if 1
event movie (t += 0.05) {
  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = 0, max = LEVEL, 
	      n = 400, box = {{0,0},{Lz,Lz}});

  foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
  boundary ({l});
  static FILE * fp2 = popen ("ppm2mpeg > velo.mpg", "w");
  output_ppm (l, fp2, min = 0, max = 2., linear = true, 
	      n = 400, box = {{0,0},{Lz,4*H0}});

  static FILE * fp3 = popen ("ppm2mpeg > f.mpg", "w");
  foreach()
    l[] = f[]*p[];
  output_ppm (l, fp3, min = 0, linear = true,
	      n = 400, box = {{0,0},{Lz,2*H0}});
}
 event pictures (t==3) {
  output_ppm (f, file = "f.png", min=0, max = 2,  spread = 2, n= 128, linear = true,
  box = {{0,0},{2,2}});
}

#endif


/**
mesh adaptation
*/
#if QUADTREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){5e-3,0.01,0.01}, LEVEL , 
		 list = {p,u,pf,uf,g,f});
}
#endif

#if gfsv
event gfsview (i += 10){
    static FILE * fp = popen ("gfsview2D -s v.gfv", "w");
    output_gfs (fp);
}
#endif


#if QUADTREE
// if no #include "grid/multigrid.h" then adapt
event adapt(i+=1){
 scalar K1[],K2[]; 

 foreach()
   K1[]=(f[0,1]+f[0,-1]+f[1,0]+f[-1,0])/4*noise();
 boundary({K1});
   
 
 for(int k=1;k<10;k++)  
 {
 foreach()
   K2[]=(K1[0,1]+K1[0,-1]+K1[1,0]+K1[-1,0])/4;
 boundary({K2});

 foreach()
   K1[]=K2[]; 
}

 foreach()
   K1[]=(K2[0,1]+K2[0,-1]+K2[1,0]+K2[-1,0])/4*noise();;
 boundary({K1});

 adapt_wavelet({K1,f},(double[]){0.001,0.01}, maxlevel = LEVELmax, minlevel = LEVELmin);
}
#endif



/**

# Run

to run

~~~bash
qcc -g -O2    -Dgfsv=1 -o granular_front granular_front.c -lm
./granular_front > out
~~~


~~~bash
 make granular_front.tst; make granular_front/plots; make granular_front.c.html  
~~~


# Results

Exemples  
 
~~~gnuplot front at several times
set xlabel "x"
set ylabel "h(x,t)" 
 p[][:]'out_f' w l,10./2.**8
~~~


front


~~~gnuplot front at several times
 p[0:20][0:2]'out_f' w l,25/2.**9,'u1.txt' u (1+$3):1 w l,'u5.txt' u (5+$3):1 w lp,'u10.txt' u  (10+$3):1 w lp
~~~



Pressure and velocity 
 
~~~gnuplot velocity and pressure
set pm3d map
set palette rgbformulae 22,13,-31;
unset colorbox
set multiplot layout 2,1
set xlabel "x  iso u"
set ylabel "y"
splot [][:1] 'field.txt' u 1:2:($5*$3)   not
set xlabel "x  iso p"
splot [][:1] 'field.txt' u 1:2:($4*$3)   not
unset multiplot
reset

~~~
 



Front  

~~~gnuplot front and analytical front
p[0:]'out_f'  w l,'out_ex' u ($1+16):2
~~~

some velocity profiles
 
~~~gnuplot velocity profile
p[][:2]'out_f' w lp,25/2.**9,'u1.txt' u (1+$2*$4):1 w l,'u5.txt' u (5+$2*$4):1 w lp,'u10.txt' u (10+$2*$4):1 w lp
~~~


velocity field


~~~gnuplot velocity arrows field 
 plot [][:1.5] 'field.txt' every 1:50 u 1:2:($5*$3):($6*$3)  w vector
~~~

vitesse du front
 
~~~gnuplot
 
 f(x) = a*x +b
 fit [3:] f(x) 'hmax' via a,b
 set key left
 p[0:]'hmax' t'debit' w lp,f(x) t'fit'
~~~
 
fronts
 
~~~gnuplot
 p[:0][0:1] \
     'field-7.00.txt' u ($1-f(7)):(($7<1 && $7 > .9 )? $2: NaN ) ,\
    'field-10.00.txt' u ($1-f(10)):(($7<1 && $7 > .9 )? $2: NaN ),\
    'field-12.00.txt' u ($1-f(12)):(($7<1 && $7 > .9 )? $2: NaN ),\
    'field-15.00.txt' u ($1-f(15)):(($7<1 && $7 > .9 )? $2: NaN ),\
    'field-20.00.txt' u ($1-f(20)):(($7<1 && $7 > .9 )? $2: NaN ),'out_ex' u ($1):2
~~~
 

[![](f.png)](./f.mpg)  
  velocity  (click on image for animation of the density)

 
  
[![](f.png)](./granular_front/velo.mpg)  
  velocity  (click on image for animation of velocity)



*/ 



/**
# Links 
* granular collapse

* granular sandglass

* bagnold



# Bibliography
* Lagrée, Staron, Popinet 
["The granular column collapse as a
continuum: validity of a two–dimensional
Navier–Stokes model with a μ(I)-rheology"](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/JFMcollapsePYLLSSP11.pdf) J. Fluid Mech. 2011 



Madeire 16 février 2015 
*/
