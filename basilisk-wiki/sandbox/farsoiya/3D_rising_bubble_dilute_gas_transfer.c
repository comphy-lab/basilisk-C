/**
# Soluble gas diffusing from a rising bubble

This is the example discussed in section 2.4.2 of [Farsoiya et al.,
2020](#farsoiya2020).

The four cases illustrated below are considered. The left-half of each
figure is the vertical velocity component and the right-half the
tracer concentration when close to the stationary regime.

<center>
<table>
<tr>
<td>![](axi_rising_bubble/final-1.png){ width="200px" }</td>
<td>![](axi_rising_bubble/final-2.png){ width="200px" }</td>
<td>![](axi_rising_bubble/final-3.png){ width="200px" }</td>
<td>![](axi_rising_bubble/final-4.png){ width="200px" }</td>
</tr>
<tr>
<td><center>Case 1</center></td> 
<td><center>Case 2</center></td> 
<td><center>Case 3</center></td> 
<td><center>Case 4</center></td> 
</tr>
</table>
</center>

The Sherwood number characterises the gas diffusion from the bubble to
the liquid. The computed Sherwood number is compared to the theory of
[Levich, 1959](#levich1959).

~~~gnuplot Sherwood number vs time
alpha_c = 0.01
r0 = 1.

sherwood(D) = (ci = $2/$5,  \
               co = $3/$6, \
               area = $7,   \
               dtr = (tr - $2)/0.01, \
               tr = $2,              \
               Sh = dtr/area/(ci*alpha_c - co)*2.*r0/D)

# Levich formula to calculate transfer rate using terminal velocity
levich(ut, D) = 2.*sqrt(D*ut/2./r0/pi)*2.*r0/D


#get the rise velocity averaged for the last 50 steps


array ut[4]
do for [i in "1 2 3 4"] {
  stats '3D.dat' index 'case '.i u 1:4 every ::250 nooutput
  ut[i] = STATS_mean_y
}
array D[4] = [ 0.4472, 0.5029, 0.17416, 0.08944 ]

tr = 0.
set key center right
set xlabel 'Time'
set ylabel 'Sh'
set xrange [0:2]
set yrange [0:14]
plot for [i = 1:4]							\
     '3D.dat' index 'case '.i u ($1*ut[i]/2):(sherwood(D[i])) w l t 'Axi Case '.i lt i,  \
	'3D.dat' index 'case '.i u ($1*ut[i]/2):(sherwood(D[i])) w l t 'Axi Case '.i lt i,  \
     for [i = 1:4] levich(ut[i], D[i]) w l t '' lt i dt 2
~~~

## References

~~~bib
@article{farsoiya2020,
  title = {Bubble mediated gas transfer of diluted component in turbulence},
  author = {P. K. Farsoiya and S. Popinet and L. Deike},
  journal = {Journal of Fluid Mechanics},
  year = {2020},
  note = {submitted}
}

@article{levich1959,
  title={Physicochemical hydrodynamics},
  author={Levich, Veniamin Grigorʹevich},
  year={1959},
  publisher={Prentice-Hall Inc.}
}
~~~
*/
#define AXIS 0

#if AXIS
  #include "axi.h"
#else
  #include "grid/octree.h"
#endif	

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "henry.h"
//#include "view.h"
// #include "navier-stokes/perfs.h"

scalar c[], * stracers = {c};
double bubble_radius = 1.;
double box_size = 20.;
double conc_liq1 = 0, conc_gas1 = 1.;
double end_time[4]={7.5,3,3,2};	
int MAXLEVEL, dcase = 1;

int main (int argc, char **argv)
{
  size (box_size);
	
  MAXLEVEL = 10;
#if !TREE
  N = 1 << MAXLEVEL;
#endif

  rho1 = 1.;
  rho2 = 0.01;
  c.alpha = 1./30.;
  TOLERANCE = 1e-4;
  
//#if AXIS
//  remove("axi.dat");
//#else
//  remove("3D.dat");
//#endif
  
  for (dcase = 1; dcase <= 4; dcase++) {
    switch (dcase) {
    case 1:
      c.D2 = 0.4472;
      f.sigma = 10.0;	
      G.x = - 2.5;	
      break;
    case 2:
      c.D2 = 0.5029;
      f.sigma = 10.0;	
      G.x = - 7.8125;	
      break;
    case 3:
      c.D2 = 0.17416;
      f.sigma = 1.0;	
      G.x = - 10.0;
      break;
    case 4:
      c.D2 = 0.08944;
      f.sigma = 10.0;	
      G.x = - 7.8125;
      break;
    default:
      fprintf (stderr, "Error: must specify case\n");
      exit (1);
    }

    mu1 = c.D2;
    mu2 = c.D2/20.;
    c.D1 = c.D2*100.;
    
   
    run();
  }
}

event init (t = 0)
{
#if TREE  
	#if AXIS
		refine (sq(2.*bubble_radius) - sq(x - box_size*0.2) - sq(y) > 0 &&
	  level < MAXLEVEL);
	#else
	refine (sq(2.*bubble_radius) - sq(x - box_size*0.2) - sq(y - box_size*0.5) - sq(z - box_size*0.5) > 0 &&
	  level < MAXLEVEL);
	#endif
	
#endif
	
  #if AXIS
	fraction (f, - (sq(bubble_radius) - sq(x - box_size*0.2) - sq(y)));
#else
	fraction (f, - (sq(bubble_radius) - sq(x - box_size*0.2) - sq(y - box_size*0.5) - sq(z - box_size*0.5)));
#endif
	
  foreach()
    c[] = conc_liq1*f[] + conc_gas1*(1. - f[]);
}

#if TREE
event adapt (i++) {
double uxemax = 0.2*normf(u.x).avg;
double uyemax = 0.2*normf(u.y).avg;
double uzemax = 0.2*normf(u.z).avg;

  double tr1emax = 0.2*normf(c).avg;
        adapt_wavelet({f,c,u},(double[]){0.01,tr1emax,uxemax,uyemax,uzemax}, maxlevel = MAXLEVEL);

}
#endif
	
event extract (t = 0; t += 0.01; t <= end_time[dcase-1])
{	
  double yb = 0., vb = 0., vbx = 0., area = 0., ci = 0., co = 0.;
  foreach (reduction(+:yb) reduction(+:vb) reduction(+:vbx)
	   reduction(+:ci) reduction(+:co)
	   reduction(+:area)) {
    double dvb = (1. - f[])*dv();
    vb += dvb;          // volume of the bubble
    yb += x*dvb;	// location of the bubble
    vbx += u.x[]*dvb;	// bubble velocity

    ci += c[]*(1. - f[])/(f[]*c.alpha + (1. - f[]))*dv();
    co += c[]*f[]*c.alpha/(f[]*c.alpha + (1. - f[]))*dv();
    
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f), p;      
      double alpha = plane_alpha (f[], n);
       // area of the bubble interface
#if AXIS
      area += y*pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p); 
#else										 
	 area += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
#endif
    }
  }
#if AXIS
  static FILE * fp = fopen ("axi.dat","a");
   if (i == 0){
    fprintf (fp, "\n\n# case %d\n", dcase);
    fprintf (fp, "t ci co vbx vb vbo area dt\n");
  }
    fprintf (fp,"%e %.12e %.12e %.12e %.12e %.12e %.12e %e\n",
   t,
   ci*2.*pi, co*2.*pi,
   vbx/vb, 2.*pi*vb, 2.*pi*statsf(f).sum, 2.*pi*area,
   dt);
#else
  static FILE * fp = fopen ("3D.dat","a");
if (i == 0){
    fprintf (fp, "\n\n# case %d\n", dcase);
    fprintf (fp, "#t ci co vbx vb vbo area dt\n");
  }

   fprintf (fp,"%e %.12e %.12e %.12e %.12e %.12e %.12e %e\n",
	   t,
	   ci, co,
	   vbx/vb, vb, statsf(f).sum, area,
	   dt);
#endif

  fflush(fp);
}

event pictures (t = end)
{
  char name[80];
#if   1
  sprintf (name, "dump-%d", dcase);
  dump (name);
#endif
/*	  #if AXIS
	  double ty[] = { - 0.25, - 0.45, - 0.5, - 0.7};
	  view (fov = 9, quat = {0.707,0.707,0,0},
		ty = ty[dcase - 1],
		width = 400, height = 800);  
	  squares ("c", spread = -1, linear = true, map = cool_warm);
	  draw_vof ("f");
	  mirror ({0,1}) {
	    squares ("u.x", spread = -1, linear = true, map = cool_warm);
	    draw_vof ("f");
	  }
	  sprintf (name, "axi_final-%d.png", dcase);
	  save (name);
		#else
		 double ty[] = { - 0.3, - 0.5, - 0.5, - 0.7};
	  view (fov = 9, quat = {0.707,0.707,0,0},
		ty = ty[dcase - 1],  tx = -0.5,
		width = 400, height = 800);  
	  squares ("c", spread = -1, linear = true, map = cool_warm, alpha =10);
	  draw_vof ("f");

	  sprintf (name, "3D_final-%d.png", dcase);
	  save (name);
*/	
}
