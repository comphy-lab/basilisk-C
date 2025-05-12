/** 
# Solution of the flow in a porous media
## Problem 
 What is the steady flow in a fully saturated porous media?
 
 This is the problem of Todd and Bear 1959. Using analog electrical computations, they reproduced the flow in a fully saturted soil.
 
 "When channel level exceeds the adjacent ground surface elevation,
 water moves through and under the confining levels into adjacent land.
 If drainage facilities are inadequate the soil becomes saturated and water often ponds on the surface."
 
 
 To solve this problem, we have the darcy law with $\beta$, the permeability :
 $$ 0 = -  \overrightarrow{\nabla} p -  \frac{1}{\beta}\overrightarrow{u}  + \rho  \overrightarrow{g}
\text{  \; with } \overrightarrow{\nabla} \cdot \overrightarrow{u} = 0$$ 
 
 the water starts in a top right corner box (the river with a dam, creating an extra pressure) and goes out in the soil, the depth of the soil is $H_d$, the water goes out the left.


 ~~~gnuplot pressure
 reset
 set view map
 set size (.5*1.6),1
 set label "river" at 8.4,8.
 set label "dam" at 7.5,9.5
 set label "soil surface" at 2,7.5
 set label "porous media" at 2,5
 set label "impervious rock" at 2,0.5
 set arrow from 8,7 to 8,9 nohead lc rgb 'black'
 set arrow from 8,8.5 to 10,8.5 nohead lc rgb 'blue'
 unset key
 unset surface
 set contour base
 set cntrparam levels incremental -15,1,15
 splot [][:] 'pressure.txt' u 1:2:($3)  w l not
 ~~~
 
 
## the problem without dimension:

Solution of the "Basic problem" in 2D, potential incompressible flow
$$ \bar u =- \frac{\partial \bar p}{\partial \bar x}, \;\; \bar v = - \frac{\partial \bar p}{\partial \bar y} - 1,\;\;\; 
  \frac{\partial \bar u}{\partial \bar x} +  \frac{\partial \bar v}{\partial \bar y}=0$$
  to be solved as:
  $$\frac{\partial^2 \bar p}{\partial \bar x^2} +  \frac{\partial^2 \bar p}{\partial \bar y^2}=0$$
  ij the river the pressure is hydrostatic
 at the wall $\bar y=0$ velocity slips.
 
 The shape of the aquifer is the overal porous media as here we are in a fully saturated configuration.
  
 Boundary conditions, those are pressure BC, teh velocity is a consequence:
 
* on the left neumann condition to mimick a large domain $x=0$ and $0<y<H_d$
 $\partial u/\partial x=0$, $\partial v/\partial x=0$,
 $\partial p/\partial x=0$

* on the bottom $y=0$ no penetration, $v=0$,  $\partial u/\partial y=0$,
 $-\partial p/\partial y=1$
 
* on top, before the dam $y=H_d$ and $x<x_0$:
 reference atmospheric pressure $p=0$, slip $\partial u/\partial y = \partial v/\partial y=0$
 
* on the right, $x=L_0$ up to the bottom of the river,
 by symmetry  $\partial p/\partial x=0$

* for the river for $x>x_0$ the
 pressure along the shape of the rive is $p= \rho g (H_0 + H_d-y)$ where $\rho g H_0$ is the over pressure due to the hight of the dam.
 
# Code

includes for advection of a field running tools and Poisson solver
*/
#include "run.h"
#include "poisson.h"

#define MAXLEVEL 6 //  11
#define MINLEVEL 5  //

scalar p[], source[];
double H0=10;
double Hd=7;
double tmax;
face vector beta[],rhop[];
mgstats mgp;
face vector u[];
/** 

## boundary conditions

Domain of size $2L_0$, 
 at the bottom hydrostatic pressure gradient, and zero pressure at  output at the right,  
pressure at the top of the boundary is a little bit larger tahn zero to compensate hydrostatic pressure 
*/ 
u.n[right]  =  neumann(0);
u.t[right]  =  neumann(0);
u.n[left]   =  neumann(0);
u.t[left]   =  neumann(0);
u.n[bottom] =  dirichlet(0);
u.t[bottom] =  neumann(0);
u.n[top]   =   neumann(0);
u.t[top]   =   neumann(0);

p[left]  = neumann(0);
p[right] = neumann(0);
p[top]   = dirichlet(0) ;
p[bottom]= neumann(1);

/**
 domain is `L0 x L0`
*/
bid river;
u.t[river] = neumann(0);
p[river] =dirichlet(H0 + Hd-y) ;

int main()
{    
    L0=10.;
    Y0=0;
    X0=0;
    init_grid (1 << MAXLEVEL);
//    DT = .05;
    tmax = 50;
   // run();
}
/**
gravity and initial values
*/
 const face vector g[] = {0.,-1.};

 event init (i = 0) {
     mask (y> Hd ? top: none);
     mask ( x>8 && y>5 ? river : none);
     
   foreach_face() {
        beta.x[] = 1;
        rhop.x[] = 1;
    }  
    boundary  ((scalar *){beta}); 
    boundary  ((scalar *){rhop}); 
 
  foreach(){
    p[] = (y<L0? Hd-y :0 );
    source[] = 0.;
  }
  boundary ({p});

// We may change the default gradient function (used for advection) to minmod-limited (rather than the centered default).
//  gradient = minmod2;
}
 
event pressurevelocities (i++; t<tmax)
{
/**
  evaluate physcical quantities $\beta$ permeability and $\rho$ density as function of
pressure
$\beta =1$
and the density is 1 by adimensionalisation
$\rho =1$
*/    
  foreach_face() {
      beta.x[] = (y>4? 1 : 1) ;
      rhop.x[] = 1;
    } 
    boundary  ((scalar *){beta}); 
    boundary  ((scalar *){rhop});  
/**

We have to solve 
$$ 0 = -  \overrightarrow{\nabla}   p -   \frac{1}{\beta} \overrightarrow{  u} - \overrightarrow{\rho e_y}
\text{  with } \overrightarrow{  \nabla} \cdot \overrightarrow{  u} = 0$$ 
 with $\beta = \beta_m$ (with $\beta_m \gg 1$) outside the water, and $\beta =1$ in the water,
 and by choice of scales $\overrightarrow g=- \overrightarrow{  e_y}$.      

this gives : 
$$0 = - \nabla \cdot (\beta \nabla p  ) +  \nabla \cdot (\beta \rho \overrightarrow g ) $$ 

We solve Poisson equation 
$\nabla \cdot (\beta \nabla p  )= s$ with the source term of the Poisson equation
$$ 
s= \nabla \cdot (\beta \rho \overrightarrow g )
$$
*/
    foreach(){
        source[]=0.0;
        } 
  boundary  ({source});
/**
solve Poisson equation 
$$\nabla \cdot (\beta \nabla p  )= s$$ 
with [http://basilisk.fr/src/poisson.h](http://basilisk.fr/src/poisson.h)
*/

  mgp = poisson (p, source, beta);

/** the velocity is then computed from the gradient (note the `face_gradient_x p=  ((p[i] - p[i-1])/Delta)`) see 
[http://basilisk.fr/src/poisson.h](http://basilisk.fr/src/poisson.h)

*/
   foreach_face()
    u.x[] = (- beta.x[]*face_gradient_x (p, 0) + beta.x[]*rhop.x[]*g.x[]);
  boundary ((scalar *){u});
}  
/**
error
*/
event logfile (i++)
{
    stats s = statsf (p);
    fprintf (stderr, "%d %g %d %g %g %g\n",
             i, t, mgp.i, s.sum, s.min, s.max);
     fprintf (stderr,"%g \n",dt);
}

/**
Save in a files...
*/
#if 0
event gfsview (i += 100) {
  static FILE * fp = popen ("gfsview2D gfsview.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

event sauve (t+=.1 )
{
FILE *  fpc = fopen("pressure.txt", "w");
output_field ({p,u.x,u.y}, fpc, linear = true);
fclose(fpc);
    
FILE *  fpq = fopen("qh.txt", "w");
    double dx=.1;
    for (double x = 0 ; x < L0; x += dx){
        for (double y = 0 ; y < L0; y += dx){
            fprintf (fpq, "%g %g %g \n",
                     x,    interpolate (p, x, 0), interpolate (u.x, x, 0));}
        fprintf (fpq,"\n");}
    fclose(fpq);
}


 event photo ( i=0 ){
 static FILE * fzb = fopen ("field.png", "w");
 output_ppm (p, fzb, min = 0, max = 15, n = 1 << MAXLEVEL, linear = true);
 }


/**
# Results
## Run
To compile and run:

 
~~~bash
 make toddbear59.tst; make toddbear59/plots ; make toddbear59.c.html;
~~~


## Plots

 
  Figure of iso pressure
 
~~~gnuplot
 reset
 set pm3d map
 set palette rgbformulae 22,13,-31;
 unset colorbox
 set xlabel "x  iso p"
 set ylabel "y"
 set size (.5*1.6),1
 splot [][ :][-10:15] 'pressure.txt' u 1:2:3   not
 reset
~~~
 

  Figure of iso pressure

~~~gnuplot pressure
 reset
 set view map
 set size (.5*1.6),1
 set label "river" at 8,8.5
 set label "soil surface" at 2,7.5
 set label "porous media" at 2,5
 set label "impervious rock" at 2,0.5
 set arrow from 8,7 to 8,9 nohead lc rgb 'black'
 set arrow from 8,8.5 to 10,8.5 nohead lc rgb 'blue'
 unset key
 unset surface
 set contour base
 set cntrparam levels incremental -15,1,15
 splot [][:] 'pressure.txt' u 1:2:($3)  w l not
~~~

 plot of Hydraulic Head  $H=p/(\rho  g)+z$:
 
~~~gnuplot Hydraulic Head
 reset
 set view map
 set size (.5*1.6),1
 unset key
 unset surface
 set contour base
 set cntrparam levels incremental 0,.5,15
 splot [][:] 'pressure.txt' u 1:2:($3+$2)  w l not
~~~


 

# Links
 
* see [http://basilisk.fr/src/hele-shaw.h](http://basilisk.fr/src/hele-shaw.h)
 
* see [http://basilisk.fr/sandbox/M1EMN/Exemples/darcyLambSneddon.c]()
 
* see [http://basilisk.fr/sandbox/M1EMN/HYDGEO/dupuit2D.c]()
 
* see [http://basilisk.fr/sandbox/M1EMN/HYDGEO/toddbear59.c]()

* see [http://basilisk.fr/sandbox/M1EMN/HYDGEO/vauclin.c]()

# Bibliography

* Joseph Bear "Dynamic of fluids in porous media" p 238
 
* [Todd Bear 1959](https://escholarship.org/uc/item/1nx0q3dd) River Seepage investigation

* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv_aquifere.pdf)
 "Ecoulements en milieux naturels:
écoulements en milieux souterrains" Cours MSF12, M1 UPMC

 
 
*/
