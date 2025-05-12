/**
# Shear cell for a granular media with gravity

We simulate a simple shear cell for a granular media of glass beads in the presence of gravity with the top plate moving with velocity $U$ and the bottom plate is kept fixed. We compare it with an Analytical solution ([Cawthorn 2011](https://www.repository.cam.ac.uk/items/6aff8476-2979-42b8-97d7-bbb5ad9922ce)) and discrete element sinulation ([Guillard et al.](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5BE7B7A60DE924ACC00D7DC906E98B41/S0022112016006054a_hi.pdf/_div_class__title__Scaling_laws_for_segregation_forces_in_dense_sheared_granular_flows__div_.pdf)).

## Equations

Constitutive equation for granular media is included in Navier-Stokes and available in [http://basilisk.fr/sandbox/M1EMN/Exemples/granular.h]()

Analytical solution can be evaluated by integrating the ODE (see [Cawthorn 2011](https://www.repository.cam.ac.uk/items/6aff8476-2979-42b8-97d7-bbb5ad9922ce)):
$$\frac{\text{d} u}{\text{d} z} =  \frac{I_0\sqrt{P_c/\rho-gz}}{d}\left(\frac{\theta-P_c/\rho+z}{\mu_2/\mu_1(P_c/\rho g - z)-\theta}\right) $$ 

where, $\theta = \tau_0/\mu_1 P_c$; $P_c,\tau_0$ being the prescribed confining pressure and shear stress at the top respectively. The vertical coordinate is defined in the range $z_1<z<0$ where the material is sub-yield. Thus $z_1=P_c/\rho g(1-\theta)$. All other symbols denote the usual parameters of granular media.

## Code

Includes and definitions, note `granular.h`

 */


#include "grid/quadtree.h"
#include "granular.h"
#include "utils.h"
#include "view.h"

// For post-processing files in Paraview
//#include "output_vtu_foreach.h"

//Media properties in SI units.
#define RHO1 1000
#define G 7.5 

//Reference scales in SI units.
#define Lref 7.5e-2
#define Pref 1600
#define Uref (sqrt(Pref/RHO1))


// Dimensional input parameters in SI units
#define Up 2.5
#define Pc Pref 


// Dimensionless input parameters
#define us (Up/Uref)
//#define us 1.9764
#define cp (Pc/Pref)
#define FrSqInv (RHO1*G*Lref/Pc)

#define tend 15


int numFOut=0;

int main()
{
  periodic(right);			// Periodic BC left-right
  L0=1;				// Domain length (x-direction)				
  origin(-0.5,0,0);				
  init_grid(128);			// Designates the number of (square) cells per L0.
   
// Gravity as body force    
  const face vector g[] = {0.,-1.*FrSqInv};
  a = g;

  DT = 1e-1;
  TOLERANCE = 1e-3;
  run(); 
}

/**
Boundary conditions
*/

//Velocity at the top
  u.n[top] = dirichlet(0);
  u.t[top] = dirichlet(us);
  
//Confining pressure at the top
  p[top] = dirichlet(cp);
  pf[top] = dirichlet(cp);

  u.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(0);
  
/**
Initial conditions
 */
event init (t=0)
  {		
  //Granular media																				//This updates ghost cells for applying BCs.
  scalar phi[];
  foreach_vertex()
    phi[] = 1;
  fractions (phi, f);
    
  //velocity field
    foreach(){
      u.x[] = 0.;
  	u.y[] = 0.;
    }
    boundary({u.x,u.y});
  }

/**
Standard logfile
*/
  
void mg_print (mgstats mg)
{
   if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0);
}
event logfile (i++) {
  stats s = statsf (f);
  fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  fflush (stderr);
}

/**
Output
 */
// ASCII output file for velocity and shear stresss profile
event profiles (t += .5)
{
    FILE * fp = fopen("zprof.txt", "w");
    scalar shear[];
    foreach()
    shear[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    boundary ({shear});
    for (double y = 0.; y < 1.0; y += 1./32)
        fprintf (fp, "%g %g %g\n", y, interpolate (u.x, 0.0, y), interpolate (shear, 0, y));
    fclose (fp);
}


event stop (t = tend){
//char name[80];
//sprintf(name, "./sol_FIN_%g",t); output_vtu((scalar *) {f,p}, (vector *) {u}, name);
return 1; 
}


/**

## Results
~~~gnuplot Velocity profile comparison with Analytical and DEM 
set xlabel "u/U" font "Times-Roman,18"
set ylabel "z/H" font "Times-Roman,18"
set key right bottom
Up=2.5
z0=0.2133
p[0:1][0:1] 'zprof.txt' u ($2/Up):1 t 'Basilisk' w p pt 4 lw 2 ps .75,\
'../Guillard.txt' u 1:2 t 'DEM(Guillard-JFM-2016)' w p pt 3 lw 2 ps .75 lc 'red',\
'../Cawthorn.txt' u 2:(($1*z0/7.5e-2)+1) t 'Analytical(Cawthorn-thesis-2011)' w l lw 2 lc 'black'

~~~

## Bibliography

* [Cawthorn thesis 2011](https://www.repository.cam.ac.uk/items/6aff8476-2979-42b8-97d7-bbb5ad9922ce)
* [Guillard, F., Forterre, Y., & Pouliquen, O. (2016). Scaling laws for segregation forces in dense sheared granular flows. Journal of Fluid Mechanics, 807, R1.](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5BE7B7A60DE924ACC00D7DC906E98B41/S0022112016006054a_hi.pdf/_div_class__title__Scaling_laws_for_segregation_forces_in_dense_sheared_granular_flows__div_.pdf)

*/
