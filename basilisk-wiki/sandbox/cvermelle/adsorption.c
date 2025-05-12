/**
# Growth of a boundary due to precipitation 

We want to investigate the growth of a boundary due to a chemical reaction.
We consider the multilayer Saint-Venant and the mass conservation.

# Equations:

## Saint-Venant:

$$
\begin{aligned}
  \partial_t h_k + \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k & =
  S,\\
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = - gh_k  \mathbf{{\nabla}} (\eta)
\end{aligned}
$$  

The fluid is at rest: $\mathbf{u} = \mathbf{0}$. 
Considering that there is a reaction at the boundary and mass conservation,
the continuity equation writes

$$
\begin{aligned}
  \partial_t h_k & = 0, k>0 \\
  \partial_t h_0 & = -\partial_t z_b 
\end{aligned}
$$

Similarly, for the momentum equation we get:
$$\mathbf{{\nabla}} (\eta) = \mathbf 0 $$


Thus, the free surface should always be equal to a constant.


## Chemical equations:
We consider chemicals that diffuse and react such as:

$$
\partial_t c_i = D_i \nabla^2 c_i - \alpha_i c_i \sum_j k_j c_j \frac{1 - \delta_{ij}}{2}
$$

with:

* $D_i$ the diffusivity
* $\alpha_i$ a stoechiometic coefficient
* $k_j$ the rate constant 

For our simple case we have at the boundary:
$$ X_0 + X_1 \rightarrow X_2 $$

It translates into

$$
\begin{aligned}
  \partial_t c_0 & = D_0 \nabla^2 c_0 -k c_0 c_1\\
  \partial_t c_1 & = D_1 \nabla^2 c_1 -k c_0 c_1
\end{aligned}
$$

## Topography equation
Dropped from a parachute we have:

$$ \partial_t zb = \frac{1}{\rho_s} k c_0 c_1 $$

# Results
Solving for these equations we get for the evolution of the bottom:


~~~gnuplot
set multiplot

#settings for main plots
set ytics auto
set xtics auto

#first plot
s1x = 0.5
s1y = 1.0
set size s1x, s1y
o1x = 0.0
o1y = 0.0
set origin o1x,o1y
plot  'hydro.txt' u 2:4 t 'zb'

#second plot
s2x = 0.5
s2y = 0.5
set size s2x, s2y
o2x = 0.5
o2y = 0.5
set origin o2x,o2y
plot  'chimney.txt' u 1:2 t 'h'


#third plot
s3x = 0.5
s3y = 0.5
set size s3x, s3y
o3x = 0.5
o3y = 0.0
set origin o3x,o3y
plot  'c_wall.txt' u 2:4 t 'c0', 'c_wall.txt' u 2:5 t 'c1'

unset multiplot 
~~~
*/


/**
# Code
## Libraries and Definitions
*/

//Librairies
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "./solvers/diffusion_neumann.h"


//Physical parameters
#define HINIT 1.
#define RATECONSTANT 1.e3
#define DIFFUSIVITY 1.
#define TMAX .1

//usefull variables and files
double altitude;
FILE * ic, * bulk, * hydro;
FILE * ch;

char name_init[30];
char name_bulk[30];
char name_hydro[30];

/**
## Definition of concentration fields in the layered framework
*/

scalar c0;
scalar c1;

event defaults(i=0){
  c0 = new scalar[nl];
  tracers = list_append (tracers, c0);
  c1= new scalar[nl];
  tracers = list_append (tracers, c1);
}
event cleanup(t= end, last){
  delete({c0,c1});
}


/**
## Main function
*/
int main(){
  L0 = 1.;
  X0 = -0.5;
  G=1.;
  nl = 64;
  N = 256;
  DT = 1e-6; // To be set manually, for stability of horizontal diffusion

  sprintf(name_init,"initialcondition.txt");
  ic = fopen(name_init,"w");
  fclose(ic);
  sprintf(name_bulk,"c_wall.txt");
  bulk = fopen(name_bulk,"w");
  fclose(bulk);
  sprintf(name_hydro,"hydro.txt");
  hydro = fopen(name_hydro,"w");
  fclose(hydro);

  ch = fopen("chimney.txt","w");
  fclose(ch);
  run();


}

/**
## Diffusion step: 
Definition of the values of the fluxes at the top and the bottom of the domain
 */
double dcb=0., dct= 0.;

event diffusionOfC (i++){
  dt = dtnext(DT);
  horizontal_diffusion({c0} ,DIFFUSIVITY, dt);
  horizontal_diffusion({c1},DIFFUSIVITY, dt);
  foreach(){
    vertical_diffusion_neumann(point, h, c0, dt, DIFFUSIVITY,
                               dct, dcb);
    vertical_diffusion_neumann(point, h, c1,dt, DIFFUSIVITY,
                               dct, dcb);
  }

  // Computation of dz/dt and dc/dt dc1/dt
  scalar dz[], dc0[], dc1[];//does not work without [] ?!?! 
  foreach(){
    int l=0;
    foreach_layer(){
      if(l==0){
        dz[]  =  RATECONSTANT*c0[]*c1[];
        dc0[] = -RATECONSTANT*c0[]*c1[];
        dc1[] = -RATECONSTANT*c0[]*c1[];
      }
      l++;
    }
    c0[] += dt*dc0[];
    c1[] += dt*dc1[];
    zb[] += dt*dz[];

    foreach_layer()
      h[] = (eta[] - zb[])/nl;

  }
  boundary({zb,h});
}

/**
## Initial condition: step functions
*/

event init(i=0){
  ic = fopen(name_init,"a");

  foreach(){
    zb[]=0.;
    altitude =0.;
    foreach_layer(){
      h[] = HINIT/nl;
      altitude += h[]/2;
      c1[] = (x < 0);
      c0[] = (x > 0);
      fprintf(ic,"%g %g %g %g\n",x, altitude, c0[], c1[]);
      altitude += h[]/2;
    }
    eta[]=zb[]+altitude;
  }
  fclose(ic);
}


/**
## Outputs
*/
event output_data(t += .001, t<= TMAX){
  bulk = fopen(name_bulk,"a");
  hydro= fopen(name_hydro,"a");

  foreach(){
    altitude =0.;
    int l = 0;
    foreach_layer(){
      altitude += h[]/2;
      if(l == 0)
        fprintf(bulk,"%g %g %g %g %g\n",t, x, altitude, c0[], c1[]);

      altitude += h[]/2;
      l++;
    }
    fprintf(hydro,"%g %g %g %g\n", t, x, altitude, zb[]);
  }
  fclose(bulk);
  fclose(hydro);
}

event growth(i += 50){
  ch = fopen("chimney.txt","a");

  double h_ch = 0.;
  foreach()
    h_ch = h_ch*(h_ch >= zb[]) + zb[]*(h_ch < zb[]);
  fprintf(ch,"%g %g\n", t, h_ch);
  fclose(ch);
}

event output_progress(i += 50){
  fprintf(stderr, "i=%d, t= %g\n",i,t);

}

/**
## Stop condition
As a safety, we want to stop the simulation when max(zb[])= 0.9
*/
event stop(i++){
  stats s = statsf(zb);

  if(s.max >= 0.9 )
    return 1;
}