/**

We want to perform mesh adaptation when solving the laplace equation for the pressure on a domain corresponding to a "liquid phase" around a "steady spherical bubble".
We only consider the liquid phase bounded by embed boundaries for the bubble interface: in other words, we modelize the bubble as a rigid cylinder.
The pressure at infinity is $p_{inf}$, and $p_0$ in the bubble.

The objective is to show that the user-interface we provide with our new adaptation method is simple to use.
An additionnal note is that even in a "simple" case, the adapted mesh is not necessarily trivial.

 */


#undef dv
#define dv() (y*sq(Delta)*cm[])  // comment for no axi case.  axi and embed not-compatible when this case was created, thus we manually implement the axi-symmetry

#include "embed.h"
#include "poisson.h" // solver

#include "./../../AMR_tools/amr.h"  // AMR

double pinf = 1.;
double p0 = 0.;

scalar psi[];
psi[top] = dirichlet (pinf);
psi[right] = dirichlet (pinf);  
psi[left] = neumann(0.);
psi[bottom] = neumann(0.);
psi[embed] = dirichlet (p0);


/** 

Bubble (= cylinder) embedded boundary.

*/

void circle_bndry(){  

  double x_c=0.;
  double R=0.1;
  
   vertex scalar phi[];
   foreach_vertex() {
      phi[] = sq(x - x_c) + sq(y) - sq(R);
   }
   boundary ({phi});
   fractions (phi, cs, fs);

   cs.refine = embed_fraction_refine;
   cs.prolongation = fraction_refine;
   foreach_dimension()
      fs.x.prolongation = embed_face_fraction_refine_x;
   boundary ({cs,fs});
   restriction ({cs,fs});
   cm = cs;
   fm = fs;
}


/**

Numerical solution

*/


void compute_soluce(){   // compute numerical solution
   foreach()
      psi[]=p0;
   psi.third = true;
   boundary({psi});  


/**

In axisymmetric (3D cylindrical + axisymmetrical hypothesis), the poisson equation $\Delta p = 0$ writes

$$ \frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial p}{\partial r} \right) + \frac{\partial^2 p}{\partial z^2} = 0$$

which may be rewritten

$$ \frac{\partial}{\partial r}\left(r\frac{\partial p}{\partial r}
\right) + \frac{\partial}{\partial z}\left(r\frac{\partial p}{\partial z}
\right) = 0$$

Note: I prefer see that using Finite Volume integration, but the result is the same.
*/
  
   face vector alphav[];

   foreach_face() {
      if (fm.x[]==0.)
        alphav.x[] = y;   
      else
        alphav.x[] = fm.x[]*y; // =fm.x[] for no axi case
   }

   boundary((scalar *){alphav});

   scalar src[];
   foreach ()
      src[] = 0.;
   boundary({src});

   poisson (psi, src, alphav, tolerance = 1e-5);   // poisson solver for pressure
}








int main() {

  L0=1.;
  
  init_grid (1 << 6);


/**

# Mesh adaptation

We perform the adaptation loop.

*/
  
  for (int k=0 ; k<= 20 ; k++){
    circle_bndry();      
    compute_soluce();

/**
The mesh adaptation user-interface is very simple to use
*/
    AMReps=6e-6;   // epsilon criterion to decide which element needs refining/coarsening
    adapt_metric( {psi} );  // adaptation procedure
  }

    
/**
We ouput the elements of the mesh inside of the considered computationnal domain.
*/

  FILE*fp=fopen("mesh","w");
  foreach(){
    if (cs[]){
    fprintf(fp,"%g %g\n", x-Delta/2., y-Delta/2.);
    fprintf(fp,"%g %g\n", x-Delta/2., y+Delta/2.);
    fprintf(fp,"%g %g\n", x+Delta/2., y+Delta/2.);
    fprintf(fp,"%g %g\n", x+Delta/2., y-Delta/2.);
    fprintf(fp,"%g %g\n\n", x-Delta/2., y-Delta/2.);
    }
  }
  fclose(fp);
    

  free_grid();
 
}

/**

The optimal mesh for the axi-symmetric case is not trivial, we clearly see that the element are more refined for increaseing $y$ than for increasing $x$.

~~~gnuplot resulting adapted mesh for the axi-symmetric case
set terminal pngcairo size 700,700

set output "mesh.png"

p "mesh" w l not
~~~

As a comparison, the mesh obtained without axi-symmetry is closer to a mesh a user would create manually.

~~~gnuplot resulting adapted mesh
set terminal pngcairo size 700,700

set output "mesh_non_axi.png"

p "mesh_non_axi" w l not
~~~

*/
