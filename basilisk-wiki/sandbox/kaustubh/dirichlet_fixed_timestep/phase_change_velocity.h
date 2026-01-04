/**
# Phase change velocity calculation

The following set of functions modifies the phase change velocity field 
$\mathbf{v}_{pc} = v_{pc} \mathbf{n}_{2\rightarrow 1}$ in interfacial cells 
and sets it to :
$$
\mathbf{v}_{pc} = \frac{1}{L_H}\left(\lambda_{1}\left.\nabla tr_{1}\right|_
{\Gamma} - \lambda_{2}\left.\nabla tr_{2}\right|_{\Gamma}\right) 
$$
which is the Stefan condition, where $L_H$ is the latent heat, $tr_{1}$ and $tr_{2}$ 
are scalars defined in phase 1 and phase 2 respectively and 
$\lambda_i$ are the thermal diffusivities.

## Anisotropy of the Gibbs-Thomson relation

We define 2 functions to take into account the anisotropy in the change velocity
where we set the following condition :
$$
\epsilon = \overline{\epsilon} \left( 1. - \alpha \cos(n \theta +
\theta_0)\right)
$$
where $\theta$ is the angle between the interface and the x-axis, $\alpha =0.5$
and $\theta_0 = 0$ are hardcoded. This will be used in the Gibbs-Thomson
formulation.
*/

scalar Teq[];


#ifndef ANISO // if the user defines ANISO, then he wants to use standard
// anisotropic function, else he has to define its own fac1() function before
// the <<#include "phase_change_velocity.h">>
double fac1(coord n, double eps4){
  if(eps4==0.)return 1.;
  double sum = 0;

  foreach_dimension()
    sum += eps4*powf(n.x,4.); // taken from doi.org/10.1016/j.jcrysgro.2010.11.013
                              // simple fourfold anisotropy.
#if dimension ==2
    double theta = atan2(n.y, n.x);
  return 1.-15.*eps4*cos(4.*theta);
#else // DIMENSION == 3 && ANISO != 0
  return (1.-3.*eps4+4*eps4*sum);
#endif
}
#endif // ifndef ANISO
/**
For the gradient calculation on the interface, we need the temperature on the
interface $T_{\Gamma}$, which is defined either by :
$$
 T_{\Gamma} = T_{m} - \epsilon_{\kappa} * \kappa - \epsilon_{v} v_{pc} = T_m -
  \left(\overline{\epsilon_{\kappa}} * \kappa - \overline{\epsilon_{v}} 
  v_{pc}\right)(1-0.5 cos(n \theta))
$$
which is the classical Gibbs-Thomson equation.
*/

double Temp_GT(Point point, double epsK, double epsV, vector vpc,
  scalar curve, face vector fs, scalar cs, double eps4){
  if(epsK==0. && epsV ==0.)return 0.;
  
  double pref = 1.;
  if(eps4 > 0.)pref = fac1(facet_normal( point, cs ,fs),eps4);
  if(epsV==0.)return epsK*curve[]*pref;
  else
  {
    double temp=0.;
    foreach_dimension()
      temp += sq(vpc.x[]);
    return  (epsK*curve[] - epsV*sqrt(temp))*pref;
  }
}

/**
## Velocity in interfacial cells

This function modifies the `vector v_pc[]` in interfacial cells using the Stefan
condition. Note that this procedure is not sufficient to obtain a proper phase
change velocity field used for a level-set, it must be complemented with a
reconstruction procedure (see [this page](LS_recons.h))*/

void phase_change_velocity_LS_embed (
 scalar cs, face vector fs, scalar tr,
 scalar tr2, double T_eq, vector v_pc, 
 double L_H, double lambda[2], double epsK, 
 double epsV, double eps4, scalar curve) {


/**
We store in `scalar Teq[]` the temperature on the interface.
*/
#if Gibbs_Thomson

  foreach(){
    // here we suppose that Tm = 0
    Teq[] =  (interfacial(point, cs) ? T_eq + Temp_GT(point, epsK, epsV, v_pc,
         curve, fs, cs, eps4) : nodata);
  }

#else
  foreach(){
    Teq[] = (interfacial(point, cs) ? T_eq : nodata);
  }
#endif

  boundary({Teq});
  restriction({Teq});

/**
To calculate the gradient $\left. \nabla tr \right|_{\Gamma}$, we use the
embed_gradient_face_x defined in embed that gives an accurate definition of the
gradients with embedded boundaries. */
  vector gtr[], gtr2[];
  boundary({tr});

  foreach(){
    foreach_dimension(){
      gtr.x[] = 0.;
    }
    if(interfacial(point,cs)){
      coord n, p;
      embed_geometry(point, &p, &n);
      double c    = 0.;
      double temp = Teq[];
      double grad = dirichlet_gradient(point, tr, cs , n, p, 
        temp, &c);

/**
For degenerate cases (stencil is too small), we add the term $tr[]*c$ to the 
gradient.
*/
      foreach_dimension(){
        gtr.x[] += grad*n.x+tr[]*c;
      }
    }
  }

  boundary((scalar*){gtr});

/**
Now we need to change the volume fraction and face fractions to calculate the
gradient in the other phase.
*/
  foreach(){
    cs[]      = 1.-cs[];
  }
  foreach_face()
  fs.x[]      = 1.-fs.x[];

  boundary({cs});
  restriction({cs});

  boundary({tr2});
  vector p_sauv[], n_sauv[];

  foreach(){
    foreach_dimension(){
      gtr2.x[] = 0.;
    }
    if(interfacial(point,cs)){
      coord n, p;
      embed_geometry(point, &p, &n);
      double c=0.;
      double temp = Teq[];
      double grad = dirichlet_gradient(point, tr2, cs , n, p, 
        temp, &c);
      foreach_dimension(){
        gtr2.x[] += grad*n.x+tr2[]*c;
      }
    }
  }
  boundary((scalar*){gtr2});
  foreach(){
    cs[]      = 1.-cs[];
  }
  foreach_face()
  fs.x[]      = 1.-fs.x[];
  boundary({cs,fs});
  restriction({cs,fs});

  /**
  With the the normal vector and the gradients of the tracers we can now 
  compute the phase change velocity $\mathbf{v}_{pc}$. */

  foreach()
    foreach_dimension()
      v_pc.x[]  = 0.; 

  foreach()
    if(interfacial(point, cs))
      foreach_dimension()
        v_pc.x[]     =  (lambda[1]*gtr2.x[] - lambda[0]*gtr.x[])/L_H;
  boundary((scalar *){v_pc});
}
