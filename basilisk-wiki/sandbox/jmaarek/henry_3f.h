/**
# Advection/diffusion of a soluble tracer

We consider the transport and diffusion of a tracer $c$ with different
solubilities in the two-phases described by
[two-phase.h](/src/two-phase.h).

The diffusion coefficients in the two phases are $D_1$ and $D_2$ and
the jump in concentration at the interface is given by
$$
c_1 = \alpha c_2
$$

The one fluid formulation for the concentration can be written as

$$
c = \chi c_1 + (1-\chi)c_2
$$


The advection/diffusion equation for $c$ can then be written
$$
\partial_t c + \nabla\cdot(\mathbf{u} c) = 
   \nabla\cdot\left(D\nabla c - D \frac{c(\alpha - 1)}
                     {\alpha \chi + (1 - \chi)}\nabla \chi\right)
$$
with $\chi$ the volume fraction and
$$
D = \frac{D_1 D_2}{D_2 f + D_1 (1 - f)}
$$
the harmonic mean of the diffusion coefficients (see [Haroun et al.,
2010](#haroun2010)). 

The diffusion coefficients and solubility are attributes of each tracer.

The *stracers* list of soluble tracers must be defined by the calling code. */

attribute {
  double D1, D2, D3, alpha1, alpha2;  //D1 for f = 1, D2 for f = 0
  scalar phi1, phi2, phi3; // private
}

extern scalar * stracers;
extern scalar f10;

/**
## Defaults

On trees we need to ensure conservation of the tracer when
refining/coarsening. */

#if TREE
event defaults (i = 0)
{
  for (scalar s in stracers) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
  }
}
#endif // TREE

/**
## Advection

To avoid numerical diffusion through the interface we use the [VOF
tracer transport scheme](/src/vof.h) for the temporary fields
$\phi_1$ and $\phi_2$, see section 2.3 of [Farsoiya et al.,
2020](#farsoiya2020). */

static scalar * phi1_tracers = NULL;
static scalar * phi2_tracers = NULL;
static scalar * phi3_tracers = NULL;

double alpha1_temp;
double alpha2_temp;


event vof (i++)
{
//fprintf(fout, "test 0"); 

  phi1_tracers = f1.tracers;
  phi2_tracers = f2.tracers;
  phi3_tracers = f3.tracers;
  
  for (scalar c in stracers) {
    scalar phi1 = new scalar, phi2 = new scalar, phi3 = new scalar;
    c.phi1 = phi1, c.phi2 = phi2, c.phi3 = phi3;
    scalar_clone (phi1, c);
    scalar_clone (phi2, c);
    scalar_clone (phi3, c);
    
    f1.tracers = list_append (f1.tracers, phi1);
    f2.tracers = list_append (f2.tracers, phi2);
    f3.tracers = list_append (f3.tracers, phi3);

    /**
    $\phi_1$ and $\phi_2$ are computed from $c$ as
    $$
    \phi_1 = c \frac{\alpha f}{\alpha f + (1 - f)}
    $$
    $$
    \phi_2 = c \frac{1 - f}{\alpha f + (1 - f)}
    $$
    */
		  
    foreach() {
      phi1[] = c[]*f1[]/(f1[]+f2[]/c.alpha1+f3[]/c.alpha2);
      phi2[] = c[]*f2[]/(f1[]*c.alpha1 + f2[] + f3[]*c.alpha1/c.alpha2);
      phi3[] = c[]*f3[]/(f1[]*c.alpha2 + f2[]*c.alpha2/c.alpha1 + f3[]);
    }
    
    alpha1_temp = c.alpha1;
    alpha2_temp = c.alpha2;
    
    
    phi1[bottom] = (1.-f10[])*f1[]/(f1[]+f2[]/alpha1_temp+f3[]/alpha2_temp);
    phi2[bottom] = (1.-f10[])*f2[]/(f1[]*alpha1_temp + f2[] + f3[]*alpha1_temp/alpha2_temp);
    phi3[bottom] = (1.-f10[])*f3[]/(f1[]*alpha2_temp + f2[]*alpha2_temp/alpha1_temp + f3[]);
    
    boundary ({phi1, phi2, phi3});
        
  } 
}

/**
## Diffusion

The advected concentration is computed from $\phi_1$ and $\phi_2$ as
$$
c = \phi_1 + \phi_2 + \phi_3
$$
and these fields are then discarded. */

#include "diff_implicit_tracer.h"

event tracer_diffusion (i++)
{
  free (f1.tracers);
  free (f2.tracers);
  free (f3.tracers);
  f1.tracers = phi1_tracers;
  f2.tracers = phi2_tracers;
  f3.tracers = phi3_tracers;
  double summ;
  for (scalar c in stracers) {	
    scalar phi1 = c.phi1, phi2 = c.phi2, phi3 = c.phi3;
    foreach(){
      summ = (f1[]+f2[]+f3[]);
      //as phi represents amount of tracer in phase it must be positive
      phi1[] = (phi1[] > 0.0 ? phi1[] : 0.);
      phi2[] = (phi2[] > 0.0 ? phi2[] : 0.);
      phi3[] = (phi3[] > 0.0 ? phi3[] : 0.);
      //concentration of tracer in phase remains constant but mass is adjusted by change in mass of tracer
      if (summ > 0.){ 
      	phi1[] = (f1[] > 0.0 ? phi1[]*clamp(f1[]/summ,0.,1.)/f1[] : 0.);
      	phi2[] = (f2[] > 0.0 ? phi2[]*clamp(f2[]/summ,0.,1.)/f2[] : 0.);
      	phi3[] = (f3[] > 0.0 ? phi3[]*clamp(f3[]/summ,0.,1.)/f3[] : 0.);}
        	      
    c[] = phi1[] + phi2[] + phi3[];}
    delete ({phi1, phi2, phi3});

    /*
    The diffusion equation for $c$ is then solved using the implicit
    discretisation
    $$
    \frac{c^{n + 1} - c^n}{\Delta t} = 
    \nabla\cdot (D \nabla c^{n + 1} + \beta c^{n + 1})
    $$
    see section 2.2 of [Farsoiya et al., 2020](#farsoiya2020). */
    
    face vector D[], beta[];
    foreach_face() {
      
      double summ_0 = (f1[-1]+f2[-1]+f3[-1]);
      double summ_1 = (f1[]+f2[]+f3[]);
      
      double f1_0 = (summ_0 > 0.0 ? clamp(f1[-1]/summ_0, 0.,1.) : f1[-1]);
      double f1_1 = (summ_1 > 0.0 ? clamp(f1[]/summ_1, 0.,1.) : f1[]);
      double f2_0 = (summ_0 > 0.0 ? clamp(f2[-1]/summ_0, 0.,1.) : f2[-1]);
      double f2_1 = (summ_1 > 0.0 ? clamp(f2[]/summ_1, 0.,1.) : f2[]);
      double f3_0 = (summ_0 > 0.0 ? clamp(f3[-1]/summ_0, 0.,1.) : f3[-1]);
      double f3_1 = (summ_1 > 0.0 ? clamp(f3[]/summ_1, 0.,1.) : f3[]);
      
    
      double ff1 = (f1_0+f1_1)/2;
      double ff2 = (f2_0+f2_1)/2;
      double ff3 = (f3_0+f3_1)/2;
      
      ff1 /= (ff1 + ff2 + ff3);
      ff2 /= (ff1 + ff2 + ff3);
      ff3 /= (ff1 + ff2 + ff3);
      
      D.x[] = fm.x[]/(ff1/c.D1 + ff2/c.D2 + ff3/c.D3);
      
      double temp = 0;
      if (ff1 > 0.0)
      	temp += (1/c.alpha1 - 1.)*ff1/(ff2/c.alpha1 + ff1)*(f2_1 - f2_0) + (1/c.alpha2 - 1.)*ff1/(ff3/c.alpha2 + ff1)*(f3_1 - f3_0);
      if (ff2 > 0.0)
      	temp += (c.alpha1 - 1.)*ff2/(ff1*c.alpha1 + ff2)*(f1_1 - f1_0) + (c.alpha1/c.alpha2 - 1.)*ff2/(ff3*c.alpha1/c.alpha2 + ff2)*(f3_1 - f3_0);
      if (ff3 > 0.0)
      	temp += (c.alpha2 - 1.)*ff3/(ff1*c.alpha2 + ff3)*(f1_1 - f1_0) + (c.alpha2/c.alpha1 - 1.)*ff3/(ff2*c.alpha2/c.alpha1 + ff3)*(f2_1 - f2_0);
      
      beta.x[] = - D.x[]/Delta*temp;

      
    }
    boundary ({c, D, beta});
    
    tr_diffusion (c, dt, D, beta);
  }
}

/**
## References

~~~bib
@article{haroun2010,
  title = {Volume of fluid method for interfacial reactive mass transfer: 
           application to stable liquid film},
  author = {Haroun, Y and Legendre, D and Raynal, L},
  journal = {Chemical Engineering Science},
  volume = {65},
  number = {10},
  pages = {2896--2909},
  year = {2010},
  doi = {10.1016/j.ces.2010.01.012}
}

@article{farsoiya2020,
  title = {Bubble mediated gas transfer of diluted component in turbulence},
  author = {P. K. Farsoiya and S. Popinet and L. Deike},
  journal = {Journal of Fluid Mechanics},
  year = {2020},
  note = {submitted}
}
~~~
*/


