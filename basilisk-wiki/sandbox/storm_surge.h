/**
# A solver for the Saint-Venant equations

The
[Saint-Venant equations](http://en.wikipedia.org/wiki/Shallow_water_equations)
can be written in integral form as the hyperbolic system of
conservation laws
$$
  \partial_t \int_{\Omega} \mathbf{q} d \Omega =
  \int_{\partial \Omega} \mathbf{f} (
  \mathbf{q}) \cdot \mathbf{n}d \partial
  \Omega - \int_{\Omega} hg \nabla z_b
$$
where $\Omega$ is a given subset of space, $\partial \Omega$ its boundary and
$\mathbf{n}$ the unit normal vector on this boundary. For
conservation of mass and momentum in the shallow-water context, $\Omega$ is a
subset of bidimensional space and $\mathbf{q}$ and
$\mathbf{f}$ are written
$$
  \mathbf{q} = \left(\begin{array}{c}
    h\\
    hu_x\\
    hu_y
  \end{array}\right), 
  \;\;\;\;\;\;
  \mathbf{f} (\mathbf{q}) = \left(\begin{array}{cc}
    hu_x & hu_y\\
    hu_x^2 + \frac{1}{2} gh^2 & hu_xu_y\\
    hu_xu_y & hu_y^2 + \frac{1}{2} gh^2
  \end{array}\right)
$$
where $\mathbf{u}$ is the velocity vector, $h$ the water depth and
$z_b$ the height of the topography. See also [Popinet, 
2011](/src/references.bib#popinet2011) for a more detailed
introduction.

## User variables and parameters

The primary fields are the water depth $h$, the bathymetry $z_b$ and
the flow speed $\mathbf{u}$. $\eta$ is the water level i.e. $z_b +
h$. Note that the order of the declarations is important as $z_b$
needs to be refined before $h$ and $h$ before $\eta$. */

#include "predictor-corrector.h"

static double update_ss (scalar * evolving, scalar * updates, double dtmax);

event defaults (i = 0){
  update = update_ss;
}
#include "saint-venant.h"

scalar pa[],fcor[];
vector ts[];
double Pa2m = 0.00009916;

/**
### Computing fluxes

Various approximate Riemann solvers are defined in [riemann.h](). */

#include "riemann.h"

double update_ss (scalar * evolving, scalar * updates, double dtmax)
{

/**
We first recover the currently evolving fields (as set by the
predictor-corrector scheme). */

  scalar h = evolving[0];
  vector u = { evolving[1], evolving[2] };

/**
The gradients are stored in locally-allocated fields. First-order
reconstruction is used for the gradient fields. */

  vector gh[], geta[], gpa[];
  tensor gu[];
  for (scalar s in {gh, geta, gpa, gu}){
    s.gradient = none;
    #if QUADTREE
      s.prolongation = refine_linear;
    #endif
  }
  gradients ({h, eta, pa, u}, {gh, geta, gpa, gu});

/**
`Fh` and `Fq` will contain the fluxes for $h$ and $h\mathbf{u}$
respectively and `S` is necessary to store the asymmetric topographic
source term. */

  vector Fh[], S[];
  tensor Fq[];

/**
The faces which are "wet" on at least one side are traversed. */

  foreach_face (reduction (min:dtmax)) {
    double hi = h[], hn = h[-1,0];
    if (hi > dry || hn > dry) {

/**
#### Left/right state reconstruction

The gradients computed above are used to reconstruct the left and
right states of the primary fields $h$, $\mathbf{u}$, $z_b$. The
"interface" topography $z_{lr}$ is reconstructed using the hydrostatic
reconstruction of [Audusse et al, 2004](/src/references.bib#audusse2004) */
    
      double dx = Delta/2.;
      double zi = eta[] - hi + Pa2m * pa[];
      double zl = zi - dx*(geta.x[] - gh.x[] + Pa2m *  gpa.x[]); 
      double zn = eta[-1,0] - hn + Pa2m *  pa[-1,0];
      double zr = zn + dx*(geta.x[-1,0] - gh.x[-1,0] + Pa2m *  gpa.x[-1,0]);
      double zlr = max(zl, zr);
      
      double hl = hi - dx*gh.x[];
      double up = u.x[] - dx*gu.x.x[];
      double hp = max(0., hl + zl - zlr);
      
      double hr = hn + dx*gh.x[-1,0];
      double um = u.x[-1,0] + dx*gu.x.x[-1,0];
      double hm = max(0., hr + zr - zlr);

/**
#### Riemann solver

We can now call one of the approximate Riemann solvers to get the fluxes. */

      double fh, fu, fv;
      kurganov (hm, hp, um, up, Delta*cm[]/fm.x[], &fh, &fu, &dtmax);
      fv = (fh > 0. ? u.y[-1,0] + dx*gu.y.x[-1,0] : u.y[] - dx*gu.y.x[])*fh;

/**
#### Topographic source term

In the case of adaptive refinement, care must be taken to ensure
well-balancing at coarse/fine faces (see [notes/balanced.tm]()). */

#if QUADTREE
      if (!is_active(cell) && is_active(aparent(0,0))) {
	hi = coarse(h,0,0);
	zi = coarse(zb,0,0) + Pa2m * coarse(pa,0,0);
      }
      if (!is_active(neighbor(-1,0)) && is_active(aparent(-1,0))) {
	hn = coarse(h,-1,0);
	zn = coarse(zb,-1,0) + Pa2m * coarse(pa,-1,0);
      }
#endif
      double sl = G/2.*(sq(hp) - sq(hl) + (hl + hi)*(zi - zl));
      double sr = G/2.*(sq(hm) - sq(hr) + (hr + hn)*(zn - zr));

/**
#### Flux update */

      Fh.x[]   = fm.x[]*fh;
      Fq.x.x[] = fm.x[]*(fu - sl);
      S.x[]    = fm.x[]*(fu - sr);
      Fq.y.x[] = fm.x[]*fv;
    }
    else // dry
      Fh.x[] = Fq.x.x[] = S.x[] = Fq.y.x[] = 0.;
  }
  
  boundary_flux ({Fh, S, Fq});

/**
#### Updates for evolving quantities

We store the divergence of the fluxes in the update fields. Note that
these are updates for $h$ and $h\mathbf{u}$ (not $\mathbf{u}$). */

  scalar dh = updates[0];
  vector dhu = { updates[1], updates[2] };
  scalar wet[];
  vector cross={ 1.,-1.};
  foreach()
    wet[] = h[] > dry;
  boundary ({wet});
  
  foreach() {
    dh[] = (Fh.x[] + Fh.y[] - Fh.x[1,0] - Fh.y[0,1])/(cm[]*Delta);
    foreach_dimension() {
        if (wet[-1,0] == 1 && wet[] == 1 && wet[1,0] == 1)
	  dhu.x[] = (Fq.x.x[] + Fq.x.y[] - S.x[1,0] - Fq.x.y[0,1])/(cm[]*Delta)+ts.x[]+cross.x*fcor[]*h[]*u.y[];
        else
	  dhu.x[] = (Fq.x.x[] + Fq.x.y[] - S.x[1,0] - Fq.x.y[0,1])/(cm[]*Delta);
	}
    //dhu.x[]+=fcor[]*h[]*u.y[];
    //dhu.y[]+=-fcor[]*h[]*u.x[];

    double dmdl = (fm.x[1,0] - fm.x[])/(cm[]*Delta);
    double dmdt = (fm.y[0,1] - fm.y[])/(cm[]*Delta);
    double fG = u.y[]*dmdl - u.x[]*dmdt;
    dhu.x[] += h[]*(G*h[]/2.*dmdl + fG*u.y[]);
    dhu.y[] += h[]*(G*h[]/2.*dmdt - fG*u.x[]);
  }

  return dtmax;
}
