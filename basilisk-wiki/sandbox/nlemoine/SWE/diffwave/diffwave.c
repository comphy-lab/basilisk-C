/**
[Return to my homepage](http://basilisk.fr/sandbox/nlemoine/README)

# Diffusive wave approximation to the Shallow Water Equations

## General formulation

The diffusive wave approximation stems from the Saint-Venant Shallow Water Equations, which can be written in full form as:

$$
\displaystyle\frac{\partial}{\partial t} \left[\begin{array}{c} h \\ h\,u_x \\ h\,u_y \end{array}\right] = -\boldsymbol{\nabla}\cdot
\left[\begin{array}{cc} h\,u_x & h\,u_y \\ h\,u_x^2+ \frac{1}{2} g h^2 & h\,u_x\,u_y \\ h\,u_x\,u_y & h\,u_y^2 + \frac{1}{2} g h^2\end{array}\right]
\ -\ g\,h\left[\begin{array}{c} 0 \\ \partial_x z_b \\ \partial_y z_b \end{array}\right]
\ -\ \frac{1}{\rho}\left[\begin{array}{c} 0 \\ \tau_x \\ \tau_y \end{array}\right]$$

The first line is the continuity equation, while the second and third are the projections of the dynamic (momentum) equation. The diffusive wave approximation is a zero-inertia approximation, obtained by neglecting the inertial acceleration terms in the left-hand side and the convective acceleration terms in the right-hand side:

$$\displaystyle\left[\begin{array}{c} \frac{\partial h}{\partial t} \\
\cancel{\frac{\partial(h\,u_x)}{\partial t}} \\
\cancel{\frac{\partial(h\,u_y)}{\partial t}} 
\end{array}\right] = -\boldsymbol{\nabla}\cdot
\left[\begin{array}{cc} h\,u_x & h\,u_y \\ \frac{1}{2} g h^2 & 0 \\ 0 & \frac{1}{2} g h^2\end{array}\right]
-\quad\cancel{
\boldsymbol{\nabla}\cdot
\left[\begin{array}{cc} 0 & 0 \\ h\,u_x^2 & h\,u_x\,u_y \\ h\,u_x\,u_y & h\,u_y^2 \end{array}\right]
}
\ -\ g\,h\left[\begin{array}{c} 0 \\ \partial_x z_b \\ \partial_y z_b \end{array}\right]
\ -\ \frac{1}{\rho}\left[\begin{array}{c} 0 \\ \tau_x \\ \tau_y \end{array}\right]$$

Denoting $q_x = h\,u_x$ and $q_y = h\,u_y$ the components of the flux density vector, the system boils down to:

$$\begin{cases}
   \frac{\partial h}{\partial t} = -\boldsymbol{\nabla}\cdot\mathbf{q} & \\
   0 = -g h \frac{\partial h}{\partial x} -g h \frac{\partial z_b}{\partial x} - \frac{1}{\rho}\tau_x & \\
   0 = -g h \frac{\partial h}{\partial y} -g h \frac{\partial z_b}{\partial y} - \frac{1}{\rho}\tau_y & \\
   \end{cases}
$$

The dynamic equation (last two scalar equations) is limited to a quasi-equilibrium between friction and the resultant force of the pressure gradient. Indeed, denoting $\eta = z_b + h$ the water surface elevation, we have:

$$\begin{cases}
   \frac{\partial h}{\partial t} = -\boldsymbol{\nabla}\cdot\mathbf{q} & \\
   0 = -g h \frac{\partial \eta}{\partial x} -\frac{1}{\rho}\tau_x & \\
   0 = -g h \frac{\partial \eta}{\partial y} - \frac{1}{\rho}\tau_y & \\
   \end{cases}
$$

In the following we chose the Manning friction model, but any other model can be chosen. We then have $\quad\boldsymbol{\tau} = \rho\,g\,n^2\,h^{-\frac{7}{3}}\|\mathbf{q}\|\mathbf{q}\quad$ and then:

$$\begin{cases}
   q_x = -\displaystyle\frac{1}{n^2\|\mathbf{q}\|}h^{\frac{10}{3}}\frac{\partial\eta}{\partial x} & \\
   & \\
   q_y = -\displaystyle\frac{1}{n^2\|\mathbf{q}\|}h^{\frac{10}{3}}\frac{\partial\eta}{\partial y} & \\
   \end{cases}
$$
Hence
$$\|\mathbf{q}\| = \displaystyle\frac{1}{n^2\|\mathbf{q}\|}h^{\frac{10}{3}}\|\boldsymbol{\nabla}\eta\|$$

We see that we now have an explicit expression for the norm of the flux density vector, and hence for the vector itself:

$$\|\mathbf{q}\| = \displaystyle\frac{1}{n}h^{\frac{5}{3}}\|\boldsymbol{\nabla}\eta\|^{\frac{1}{2}}$$

$$\mathbf{q} = \displaystyle-\frac{1}{n}\big(\eta-z_b\big)^{\frac{5}{3}}\|\boldsymbol{\nabla}\eta\|^{-\frac{1}{2}}\boldsymbol{\nabla}\eta$$

Finally we are left with a single equation (the continuity equation) which can be written in the form of a diffusion equation:

$$\frac{\partial\eta}{\partial t} = +\boldsymbol{\nabla}\cdot\left(
\underbrace{\frac{1}{n}\big(\eta-z_b\big)^{\frac{5}{3}}\|\boldsymbol{\nabla}\eta\|^{-\frac{1}{2}}}_{
\normalsize D\left(\eta,\|\boldsymbol{\nabla}\eta\|\right)
}\boldsymbol{\nabla}\eta
\right)$$

The problem is nonlinear since the diffusivity $D$ depends both on $\eta$ and on the norm of the gradient of $\eta$. However it can be easily implemented in Basilisk thanks to the [Poisson solver](basilisk.fr/src/diffusion.h), using an implicit scheme with a [staggered grid](http://basilisk.fr/Basilisk%20C#face-and-vertex-fields) diffusivity:

![Centered, face and vertex staggering. The diffusivity field uses the second type.](http://basilisk.fr/src/figures/staggering.svg)

It is necessary to desingularize the norm of the gradient of $\eta$, which has exponent $-\frac{1}{2}$ in the diffusivity. We use:

$$ D\left(\eta,\|\boldsymbol{\nabla}\eta\|\right) = \frac{1}{n}\big(\eta-z_b\big)^{\frac{5}{3}}\Big(\|\boldsymbol{\nabla}\eta\|+\epsilon\Big)^{-\frac{1}{2}} $$

## Formulation with detrended topography

Just like for the full Saint-Venant system, solving the diffusive wave equation using a [detrended topography](../manning-tilt.h) requires some rewriting of the problem. Denoting $\mathbf{I}_\mathrm{reg} = (I_x,I_y)^T$ the regional ``tilt'' driving the flow, we have:
$$\begin{cases}
   \boldsymbol{\nabla}\eta = \boldsymbol{\nabla}\eta' - \mathbf{I}_\mathrm{reg} & \\
   & \\
   h = (\eta-z_b) = (\eta'-z_b') & \\
   \end{cases}
$$

Hence:
$$\frac{\partial\eta'}{\partial t} = +\boldsymbol{\nabla}\cdot\left(\underbrace{\frac{1}{n}\big(\eta'-z_b'\big)^{\frac{5}{3}}\|\boldsymbol{\nabla}\eta'-\mathbf{I}_\mathrm{reg}\|^{-\frac{1}{2}}}_{\normalsize D\left(\eta',\|\boldsymbol{\nabla}\eta'\|\right)}\left({\boldsymbol{\nabla}\eta'-\mathbf{I}_\mathrm{reg}}\right)\right)=\boldsymbol{\nabla}\cdot\big(D\boldsymbol{\nabla}\eta'\big) - \mathbf{I}_\mathrm{reg}\cdot\boldsymbol{\nabla}D$$

We see that we can keep using the Poisson solver to compute the evolution of the new variable $\eta'$ (detrended water surface elevation) with the detrended topography $z_b'$ as input, provided that we add a ficticious source term given by $-\mathbf{I}_\mathrm{reg}\cdot\boldsymbol{\nabla}D$. The desingularization of the norm of the gradient is of course still required:

$$D\left(\eta',\|\boldsymbol{\nabla}\eta'\|\right) = \frac{1}{n}\big(\eta'-z'_b\big)^{\frac{5}{3}}\Big(\|\boldsymbol{\nabla}\eta'-\mathbf{I}_\mathrm{reg}\|+\epsilon\Big)^{-\frac{1}{2}}$$
and we set $\epsilon = 0.001\ \|\mathbf{I}_\mathrm{reg}\|$.

## Issues with mass conservation

One of the specificities of the diffusive wave formulation is that we keep setting the water surface elevation $\eta$ equal to the topographic elevation $z_b$ even where the surface is dry. It is therefore possible that, over a given time step, the value of $\eta$ changes from $\eta=z_b$ to $\eta<z_b$ (or, equivalenty, from $\eta'=z'_b$ to $\eta'<z'_b$), typically
with a non-zero flux between a wet cell and a neighboring dry cell if the interface diffusivity is non-zero. At the end of the iteration, we reset $\eta\leftarrow\textrm{max}\left\{\eta,z_b\right\}$ but this amounts to artificially create a small volume (per unit area) $\Delta\eta^{(i)} = \left(z_b-\eta^{(i)}\right)$ on the corresponding, dry cell. In order to ensure mass conservation, we store the quantity $\Delta\eta^{(i)}$ in memory (using the scalar field $\texttt{excess[]}$): it will be compensated at the next iteration using a sink term $-\frac{\Delta\eta^{(i)}}{\Delta t}$. We can even benefit from the implicit scheme of the Poisson solver to improve the estimate of the required balancing / compensating term; indeed, if resetting $\eta = z_b$ was necessary for a given cell at a given time step, we can bet that it will again be necessary at the next iteration. The corresponding balancing sink term would then have the value $-\left(\frac{z_b-\eta^{(i+1)}}{\Delta t}\right)$. We can then build a semi-implicit estimate of the sink term, using a parameter $\alpha$ and the boolean $\mathbb{1}_{(\Delta\eta^{(i)}>0)}$ (which has value 1 if $\Delta\eta^{(i)}>0$ on the cell at the previous iteration, and 0 otherwise):

$$r=-\left[(1-\alpha)\underbrace{\frac{\Delta\eta^{(i)}}{\Delta t}}_{\textrm{explicit}}\ +\ \alpha\,\underbrace{\mathbb{1}_{(\Delta\eta^{(i)}>0)}\left(\frac{z_b-\eta^{(i+1)}}{\Delta t}\right)}_{\textrm{implicit}}\right]$$

$$r=-\left[(1-\alpha)\frac{\Delta\eta^{(i)}}{\Delta t}+\alpha\,\mathbb{1}_{(\Delta\eta^{(i)}>0)}\frac{z_b}{\Delta t}\right]\quad + \quad \left[\alpha\,\mathbb{1}_{(\Delta\eta^{(i)}>0)}\frac{1}{\Delta t}\right]\eta^{(i+1)}$$

It yields a reaction term with an affine form, directly handled by the Poisson solver in Basilisk:
$$r = r_0 + \beta\ \eta^{(i+1)}$$
This correction is stable for $\alpha<1$ (a more thorough study of this property would be necessary).

It is worth noting that the situation $\eta<z_b$ can result from the diffusion equation itself, but also from the adaptivity of the grid. It is indeed possible that, after the adaptation step, the surface $\texttt{eta[]}$ be locally below the topographic surface $\texttt{zb[]}$; in the latter case, if the surface is locally dry, the compensating sink term is disabled.

## Basilisk implementation

We test the scheme on a short reach (2 km) of the Garonne River. It is extracted from the larger dataset described in [Le Moine & Mahdade (2021)](https://doi.org/10.3390/rs13214435). The boundary conditions consist in a prescribed stage hydrograph at the inflow section, and a uniform flow condition at the outflow section. The results of the diffusive wave scheme are compared with those of the [full Saint-Venant system](fullsaintvenant.c) with the exact same setting.

![Detrended topo-bathymetry $z_b'$ of the test reach. Color map ranges from -7 m to +4 m relative to bankfull level. The square domain is 2048 x 2048 m.](diffwave/topo.png)
*/

#include "grid/quadtree.h"
#include "run.h"
#include "diffusion.h"
#include "terrain.h"
#include "utils.h"
#include "output.h"

#define ETABF 10.6
#define n_ch 0.03
#define n_fp 0.06

#define sec_per_day 86400.
#define M_PI acos(-1.0)

#define ETAE     1e-2 // error on water elevation (1 cm)
#define HE       1e-2 // error on water depth (1 cm)
#define UE       1e-2 // 0.01 m/s
#define SE       1.0e-6

double G = 9.81, dry = 1e-3;
double prec,facExplicit;
scalar zb[],eta[],h[],nmanning[],l[],r[],beta[],excess[];
face vector D[];
int nsteps,tag;

double dt, tini;
mgstats mgd;

coord tilt;

#define LEVEL 7
#define MINLEVEL 5

/**
### Inflow section: sine-shaped stage hydrograph
*/
double stage_hydrograph (double tt)
{
   return( ETABF-2.0*cos(2.0*M_PI*tt/sec_per_day) );
}

double segment_eta (coord segment[2], scalar h, scalar eta)
{
  coord m = {segment[0].y - segment[1].y, segment[1].x - segment[0].x};
  normalize (&m);
  double tot = 0., wtot = 0.;

  foreach_segment (segment, p) {
    double dl = 0.;
    foreach_dimension() {
      double dp = (p[1].x - p[0].x)*Delta/Delta_x*(fm.y[] + fm.y[0,1])/2.;
      dl += sq(dp);
    }
    dl = sqrt (dl);    
    for (int i = 0; i < 2; i++) {
      coord a = p[i];
      tot += dl/2.*
	interpolate_linear (point, (struct _interpolate)
			    {h, a.x, a.y, 0.})*
	interpolate_linear (point, (struct _interpolate)
				 {eta, a.x, a.y, 0.}) ;

      wtot += dl/2.* interpolate_linear (point, (struct _interpolate)
			    {h, a.x, a.y, 0.});
    }
  }
  // reduction
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &wtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return (tot/wtot);
}

/**
### Diffusivity function */

int update_diffusivity(scalar zb, scalar eta, face vector D, bool write)
{
  double eps_slope = sqrt(sq(tilt.x)+sq(tilt.y))/1000.;

  //    /!\ reminder : D.x and D.y are not co-located
  foreach_face()
  {
    double hl = eta[-1,0]-zb[-1,0] ;
    double hm = (eta[]+eta[-1,0]-zb[]-zb[-1,0])/2.; // depth interpolated at face center
    hm = max(hm,0.);
    double sn = (eta[]-eta[-1,0])/Delta; // slope component along face normal
    double st = (eta[0,-1]+eta[-1,-1]-eta[0,1]-eta[-1,1])/4./Delta;	// slope component transverse to face normal
    sn -= tilt.x ;
    st -= tilt.y ;
    double S = sqrt(sq(sn)+sq(st)) + eps_slope ; // desingularize slope
    D.x[] = 2.*pow(hm,5./3.)/sqrt(S)/(nmanning[]+nmanning[-1,0]) ;
  }

  boundary((scalar *){D});
  return(0);
}

int update_sink(double Dt)
{
  double a = 0.5;
  foreach()  // cell-centered source term compensating for detrending and mass unbalance
  {
    beta[] = excess[]>0. ? a/Dt : 0.;
    r[] = -tilt.x * (D.x[1,0]-D.x[])/Delta - tilt.y * (D.y[0,1]-D.y[])/Delta - (1.-a)*excess[]/Dt;
    r[] += (excess[]>0.) ? -a*zb[]/Dt : 0.;
  }
  boundary({r});
  return(0);
}

FILE * fpcontrol;

int main (int argc, char * argv[])
{

  size (2048.);
  origin (-424.,-2332.);
  init_grid (1 << LEVEL);

  // Get .kdt file
  system("wget https://dropsu.sorbonne-universite.fr/s/sT5GkqPqzsaMYab/download");
  system("mv download garonne_local.kdt");
  // Get .pts file
  system("wget https://dropsu.sorbonne-universite.fr/s/5TWiFJaAXGKcb6L/download");
  system("mv download garonne_local.pts");
  // Get .sum file
  system("wget https://dropsu.sorbonne-universite.fr/s/sz5qF6fc7MwkPzz/download");
  system("mv download garonne_local.sum");
  
  G = 9.81;
  tilt.x = 0.;
  tilt.y = 1.103e-3;

  facExplicit = 10.;

  run();

  if(pid()==0.)
   fclose(fpcontrol);

}

event initialize (i=0)
{
    terrain(zb,"garonne_local", NULL);
    scalar zmin = zb.dmin;
  
    output_ppm (zb, min = 3.6, max = 14.6, file = "topo.png",
		n = 512, linear = true);
    
   foreach()
   {
     eta[] = fmax(zb[],stage_hydrograph(0.));
     h[] = eta[]-zb[];
     double zzmin = zmin[]<nodata ? zmin[] : zb[];
     nmanning[] = zzmin<ETABF ? n_ch : n_fp;
     excess[] = 0.;
   } 
   boundary({zb,eta,excess});   
   (void) update_diffusivity(zb,eta,D,true);
   (void) update_sink(1.0);   

  if(pid()==0)
   fpcontrol = fopen("outlet_info.txt","w");

}

int adapt();

event integration (i++;t<86400.)
{
  foreach_boundary(bottom)
  { 
    eta[] = fmax(zb[],stage_hydrograph(t));
    h[] = eta[]-zb[];
  }

  eta[bottom] = neumann(0.);
  h[bottom] = neumann(0.);
 
  boundary({zb,eta,h});

  (void) update_diffusivity(zb,eta,D,false);

  stats sx = statsf (D.x);
  stats sy = statsf (D.y);

  double dtExplicit = 0.25*L0*L0/N/N/fmax(sx.max,sy.max);
//  dt = dtnext(facExplicit*dtExplicit);
  dt = dtnext(1.7);

  fprintf(stdout,"t = %g, Dmax : %g, dt = %g\n",t,fmax(sx.max,sy.max),dt);
  fflush(stdout);

  // recompute sink term at each time step (but not diffusivity) 
  (void) update_sink(dt);
  
  mgd = diffusion(eta,dt,D,r,beta);

  (void) adapt();

  foreach()
  {
    excess[] = (eta[]<zb[]) ? zb[]-eta[] : 0.; // extra water needed to maintain eta >= zb
    eta[] = fmax(eta[],zb[]);
  }
  
  boundary({zb,eta,excess});
}

// Adaptivity

int adapt() {

#if TREE
  scalar etamask[], slope[];

  foreach()
  {
    h[] = eta[]-zb[] > dry ? eta[]-zb[] : 0.;
    etamask[] = h[] > dry ? eta[] : 0.;
    slope[] = h[] > dry ? sqrt( sq(eta[1,0]-eta[-1,0]-2.*tilt.x*Delta) + sq(eta[0,1]-eta[0,-1]-2.*tilt.y*Delta) )/2./Delta : 0.; 
  }
  boundary ({h,etamask,slope});


  astats s = adapt_wavelet ({h,etamask,slope}, (double[]){HE,ETAE,SE},
			    LEVEL, MINLEVEL);
  boundary({h});

  scalar zmin = zb.dmin;
  
  foreach(){
    double zzmin = zmin[]<nodata ? zmin[] : zb[];
    nmanning[] = zzmin<ETABF ? n_ch : n_fp;
    l[] = level;
    bool wet_stencil = false;
    foreach_neighbor(1)
      wet_stencil = wet_stencil | (h[]>dry);
    // If 3x3 neighborhood is water-free according to newly adapted h[], just ignore adapted eta[], set it to newly adapted zb[].
    // It will disable sink term.
    eta[] = wet_stencil ? eta[] : zb[];
  }

  boundary(all);

//  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
#else // Cartesian
  return 0;
#endif
}

event snapshot (t+=300. ; t<=86400.)
{
  scalar m[], etam[];
  foreach() {
    m[] = eta[]*(h[] > dry) - zb[];
    etam[] = h[] < dry ? 0. : (zb[] > 0. ? eta[] - zb[] : eta[]);
  }
  output_ppm (h, mask = m, min = 0, max = 10, 
	      n = 2048, linear = true, file = "h.mp4");
}

event info_conservation (t+=30. ; t<=86400.)
{ 
    coord Section[2];
    Section[0] = (coord) { 0. , Y0+L0-L0/N };
    Section[1] = (coord) { 600. , Y0+L0-L0/N };

    // "Measure" mean elevation at top (outlet) section

    double etam = segment_eta (Section,h,eta);

    // "Measure" discharge at top (outlet) section by looping on faces

    double Q = 0.;

    foreach_boundary(top)
       Q += abs(x-300.)<300 ? Delta * D.y[] * ((eta[0,-1]-eta[])/Delta + tilt.y) : 0.; 

    // Total water volume in domain
    
    double Vol = 0.;
    foreach()
      Vol += h[] * Delta * Delta ;

    if(pid()==0){
      fprintf(fpcontrol,"%g %g %g %g\n",t/3600.,Q,etam,Vol);
      fflush(fpcontrol);
    }
}

/**
## Results 

![Simulated water depth (flow direction is from bottom to top). Left: full Saint-Venant solution; Right: diffusive wave approximation.](hstack/h_full_diffwave.mp4)(width=80% )

~~~gnuplot Comparison with full Saint-Venant: domain volume
   reset
   set xlabel "t (h)"
   set ylabel "Total volume in domain (m3)"
   ydim = 800 
   xdim = 480
   plot '../fullsaintvenant/outlet_info.txt' u 1:4 w lines lc rgb "black" lw 3 title 'full Saint-Venant', \
        'outlet_info.txt' u 1:4 w lines lc rgb "black" lw 1 title 'diffusive wave'
~~~

~~~gnuplot Comparison with full Saint-Venant: outflow discharge
   reset
   set xlabel "t (h)"
   set ylabel "Outflow discharge (m3/s)"
   ydim = 800 
   xdim = 480
   plot '../fullsaintvenant/outlet_info.txt' u 1:2 w lines lc rgb "black" lw 3 title 'full Saint-Venant', \
        'outlet_info.txt' u 1:2 w lines lc rgb "black" lw 1 title 'diffusive wave'
~~~

![outlet_info.txt](diffwave/outlet_info.txt)
*/
