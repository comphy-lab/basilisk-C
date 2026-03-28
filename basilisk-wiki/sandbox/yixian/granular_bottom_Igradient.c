/**
# Flow in a Sandglass / Hourglass / Silo with a bottom orifice
 We propose an implementation of the I-gradient nonlocal rheology proposed by Bouzid et al. for the flow in a silo with a bottom orifice.
 
# Code
 Includes and definitions
 */
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
// Domain extent
#define LDOMAIN 1.
// heap definition
double  H0,R0,D,W,tmax,zetaA;
/**
 passive fluid small density to preserve 0 pressure
 and small viscocity
 */
#define RHOF 1e-4
#define mug  1e-5
// Maximum refinement
#define LEVEL 6
char s[80];
FILE * fpf,*fwq;
scalar f[];
scalar * interfaces = {f};
face vector alphav[];
face vector muv[];
scalar rhov[];
scalar duxx[],duxy[],duyx[],duyy[],dpdx[],dpdy[];
scalar muI[],muInonloc[];
scalar Ifield[],I_new[];
scalar eta[];
double rhos     = 1.;  // density of particle
/**
 Boundary conditions for granular flow, pressure must be zero at the surface.
 The pressure is zero in the hole
 No slip boundary conditions on the other walls.
 */

p[top] = dirichlet(0);
pf[top] = dirichlet(0);
u.n[top] = neumann(0);
u.t[bottom] =  fabs(x-LDOMAIN/2)<= W/2 ? neumann(0):  dirichlet(0);
u.n[bottom] =  fabs(x-LDOMAIN/2)<= W/2 ? neumann(0):  dirichlet(0);
p[bottom]   =  fabs(x-LDOMAIN/2)<= W/2 ? dirichlet(0): neumann(0);
pf[bottom]   =  fabs(x-LDOMAIN/2)<= W/2 ? dirichlet(0): neumann(0);
u.n[right] = dirichlet(0);
u.n[left] = dirichlet(0);
u.t[right] = dirichlet(0);
u.t[left] = dirichlet(0);


Ifield[left]   = dirichlet(0);
Ifield[right]  = dirichlet(0);
Ifield[top]    = neumann(0);
Ifield[bottom] = fabs(x-LDOMAIN/2)<= W/2 ? neumann(0):  dirichlet(0);
/**
 the three cases in the main */

int main(int argc, char ** argv) {
    
    disable_fpe (FE_DIVBYZERO|FE_INVALID);
    L0 = LDOMAIN;
    // number of grid points
    N = 1 << LEVEL;
    TOLERANCE = 1e-3;
    H0=0.9;
    R0=20.000;
    
    const face vector g[] = {0.,-1.,0};
    a = g;
    alpha = alphav;
    mu = muv;
    rho=rhov;
    // maximum timestep
    DT = 0.001;
    tmax = 6.;
    W = 0.25;
    zetaA = 0.5;
    D = 0.04375;
    fpf = fopen ("interface.txt", "w");
    run();
    fclose (fpf);
    fprintf (stdout,"\n");
    
}
/**
 initial heap, a rectangle
 */
event init (t = 0) {
    scalar phi[];
    foreach_vertex()
    phi[] = min(H0 - y, R0 - x);
    fractions (phi, f);
    /**
     lithostatic pressure, with a zero pressure near the hole
     to help
     */
    foreach() {
        p[] = (fabs(x-LDOMAIN/2)<= W/2 && fabs(y)<= .1) ?  0 : max(H0 - y,0) ;
        muI[] = 0.;
    }
    boundary ({muI});
}
/**
 total density
 */
#define rho(f) ((f) + RHOF*(1. - (f)))
/**
 Viscosity computing $D_2=D_{ij}D_{ji}$;
 
 In the pure shear flow
 $D_{11}=D_{22}=0$ et $D_{12}=D_{21}=(1/2)(\partial u/\partial y)$,
 so that
 $D_2=\sqrt{D_{ij}D_{ij}} =\sqrt{ 2 D_{12}^2} = \frac{\partial u}{ \sqrt{2}\partial y}$.
 In a pure shear flow, $\partial u/\partial y= \sqrt{2} D_2$.
 
 The inertial number $I$ is $D \sqrt{2} D_2/\sqrt(p)$
 and $\mu = \mu_s+ \frac{\Delta \mu}{1+I/I_0}$
 
 $\mu_{nl} = \mu(I) (1-\nu d^2 )\frac{\nabla^2 I}{I}$
 
 $\eta_{l} = \frac{\mu_{s} p}{\sqrt{|\dot{\gamma}|^{2} + \lambda_{r}^{2}}} +\frac{(\mu_{2} - \mu_{s})p}{(I_{0} / d)\sqrt{p / \rho_{s}} + |\dot{\gamma}|}$
 
 $\eta_{nl} = \eta_l \left( 1- \frac{\nu d \sqrt{p/\rho_s}}{\sqrt{|\dot{\gamma}|+\lambda_r^2}}\right)$
 
 here $\nu = (zetaA)^2$ where zetaA is the nonlocal constant A in Kamrin et al.
 
 
 note that if $\eta$ is too small an artificial small viscosity $\rho D \sqrt{gD}$
 is taken see Lagrée et al. 11 § 2.3.1
 
 
 */
event properties (i++) {
    trash ({alphav});
    double mus = 0.4,mu2 = 0.68,I0 = 0.4;
    double lambda = 0.01;
    foreach() {
        eta[] = mug;
        Ifield[]    = 0.;
        
        if (p[] > 0. && f[] > 0.) {
            double D2 = 0.;
            foreach_dimension() {
                double dxx = u.x[1,0] - u.x[-1,0];
                double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
                D2 += sq(dxx) + sq(dxy);
            }
            
            if (D2 > 0.) {
                D2 = sqrt(D2)/(2.*Delta);  // this is D2
                double sD2 = sqrt(2.)*D2; // this sD2 is (sqrt(2) D2)
                double In = sD2*D/sqrt(p[]);
                Ifield[] = In;
            }
        }
    }
    boundary ({Ifield});
    
    
    foreach() {
        if (p[] > 0. && f[] > 0.) {
            I_new[] = (Ifield[1,0] + Ifield[-1,0] + Ifield[0,1] + Ifield[0,-1]- 4.*Ifield[])/sq(Delta);
            double D2 = 0.;
            foreach_dimension() {
                double dxx = u.x[1,0] - u.x[-1,0];
                double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
                D2 += sq(dxx) + sq(dxy);
                
            }
            if (D2 > 0.) {
                D2 = sqrt(D2)/(2.*Delta);  // this is D2
                double sD2 = sqrt(2.)*D2; // this sD2 is (sqrt(2) D2)
                double eta1 = mus*p[]/sqrt(sq(sD2)+sq(lambda)) + (mu2-mus)*p[]/((I0/D)*sqrt(p[]/rhos)+sD2);
                double eta2 = eta1*(1-sq(zetaA)*D*sqrt(p[]/rhos)*I_new[]/sqrt(sq(sD2)+sq(lambda)));
                double etamin = sqrt(D*D*D);
                eta[] = max(eta2, etamin);
                eta[] = min(eta[],100);
            }
            else {
                eta[]     = 100.;
            }
        }
        else {
            eta[]     = mug;
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
        muv.x[] = (fm*(eta[] + eta[-1,0])/2. + (1. - fm)*mug);
        alphav.x[] = 1./rho(fm);
    }
    foreach()
    rhov[] = rho(fa[]);
    boundary ({muv,alphav,rhov});
    foreach(){
        duxy[] = (u.x[0,1] - u.x[0,-1])/2./Delta + (u.y[1,0] - u.y[-1,0])/2./Delta ;
        duyx[] = (u.y[1,0] - u.y[-1,0])/2./Delta + (u.x[0,1] - u.x[0,-1])/2./Delta ;
        duxx[] = 2 * (u.x[1,0] - u.x[-1,0])/2./Delta ;
        duyy[] = 2 * (u.y[0,1] - u.y[0,-1])/2./Delta;
        dpdx[] = (p[1,0] - p[-1,0])/2./Delta;
        dpdy[] = (p[0,1] - p[0,-1])/2./Delta;
    }
}
/**
 convergence outputs
 */
void mg_print (mgstats mg)
{
    if (mg.i > 0 && mg.resa > 0.)
    fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
             exp (log (mg.resb/mg.resa)/mg.i));
}
/**
 convergence stats
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
 save some interfaces
 */
event interface (t = 2. ; t+=1. ; t<= 4.) {
#if dimension == 2
    output_facets (f, fpf);
#endif
    char s[80];
    sprintf (s, "field-%g.txt", t);
    FILE * fp = fopen (s, "w");
    output_field ({f,p,u,uf,pf,muI,muInonloc,Ifield,I_new,eta,duxx,duxy,duyx,duyy,dpdx,dpdy}, fp, linear = true);
    fclose (fp);
}



/**
 Rate of flowing materials across the hole
 */
event debit (t += 0.05; t <= tmax ) {
    static double Vold,V=1,Qinst=0;
    Vold=V;
    V=0;
    double dx = L0 / N;
    foreach()
    V = V + f[] * dx * dx;
    Qinst = -(V-Vold)/0.05;
    if(t>=.1) fprintf (stdout,"%lf %lf %lf \n",t,V,Qinst);
    fflush (stdout);
}



/**
# to run
 
~~~bash
 qcc  -g -O2 -o granular_bottom_Igradient granular_bottom_Igradient.c -lm
 ./granular_bottom_Igradient > vol
~~~
 
Plot of vertical velocity at time 2
 
~~~gnuplot 
 reset
 set pm3d map
 set pm3d; set palette rgbformulae 22,13,-31;
 splot 'field-2.txt' u 1:2:(abs($6)*$3) with pm3d notitle
 set label "Vy" at graph 0.45, 0.85 font ",24" front
# set terminal postscript color enhanced
# set output "t2vy.eps"
# replot
~~~  
Plot of pressure at time 2
 
~~~gnuplot 
 
 reset
 set pm3d map
 set pm3d; set palette rgbformulae 22,13,-31;
 splot 'field-2.txt' u 1:2:(abs($4)**$3) with pm3d notitle
 set label "pressure" at graph 0.45, 0.85 font ",24" front
 set terminal postscript color enhanced
# set output "t2p.eps"
# replot
 ~~~
 
 Plot of viscosity at time 2
 
~~~gnuplot 

 reset
 set cbrange [0:2]
 set pm3d map
 set pm3d; set palette rgbformulae 22,13,-31;
 set label "Viscosity" at graph 0.45, 0.85 font ",24" front
 splot 'field-2.txt' u 1:2:(abs($14)**$3) with pm3d notitle
 set terminal postscript color enhanced
 set output "t2eta.eps"
 replot
 ~~~
 

# Bibliography
 
 * Y. Zhou PhD
 
 * L. Staron, P.-Y. Lagrée, & S. Popinet (2014)
 "Continuum simulation of the discharge of the granular silo, A validation test for the μ(I) visco-plastic flow law" 
 Eur. Phys. J. E (2014) 37: 5 DOI 10.1140/epje/i2014-14005-6
 
 * L. Staron, P.-Y. Lagrée & S. Popinet (2012)
 "The granular silo as a continuum plastic flow: the hour-glass vs the clepsydra" 
 Phys. Fluids 24, 103301 (2012); doi: 10.1063/1.4757390  

 * Bouzid, M., Trulsson, M., Claudin, P., Cl ́ement, E. & Andreotti, B. (2013), ‘Nonlocal rheology of granular flows across yield conditions’, Phys. Rev. Lett. 111, 238301. URL: https://link.aps.org/doi/10.1103/PhysRevLett.111.238301
 
 * Lin, C. & Yang, F. 2021 Continuum simulation of non-local effects in a granular silo discharge flow using a regularized μ(I) rheology model. Physics of Fluids 33 (9), 093302.
 
 * K. Kamrin and G. Koval, Phys. Rev. Lett.,Nonlocal constitutive Relation for Steady Granular Flow 2012, 108, 178301.

*/
