/**
# Flow in a Sandglass / Hourglass / Silo with a bottom orifice
 We propose an implementation of the non-local granular fluidity (NGF) model proposed by Kamrin and Koval (2012) for the flow in a silo with a bottom orifice.
 
# Code
 Includes and definitions
 */


#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "poisson.h"
#define LDOMAIN 1.
double  H0,R0,D,W,tmax,zetaA;
#define RHOF 1e-4
#define mug  1e-5
#define LEVEL 6

char s[80];
FILE * fpf,*fwq;
scalar f[];
scalar * interfaces = {f};
face vector alphav[];
face vector muv[];
scalar rhov[];
scalar duxx[],duxy[],duyx[],duyy[],dpdx[],dpdy[];

/* Nonlocal fields and params */
scalar gloc[];      /* local fluidity */
scalar gnonloc[];   /* nonlocal fluidity */
scalar muIloc[];    /* local mu(I) */
scalar muI[];
scalar eta[];
scalar zeta_loc[];
scalar lambda[],rhs[];

face vector alphaf[];
double mus     = 0.4;      //
double mu2     = 0.68;      //
double I0       = 0.4;       //
double t_switch = 0.1;
// ==== Boundary conditions for g ====
gnonloc[top] = neumann(0);
gnonloc[left]  = dirichlet(0.);
gnonloc[right] = dirichlet(0.);
gnonloc[bottom] = fabs(x-LDOMAIN/2)<= W/2 ? neumann(0):  dirichlet(0);


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

int main(int argc, char ** argv) {

    disable_fpe (FE_DIVBYZERO|FE_INVALID);
    L0 = LDOMAIN;
    N = 1 << LEVEL;
    TOLERANCE = 0.001;
    H0 = 0.9;
    R0 = 20.000;
    const face vector g[] = {0.,-1.,0};
    a = g;
    alpha = alphav;
    mu = muv;
    rho = rhov;
    DT = 0.001;
    tmax = 4.;
    W = 0.25;
    zetaA = 0.5;
    D = 0.04375;
    fpf = fopen ("interface.txt", "w");
    run();
    fclose (fpf);
    fprintf (stdout,"\n");
    return 0;

}

/* initial heap, a rectangle */
event init (t = 0) {
    scalar phi[];
    foreach_vertex()
        phi[] = min(H0 - y, R0 - x);
    fractions (phi, f);
    foreach() {
        p[] = (fabs(x-LDOMAIN/2)<= W/2 && fabs(y)<= .1) ?  0 : max(H0 - y,0) ;
        gloc[] = 0.;
        gnonloc[] = 0.;
        muIloc[] = 0.;
        muI[]  = 0.;
    }
    boundary({gloc, gnonloc, muIloc, muI});
}

#define rho(f) ((f) + RHOF*(1. - (f)))

/**
We solve for each time step the following additional set of equations using values obtained from the previous local formulation

$$g_{\mathrm{loc}} = \frac{|\dot{\gamma} |}{ \mu_{\mathrm{loc}}(I)}$$
$$ -[\xi(\mu_\mathrm{loc})]^2 \nabla^{2}g  +g = g_{\mathrm{loc}}$$
$$\mu =\frac{ |\dot{\gamma} | }{ g}$$
$$ \eta_{\mathrm{eff}} = \frac{\mu p}{  |\dot{\gamma} |}$$

where $\xi(\mu)$ is the stress-dependent cooperativity length given by:

$xi (\mu) = A\sqrt{\frac{\mu_2 - \mu}{\Delta\mu(\mu - \mu_s)}} d$

with $A$ a dimensionless non-local amplitude that quantifies the spatial extent of cooperativity in the flow. To avoid division by zero when computing the friction coefficient from the fluidity, we add a small regularisation term $\varepsilon = 10^{-16}$ to the denominator, i.e., $\mu = |\dot{\gamma}| / (g + \varepsilon)$.

 */

event properties (i++) {
    trash ({alphav});
    trash ((scalar *){alphaf});

    foreach() {
        muIloc[]   = 0.;
        gloc[]     = 0.;
        rhs[]      = 0.;
        lambda[]   = 1.0;
        zeta_loc[] = 0.;

        if (p[] > 0. && f[] > 0.) {
            double D2 = 0.;
            foreach_dimension() {
                double dxx = u.x[1,0] - u.x[-1,0];
                double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]) / 2.;
                D2 += sq(dxx) + sq(dxy);
            }

            if (D2 > 0.) {
                D2 = sqrt(D2) / (2.0 * Delta);
                double sD2 = sqrt(2.) * D2;

                double p_eff  = p[];
                double In_loc = sD2 * D / sqrt(p_eff);
                double Delta_mu = mu2 - mus;
                double mu_loc = mus + Delta_mu * In_loc / (I0 + In_loc);
                muIloc[] = mu_loc;
                gloc[]   = sD2 / mu_loc;
                double num = mu2 - mu_loc;
                double den = Delta_mu * (mu_loc - mus);
                double xi = zetaA * D * sqrt(num / den);
 
                zeta_loc[] = xi;
            }
        }
        rhs[] = gloc[];
    }
    boundary ({muIloc, gloc, rhs, lambda, zeta_loc});


    if (t < t_switch) {

        foreach() {
            eta[] = mug;
            muI[] = 0.;

            if (p[] > 0. && f[] > 0.) {
                double D2 = 0.;
                foreach_dimension() {
                    double dxx = u.x[1,0] - u.x[-1,0];
                    double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]) / 2.;
                    D2 += sq(dxx) + sq(dxy);
                }
                if (D2 > 0.) {
                    D2 = sqrt(D2) / (2.0 * Delta);
                    double sD2 = sqrt(2.) * D2;

                    double muI_val = muIloc[];
                    muI[] = muI_val;

                    double etamin = sqrt(D*D*D);
                    eta[] = max((muI_val * p[]) / sD2, etamin);
                    eta[] = min(eta[], 100.);
                }else {
                    eta[]     = 100.;
                     muI[]     = 0.;
                }
            }
            else {
                eta[]     = mug;
                 muI[]     = 0.;
                 }
        }
        boundary ({eta, muI});

    }

    else {

        foreach_face(x) {
            double zf = 0.5*(zeta_loc[] + zeta_loc[-1,0]);
            alphaf.x[] = -sq(zf);
        }
        foreach_face(y) {
            double zg = 0.5*(zeta_loc[] + zeta_loc[0,-1]);
            alphaf.y[] = -sq(zg);
        }
        boundary ((scalar *){alphaf});

        mgstats mg = poisson (gnonloc, rhs, alphaf, lambda);

        if (mg.i > 0 && mg.resa > 0.)
            fprintf(stderr, "poisson mg: iter=%d resb=%g resa=%g\n", mg.i, mg.resb, mg.resa);

        foreach() {
            if (gnonloc[] < 0.) {
                fprintf(stderr, "Non-physical gnonloc < 0 at (x=%g,y=%g): %g\n", x, y, gnonloc[]);
                gnonloc[] = 0.;
            }
        }
        boundary ({gnonloc});

        foreach() {
            eta[] = mug;
            muI[] = 0.;

            if (p[] > 0. && f[] > 0.) {
                double D2 = 0.;
                foreach_dimension() {
                    double dxx = u.x[1,0] - u.x[-1,0];
                    double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]) / 2.;
                    D2 += sq(dxx) + sq(dxy);
                }
                if (D2 > 0.) {
                    D2 = sqrt(D2) / (2.0 * Delta);
                    double sD2  = sqrt(2.) * D2;

                    double g_eff = gnonloc[];
                    
                    double muI_val = sD2 / (g_eff + 1e-16);
                    muI[] = muI_val;

                    double etamin = sqrt(D*D*D);
                    eta[] = max((muI_val * p[]) / sD2, etamin);
                    eta[] = min(eta[], 100.);
                }
                else {
                    eta[]     = 100.;
                    gnonloc[] = 0.;
                    muI[]     = 0.;
                }
            }
            else {
                eta[]     = mug;
                gnonloc[] = 0.;
                muI[]     = 0.;
            }
        }
        boundary ({eta, gnonloc, muI});
    }

    scalar fa[];
    foreach()
        fa[] = (4.*f[] +
                2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
                f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1]) / 16.;
    boundary ({fa});

    foreach_face() {
        double fm = (fa[] + fa[-1,0]) / 2.;
        muv.x[]    = (fm * (eta[] + eta[-1,0]) / 2. + (1. - fm) * mug);
        alphav.x[] = 1.0 / rho(fm);
    }
    foreach()
        rhov[] = rho(fa[]);
    boundary ({muv, alphav, rhov});
    foreach(){
        duxy[] = (u.x[0,1] - u.x[0,-1])/2./Delta + (u.y[1,0] - u.y[-1,0])/2./Delta ;
        duyx[] = (u.y[1,0] - u.y[-1,0])/2./Delta + (u.x[0,1] - u.x[0,-1])/2./Delta ;
        duxx[] = 2 * (u.x[1,0] - u.x[-1,0])/2./Delta ;
        duyy[] = 2 * (u.y[0,1] - u.y[0,-1])/2./Delta;
        dpdx[] = (p[1,0] - p[-1,0])/2./Delta;
        dpdy[] = (p[0,1] - p[0,-1])/2./Delta;
    }
}

/* convergence outputs */
void mg_print (mgstats mg)
{
    if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 exp (log (mg.resb/mg.resa)/mg.i));
}

event logfile (i++) {
    stats s = statsf (f);
    fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
    mg_print (mgp);
    mg_print (mgpf);
    mg_print (mgu);
    fflush (stderr);
}

/* output */
event interface (t = 2. ; t+=1. ; t<= 4.) {
#if dimension == 2
  output_facets (f, fpf);
#endif
  char s[80];
  sprintf (s, "field-%g.txt", t);
  FILE * fp = fopen (s, "w");
  output_field ({f,p,u,uf,pf,muIloc,muI,gloc,gnonloc,eta,duxx,duxy,duyx,duyy,dpdx,dpdy}, fp, linear = true);
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



/* End of file */


/**
# to run
 
~~~bash
 qcc  -g -O2 -o granular_bottom_NGF granular_bottom_NGF.c -lm
 ./granular_bottom_NGF > vol
~~~
 
 
# plots

 

 

~~~gnuplot
 reset
 set pm3d map
 set pm3d; set palette rgbformulae 22,13,-31;
 splot 'field-2.txt' u 1:2:(abs($6)*$3) with pm3d notitle
 set label "Vy" at graph 0.45, 0.85 font ",24" front 
~~~

 set terminal postscript color enhanced
 set output "t2vy.eps"
 replot
 
 Plot of pressure at time 2
 
~~~gnuplot
 
 reset
 set pm3d map
 set pm3d; set palette rgbformulae 22,13,-31;
 splot 'field-2.txt' u 1:2:(abs($4)**$3) with pm3d notitle
 set label "pressure" at graph 0.45, 0.85 font ",24" front
 set terminal postscript color enhanced
~~~
 
 set output "t2p.eps"
 replot
 
 
 Plot of viscosity at time 2
 
~~~gnuplot
 
 reset
 set cbrange [0:2]
 set pm3d map
 set pm3d; set palette rgbformulae 22,13,-31;
 set label "Viscosity" at graph 0.45, 0.85 font ",24" front
 splot'field-2.txt' u 1:2:(abs($14)**$3) with pm3d notitle
~~~
 
 set terminal postscript color enhanced
 set output "t2eta.eps"
 replot
 
 
 

# Links
[https://basilisk.fr/sandbox/yixian/granular_bottom_NGF.c]()

[https://basilisk.fr/sandbox/yixian/granular_bottom_constantNGF.c]()

[https://basilisk.fr/sandbox/yixian/granular_bottom_dynamicNGF.c]()

[https://basilisk.fr/sandbox/yixian/granular_bottom_linearisedconstantNGF.c]()

[https://basilisk.fr/sandbox/yixian/granular_bottom_linearisedNGF.c]()

[https://basilisk.fr/sandbox/yixian/granular_bottom_Local.c]()

[https://basilisk.fr/sandbox/yixian/granular_bottom_muItheta.c]()

[https://basilisk.fr/sandbox/yixian/granular_bottom_Igradient.c]()
 
# Bibliography
 
 * Y. Zhou PhD [https://theses.fr/2016AIXM4731]()
 
  * L. Staron,  <a href="http://www.lmm.jussieu.fr/%7Elagree/TEXTES/PDF/epje130141.pdf">P.-Y. Lagr&eacute;e</a>,   &amp; S. Popinet  (2014)<br>
"Continuum simulation of the discharge of the granular silo,
 A validation test for the &#956;(I) visco-plastic flow law"
<br>
 Eur. Phys. J. E (2014) 37: 5 DOI 10.1140/epje/i2014-14005-6<br>
 
 * L. Staron,   <a href="http://www.lmm.jussieu.fr/%7Elagree/TEXTES/PDF/PhysFluids_24_103301.pdf">P.-Y. Lagr&eacute;e</a>  &amp; S. Popinet (2012)<br>
"The granular silo as a continuum plastic flow: the hour-glass vs the clepsydra"
<br>Phys. Fluids 24, 103301 (2012); doi: 10.1063/1.4757390 <br>
 
 
 * Kamrin, K. and Koval, G. (2012). Nonlocal constitutive relation for steady granular flow. Phys. Rev. Lett.,
 108:178301.
 
 */




