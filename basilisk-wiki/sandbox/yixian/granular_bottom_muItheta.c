/**
# Flow in a Sandglass / Hourglass / Silo with a bottom orifice
 We propose an implementation of the $\mu(I,\Theta)$ granular temperature nonlocal rheology proposed by Kim and Kamrin for the flow in a silo with a bottom orifice.
 
# Code
 Includes and definitions
 */
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "poisson.h"
#include "diffusion.h"   
#define LDOMAIN 1.

double  H0,R0,D,W,tmax;
double zetaA;      // nonlocal amplitude A
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

/* Nonlocal fields and params*/
scalar eta[];
scalar zeta_loc[];
face vector alphaf[];

scalar muIloc[];    /* local mu(I) */
scalar muI[];
scalar Inn[];
scalar thetaloc[];      /* local theta temperature */
scalar thetanonloc[];   /* nonlocal theta temperature (solution) */
double mu_s     = 0.4;
double mu_2     = 0.68;      //
double I0       = 0.4;       //
double t0_g     = 0.001;      //
double rhos     = 1.;        // density of granular
double t_switch = 0.1;
double B_new = 1;    //
double A_new = 0.15;    //

/* Boundary conditions */
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


// ==== boundary condition for theta  ====
thetanonloc[left]  = neumann(0);
thetanonloc[right] = neumann(0);
thetanonloc[top] = neumann(0);
thetanonloc[bottom] = neumann(0);

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
    tmax = 6.;
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
    foreach(){
        p[] = (fabs(x-LDOMAIN/2)<= W/2 && fabs(y)<= .1) ?  0 : max(H0 - y,0) ;
    }
}

#define rho(f) ((f) + RHOF*(1. - (f)))

event properties (i++) {
    trash ({alphav});
    foreach()
        eta[] = mug;
    
    // ===== 1. Compute local fields：μ and theta_loc =====
    foreach () {
      thetaloc[]     = 0.;
      muIloc[]   = 0.;
      zeta_loc[] = 0.;
        
        if (p[] > 0. && f[] > 0.) {
            double D2 = 0.;
            foreach_dimension() {
                double dxx = u.x[1,0] - u.x[-1,0];
                double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]) / 2.;
                D2 += sq(dxx) + sq(dxy);
            }
            
            if (D2 > 0.) {
                D2 = sqrt(D2)/(2.*Delta);
                double sD2 = sqrt(2.)*D2;        // shear rate
                double In  = sD2*D/sqrt(p[]);    // I
                double thetalocvalue = A_new/B_new* pow(In, 1.5);
                double mu_loc = mu_s + (mu_2 - mu_s)*In/(I0 + In);
                muIloc[] = mu_loc;
                Inn[] = In;
                thetaloc[]   = thetalocvalue;
                double xi  = zetaA * D;
                zeta_loc[] = xi;
                
            }
        }
        else {
            thetaloc[]     = 0.;
            muIloc[]   = 0.;
            muI[] = 0.;
            zeta_loc[] = 0.;
        }
    }
    boundary ({thetaloc, muIloc, zeta_loc});

    // ===== 2. t < t_switch：Local rheology =====
    if (t < t_switch) {
        foreach() {
          if (p[] > 0. && f[] > 0.) {
            double D2 = 0.;
            foreach_dimension() {
              double dxx = u.x[1,0] - u.x[-1,0];
              double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]) / 2.;
              D2 += sq(dxx) + sq(dxy);
            }

            if (D2 > 0.) {
              D2 = sqrt(D2)/(2.*Delta);
              double sD2 = sqrt(2.)*D2;
              double mu_loc = muIloc[];

              double etamin = sqrt(D*D*D);
              double eta_loc = (mu_loc * p[]) / sD2;
              eta_loc = max(eta_loc, etamin);
              eta_loc = min(eta_loc, 100.);
              eta[]   = eta_loc;

              thetanonloc[] = thetaloc[];
              muI[]     = mu_loc;
            } else {
              eta[]     = 100.;
              thetanonloc[] = 0.;
              muI[]     = 0.;
            }
          }
          else {
            eta[]     = mug;
            thetanonloc[] = 0.;
            muI[]     = 0.;
          }
        }
        boundary ({eta, thetanonloc, muI});
      }
    // ============ 3. t >= t_switch：Extended mu(I) granular temperature rheology============
    if (t >= t_switch) {
        foreach_face() {
            double zL = zeta_loc[];
            double zR = zeta_loc[-1,0];
            double zf = 0.5*(zL + zR);
            alphaf.x[] = sq(zf);   //  alphaf = A^2 d^2
        }
        boundary ((scalar *){alphaf});

        
        
        // 3. beta = - B
        scalar beta[];
        foreach() {
            beta[] = -B_new;
        }
        boundary({beta});
        // 4. theta = t0
        scalar theta[];
        foreach() {
            theta[] = t0_g;
        }
        boundary({theta});
        
        //rhs = A * I^1.5
        scalar rhs[];
        foreach() {
            rhs[] = 0.;
            thetanonloc[] = 0.;
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
                    double p_eff = p[];
                    double I_val = sD2 * D / sqrt(p_eff);
                    rhs[] = A_new * pow(I_val, 1.5);
                }
            }
        }
        boundary({rhs});
        
        foreach()
          thetanonloc[] = thetaloc[];
        boundary ({thetanonloc});
        
        // 3.1 solve： t0 \dot \theta = \nabla(alphaf \nabla \theta) + \beta \theta + rhs
        // alphaf= A^2 d^2, \beta = -B, rhs = A*I^1.5
        mgstats mg_g = diffusion (thetanonloc, dt, alphaf, rhs, beta, theta);

        if (mg_g.i > 0 && mg_g.resa > 0.)
            fprintf (stderr,
                     "#   thetanonloc diffusion: i=%d resa=%g resb=%g\n",
                     mg_g.i, mg_g.resa, mg_g.resb);

    
        boundary ({thetanonloc});
            
        
        // --- 3.2 solve g^{n+1} ---
        foreach() {
            if (p[] > 0. && f[] > 0.) {
                
                double D2 = 0.;
                foreach_dimension() {
                    double dxx = u.x[1,0] - u.x[-1,0];
                    double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]) / 2.;
                    D2 += sq(dxx) + sq(dxy);
                }
                
                if (D2 > 0.) {
                    D2 = sqrt(D2)/(2.*Delta);
                    double sD2 = sqrt(2.)*D2;
                    double mu_loc = muIloc[];
                    double thetalocvalue = thetaloc[];
                    double thetanonlocvalue = thetanonloc[];
                    double mu_eff =  mu_loc*pow(thetalocvalue/(thetanonlocvalue), 1./8.);//p=1/8
                    muI[] = mu_eff;
                    double etamin = sqrt(D*D*D);
                    double eta_nonloc = (mu_eff * p[]) / sD2;
                    eta_nonloc = max(eta_nonloc, etamin);
                    eta_nonloc = min(eta_nonloc, 100.);
                    eta[] = eta_nonloc;
                    
                } else {
                    eta[]     = 100.;
                    thetanonloc[] = 0.;
                    muI[]     = 0.;
                }
            }
            else {
                eta[]     = mug;
                thetanonloc[] = 0.;
                muI[]     = 0.;
            }
        }
        boundary ({thetanonloc, eta, muI});
    }
    
    scalar fa[];
    foreach()
        fa[] = (4.*f[] +
                2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
                f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1]) / 16.;
    boundary ({fa});

    foreach_face() {
        double fm = (fa[] + fa[-1,0]) / 2.;
        muv.x[]   = (fm * (eta[] + eta[-1,0]) / 2. + (1. - fm) * mug);
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


/* output including thetanonloc for diagnostics */
event interface (t = 2 ; t+=1. ; t<= 4.) {
#if dimension == 2
  output_facets (f, fpf);
#endif
  char s[80];
  sprintf (s, "field-%g.txt", t);
  FILE * fp = fopen (s, "w");
  output_field ({f,p,u,uf,pf,muIloc,muI,thetaloc,thetanonloc,eta,duxx,duxy,duyx,duyy,dpdx,dpdy}, fp, linear = true);
  fclose (fp);
}

/* debit: use explicit dx */
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
 qcc  -g -O2 -o granular_bottom_muItheta granular_bottom_muItheta.c -lm
 ./granular_bottom_muItheta > vol
 ~~~
 
# Plot

 Plot of vertical velocity at time 2
 
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
 splot 'field-2.txt' u 1:2:(abs($14)**$3) with pm3d notitle
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
 
 * Kim, S. and Kamrin, K. (2020). Power-law scaling in granular rheology across flow geometries. Phys. Rev.
 Lett., 125:088002.
 
 * Irmer, M. G., Brodsky, E. E., and Clark, A. H. (2025). Granular temperature controls local rheology of vibrated granular flows. Phys. Rev. Lett., 134:048202.
 
 * Yuan, Z., Zhao, H., and Wang, D. (2026). Toward a unified continuum framework for dense granular flows:
 μ(I) rheology extended with granular temperature. J. Fluid Mech., 1030.
 
 
 */

