#include "grid/octree.h" 		// For 3D
#include "view.h"			// For bview
#include "navier-stokes/centered.h"     // Navier stokes 
// #include "tracer.h"			// Tracers
// #include "diffusion.h"			// Diffusion 

#include "embed.h"
#include "PointTriangle.h"

/** Global variables */
int minlevel, maxlevel;         	// Grid depths
double meps, eps;			// Maximum error and error in u fields

double TEND = 15;
double H = 5, W = 80, D = 2; //Height (y), width (z) and `depth`(x)
double xp = 200, zp = 200;

#define CP 1005.	// C_p for air 
#define gCONST 9.81	// Gravitational constant
#define TREF 273.	// Kelvin
#define INVERSION .3 	// Kelvin per meter
#define karman 0.4      // von Karman constant 

#define roughY0u 0.1    // roughness wind length 
#define roughY0h 0.1     // roughness heat length

#define WIND(s) (max(0.1*log(s/roughY0u),0.))   // log Wind profile 
// #define QFLX 0. 	// 0 (0.001 = 20wm2)
// #define BSURF (1.5*b[] - 0.5*b[0, 1])
// #define GFLX (-Lambda*(BSURF - bd))
// double Lambda = 0.005, bd = 0.;   // Grass coupling
// #define STRAT(s) gCONST/TREF*s*INVERSION

// scalar b[];
// scalar * tracers = {b};

double crho = 1.;	
// face vector av[]; 

// #include "physics.h"			// Physics of the simulation 
// #include "fan.h"			// Include a fan
#include "lambda2.h"

face vector muc[];
double nu;

int main() {	
    nu = 1./300.;
    minlevel = 4;
    maxlevel = 9;

    L0 = 400.;
    X0 = Y0 = Z0 = 0.;
    mu = muc;

    init_grid(1<<8);
    // a = av; 
    // meps = 10.;					// Maximum adaptivity criterion
    // DT = 0.5;					// For poisson solver 
    // TOLERANCE=10E-4;				// For poisson solver 
    // CFL = 0.8;					// CFL condition
    run();						// Start simulation 
}

event init(t=0) {
    eps = .05;
    
	  // u.n[bottom] = dirichlet(0.);
	  // u.t[bottom] = dirichlet(0.);
	  // u.n[top] = dirichlet(0.);
	  // u.t[top] = dirichlet(WIND(y));
    // u.n[embed] = dirichlet(0.);
    // u.t[embed] = dirichlet(0.);
    // periodic (left);
    // #if dimension == 3
    //     u.r[embed] = dirichlet(0.);
    //     periodic(front);
    // #endif 
        
    u.n[left]  = dirichlet(WIND(y));
    p[left]    = neumann(0.);
    pf[left]   = neumann(0.);

    u.n[right] = neumann(0.);
    p[right]   = dirichlet(0.);
    pf[right]  = dirichlet(0.);

    u.n[embed] = dirichlet(0.);
    u.t[embed] = dirichlet(0.);
    #if dimension == 3
        u.r[embed] = dirichlet(0.);
        // periodic(back);
    #endif


    vertex scalar phw[];
    D /= 2.;
    H -= D/2;
    W -= W/2;
    foreach_vertex() {
        coord cc = {x, y, z};
        if ((z - zp) > W/2) { // distance to a side edge
        coord p1 = {xp, Y0, zp + W/2.}, p2 = {xp, Y0 + H, zp + W/2.};
        coord tmp1[1]; double tmp2[1];
        phw[] = sqrt (PointSegmentDistance (&cc, &p1, &p2, tmp1, tmp2));
        }
        else if ((z - zp) < -W/2) {
        coord p1 = {xp, Y0, zp - W/2.}, p2 = {xp, Y0 + H, zp - W/2.};
        coord tmp1[1]; double tmp2[1];
        phw[] = sqrt (PointSegmentDistance (&cc, &p1, &p2, tmp1, tmp2));
        }
        else if (y - Y0 > H) {
        coord p1 = {xp, Y0 + H, zp - W/2.}, p2 = {xp, Y0 + H, zp + W/2.};
        coord tmp1[1]; double tmp2[1];
        phw[] = sqrt (PointSegmentDistance (&cc, &p1, &p2, tmp1, tmp2));
        }
        else
        phw[] = fabs(cc.x - xp);
    }
    fractions (phw, cs, fs, D);
    boundary ({cs, fs});
    foreach()
      u.x[] = cs[] ? WIND(y) : 0.;
    // while(adapt_wavelet((scalar *){cs, u},(double []){0.01,eps,eps,eps},maxlevel,minlevel).nf) {
    //   foreach() {
    //       u.x[] = cs[]? WIND(y): 0.;
    //   }
    // }
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*nu; 
  boundary ((scalar*){muc});
}

// event inflow(i++){
//     double sides = 50;
//     double relaxtime = dt/50;
//     foreach(){
// 	if((x < sides || x > L0-sides) ||
// 	   (z < sides || z > L0-sides) ||
// 	   (y > L0-2*sides )) {
// 	    double a = (x < sides) ? x : fabs(x-L0);
// 	    a = 1.; 
// 	    u.x[] = u.x[] + a*(WIND(y)-u.x[])*relaxtime;
//  	    // b[] = b[] + a*(STRAT(y) - b[])*relaxtime;
// 	    u.y[] = u.y[] - a*u.y[]*relaxtime;
// 		  u.z[] = u.z[] - a*u.z[]*relaxtime;
// 	}
//     }
// }

// event tracer_diffusion(i++){
//     scalar r[];
//     foreach() {
//         r[] = 0;
//         if (y < Delta)
//             r[] = (QFLX + GFLX)/sq(Delta); // div needed as normalization 
//     }
    
//     double flx = 0, bt = 0;
//     double fctr = CP*TREF/gCONST;
//     foreach_boundary(bottom reduction(+:flx) reduction(+:bt)) {
//         flx = flx + (QFLX + GFLX) * sq(Delta);
//          bt = bt + BSURF * sq(Delta);
//     }
//     bt = bt/sq(L0);
//     flx = flx/sq(L0);
//     fprintf(stderr, "soil=%g %g %g %d\n", t, fctr*flx, fctr*bt/CP, i);  
    
//     mgb = diffusion(b, dt, mu, r = r);
// }

event adapt(i++) {
    adapt_wavelet((scalar *){cs,u},(double []){0.01,eps,eps,eps},maxlevel,minlevel);
}

event progress(i++) {
    fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}

event mov (t += 0.5) {
#if (dimension == 2)
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  // draw_vof ("cs", "fs", filled = -1, fc = {0,0,0});
  // draw_vof ("cs", "fs");
  squares ("omega", linear = true, map = cool_warm);
  mirror ({0,-1})
  cells();
#elif (dimension == 3)
  scalar l2[], vyz[];
  foreach()
    vyz[] = ((u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] - u.z[0,-1]))/(2.*Delta);
  boundary ({vyz});
  lambda2 (u, l2);
  view (fov = 30, theta = 0.5, phi = 0.4, 
	tx = -0.1, ty = 0.1, bg = {65./256,157./256,217./256},
	width = 1080, height = 1080);
  box();
  draw_vof ("cs", "fs", fc = {0.5,0.1,0.2});
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.01);
  translate (y = -5){
    squares ("u.x", n = {0,1,0}, alpha = 5.,
	     min = -1.2, max = 2, map = cool_warm);}
  // translate (z = -L0/2.){
  //   squares ("b", n = {0,0,1}, alpha = L0/2.,
	//      min = -0.1, max = 0.6, map = cool_warm);}
  // translate (x = -L0/2.){
  //   squares ("b", n = {1,0,0}, alpha = L0/2.,
	//      min = -0.1, max = 0.6, map = cool_warm);}
#endif
  save ("lam_3D.mp4");
}

event end(t=TEND) {
}
