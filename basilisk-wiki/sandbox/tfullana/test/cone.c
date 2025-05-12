#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
// #include "contact.h"
#include "tag.h"
#include "navier-stokes/perfs.h"
#include "view.h"

double Bo = 10. [0];            // Bond number
double Oh = 10. [0];            // Ohneshorge number
#define l_cap 1./sqrt(Bo)                 // capillary length (reference length)
#define l_ref 1.

double v0 = 15.0 [0];            // dimensionless initial volume
double h0 = 0.0384149 [0];           // dimensionless initial film thickness
double lc = 4.5 [0];             // dimensionless cone length
double beta = 25.*pi/180. [0];  // cone angle
double rho_r = 1e-3 [0];        // density ratio
double mu_r = 1e-3 [0];         // viscosity ratio

#define V0 v0*cube(l_cap)/cube(l_ref)
#define LC lc*l_cap/l_ref //(lc/l_cap)

// #define _a (2.*LC*tan(beta) + LC)
// #define _b (2.*sq(LC)*tan(beta))
// #define _c (3.*v0/pi)
// #define _T (_a*_b/6. + _c/2. - cube(_a)/27.)
// #define _x (pow(_T + sqrt(sq(_T) + cube(_b/3. - sq(_a)/9.)), 3./2.) + pow(_T - sqrt(sq(_T) + cube(_b/3. - sq(_a)/9.)), 3./2.) - _a/3.)

#define _T (cbrt(cube(LC*sin(beta)) + (3.0*(V0)/pi)*sq(cos(beta))*sin(beta)))
#define _x (_T - LC*sin(beta))

#define H0 _x //(h0/l_cap)

// https://math.stackexchange.com/questions/786335/is-there-a-function-that-draws-a-triangle-with-rounded-edges
double Glame(double x, double y, int nn) {
  double sum = 0.;
  for (int ii = 1; ii <= 3; ii++){
    sum += pow(fabs(x*cos(2.*pi*ii/3.) - y*sin(2.*pi*ii/3.) - 1./3.), nn);
  }
  return sum - 1.;
}

#define L1 2.*LC

#define x1 (L1 - LC)
#define x2 (L1)
#define y1 0.
#define y2 (LC*tan(beta))
#define a1 ((y2 - y1) / (x2 - x1))
#define b1 (y1 - a1*x1)
#define DIAG (a1 * x + b1 - y)
// #define DOM intersection(-DIAG, -(y - 2.*(LC*tan(beta) + H0/cos(beta))))
#define DOM intersection(-DIAG, -(y - 1.*LC))
// #define _DOM intersection(-DIAG, -(y - LC))
// #define _circle (sq(x - LC - 0.5) + sq(y) - sq(0.5*tan(beta)))
// #define DOM intersection(_circle, _DOM)
// #define DOM Glame((x - LC), y, 6)
// #define DOM Glame((-(cos(beta)*x + sin(beta)*y) - LC), (-sin(beta)*x + cos(beta)*y), 6)

#define f_x1 (L1 - LC - 1.0*H0/cos(pi/2. - beta))
#define f_x2 (L1)
#define f_y1 0.
#define f_y2 (LC*tan(beta) + H0/cos(beta))
#define f_a1 ((f_y2 - f_y1) / (f_x2 - f_x1))
#define f_b1 (f_y1 - f_a1*f_x1)
// #define f_DIAG (f_a1 * x + f_b1 - y)
#define f_DIAG max(max(intersection((f_a1 * x + f_b1 - y), (x - x1)), -(sq(x - x1) + sq(y) - sq(H0/cos(pi/2. - beta)*tan(beta)))), ((x - L1 + 0.1*l_cap)))

#define fp_x1 (L1 - LC - 3.*H0/cos(pi/2. - beta))
#define fp_x2 (L1)
#define fp_y1 0.
#define fp_y2 (LC*tan(beta) + 3.*H0/cos(beta))
#define fp_a1 ((fp_y2 - fp_y1) / (fp_x2 - fp_x1))
#define fp_b1 (fp_y1 - fp_a1*fp_x1)
#define fp_DIAG (fp_a1 * x + fp_b1 - y)

#define fm_x1 (L1 - LC + 3.*(L1/pow(2,LEVEL))/cos(pi/2. - beta))
#define fm_x2 (L1)
#define fm_y1 0.
#define fm_y2 (LC*tan(beta) - 3.*(L1/pow(2,LEVEL))/cos(beta))
#define fm_a1 ((fm_y2 - fm_y1) / (fm_x2 - fm_x1))
#define fm_b1 (fm_y1 - fm_a1*fm_x1)
#define fm_DIAG (fm_a1 * x + fm_b1 - y)
// #define DOM_toref intersection(-fm_DIAG, -(y - 2.*(LC*tan(beta) + H0/cos(beta))))
#define DOM_toref intersection(-fm_DIAG, -(y - 1.*LC))

#define TOREF difference(-fm_DIAG, -max(-(y - 2.*H0), fp_DIAG))

int LEVEL = 7;

uf.n[bottom] = 0.;

// u.t[right] = neumann(0);
// u.n[right] = neumann(0);
u.t[right] = dirichlet(0);
u.n[right] = dirichlet(0);
// uf.n[right] = 0.;
// uf.t[right] = 0.;
f[right] = 1.;
// f[right] = neumann(0);
// vector h[];
// h.t[right] = contact_angle(60.*pi/180.);

p[left] = dirichlet(0);
u.t[left] = neumann(0);
u.n[left] = neumann(0);

#if EMBED
u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);
f[embed] = neumann(0);
#endif

char movie_name[360], snap_name[360];



int main(int argc, char *argv[])
{
  if (argc > 1)
    Bo = atof(argv[1]);
  if (argc > 2)
    Oh = atof(argv[2]);
  if (argc > 3)
    lc = atof(argv[3]);
  if (argc > 4)
    v0 = atof(argv[4]);
  if (argc > 5)
    beta = atof(argv[5])*pi/180.;
  if (argc > 6)
    LEVEL = atoi(argv[6]);


  // for(beta = 8.*pi/180.; beta <= 16.*pi/180.; beta *= 2.){
  //   for(LC = 2.; LC <= 4.; LC *= 2.){
  //     for(H0 = 0.1; H0 <= 0.2; H0 *= 2.){
  //       for(Bo = 4.; Bo <= 16.; Bo *= 2.){
  //         for(Oh = 4.; Oh <= 16.; Oh *= 2.){
            size(L1);

            origin(0., 0.);
          
            init_grid(1 << 6);
          
            DT = HUGE [0];
          
            rho1 = 1. [0], mu1 = Oh/sqrt(Bo);
            rho2 = rho_r, mu2 = Oh*mu_r/sqrt(Bo);
          
            stokes = true;
          
            f.sigma = 1./Bo;
            // f.height = h; 
            // TOLERANCE = 1e-4 [*];
            run();
  //         }
  //       }
  //     }
  //   }
  // }
}

event init(t = 0){
    refine (DOM_toref > 0 && level < LEVEL);
    unrefine (DOM_toref < 0 && level > 2);
#if EMBED
    solid(cs, fs, DOM);
    fractions_cleanup (cs, fs);
    cm_update (cm, cs, fs);
    fm_update (fm, cs, fs);
    restriction ({cm, fm, cs, fs});
#endif
    fraction(f, f_DIAG);
    sprintf (movie_name, "cone_Bo%.0f_Oh%.0f_LC%.0f_H%.1f_beta%.0f_LEVEL%d.mp4", Bo, Oh, LC, H0, beta*180./pi, LEVEL);
    char str[360];
    sprintf (str, " Bo=%.0f Oh=%.0f L=%.0f H=%.1f angle=%.0fdeg\n t/t_mu=%.1f", Bo, Oh, LC, H0, beta*180./pi, t);
    view(quat = {0.000, 0.000, -0.707, 0.707}, ty = -0.5,
      width = 1024, height = 1024);
    clear();
    draw_string (str, 1, lw = 3, lc = {0, 0, 0}, pos = 1, size = 40);
    box(notics = true, lw = 2.);
    draw_vof(c = "f", fc = {0.5, 0, 0.5}, lc = {0.5, 0, 0.5}, lw = 3);
    draw_vof(c = "cs", fc = {0, 0, 0}, lc = {0, 0, 0}, lw = 2, filled = -1);
    squares("u.x*cs", cbar = true, min = -0.5, max = 0.5, label = "u.z", pos = {-.75, -.75}, lscale = 1, size = 40, lw = 3);
    mirror ({0,1}) {
        box(notics = true, lw = 2.);
        draw_vof(c = "f", fc = {0.5, 0, 0.5}, lc = {0.5, 0, 0.5}, lw = 3);
        draw_vof(c = "cs", fc = {0, 0, 0}, lc = {0, 0, 0}, lw = 2, filled = -1);
        squares ("p*cs", cbar = true, min = -1., max = 1., label = "P", pos = {.6, -.75}, lscale = 1, size = 40, lw = 3);
    }
    save("test.png");
    double vd1 = 0.;
  
    foreach (reduction(+:vd1))
    {
      double dv1 = 2.*pi*f[]*dv();
      vd1 += dv1;
    }
    // fprintf (stdout, "%g %g %g", t, _x, _x*cos(pi/2. - beta));
    fprintf (stderr, "%g %g %g %g %g %g", t, v0, vd1, vd1/cube(l_cap), V0, _x);
}

event acceleration(i++)
{

  face vector av = a;
  foreach_face(x)
      av.x[] -= 1.; //Bo;
  remove_droplets (f, 3, 1e-3);
}

event savemovie (t+=1e-1){
    // char str[360];
    // sprintf (str, " Bo = %.1f Oh = %.1f\n LC = %.1f H0 = %.2f\n beta = %.0f deg \n t/t_mu = %.2f", Bo, Oh, LC, H0, beta*180./pi, t);
    // view (quat = {0.000, 0.000, -0.707, 0.707},
    //     fov = 35, near = 0.01, far = 1000,
    //     tx = 0.0, ty = -0.575, tz = -1.743,
    //     width = 512, height = 1024);
    // clear();
    // // squares("sqrt(u.x*u.x+u.y*u.y)*cs", cbar = true, min = 0., max = 1., label = "u mag.", pos = {-.95, -.95});
    // draw_vof(c = "f", fc = {0, 0, 0}, lc = {0, 0, 0}, lw = 2);
    // draw_string (str, 1, lw = 3, lc = {0, 0, 0}, pos = 1, size = 30);
    // draw_vof(c = "cs", fc = {0, 0, 0}, lc = {0, 0, 0}, lw = 2, filled = -1);
    // // scalar psi[], stream[], omega[];
    // // psi[bottom] = dirichlet (0);
    // // psi[top]    = dirichlet (0);   //Also a consistent BC for psi.
    // // psi[left]   = neumann (0);
    // // psi[right]  = neumann (0);
    // // // psi[embed]  = dirichlet (0);
    // // vorticity (u, omega);
    // // poisson (psi, omega);
    // // boundary ({psi});
    // // foreach_vertex()
    // //   stream[] = cs[]*(psi[0,-1] + psi[-1,-1]
    // //       + psi[] + psi[-1])/4; 
    // // isoline ("stream", n = 32);
    // vectors("u", scale = 0.1, level = 6);
    // mirror ({0,1}) {
    //     squares ("p*cs", cbar = true, min = -1., max = 1., label = "P", pos = {.65, -.9}, lscale = 0.5, size = 30);
    //     draw_vof(c = "f", fc = {0, 0, 0}, lc = {0, 0, 0}, lw = 2);
    //     draw_vof(c = "cs", fc = {0, 0, 0}, lc = {0, 0, 0}, lw = 2);
    // }
    // view (quat = {0.000, 0.000, -0.707, 0.707},
    //   fov = 35, near = 0.01, far = 1000,
    //   tx = 0.0, ty = -0.575, tz = -1.743,
    //   width = 512, height = 1024);
    // clear();
    // draw_string (str, 1, lw = 3, lc = {0, 0, 0}, pos = 1, size = 30);
    // box(notics = true, lw = 2.);
    // draw_vof(c = "f", fc = {0.5, 0, 0.5}, lc = {0.5, 0, 0.5}, lw = 3);
    // draw_vof(c = "cs", fc = {0, 0, 0}, lc = {0, 0, 0}, lw = 2, filled = -1);
    // squares("u.x*cs", cbar = true, min = -0.5, max = 0.5, label = "u.z", pos = {-.95, -.9}, lscale = 2, size = 15);
    // mirror ({0,1}) {
    //     box(notics = true, lw = 2.);
    //     draw_vof(c = "f", fc = {0.5, 0, 0.5}, lc = {0.5, 0, 0.5}, lw = 3);
    //     draw_vof(c = "cs", fc = {0, 0, 0}, lc = {0, 0, 0}, lw = 2, filled = -1);
    //     squares ("p*cs", cbar = true, min = -1., max = 1., label = "Pressure", pos = {0.35, -.9}, lscale = 2, size = 15);
    // }
    char str[360];
    sprintf (str, " Bo=%.0f Oh=%.0f L=%.0f H=%.1f angle=%.0fdeg\n t/t_mu=%.1f", Bo, Oh, LC, H0, beta*180./pi, t/sqrt(Bo));
    view(quat = {0.000, 0.000, -0.707, 0.707}, ty = -0.5,
      width = 1024, height = 1024);
    clear();
    draw_string (str, 1, lw = 3, lc = {0, 0, 0}, pos = 1, size = 40);
    box(notics = true, lw = 2.);
    draw_vof(c = "f", fc = {0.5, 0, 0.5}, lc = {0.5, 0, 0.5}, lw = 3);
    draw_vof(c = "cs", fc = {0, 0, 0}, lc = {0, 0, 0}, lw = 2, filled = -1);
    squares("u.x*cs", cbar = true, min = -0.1, max = 0.1, label = "u.z", pos = {-.75, -.75}, lscale = 1, size = 40, lw = 3);
    mirror ({0,1}) {
        box(notics = true, lw = 2.);
        draw_vof(c = "f", fc = {0.5, 0, 0.5}, lc = {0.5, 0, 0.5}, lw = 3);
        draw_vof(c = "cs", fc = {0, 0, 0}, lc = {0, 0, 0}, lw = 2, filled = -1);
        squares ("p*cs", cbar = true, min = -1., max = 1., label = "P", pos = {.6, -.75}, lscale = 1, size = 40, lw = 3);
    }
    save(movie_name);
}

event output (t = 100/sqrt(Bo)){
    sprintf (snap_name, "snap_cone_Bo%.0f_Oh%.0f_LC%.0f_H%.1f_beta%.0f_t%g_LEVEL%d", Bo, Oh, LC, H0, beta*180./pi, t, LEVEL);
    dump(snap_name);
}