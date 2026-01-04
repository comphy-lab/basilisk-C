/**
<center>
Benchmark testcase "Single Bubble Rising" originates from [Rudman (1998)](https://onlinelibrary.wiley.com/doi/abs/10.1002/%28SICI%291097-0363%2819980815%2928%3A2%3C357%3A%3AAID-FLD750%3E3.0.CO%3B2-D)
</center>


# Benchmark Test Condition

* Dimensionless Numbers
* Geometric Environment
* Boundary Condition

* Fixed Timestep
* Fixed Scale (Resolution)

<div class="figure"> <video controls="" preload="metadata"
width="80%"> <source
src="https://antoonvanhooft.nl/media/rings.mp4"
type="video/mp4"> Your browser does not support the video
tag. </video> <p class="caption"> 400000 dye particles (Video via
antoonvanhooft.nl) </p> </div>

## Set-up
*/

#define SET_EMBED
#define SET_TREE
// #define SET_AMR
#define SET_TIMESTEP
#define SET_USER_MU

#include "grid/quadtree.h"

#ifdef SET_EMBED
#include "embed.h"
#endif
#ifdef SET_USER_MU
#define mu(f) (1. / (clamp(f, 0, 1) * (1. / mu1 - 1. / mu2) + 1. / mu2))
#endif
#ifdef SET_AMR
#define MAX_SCALE 10
#define MIN_SCALE 5
#endif

#include "tag.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "maxruntime.h"
#include <time.h>


// time range config
#define MAXTIME 12

// residual goal
#define ERR_U
#define ERR_DIV

#ifndef EPS_GAS
#define EPS_GAS 1e-12 
#endif

// dimensionless numb config
const double RHOR = 1000.;
const double MUR = 1;
const double Re = 500.;
const double We = 25.;
// geometric constraints
const double H = 9;
const double R = 0.5;
const double W = 3.;
// bubble & free surface
const double Hf = 6;
const double xc = 0.;
const double yc = 2.;

#define RUD_64X192
// #define RUD_128X384
// #define RUD_192X576
// #define RUD_256X768

// X_SCLALE = (NSCALE * W) / DOMAIN_SIZE
// Y_SCALE = (NSCALE * H) / DOMAIN_SIZE
#ifdef RUD_64X192
const double DOMAIN_SIZE = 12;
const int NSCALE = (1 << 8);
const int X_SCALE = 64;
const int Y_SCALE = 192;
const double MAX_TIMESTEP = 1e-4;
#endif
#ifdef RUD_128X384
const double DOMAIN_SIZE = 12;
const int NSCALE = (1 << 9);
const int X_SCALE = 128;
const int Y_SCALE = 384;
const double MAX_TIMESTEP = 1e-4;
#endif
#ifdef RUD_192X576
const double DOMAIN_SIZE = 16;
const int NSCALE = (1 << 10);
const int X_SCALE = 192;
const int Y_SCALE = 576;
const double MAX_TIMESTEP = 1e-4;
#endif
#ifdef RUD_256X768
const double DOMAIN_SIZE = 12;
const int NSCALE = (1 << 10);
const int X_SCALE = 256;
const int Y_SCALE = 768;
const double MAX_TIMESTEP = 1e-4;
#endif

// output configs
const int PX = 5 * 1024;
const int PY = 5 * 1024;
const double FRAME_PER_SEC = 0.05;
const int FRAME_PER_ITER = (int) (FRAME_PER_SEC / MAX_TIMESTEP);
const int SNAP_PER_ITER = (int) (0.5 / MAX_TIMESTEP);
const float white[3] = {1, 1, 1};

// translation
double Hb;
double Hft;
double yt;

u.n[top] = neumann(0);
u.t[top] = neumann(0);
p[top] = dirichlet(0);
pf[top] = dirichlet(0);
f[top] = neumann(0);

int main()
{
   size(DOMAIN_SIZE);
   origin(-L0 / 2, 0);
   init_grid(NSCALE);

   // two-phase.h inputs
   rho1 = 1. [0];
   rho2 = rho1 / RHOR;
   mu1 = 1 / Re;
   mu2 = mu1 / MUR;
   f.sigma = 1 / We; // Fr=1, g=1 => U=1

   Hb = DOMAIN_SIZE - H;
   Hft = Hb + Hf;
   yt = Hb + yc;

   TOLERANCE = 1e-6 [*]; // just fit
   NITERMAX = 300;
   CFL = 0.2;

   char env_name[256] = "UNKNOWN";
   const char *im = getenv("IM");
   if (im && *im)
   {
      strncpy(env_name, im, sizeof(env_name) - 1);
      env_name[sizeof(env_name) - 1] = '\0';
   }

   char sol_method[64] = "";
#ifndef SET_TREE
   sprintf(sol_method, "MULTIGRID");
#else
   sprintf(sol_method, "QUADTREE");
#endif

   char mesh_info[512];
#ifndef SET_AMR
   sprintf(mesh_info,
           "NX = %d\n"
           "NY = %d\n",
           X_SCALE, Y_SCALE);
#else
   sprintf(mesh_info,
           "FINE = %d\n"
           "COARSE = %d\n",
           MAX_SCALE, MIN_SCALE);
#endif

   time_t t_start = time(NULL);
   struct tm st = *localtime(&t_start);

   fprintf(stderr,
           "--------------------summary--------------------\n"
           "Re = %.4f\n"
           "We = %.4f\n"
           "RHOR = %.4f\n"
           "MUR = %.4f\n"
           "R = %.4f\n"
           // mesh info begins
           "SOLUTION = %s\n"
           "DOMAIN = %.4f\n"
           "%s"
           // mesh info ends
           "PX = %d\n"
           "PY = %d\n"
           "T = %.2f\n"
           "CFL = %.4f\n"
           "MAX_TIMESTEP = %.6f\n"
           "MAX_ITERATION = %d\n"
           "TOLERANCE = %.6f\n"
           "FRAME_PER_SEC = %.3f\n"
           "ENVIRONMENT = %s\n"
           "START_TIME = %04d-%02d-%02d %02d:%02d:%02d\n"
           "--------------------summary--------------------\n"
           "------------------case begins------------------\n",
           Re, We, RHOR, MUR, R, sol_method, DOMAIN_SIZE, mesh_info, PX, PY,
           (double)MAXTIME, CFL, MAX_TIMESTEP, NITERMAX, TOLERANCE, FRAME_PER_SEC, env_name,
           st.tm_year + 1900, st.tm_mon + 1, st.tm_mday, st.tm_hour, st.tm_min, st.tm_sec);

#ifndef SET_AMR
   assert(X_SCALE == (NSCALE * W) / DOMAIN_SIZE);
   assert(Y_SCALE == (NSCALE * H) / DOMAIN_SIZE);
#endif

   clock_t tic = clock();
   run();
   clock_t tok = clock();
   double elapsed = (double)(tok - tic) / CLOCKS_PER_SEC;

   time_t t_end = time(NULL);
   struct tm et = *localtime(&t_end);

   fprintf(stderr,
           "-------------------case ends-------------------\n"
           "--------------------summary--------------------\n"
           "Re = %.4f\n"
           "We = %.4f\n"
           "RHOR = %.4f\n"
           "MUR = %.4f\n"
           "R = %.4f\n"
           // mesh info begins
           "SOLUTION = %s\n"
           "DOMAIN = %.4f\n"
           "%s"
           // mesh info ends
           "PX = %d\n"
           "PY = %d\n"
           "T = %.2f\n"
           "CFL = %.4f\n"
           "MAX_TIMESTEP = %.6f\n"
           "MAX_ITERATION = %d\n"
           "TOLERANCE = %.6f\n"
           "FRAME_PER_SEC = %.3f\n"
           "ENVIRONMENT = %s\n"
           "START_TIME = %04d-%02d-%02d %02d:%02d:%02d\n"
           "END_TIME = %04d-%02d-%02d %02d:%02d:%02d\n"
           "ELAPSED_TIME = %g\n"
           "--------------------summary--------------------\n",
           Re, We, RHOR, MUR, R, sol_method, DOMAIN_SIZE, mesh_info, PX, PY,
           (double)MAXTIME, CFL, MAX_TIMESTEP, NITERMAX, TOLERANCE, FRAME_PER_SEC, env_name,
           st.tm_year + 1900, st.tm_mon + 1, st.tm_mday, st.tm_hour, st.tm_min, st.tm_sec,
           et.tm_year + 1900, et.tm_mon + 1, et.tm_mday, et.tm_hour, et.tm_min, et.tm_sec, elapsed);
}

event init(i = 0)
{
   // return if already defined
   if (restore(file = "restart"))
   {
      return 0;
   }
   system("mkdir -p dumps");
   system("mkdir -p outputs");

#ifndef SET_TREE
   // wall geometry
   scalar phi_w[];
   foreach ()
   {
      double width = W / 2. - fabs(x);
      double height = y - Hb;
      phi_w[] = fmin(width, height);
   }
   boundary({phi_w});
   fractions(phi_w, cs, fs);
#else
   mask(x > +W / 2 ? right : none);
   mask(x < -W / 2 ? left : none);
   mask(y < Hb ? bottom : none);
#endif

   // bubble geometry
   vertex scalar phi_v[];
   foreach_vertex()
   {
      double p_liq = Hft - y;
      double p_bub = sq(x - xc) + sq(y - yt) - sq(R);
      phi_v[] = fmin(p_liq, p_bub);
   }
   boundary({phi_v});
   fractions(phi_v, f);
   boundary({f});

#ifdef SET_TIMESTEP
   DT = MAX_TIMESTEP;
   dtnext(DT);
#endif
}

event check_init(i = 0)
{
   clear();
   view(fov = 0, ty = -0.5, width = PX, height = PY);
   squares("f", min = 0, max = 1);
   draw_vof("f");
   draw_vof("cs");
   cells((coord){0, 0, 1}, 0., white, 1.0f);
   save("initial_state.png");
}

event acceleration(i++)
{
   face vector av = a;
   foreach_face(y)
       av.y[] -= 1.;
}

#ifdef SET_AMR
event adapt(i++)
{
   adapt_wavelet({f, u}, (double[]){1e-3, 1e-3, 1e-3}, MAX_LEVEL, MIN_LEVEL);
}
#endif

static inline double local_ymax_on_interface(coord n, double alpha)
{
   coord cand[4];
   int k = 0;
   

   if (fabs(n.x) > EPS_GAS)
   {
      double x = (alpha - n.y * 0.5) / n.x; // y=+1/2
      if (x >= -0.5 && x <= 0.5)
         cand[k++] = (coord){x, 0.5};
      x = (alpha + n.y * 0.5) / n.x; // y=-1/2
      if (x >= -0.5 && x <= 0.5)
         cand[k++] = (coord){x, -0.5};
   }
   if (fabs(n.y) > EPS_GAS)
   {
      double y = (alpha - n.x * 0.5) / n.y; // x=+1/2
      if (y >= -0.5 && y <= 0.5)
         cand[k++] = (coord){0.5, y};
      y = (alpha + n.x * 0.5) / n.y; // x=-1/2
      if (y >= -0.5 && y <= 0.5)
         cand[k++] = (coord){-0.5, y};
   }

   if (k == 0)
      return -HUGE;
   double ymax = -HUGE;
   for (int i = 0; i < k; i++)
      ymax = fmax(ymax, cand[i].y);
   return ymax;
}

event logfile(i += 10)
{
   scalar gasmask[];
   foreach ()
      gasmask[] = (1 - f[] > EPS_GAS);

   int ncomp = tag(gasmask);
   static int atmos_id = -1;
   foreach_boundary(top)
   {
      if (gasmask[] > 0)
         atmos_id = (int)gasmask[];
   }

   double xb = 0., yb = 0., sb = 0.;
   double vbx = 0., vby = 0.;
   double max_h = -HUGE;
   foreach (
       reduction(+ : xb) reduction(+ : yb)
           reduction(+ : vbx) reduction(+ : vby)
               reduction(+ : sb) reduction(max : max_h))
   {
      int id = (int)gasmask[];
      if (id <= 0 || id == atmos_id)
         continue;
      double w = (1. - f[]);
      w *= cs[];
      w *= dv();
      xb += x * w;
      yb += y * w;
      vbx += u.x[] * w;
      vby += u.y[] * w;
      sb += w;

      if (f[] > 0. && f[] < 1.)
      {                                        
         coord n = interface_normal(point, f); 
         double nn = sqrt(sq(n.x) + sq(n.y));  
         if (nn > 0.)
         {
            n.x /= nn;
            n.y /= nn;
         }
         double alpha = plane_alpha(f[], n); 

         double yloc = local_ymax_on_interface(n, alpha); 
         if (yloc != -HUGE)
         {
            double yabs = y + yloc * Delta;
            max_h = fmax(max_h, yabs - Hb);
         }
      }
   }
   static double sb0 = -1;
   if (sb0 < 0)
      sb0 = sb;
   double ratio = sb0 > 0. ? sb / sb0 : 0.;
   double rel_err_gas = (sb0 > 0) ? fabs(sb - sb0) / sb0 : 0.0;

   fprintf(stderr,
           "%.8f %g %.8f %.8f %.8f %.8f %.8f ratio=%.8f err=%g y_max=%.8f\n",
           t, dt, sb, xb / sb, yb / sb, vbx / sb, vby / sb, ratio, rel_err_gas, max_h);
   fflush(stderr);
}

event snapshot(i=0; i += SNAP_PER_ITER)
{
   scalar l2[], omegay[];
   char name[256];
   sprintf(name, "dumps/dump-%03d", (int)t);
   dump(file = name);

   // set latest dump as restart point
   unlink("restart");
   symlink(name, "restart");
}

event movie(t = 0; t<=MAXTIME; i += FRAME_PER_ITER)
{
   static int frame = 0;
   scalar kappa[];
   curvature(f, kappa);

   clear();
   view(fov = 0, ty = -0.5, width = PX, height = PY);
   squares("f", min = 0, max = 1);

   draw_vof("f");
   draw_vof("cs");
   // vectors("u", scale = 0.005, lc = {0, 0, 0});

   box();
   char filepath[256];
   sprintf(filepath, "outputs/frame-%04d.png", frame++);
   save(filepath);
}