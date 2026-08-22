/**
# Flow past a cylinder in 3D

![](https://antoonvanhooft.nl/media/mov20003D.mp4)

 */

/**
This code is a modified version of [Antoon's 3D cylinder](https://basilisk.fr/sandbox/Antoonvh/tree/cyl3d.c). It includes a mask and viscous forces. I think it's a bit more stable.

It lacks of precision... Maybe it needs an extreme mesh to be precise, I don't know :(*/





#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
//#include "navier-stokes/double-projection.h"
//#include "navier-stokes/perfs.h"
/** The code is unstable, so it could be usefull for some applications to add a small shift (second line) */
#define CYLINDER (sqrt(sq(x) + sq(y)) - R)
//#define CYLINDER (sqrt(sq(x) + sq(y - 0.01*R)) - R)

double R = 0.5, U = 1, Re = 40., c = 30., ue, nu;
int maxlevel = 11;

/** Mask size (spanwise); Lz=4 equal to a cylinder of length 4*D :*/
double Lz = 12.;

u.n[left]  = dirichlet (U);
uf.n[left] = U;
u.t[left]  = dirichlet (0);
p[left]    = dirichlet (0.);
pf[left]   = dirichlet (0.);

u.n[right] = neumann (0.);
p[right]   = neumann (0.);
pf[right]  = neumann (0.);

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);
#if (dimension == 3)
u.r[embed] = dirichlet (0.);
u.r[left]  = dirichlet (0);
#endif
FILE * fp;
face vector muc[];
int main() {
  nu = U*R/Re;
  periodic (top);
#if (dimension == 3)
  periodic (back);
#endif
  L0 = 50;
  X0 = -10;
  Y0 = Z0 = -L0/2.;
  mu = muc;
  char logname[99];
  ue = U/c;
  sprintf (logname, "log3rd%g3D%d-%g_Lz%g",Re, maxlevel, c, Lz);
  fp = fopen (logname, "w");

    N = 64;
    run();

}

event init (t = 0) {
  /**
        I partially mask the spanwise (z) direction of the domain.
        I do so in order to reduce the number of computational points.
        Hoping to make my computation faster...*/
  mask (z > Lz *R ? front : z < - Lz *R ? back : none);
        /**
        Please note that, doing so, the computational cells are always there.
        The only thing is that our foreach() iterators will not explore these cells,
        possibly making the code faster.*/
      /**
      I give an approximate refinement in the zone that will get the
      embedded cylinder. I do not refine too much here, just enough
      to "see" it .*/
  refine (CYLINDER < 0.2*R && CYLINDER > -0.2*R && level < maxlevel-2);
      /**
      I cut-out the embedded boundary, at lower refinement...*/
        solid(cs,fs,CYLINDER);
      /**
      I now refine the mesh at the maximum level I'll employ during simulation.*/
  refine (CYLINDER < 0.2*R && CYLINDER > -0.2*R && level < maxlevel);
      /**
      Finally, I re-call the solid function so to be sure that
      embed has been defined on the finest possible mesh...*/
        solid(cs,fs,CYLINDER);

  foreach() {
    u.x[] = cs[] > 0;
    u.y[] = noise()*U/200.;
#if (dimension == 3)
    u.z[] = noise()*U/200.;
#endif
  }
}

event damp (i++) {
  coord Uinf = {U, 0, 0};
  foreach() {
    if (fabs(x - (X0 + L0/2.)) > 4*L0/10.)
      foreach_dimension()
        u.x[] += dt*(Uinf.x - u.x[])/2.;
  }
  boundary ((scalar*){u});
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*nu;
  boundary ((scalar*){muc});
}

void prolongate_ratio (Point point, scalar s) {
  foreach_child() {
    if (s[] != nodata)
      s[] += s[]*Delta;
  }
}

event adapt (i++) {
  scalar res[];
  foreach() {
    res[] = nodata;
    if (cs[] > 0 && cs[] < 1) {
      res[] = U/sqrt(R*nu);
    }
  }
  res.prolongation = prolongate_ratio;
  adapt_wavelet ((scalar*){res, u}, (double[]){0.01, ue, ue, ue}, maxlevel, 4);
//  unrefine (level > 4 && (x - X0) < L0/10. && (X0 + L0 - x) < L0/10.);
    unrefine (level > 4 && ((x - X0) < L0/10. || (X0 + L0 - x) < L0/10.));

}

/**
To ensure that the embedded geometry is well defined (regardless on the
refinement level), I re-run the embed definition, here every 100th iteration.*/

//event re_define_cylinder_geometry(i+=100){
//      solid(cs,fs,CYLINDER);
//
//}



/** my modified embed (that calculates pressure drag AND viscous drag):*/

double embed_interpolate_3D (Point point, scalar s, coord p) {
  int i = sign(p.x), j = sign(p.y);
#if dimension == 2
  if (cs[i] && cs[0,j] && cs[i,j])
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) +
            (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
#else // dimension == 3
  int k = sign(p.z);
  // CORRECTION : variables locales au lieu d'ecraser les macros x, y, z
  double lx = fabs(p.x), ly = fabs(p.y), lz = fabs(p.z);
  if (cs[i] && cs[0,j] && cs[i,j] && cs[0,0,k] &&
      cs[i,0,k] && cs[0,j,k] && cs[i,j,k]) {
    return (((s[]*(1. - lx) + s[i]*lx)*(1. - ly) +
             (s[0,j]*(1. - lx) + s[i,j]*lx)*ly)*(1. - lz) +
            ((s[0,0,k]*(1. - lx) + s[i,0,k]*lx)*(1. - ly) +
             (s[0,j,k]*(1. - lx) + s[i,j,k]*lx)*ly)*lz);
  }
#endif
  else {
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (cs[i])
        val += fabs(p.x)*(s[i] - s[]);
      else if (cs[-i])
        val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
}

void embed_force_3D (scalar p, vector u, face vector mu, coord * Fp, coord * Fmu)
{
  double Fp_x = 0., Fp_y = 0., Fp_z = 0.;
  double Fmu_x = 0., Fmu_y = 0., Fmu_z = 0.;
  foreach (reduction(+:Fp_x) reduction(+:Fp_y) reduction(+:Fp_z)
           reduction(+:Fmu_x) reduction(+:Fmu_y) reduction(+:Fmu_z)) {
    if (cs[] > 0. && cs[] < 1.) {

      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);

      double Fn = area*embed_interpolate_3D (point, p, b);
      foreach_dimension()
        Fp_x += Fn*n.x;

      if (constant(mu.x) != 0.) {
        double mua = 0., fa = 0.;
        foreach_dimension() {
          mua += mu.x[] + mu.x[1];
          fa  += fs.x[] + fs.x[1];
        }
        mua /= fa;

        coord dudn = embed_gradient (point, u, b, n);
#if dimension == 2
        foreach_dimension()
          Fmu_x -= area*mua*(dudn.x*(sq(n.x) + 1.) +
                             dudn.y*n.x*n.y);
#else // dimension == 3
        foreach_dimension()
          Fmu_x -= area*mua*(dudn.x*(sq(n.x) + 1.) +
                             dudn.y*n.x*n.y +
                             dudn.z*n.x*n.z);
#endif
      }
    }
  }
  foreach_dimension() {
    Fp->x  = Fp_x;
    Fmu->x = Fmu_x;
  }
}




event logger (i += 5) {
#if (dimension == 2)
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  fprintf (fp, "%d %g %g %g %g %g %ld\n",
           i, t, Fp.x, Fp.y, Fmu.x, Fmu.y, grid->n);
#elif (dimension == 3)
  /** My way to compute forces :*/
  coord Fp, Fmu;
  embed_force_3D (p, u, mu, &Fp, &Fmu);
  double Fx = 2.*(Fp.x + Fmu.x) / (Lz * 2*R);
  double Fy = 2.*(Fp.y + Fmu.y) / (Lz * 2*R);
  double Fz = 2.*(Fp.z + Fmu.z) / (Lz * 2*R);

  if (pid() == 0 && fp != NULL) {
    fprintf(fp, "%+3.2e %06d %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e %+6.5e\n",
            Re, i, t, dt, Fx, Fp.x, Fmu.x, Fy, Fz);
    fflush(fp);
  }

#endif
}

//#include "view.h"
//#include "lambda2.h"
//event movies (t += 0.2) {
//#if (dimension == 2)
//  view (fov = 7, width = 1500, height = 400, tx = -0.25);
//  scalar omega[];
//  vorticity (u, omega);
//  boundary ({omega});
//  translate (z = 0.05) {
//    draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
//    draw_vof ("cs", "fs", lw = 2);
//  }
//  squares ("omega", min = -2, max = 2, linear = true, map = cool_warm);
//  cells();
//#elif (dimension == 3)
//  scalar l2[];
//  lambda2 (u, l2);
//  view (fov = 18.1854, quat = {0.431384,-0.216693,-0.317091,0.816338},
//      tx = 0, ty = 0, bg = {0.3,0.4,0.6}, width = 1080,
//      height = 1080, samples = 3);
//  isosurface ("l2", -0.01);
//  cells (alpha = -L0/2);
//  draw_vof ("cs", "fs", fc = {0.5,0.1,0.2});
//#endif
//  char str[99];
//  sprintf (str, "Re = %g, C = %g, ML= %d", Re, c, maxlevel);
//  draw_string (str, 1, lw = 3, lc = {1, 0, 1});
//  save ("mov20003D.mp4");
//}




//event dumper (t += 10) {
//  
//  p.nodump = false;
//  char str[99];
//  sprintf (str, "dump3D%g", t);
//  dump(str);
//}

 event stop (t = 200) {
  fclose (fp);
}                      