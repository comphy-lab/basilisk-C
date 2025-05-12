/**
# Stokes flow past a periodic array of cylinders
It was compared the numerical results with the solution given by the
multipole expansion of [Sangani and Acrivos, 1982](#sangani1982). [This example](http://basilisk.fr/src/test/cylinders.c) is inspired me checking the module. There is a cylinder with radius r and periodic boundary conditions. Initial velocity at the whole domain is equal to 0. The external force $\mathbf{g}$ drives the fluid. 

![Sangani and Acrivos problem](https://i.ibb.co/hF1tmXD/Sangani.png)

Here we switch on the Brinkman Penalization Module with debugging mode, where defined value of Brinkman Penalization term dbp[], total Right Hand Side of Navier--Stokes equation total_rhs[], tangential component of the velocity utau[], gradient of tangential velocity grad_utau_n[].

*/

#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1

#undef SEPS
#define SEPS 1e-30

/*
Debug mode for MinMaxValues: output min and max value of chosen fields. 
*/
//#define DEBUG_MINMAXVALUES
/*
After exporting vtu files, write to stderr some information: iteration of output iter_fp, time t, step of time dt
*/
//#define DEBUG_OUTPUT_VTU_MPI

//fs[] is a mask field for solids, omega[] is a vorticity field, Us[] is a solid velocity field.
scalar fs[], omega[];
vector Us[];

#include "centered-weugene.h"
#include "view.h"
#include "output_vtu_foreach.h"

/**
This is Table 1 of [Sangani and Acrivos, 1982](#sangani1982), where
the first column is the volume fraction $\Phi$ of the cylinders and
the second column is the non-dimensional drag force per unit length of
cylinder $F/(\mu U)$ with $\mu$ the dynamic vicosity and $U$ the
average fluid velocity. */

static double sangani[9][2] = {
        {0.05, 15.56},
        {0.10, 24.83},
        {0.20, 51.53},
        {0.30, 102.90},
        {0.40, 217.89},
        {0.50, 532.55},
        {0.60, 1.763e3},
        {0.70, 1.352e4},
        {0.75, 1.263e5}
};

/**
We will vary the maximum level of refinement, *nc* is the index of the
case in the table above, the radius of the cylinder will be computed
using the volume fraction $\Phi$. */

int maxlevel = 10, minlevel = 4, nc;
double radius;

void calc_solid(scalar fs){
  vertex scalar phi[];
  face vector face_fs[];
  foreach_vertex() {
    phi[] = (sq(x) + sq(y) < sq(radius) ) ? 1 : -1;
  }
  boundary ({phi});
  fractions (phi, fs, face_fs);
  boundary ({fs});
}

int main(int argc, char * argv[]){
    eta_s =1e-6;
    if (argc > 1) {
      eta_s = atof(argv[1]); //convert from string to float
    }
    if (argc > 2) {
      maxlevel = atoi(argv[2]); //convert from string to float
    }
    fprintf(fout, "eta_s=%g maxlevel=%d", eta_s, maxlevel);
    /**
    The domain is the periodic unit square, centered on the origin. */
    L0=1;
    origin (-L0/2., -L0/2.);
    periodic (right);
    periodic (top);

    /**
    We turn off the advection term.*/

    stokes = true;
    DT = 1e-3;
    TOLERANCE = 1e-8;
    NITERMIN = 5;
    /**
    We do the 9 cases computed by Sangani & Acrivos. The radius is
    computed from the volume fraction. */

    for (nc = 0; nc < 9; nc++) {
        N = 1 << maxlevel;
        radius = sqrt(sq(L0) * sangani[nc][0] / pi);
        run();
    }
}

/**
We need an extra field to track convergence. */

scalar un[];

event init (t = 0){
    int it = 0;
    do {
        calc_solid(fs);
    }while (adapt_wavelet({fs}, (double []){0.001},
                          maxlevel = maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
    /**
    And set acceleration and viscosity to unity. */

    const face vector g[] = {1.,0.};
    a = g;
    mu = fm;

    /**
    We initialize the reference velocity. */
    foreach() un[] = u.y[];
    event("vtk_file");
}

/**
We check for a stationary solution. */

event logfile (i++; i <= 5000){
    double avg = normf_weugene(u.x, fs).avg;
    double du = change_weugene (u.x, un, fs)/(avg + SEPS); //change 1) Linf  2) un = u
    fprintf (fout, "%d %d %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
    maxlevel, i,
    mgp.i, mgp.nrelax, mgp.minlevel,
    mgu.i, mgu.nrelax, mgu.minlevel,
    du, mgp.resa*dt, mgu.resa, statsf_weugene(u.x, fs).sum, normf_weugene(p, fs).max);
    fflush(fout);

    if (((i > 1) && (du/dt < 1e-2))|| (i == 5000)) {
        /**
        We output the non-dimensional force per unit length on the
        cylinder $F/(\mu U)$, together with the corresponding value from
        Sangani & Acrivos and the relative error. */
        stats s = statsf_weugene(u.x, fs);
        double Phi = 1. - s.volume/sq(L0);
        double Phia = pi*sq(radius)/sq(L0);
        double U = s.sum/s.volume;
        double F = sq(L0)/(1. - Phi);
        fprintf (ferr, "%d %g %g %g %g i=%d ifp=%d| %g=%g? F:%g dt:%g t:%g U:%g Uw:%g Ua:%g\n", maxlevel, sangani[nc][0], F/U, sangani[nc][1], fabs(F/U - sangani[nc][1])/sangani[nc][1], i, iter_fp, Phi, Phia, F, dt, t, U);
        fprintf(fout, "stationary flow nc = %d i = %d du = %g", nc, i, du);
        fflush(fout);
        return 9;
    }
}

//Output
event vtk_file (t += 0.01){
    char subname[80]; sprintf(subname, "br");//it sets prefix of vtu files
    scalar l[];
    vorticity (u, omega);
    foreach() {l[] = level; omega[] *= 1 - fs[];}
//    output_vtu_MPI( (scalar *) {l, omega, fs, p}, (vector *) {u, uf, dbp, utau, grad_utau_n, n_sol, target_U}, subname, L0/pow(2, minlevel));
    output_vtu_MPI( (scalar *) {l, omega, fs, p}, (vector *) {u, uf, dbp, utau}, subname, L0/pow(2, minlevel));
}

//we defined an array of adaptation fields and their threshold values.
#define ADAPT_SCALARS {fs, omega}
#define ADAPT_EPS_SCALARS {1e-3, 1e-2}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
    //calc_solid(fs);// we can recalculate the value of solid value. (It is necessary only at the beggining, when refinement level at the interface changes and not necessary for subsequent time steps)
}

/**
The non-dimensional drag force per unit length closely matches the
results of Sangani & Acrivos. For $\Phi=0.75$ and level 8 there is
only about 6 grid points in the width of the gap between cylinders.

![Popinet's trick](https://www.docdroid.net/esTiv7K/popinet-trick-l8-10.pdf)
![L=8 convergence](https://www.docdroid.net/2EvBbFm/pic8.pdf)
![L=9 convergence](https://www.docdroid.net/HBjrYsE/pic9.pdf)
![L=10 convergence](https://www.docdroid.net/VMvpDMS/pic10.pdf)
![L=11 convergence](https://www.docdroid.net/vt3cveG/pic11.pdf)
![Comparison of methods](https://www.docdroid.net/XTa44Tn/1e-6.pdf)

## References

~~~bib
@article{sangani1982,
  title={Slow flow past periodic arrays of cylinders 
  with application to heat transfer},
  author={Sangani, AS and Acrivos, A},
  journal={International Journal of Multiphase Flow},
  volume={8},
  number={3},
  pages={193--206},
  year={1982},
  publisher={Elsevier}
}
~~~
*/