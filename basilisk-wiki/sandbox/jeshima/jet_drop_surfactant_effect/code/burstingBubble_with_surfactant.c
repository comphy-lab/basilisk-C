/**
 * This simulation studies the effect of surfactants on axisymmetric jet drop formation.
 *
 * Nondimensionalisation:
 * Length units: initial bubble radius
 * Mass units: density of liquid phase * (initial bubble radius)^3
 * Time units: sqrt(density of liquid phase * (initial bubble radius)^3/(clean surface tension))
 * Surfactant concentration units: initial surfactant concentration
 *
 * Based on:
 * - Alexis Berny et al. (/sandbox/aberny/bubble/burstingBubble.c)
 * - Marangoni flow package by Palas Kumar Farsoiya et al. (/sandbox/farsoiya/marangoni_surfactant/axi_rising_bubble.c)
 *
 * References: [Berny et al., 2020](#berny2020), [Farsoiya et al., 2024](#farsoiya2024)
 *
 * Note: Initialisation may not work for parallel codes. To work around this, run the first timestep in serial code, then load the dump file back to run in parallel.
 */

/*
 * To run the code:
 * Compile with: qcc -O2 -Wall -I$BASILISK/gl -L$BASILISK/gl burstingBubble_with_surfactant.c -o bub_sim -lfb_tiny -lglutils -lm
 *
 * Requirements:
 * - bubbleShape.h, findBond.h, measure.h, redistance2.h, surfactant-transport.h
 * - Create a directory called "dat" to output simulation data.
 */

#include <time.h>

/**
 * Axisymmetric Navier Stokes with surfactant
 */
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "curvature.h"
#include "integral.h"
#include "surfactant-transport.h"
#include "view.h"
#include "tag.h"
#include "distance.h"

/**
 * For performance tracking
 */
#include "navier-stokes/perfs.h"
#include "profiling.h"
#include <sys/stat.h>

/**
 * Physical parameters in the problem:
 *
 * 1. La: La = rho*gamma_c*R_0/mu^2, Laplace number
 * 2. beta: Surfactant parameter (surface tension isotherm)
 * 3. dsiginf: Non-dimensional surface tension deficit at saturation
 * 4. Bo: Bo = rho g R_0^2/gamma_c, Bond number
 * 5. E: Surfactant excess amount at the tip
 * 6. sd: Standard deviation of surfactant excess distribution
 * 7. Pe: Peclet number of the flow
 *
 * [rho: density, gamma_c: clean surface tension, R_0: initial bubble radius, mu: dynamic viscosity]
 * Surface tension isotherm (nondimensional): sigma(gamma) = 1 - dsiginf * tanh((beta/dsiginf) * gamma)
 */

// typical values
double La = 2400;
double beta = 0.3;
double dsiginf = 0.55;
double Bond = 0.16;
double E = 0.5;
double sd = 0.02;
double Pe = 100;

/**
 * The convention for nondimensionalisation is that the initial surfactant concentration = 1
 */
double gamma_0 = 1;

/**
 * The simulation ends at t = END
 */
#define END 1.4

/**
 * AdaptOn is the option for AMR
 */
#define AdaptOn 1

/**
## Numerical Parameters

LEVEL = maximum refinement of the grid

initial_level = refinement level for calculating the initial bubble profile shapeBond
[Note: the initial grid taken for the simulation is adjusted after the initial bubble profile is calculated]

L0 = domain length
*/

// Level used in the main text
// int LEVEL = 12;
// int initial_level = 11;

// Level for debugging
int LEVEL = 10;
int initial_level = 8;

// domain size
#define L0 10.

/**
 * Minimum level of refinement for AMR
 */
#define MinLevel 4

/**
 * DumpSimulation: option to dump simulation when a new liquid phase is created (Berny et al. 2020)
 */
#define DumpSimulation 1

/**
 * Boundary Condition: liquid passes through the top
 */
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

/**
 * Parameters for adaptive mesh refinement
 * intemax: refinement for f
 * umax: refinement for u
 * c1emax: refinement for c1
 * pfemax: refinement for pfield
 */
double intemax = 0.005;
double uemax = 0.01;
double c1emax = 0.001;
double pfemax = 0.001;

/**
 * Bview option to use Basilisk to output data
 */
#define Bview 1

double tipcenter_r = 0.;
double tipcenter_z = 0.;
scalar sigmaf[];
int main(int argc, char const *argv[])
{

    size(L0);
    origin(-L0 / 2., 0.);
    init_grid(1 << initial_level); // initial grid for computing initial bubble shape (note: not the actual initial grid, which is adjusted afterwards)

    /**
     Load in run time parameters
     */
    if (argc >= 2)
    {
        La = atof(argv[1]);
    }

    if (argc >= 3)
    {
        beta = atof(argv[2]);
    }

    if (argc >= 4)
    {
        dsiginf = atof(argv[3]);
    }

    if (argc >= 5)
    {
        Bond = atof(argv[4]);
    }

    if (argc >= 6)
    {
        E = atof(argv[5]);
    }

    if (argc >= 7)
    {
        sd = atof(argv[6]);
    }

    /**
    Physical parameters must be positive.*/

    assert(La > 0.);
    assert(dsiginf > 0.);

    double Oh = sqrt(1. / La);

    /**
    Physical parameters of fluid:
    Fluid 1 = liquid, fluid 2 = air
    */

    rho1 = 1., mu1 = Oh;
    rho2 = 1. / 998., mu2 = Oh / 55.;
    d.sigmaf = sigmaf;

    // turn down tolerance for Poisson
    TOLERANCE = 1e-6 [*];
    // surfactant properties, see Marangoni.h
    D_s = 1. / Pe;
    surfactant = 1;

    // since using clsvof, need to add onto existing tracers
    tracers = list_concat(tracers, {c1, d2, pfield});

    /**
    We print the characteristics of the fluid in the log file. */
    fprintf(stderr, "props %f %f %f %f %f %i\n \n",
            rho2, mu2, La, beta, dsiginf, LEVEL);

    // perfs update (from Palas Farsoiya)
    // keep old perfs files, but rename
    struct stat st = {0};
    char name[80];
    char newname[80];
    sprintf(name, "perfs");

    if (stat(name, &st) == 0)
    {
        sprintf(newname, "perfs-%ld", time(0));
        rename(name, newname);
    }

    run();
}

/**
 * Gravity is applied here
 */
#define Gravity 1
#if Gravity
event acceleration(i++)
{

    face vector av = a;

    foreach_face(x)
        av.x[] += Bond * (-1. + alphav.x[] * rho1 / max(y, 1e-20));

    boundary((scalar *){av});
}
#endif
/**

## Initialisation


Initialise the bubble shape using findBond.h.

*/

#include "findBond.h"

event init(t = 0)
{
    // If "restart" exists, load "restart" in order to restart from last checkpoint.
    // At the moment, the first timestep has to be done serially. After the first timestep, load back in to run in parallel.
    if (!restore("restart", list = all))
    {
        /**
        Initialisation

        Create the bubble shape, and then adjust grids appropriately
         */

        coord *dataShape; // will contain the coordinates of the interface

        Circle *hollow = NULL;
        Circle fillet;
        hollow = &fillet;

        Circle *cap = NULL;
        Circle topCap;
        cap = &topCap;

        // The initial bubble shape is set by Bond_eff
        double Bond_eff;
        Bond_eff = Bond / (1. - dsiginf * tanh((beta / dsiginf)));

        dataShape = shapeBond(Bond_eff, hollow, cap);
        // store tip location
        tipcenter_r = hollow->x;
        tipcenter_z = hollow->y;

        /**
        As per Berny et al.:
        We are using the distance function, applied on the list of coordinate
        dataShape. This will generate a distance function. */

        scalar d_temp[];
        distance(d_temp, dataShape);

        foreach ()
        {
            // dataShape does the opposite sign
            d[] = -d_temp[];
        }
        /**
        Create the initial grid manually.
        Take grid points refined near the interface
        */
        refine(fabs(d[]) < 0.15 && level < (LEVEL - 2));
        refine(fabs(d[]) < 0.1 && level < (LEVEL - 1));
        refine(fabs(d[]) < 0.05 && level < LEVEL);
        unrefine(fabs(x) > 0.2 && fabs(y) > 1.5 && level > 6);
        unrefine(fabs(x + 1) > 1.5 && level > 6);

        /**
        As per Berny et al., output initial profile */
        static FILE *fp = fopen("initial.png", "w");
        output_ppm(f, fp, 400, min = 0, max = 1, box = {{-2.5, 0}, {2.5, 3}}); // a png  file, without bview
        scalar *list = {d};
        clear();
        view(fov = 4.69417, quat = {0, 0, -0.707108, 0.707106}, tx = 0., ty = 0.0855609, bg = {1, 1, 1}, width = 1280, height = 720, samples = 4);
        draw_vof("f", filled = 1);
        mirror({0, 1, 0}, 0)
        {
            draw_vof("f", filled = 1);
        }
        save("initial.ppm");
        // Output initial grid
#if Bview
        // The whole grid
        static FILE *fpinit = fopen("i0grid.png", "w");
        scalar l[];
        foreach ()
            l[] = level;
        clear();
        view(fov = 12.699, quat = {0, 0, 0.707107, -0.707107}, ty = 0.00253494, bg = {1, 1, 1}, width = 1920, height = 1080, samples = 4);
        draw_vof("f", filled = 1);
        squares("l", min = 1, max = LEVEL);
        mirror({0, 1, 0}, 0)
        {
            draw_vof("f");
            squares("l", min = 1, max = LEVEL);
        }
        save(fp = fpinit);
        // A more zoomed in grid
        static FILE *fplzinit = fopen("i0gridzoom.png", "w");
        clear();
        view(fov = 4.42681, quat = {0, 0, 0, 1}, tx = 0.0579169, bg = {1, 1, 1}, width = 1920, height = 1080, samples = 4);
        draw_vof("f");
        squares("l", min = 1, max = LEVEL);
        cells();
        mirror({0, 1, 0}, 0)
        {
            draw_vof("f");
            squares("l", min = 1, max = LEVEL);
            cells();
        }
        save(fp = fplzinit);

        // An even more zoomed in grid
        static FILE *fplzmoreinit = fopen("i0gridmorezoom.png", "w");
        clear();
        view(fov = 1, quat = {0, 0, 0, 1}, ty = 0.025, bg = {1, 1, 1}, width = 1920, height = 1080, samples = 4);
        draw_vof("f");
        squares("l", min = 1, max = LEVEL);
        cells();
        mirror({0, 1, 0}, 0)
        {
            draw_vof("f");
            squares("l", min = 1, max = LEVEL);
            cells();
        }
        save(fp = fplzmoreinit);
#endif

        event("properties2"); // Populate pfield
    }
    else
    { // restart dump
        FILE *fpdump = fopen("dumpRestore.png", "w");
        output_ppm(cm, fpdump, 400);
    }
}

/**
 * Initialise c1 and sigmaf
 * (sigmaf and c1 had to be set separately from the init loop)
 */
event initialise_surfactant(i = 0)
{
    foreach ()
    {
        double deltas = (pfield[] * (1. - pfield[])) / EPSILON;
        c1[] = gamma_0 * (1 + E * (1 / (sd * sqrt(2 * pi))) * exp(-(1 / (2 * sd * sd)) * (pow((x - tipcenter_z), 2) + pow((y - tipcenter_r), 2)))) * deltas;
        double gamma = c1[] * 4. * EPSILON;
        sigmaf[] = 1. - dsiginf * tanh((beta / dsiginf) * gamma);
    }
}

/**
 * Isotherm event for Marangoni effects
 */
event stability(i++)
{
    foreach ()
    {
        double gamma = c1[] * 4. * EPSILON;
        sigmaf[] = 1. - dsiginf * tanh((beta / dsiginf) * gamma);
    }
}

/**
 * Adaptation
 */
#if AdaptOn
event adapt(i++)
{
    adapt_wavelet({f, u, pfield, c1}, (double[]){intemax, uemax, uemax, uemax, pfemax, c1emax},
                  maxlevel = LEVEL, minlevel = MinLevel);
}
#endif

/**
 * Dump after one timestep (for parallelisation workaround)
 */
event secondtimestep(i = 1)
{
    dump("restart");
}

/**
 * Output interface videos (interface, omega, level)
 * Note: omega video may be spurious at r = 0 due to postprocessing error.
 * For the actual omega field, load in the dump file and extract vorticity from the velocity directly. Included only for reference.
 */
event outputInterface(t += 0.005)
{

#if Bview
    scalar omega[];
    vorticity(u, omega);
    /**
    Interface evolution
    */
    {
        static FILE *fp1 = popen("ppm2mp4 interface.mp4", "w");
        clear();
        view(fov = 12.699, quat = {0, 0, 0.707107, -0.707107}, ty = 0.00253494, bg = {1, 1, 1}, width = 1920, height = 1080, samples = 4);
        draw_vof("f", filled = 1);
        mirror({0, 1, 0}, 0)
        {
            draw_vof("f", filled = 1);
        }
        save(fp = fp1);
    }
    /*
    Vorticity evolution
    */
    {
        static FILE *fp2 = popen("ppm2mp4 omega.mp4", "w");
        clear();
        view(fov = 12.699, quat = {0, 0, 0.707107, -0.707107}, ty = 0.00253494, bg = {1, 1, 1}, width = 1920, height = 1080, samples = 4);
        draw_vof("f");
        squares("omega", map = jet, min = -50, max = 50, linear = true);
        mirror({0, 1, 0}, 0)
        {
            draw_vof("f");
            squares("-omega", map = jet, min = -50, max = 50, linear = true);
        }
        save(fp = fp2);
    }

    /**
    Grid evolution
    */
    {
        static FILE *fp3 = popen("ppm2mp4 level.mp4", "w");
        scalar l[];
        foreach ()
            l[] = level;
        clear();
        view(fov = 12.699, quat = {0, 0, 0.707107, -0.707107}, ty = 0.00253494, bg = {1, 1, 1}, width = 1920, height = 1080, samples = 4);
        draw_vof("f", filled = 1);
        squares("l", min = 1, max = LEVEL);
        mirror({0, 1, 0}, 0)
        {
            draw_vof("f");
            squares("l", min = 1, max = LEVEL);
        }
        save(fp = fp3);
    }
    /**
    Grid evolution zoomed
    */
    {
        static FILE *fplz = popen("ppm2mp4 levelZoom.mp4", "w");
        scalar l[];
        foreach ()
            l[] = level;
        clear();
        view(fov = 4.42681, quat = {0, 0, 0, 1}, tx = 0.0579169, bg = {1, 1, 1}, width = 1920, height = 1080, samples = 4);
        draw_vof("f");
        squares("l", min = 1, max = LEVEL);
        cells();
        mirror({0, 1, 0}, 0)
        {
            draw_vof("f");
            squares("l", min = 1, max = LEVEL);
            cells();
        }
        save(fp = fplz);
    }

    /**
    Vorticity evolution zoomed*/
    {
        static FILE *fpz = popen("ppm2mp4 omegaZoom.mp4", "w");
        clear();
        view(fov = 6.69986, quat = {0, 0, 0, 1}, tx = -0.05, ty = 0, bg = {1, 1, 1}, width = 1920, height = 1080, samples = 4);

        squares("omega", min = -200, max = 200, map = jet, linear = true);
        // squares("omega", min = -200, max = 200, map = bwr, linear = true);
        draw_vof("f");
        mirror(n = {0, 1, 0}, alpha = 0)
        {
            squares("-omega", min = -200, max = 200, map = jet, linear = true);
            // squares("omega", min = -200, max = 200, map = bwr, linear = true);
            draw_vof("f");
        }
        save(fp = fpz);
    }
#endif
    dump("restart");
}

/**
 * Output cell Peclet number
 * Large cell Peclet number can cause issues with surfactant transport (see Farsoiya et al. 2024)
 */

event pecletoutput(t += 0.001; t <= END)
{
    double maxu_local = 0.;
    double dmin_local = 1.;
    foreach (reduction(max : maxu_local) reduction(min : dmin_local))
    {
        double absu = sqrt(u.x[] * u.x[] + u.y[] * u.y[]);
        if (absu > maxu_local)
        {
            maxu_local = absu;
        }
        if (Delta < dmin_local)
        {
            dmin_local = Delta;
        }
    }

#if _MPI
    double maxu_global, dmin_global;
    MPI_Allreduce(&maxu_local, &maxu_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&dmin_local, &dmin_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else
    double maxu_global = maxu_local;
    double dmin_global = dmin_local;
#endif

    if (pid() == 0)
    {
        FILE *log_pec = fopen("peclet", "a");
        double cellpec = dmin_global * maxu_global * Pe;
        fprintf(log_pec, "%g %g\n", t, cellpec);
        fclose(log_pec);
    }
}

/**
 * Dump files infrequently (dump files are large)
 */
event simuInfo(t += 0.1; t <= END)
{
    vector h2[];
    heights(f, h2);

    scalar omega[];
    vorticity(u, omega);

    char dumpFile[80];
    sprintf(dumpFile, "simulation-%07d", i);
    dump(file = dumpFile);
}

/**
 * Interface output
 *
 *
 *
 * Column 1: r coordinate
 * Column 2: ignore (for 3D extension)
 * Column 3: gamma value at the interface (see note below)
 * Column 4: z coordinate
 *
 * Note that the value of gamma is not exactly the interfacial surfactant concentration Gamma.
 * This is because the gamma field is a volumetric quantity as the code is VoF rather than interface tracking.
 * More precisely, gamma = Gamma only when phi = 1/2. See Farsoiya et al. 2024 for details.
 * In practice for sufficiently fine grid, Gamma can be obtained by taking the local maxima of gamma should suffice since in a given number of fixed grids, it is likely that phi is approximately 1/2.
 *
 * Code by Nicolo Scapin
 */
int compare(const void *a, const void *b)

{

    double x = *(double *)a;
    double y = *(double *)b;

    if (x < y)
    {
        return -1;
    }
    else if (x > y)
    {
        return +1;
    }
    else
    {
        return +0;
    }
}
void output_int_qtn(char *fname, int istep, scalar c, scalar sca, double stp_eta)
{

    // We first loop over all the interfacial points
    // and we count them (per processor)

    int int_pt = 0;
    foreach (serial)
    {
        if (interfacial(point, c))
        {
            if (point.level == LEVEL)
            {
                coord n = interface_normal(point, c), pp;
                double alpha1 = plane_alpha(c[], n);
                plane_area_center(n, alpha1, &pp);

                double xc = x + Delta * pp.x;
                double yc = y + Delta * pp.y + stp_eta;

#if dimension > 2
                double zc = z + Delta * pp.z; // we keep here to avoid warning (otherwise: unused variable)
                Point point = locate(xc, yc, zc);
#else
                Point point = locate(xc, yc);
#endif
                if (point.level > 0)
                {                    // best case
                    POINT_VARIABLES; // older Basilisk version
                    // POINT_VARIABLES(); // latest version
                    int_pt++;
                }
            }
        }
    }

    int tot_column = 4;
    double t_mat[int_pt][tot_column];

    for (int j = 0; j < tot_column; j++)
    {
        for (int i = 0; i < int_pt; i++)
        {
            t_mat[i][j] = 0;
        }
    }
    fprintf(stderr, "First pass over interfacial points\n");

    scalar pos[];
#if dimension > 2
    coord G = {0., 1., 0.}, Z = {0., 0., 0.};
#else
    coord G = {0., 1.}, Z = {0., 0.};
#endif
    position(c, pos, G, Z);

    int count = 0;
    foreach (serial)
    {
        if (interfacial(point, c))
        {
            if (point.level == LEVEL)
            {
                coord n = interface_normal(point, c), pp;
                double alpha1 = plane_alpha(c[], n);
                plane_area_center(n, alpha1, &pp);
                double eta = pos[];
                double xc = x + Delta * pp.x;
                double yc = y + Delta * pp.y + stp_eta;
                double zc = z + Delta * pp.z;
#if dimension > 2
                double zc = z + Delta * pp.z; // we keep here to avoid warning (otherwise: unused variable)
                Point point = locate(xc, yc, zc);
#else
                Point point = locate(xc, yc);
#endif
                if (point.level > 0)
                { // best case

                    POINT_VARIABLES; // older basilisk version
                    // POINT_VARIABLES(); // latest basilisk version

                    t_mat[count][0] = xc;

#if dimension > 2
                    t_mat[count][1] = zc;
#else
                    t_mat[count][1] = 0;
#endif
                    // t_mat[count][2] = my_interpolation(point, sca, xc, yc, zc);
                    t_mat[count][2] = sca[];
                    // to mesaure from x= 0
                    t_mat[count][3] = eta; // already defined at the interface

                    count++;
                }
            }
        }
    }
    fprintf(stderr, "Second pass over interfacial points\n");

    // We sort locally t_mat by the x coordinate (the first, i.e. 0-th, column)

    qsort(t_mat, int_pt, tot_column * sizeof(double), compare);
    fprintf(stderr, "First sort\n");

    double tot_row;

#if _MPI

    // On multiple cores, we gather all the int_pt to the root pid
    // and, then, we broadcast this information to all the processes

    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    int counts_it[nproc];
    if (pid() == 0)
    {
        MPI_Gather(&int_pt, 1, MPI_INT, counts_it, 1, MPI_INT, 0, MPI_COMM_WORLD); // MPI_gather gathers by rank order
    }
    else
    {
        MPI_Gather(&int_pt, 1, MPI_INT, NULL, 1, MPI_INT, 0, MPI_COMM_WORLD); // MPI_gather gathers by rank order
    }
    MPI_Bcast(counts_it, nproc, MPI_INT, 0, MPI_COMM_WORLD);

    // Each processor knows the int_pt of the others. So we can compute
    // the displacement, the total number of interfacial points and
    // the number of elements owned by each processor

    int tot_int_pt = 0;
    int tot_el_p[nproc], disp_r[nproc], disp[nproc];
    for (int i = 0; i < nproc; i++)
    {
        disp_r[i] = tot_int_pt;
        disp[i] = disp_r[i] * tot_column;
        tot_int_pt += counts_it[i];
        tot_el_p[i] = counts_it[i] * tot_column;
    }
    tot_row = tot_int_pt;

    // --> Gather to the root pid
    // --> Sort by first column

    double t_mat_tot[tot_int_pt][tot_column];

    if (pid() == 0)
    {
        int tot_el_p0[nproc], disp0[nproc];
        for (int i = 0; i < nproc; i++)
        {
            tot_el_p0[i] = tot_el_p[i];
            disp0[i] = disp[i];
        }
        MPI_Gatherv(&t_mat, int_pt * tot_column, MPI_DOUBLE, t_mat_tot, tot_el_p0, disp0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        fprintf(stderr, "GatherV 0\n");
        qsort(t_mat_tot, tot_int_pt, tot_column * sizeof(double), compare);
        fprintf(stderr, "Second sort 0\n");
    }
    else
    {
        int tot_el_p1[nproc], disp1[nproc];
        for (int i = 0; i < nproc; i++)
        {
            tot_el_p1[i] = tot_el_p[i];
            disp1[i] = disp[i];
        }
        MPI_Gatherv(&t_mat, int_pt * tot_column, MPI_DOUBLE, NULL, tot_el_p1, disp1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

#else

    tot_row = int_pt;
    double t_mat_tot[int_pt][tot_column];

    for (int j = 0; j < tot_column; j++)
    {
        for (int i = 0; i < int_pt; i++)
        {
            t_mat_tot[i][j] = t_mat[i][j];
        }
    }

#endif

    if (pid() == 0)
    {

        fflush(stderr);
        FILE *eta_loc = fopen(fname, "w");

        // Print in ASCII format

        for (int i = 0; i < tot_row; i++)
        {
            fprintf(eta_loc, "%8E %8E %8E %8E\n",
                    t_mat_tot[i][0], t_mat_tot[i][1], t_mat_tot[i][2], t_mat_tot[i][3]);
        }

        // Print in binary format

        /*
        for (int i = 0; i < tot_int_pt; i++) {
          for (int j = 0; j < tot_column; j++) {
      fwrite ( &t_mat_tot[i][j], sizeof(double), 1, eta_loc );
          }
        }
        */

        fclose(eta_loc);
        fflush(eta_loc);
    }
    fprintf(stderr, "Print\n");
}

/**
 * Save the interfacial values
 */
// event eta_loc(t += (t >= 0.35 && t <= 0.55) ? 0.001 : 0.01; t <= END) // use this for more frequent outputs close to curvature reversal
event eta_loc(t += 0.01; t <= END)

{
    foreach ()
    {
        gamma2[] = c1[] * 4. * EPSILON;
    }
    fflush(stderr);
    char etaname_in[100];
    sprintf(etaname_in, "dat/interface-vals-%g.txt", t);
    double stp_eta = 0.;
    boundary({gamma2});
    output_int_qtn(etaname_in, i, f, gamma2, stp_eta);
    // interface-vals contains the interface points and surfactant concentraion at the interface point
    // Column 1: r coordinate
    // Column 2: ignore (this option is for 3D)
    // Column 3: gamma
    // Column 4: z coordinate

    fflush(stderr);
    if (pid() == 0)
    {

        char name_1[80];
        sprintf(name_1, "dat/log_eta.out");
        FILE *log_sim = fopen(name_1, "a");
        fprintf(log_sim, "%8E %8E\n", t, 1.0 * i);
        fclose(log_sim);
    }
}

/**
 * Time step detection utility
 * Used to detect a time step for output at desired times without forcing solver to reach it precisely
 * As per Berny et al.
 */
int timeDetect(double t, double tPrevious, double target)
{
    if (tPrevious == -1)
        return -1;
    if ((tPrevious - target < 0.) && (t - target > 0.))
        return 1;
    else
        return -1;
}

/**
 * Droplet measurement data (from measure.h)
 */
#define MEASURE 1
#if MEASURE
#include "measure.h"
#endif

/**
 * Final state output
 */
event finalState(t = end)
{
    dump(file = "final");
}

/**
## References

~~~bib
@article{berny2020,
  title={Role of all jet drops in mass transfer from bursting bubbles},
  author={Berny, Alexis and Deike, Luc and S{\'e}on, Thomas and Popinet, St{\'e}phane},
  journal={Physical Review Fluids},
  volume={5},
  number={3},
  pages={033605},
  year={2020},
  publisher={APS}
}
@article{farsoiya2024,
  title={Coupled volume of fluid and phase field method for direct numerical simulation of insoluble surfactant-laden interfacial flows and application to rising bubbles},
  author={Farsoiya, Palas Kumar and Popinet, St{\'e}phane and Stone, Howard A and Deike, Luc},
  journal={Physical Review Fluids},
  volume={9},
  number={9},
  pages={094004},
  year={2024},
  publisher={APS}
}
~~~
*/
