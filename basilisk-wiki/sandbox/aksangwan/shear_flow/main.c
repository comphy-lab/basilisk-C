
// **Two-Phase Shear Flow Simulation**
// This code simulates a two-phase shear flow using the Basilisk C framework.
// The setup involves two fluids with different densities, gravity, and a moving top boundary.

// Author - Anil Kumar, IIT Bombay, India (Email: ak.sangwan.iitb@gmail.com)
#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase.h"
#include "reduced.h"
#include "tension.h"
#include "output.h"
#include "distance.h" // works only in serial
#include "output_vtu_foreach.h" // Works only in serial; For vtu file format output

#define MAXLEVEL 10 // maximum level of refinement; higher value means finer grid
#define t_final 2.5 // Final time 
double tdump = 0.005; // time step to dump output files
#define density_ratio 0.5 // Density ratio upper by lower fluid
#define gG 980.0 // Gravity
#define bB 2.54 // Delta for air (Shear height in air)
#define eEps 0.02 // Steepness parameter
#define cr 61.25600429 // Real part of wave speed
#define ci 22.55092783 // Imaginary part of wave speed
#define kK 0.2920977808 // Wave number
#define aA (eEps/kK) // Amplitude
#define hH (pi / kK) // Half domain Height
#define lL (2.0 * hH) // Domain Length
#define U0 199.5675324 // Top boundary velocity
#define gridSize 1024 // Grid size for initial condition generation
char restoreFileName[100]; // Name of the restore file

// defined at the end of the script
void make_dir(char *name); // Function to make directory
void my_init_cond(); // Function to set initial condition
void make_interface(); // Function to make interface
void print_time(FILE *ptr); // Function to print time

// Boundary conditions
u.n[top] = dirichlet(0.); // Normal velocity at top boundary
u.n[bottom] = dirichlet(0.); // Normal velocity at bottom boundary
u.t[top] = dirichlet(U0); // Tangential velocity at top boundary
u.t[bottom] = neumann(0.); // Tangential velocity at bottom boundary
int main(int argc, char *argv[])
{
    if (argc > 1)
    {
        for (int i = 1; i < argc; i++)
        {
            (strcmp(argv[i], "-r") == 0) ? sprintf(restoreFileName, "dumpfile/dump-%d", atoi(argv[i + 1])) : none;
        }
    }
    periodic(right); // Periodic boundary condition at right
    G.y = -gG; // Gravity in y-direction
    f.sigma = 72.; // Surface tension
    rho1 = 1.0; // Density of fluid 1 (lower fluid)
    rho2 = density_ratio * rho1; // Density of fluid 2 (upper fluid)
    mu1 = 0.;   mu2 = 0.; // Viscosity of both fluids
    size(lL); // Domain size
    init_grid(512); // Initial grid size
    origin(-0.5 * lL, -0.5 * lL); // Origin of the domain
    make_dir("dumpfile vtufile interface"); // Make necessary directories
    run(); // Start the simulation
}

event init(t = 0) // Initialization event at time t=0
{
    // sprintf(restoreFileName, "dumpfile/dump-%d",0); // Default restore file name
    if (restore(file = restoreFileName)) // Restore from file if available
    {
        fprintf(stderr, "\n\n*******************");
        fprintf(stderr, "Restoring from file: \"%s\"", restoreFileName);
        fprintf(stderr, "*******************\n\n");
    }
    else
    {
        refine((y > -0.1 * hH && y < 0.1 * hH) && level < 10); // Initial refinement near interface
        make_interface(); // Create the interface
        my_init_cond(); // Set initial conditions
    }
}

event logfile(i++) // Logging event at each iteration
{
    if (pid() == 0) // Only the master process logs
    {
        FILE *ftr = fopen("track.out", "a"); // Open log file
        fprintf(ftr, "t: %-20g dt: %-20g", t, dt); // Log time and time step
        print_time(ftr); // Print current time
        fclose(ftr); // Close log file
        fprintf(stderr, "\033[2K\r\033[A t: %-20g dt: %-20g\n", t, dt); // Print to stderr
        
    }
}

event adaptivity_conditions(i++) // Adaptation event at each iteration
{
    scalar KAPPA[]; // Curvature scalar
    curvature(f, KAPPA); // Compute curvature
    boundary((scalar *){KAPPA}); // Apply boundary conditions
    adapt_wavelet({f, KAPPA, u}, (double[]){0.001, 1e-4, 0.01, 0.01, 0.01}, MAXLEVEL, 5); // Adaptive mesh refinement
}

int pt;
event code_output(t = 0.0; t += tdump; t <= t_final) // Output event at specified time intervals
{
    pt = (round(t / tdump));
    sprintf(restoreFileName, "dumpfile/dump-%d", pt);
    if (pid() == 0)
    {
        FILE *fp = fopen("restoreFile.txt", "w");
        fprintf(fp, "%s", restoreFileName);
        fclose(fp);
    }
    p.nodump = false;
    dump(restoreFileName); // Dumping the restore file used by Basilisk
    // For VTU output
    char name[80];
    sprintf (name, "vtufile/snap-%d", pt);
    output_vtu ((scalar *) {f, p}, (vector *) {u}, name); // Binary VTU output

    //This function “draws” interface facets in a file. The segment endpoints are defined by pairs of coordinates. Each pair of endpoints is separated from the next pair by a newline, so that the resulting file is directly visualisable with gnuplot.
    sprintf (name, "interface/int-%d.dat", pt);
    FILE *fp = fopen (name, "w");
    output_facets (f, fp); // Outputting interface facets
    fclose (fp);
}




// Below are custom defined functions

void make_dir(char *name)
{
    char s1[100] = "mkdir -p ";
    strcat(s1, name);
    system(s1);
}

void my_init_cond()
{

    char data[200];
    char str[150];
    int count;
    FILE *ptr = fopen("data.csv", "w");
    foreach ()
        fprintf(ptr, "%lf, %lf, %f\n", x, y, f[]);
    fclose(ptr);
    sprintf(data,
            "python3 genData.py -a %f -cr %f -ci %f -k %f -U0 %f -g %f -d %f -ha %f -vel",
            aA, cr, ci, kK, U0, gG, density_ratio, bB);
    system(data);
    double *xdata = malloc((gridSize * gridSize + 10) * sizeof(double));
    double *ydata = malloc((gridSize * gridSize + 10) * sizeof(double));
    // double *pdata = malloc((gridSize * gridSize + 10) * sizeof(double));
    FILE *fpt = fopen("veldata.dat", "r");
    if (fpt == NULL)
    {
        printf("Can't open file ");
    }
    else
    {
        count = 0;
        while (fgets(str, 150, fpt))
        {
            char *ptr = strtok(str, ",");
            int column = 1;
            while (ptr != NULL)
            {
                if (column == 1)
                    xdata[count] = atof(ptr);
                else if (column == 2)
                    ydata[count] = atof(ptr);
                // else if (column == 3)
                //     pdata[count] = atof(ptr);
                else
                {
                    printf("column is large");
                    fflush(stdout);
                    exit(0);
                }
                column++;
                ptr = strtok(NULL, ",");
            }
            count++;
        }
    }
    fclose(fpt);
    count = 0;
    foreach ()
    {
        u.x[] = xdata[count];
        u.y[] = ydata[count];
        // p[] = pdata[count];
        count++;
    }
    free(xdata);
    free(ydata);

    remove("data.csv");
    remove("veldata.dat");
}

void make_interface()
{
    fprintf(stderr, "Generating interface  \n"); // \033[A\33[2K\r
    char data[200];
    FILE *ptr = fopen("etadata.csv", "w");
    foreach ()
        fprintf(ptr, "%lf\n", x);

    fclose(ptr);
    sprintf(data,
            "python3 genData.py -a %f -cr %f -ci %f -k %f  -U0 %f -g %f -ha %f -eta",
            aA, cr, ci, kK, U0, gG, bB);
    system(data);
    ptr = fopen("eta.dat", "rb");
    if (ptr == NULL)
    {
        fprintf(ferr, "There is no file eta file \n");
        exit(0);
    }

    // fprintf(stderr, "Herereeeeeee");
    coord *shapedata = input_xy(ptr);
    fclose(ptr);
    scalar a[];
    distance(a, shapedata);
    vertex scalar phi[];
    foreach_vertex()
    {
        double xx = (a[] + a[-1] + a[0, -1] + a[-1, -1]) / 4.0;
        phi[] = -max(-xx, y - 0.5 * lL);
    }
    fractions(phi, f);
    remove("etadata.csv");
    remove("eta.dat");
}


void print_time(FILE *ptr)
{
    time_t s;
    time(&s);
    fprintf(ptr, "%s", ctime(&s));
}

