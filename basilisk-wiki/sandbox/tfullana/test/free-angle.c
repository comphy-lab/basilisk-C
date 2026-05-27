#define JVIEW 0
#define ANGLE 0
#define ADAPT 0
#define GRID_CONV 1

#include "advection.h"
#include "contact.h"
#include "vof.h"

#define shift 0.05

#include <vofi.h>
#pragma autolink -L$HOME/local/lib -lvofi

typedef vofi_creal creal;
#define Get_fh vofi_Get_fh
#define Get_cc vofi_Get_cc

static double sphere (creal p[dimension]){
  return sq(p[0]) + sq(p[1] + shift) - sq(0.5);
}

static void vofi (scalar c, int levelmax)
{
  double fh = Get_fh (sphere, NULL, 1./(1 << levelmax), dimension, 0);
  foreach() {
    creal p[2] = {x - Delta/2., y - Delta/2.};
    c[] = Get_cc (sphere, p, Delta, fh, dimension);
  }
}

#if JVIEW
#include "display.h"
#endif

#define TEND 10.0
#define c1 0.5
#define c2 0.2
#define tau 1.0

#define radius 0.5

#define theta_0 acos(shift / radius)

#define f_err 1e-3
#define u_err 1e-3

int LEVEL = 8, MINLEVEL = 6, MAXLEVEL = 9;

scalar f[], cf[], kappa[];
scalar *interfaces = {f}, *tracers = NULL;
vector h[];

// Global variables for the boundary condition
double angle_imposed_left = pi/2., angle_imposed_right = pi/2.;

// We store the "measured" angle globally so we can log exactly what was used
double angle_recorded_left = pi/2., angle_recorded_right = pi/2.;

double kappa_num = 1./(0.5 );
double time_prev = 0.;

// --- Analytical Functions ---
double analytical_position_right(double t) {
    double s = tau * sin((pi * t) / tau) / pi;
    return radius * sin(theta_0) * exp(c1 * s);
}

double analytical_angle_right(double t) {
    double s = tau * sin((pi * t) / tau) / pi;
    return pi / 2. + atan((-1. / tan(theta_0)) * exp(2. * c1 * s) + (c2 / (2. * c1) * (exp(2. * c1 * s) - 1.)));
}

double analytical_angle_left(double t) {
    double s = tau * sin((pi * t) / tau) / pi;
    return pi / 2. + atan((-1. / tan(theta_0)) * exp(2. * c1 * s) - (c2 / (2. * c1) * (exp(2. * c1 * s) - 1.)));
}

double analytical_curvature_right(double t) {
    return (1. / (radius)) * exp(3. * c1 * sin(pi * t / tau) * (tau / pi));
}

double analytical_curvature_right2(double t) {
    double s = tau * sin((pi * t) / tau) / pi;
    double k0 = (1. / (radius));
    double cot = 1. / tan(theta_0);
    return k0*pow(1.+sq(cot),3./2.)*exp(3.*c1*s)/pow((1. + sq(cot)*exp(4.*c1*s)),3./2.);
}

double k_ode_func(double t_val, double k_val) {
    double ang = analytical_angle_right(t_val);
    double a = -cos(ang);
    double b = sin(ang);
    return -3. * k_val * cos(pi * t_val / tau) * (sq(a)*c1 + a*b*c2 - sq(b)*c1);
}

double velocity_x(double t, double x, double y) {
    return c1 * cos(pi * t / tau) * x + c2 * cos(pi * t / tau) * y;
}

double velocity_y(double t, double x, double y) {
    return -c1 * cos(pi * t / tau) * y;
}

// The Rule: Applies the current global variables to the boundary
#if ANGLE
// h.t[bottom] = contact_angle((x < 0.) ? (angle_imposed_left) : (angle_imposed_right));
h.t[bottom] = contact_angle((x < 0.) ? (analytical_angle_left(t)) : (analytical_angle_right(t)));
#else
h.t[bottom] = contact_angle((x < 0.) ? (angle_imposed_left) : (angle_imposed_right));
#endif
int main()
{
    if (pid() == 0) system("mkdir -p data");
    origin(-1., 0.);

    f.height = h;
    f.tracers = {cf};

#if GRID_CONV
    for (LEVEL = MINLEVEL; LEVEL <= MAXLEVEL; LEVEL += 1)
    {   size(2);
        init_grid(1 << LEVEL);
        
        // Reset angles to equilibrium at start of every level
        angle_imposed_left = theta_0;
        angle_imposed_right = theta_0;
        angle_recorded_left = theta_0;
        angle_recorded_right = theta_0;
        

        DT = 0.02*(2./pow(2, LEVEL))/(fabs(c1) + fabs(c2));

        time_prev = 0.;
        kappa_num = 1./(radius);
        run();
    }
#else
    init_grid(1 << LEVEL);
    DT = 0.1/(max(c1, c2)*(1 << LEVEL));
    run();
#endif
}

event init(t = 0)
{
    // fraction(f, -(sq(x) + sq(y + shift) - sq(radius)));
    // foreach ()
    //     cf[] = f[];
    // boundary({f});
	vertex scalar phi[];
	foreach_vertex()
		phi[] = sq(0.5) - (sq(x) + sq(y + shift));

	fractions (phi, f);
	vofi (f, LEVEL);
	boundary({f});
}

event impose_velocity(i++)
{
    foreach_face(x) u.x[] = velocity_x(t, x, y);
    foreach_face(y) u.y[] = velocity_y(t, x, y);
}

// --- NEW EVENT: Update Angle BEFORE VOF Solver (Removes Lag) ---
// "stability" runs early in the timestep, before advection.
event update_boundary (i++) {
    
    // 1. Calculate Curvature of the CURRENT shape
    curvature(f, kappa);

    double new_angle_left = angle_imposed_left; 
    double new_angle_right = angle_imposed_right;

    // 2. Extrapolate Angles (One Layer Above)
    foreach (reduction(max:new_angle_left) reduction(max:new_angle_right))
    {
        if ((y <= (2. * Delta)) & (y > (Delta)) & (h.x[] != nodata) & (f[] < (1. - 1e-6)) & (f[] > 1e-6))
        {
            if (x < 0.) {
                coord n = interface_normal(point, f);
                double theta_a = pi / 2. - asin(n.y);
                double tmp = pow(((h.x[0, 1] - h.x[]) / 1.0), 2);
                new_angle_left = theta_a + 1.5 * Delta * kappa[] * sqrt(1 + tmp) / sin(theta_a);
            }
            if (x > 0.) {
                coord n = interface_normal(point, f);
                double theta_a = pi / 2. - asin(n.y);
                double tmp = pow(((h.x[0, 1] - h.x[]) / 1.0), 2);
                new_angle_right = theta_a + 1.5 * Delta * kappa[] * sqrt(1 + tmp) / sin(theta_a);
            }
        }
    }

    // 3. Update Global Variables IMMEDIATELY
    // The VOF solver will run *after* this and see these new values.
#if !ANGLE
    angle_imposed_left = new_angle_left;
    angle_imposed_right = new_angle_right;
#endif
    
    // Store for logging later
    angle_recorded_left = new_angle_left;
    angle_recorded_right = new_angle_right;
}


// --- LOGGING EVENT (Runs at end of step) ---
event log_results(i++; t <= TEND)
{
    char name[80], name2[80], name3[80], name4[80], name5[80];
#if GRID_CONV
    sprintf(name, "data/positions_c1_%g_c2_%g_shift%g_ANGLE%d_LEVEL%d", c1, c2, shift, ANGLE, LEVEL);
    sprintf(name2, "data/angles_left_c1_%g_c2_%g_shift%g_ANGLE%d_LEVEL%d", c1, c2, shift, ANGLE, LEVEL);
    sprintf(name3, "data/angles_right_c1_%g_c2_%g_shift%g_ANGLE%d_LEVEL%d", c1, c2, shift, ANGLE, LEVEL);
    sprintf(name4, "data/curvatures_c1_%g_c2_%g_shift%g_ANGLE%d_LEVEL%d", c1, c2, shift, ANGLE, LEVEL);
    sprintf(name5, "data/errors_c1_%g_c2_%g_shift%g_ANGLE%d_LEVEL%d", c1, c2, shift, ANGLE, LEVEL);
#else
    sprintf(name, "data/positions");
    sprintf(name2, "data/angles_left");
    sprintf(name3, "data/angles_right");
    sprintf(name4, "data/curvatures");
    sprintf(name5, "data/errors");
#endif

    // We only calculate position/curvature here for Reporting.
    // The angle used was already calculated in 'update_boundary'.
    
    double position_left = 0., position_right = 0., position_ref_right = analytical_position_right(t);
    double curvature_left = 0., curvature_right = 0., curvature_ref_right = analytical_curvature_right2(t);
    double angle_ref_right = analytical_angle_right(t);

#if ADAPT
    adapt_wavelet({f, u}, (double[]){f_err, u_err, u_err, u_err}, LEVEL);
#endif

    // Re-calculate curvature of the RESULTING shape for logging
    curvature(f, kappa);

    foreach_boundary(bottom,
        reduction(min:position_left) reduction(max:position_right)
        reduction(max:curvature_left) reduction(max:curvature_right))
    {
        if ((h.x[] != nodata) & (f[] < (1. - 1e-6)) & (f[] > (1e-6))) {
            if (x < 0.) {
                position_left = x + Delta * height(h.x[]);
                curvature_left = kappa[];
            }
            if (x > 0.) {
                position_right = x + Delta * height(h.x[]);
                curvature_right = kappa[];
            }
        }
    }

    // RK4 Reference Update
    double my_dt = t - time_prev;
    if (my_dt > 0.0) {
        double k1 = k_ode_func(time_prev, kappa_num);
        double k2 = k_ode_func(time_prev + 0.5*my_dt, kappa_num + 0.5*my_dt*k1);
        double k3 = k_ode_func(time_prev + 0.5*my_dt, kappa_num + 0.5*my_dt*k2);
        double k4 = k_ode_func(time_prev + my_dt, kappa_num + my_dt*k3);
        kappa_num += (my_dt / 6.0) * (k1 + 2.*k2 + 2.*k3 + k4);
    }
    time_prev = t;

    if (pid() == 0) {
        static FILE *fp = fopen(name, "w");
        static FILE *fp2 = fopen(name2, "w");
        static FILE *fp3 = fopen(name3, "w");
        static FILE *fp4 = fopen(name4, "w");
        static FILE *fp5 = fopen(name5, "w");

        double position_err_right = fabs(position_right - position_ref_right)/fabs(position_ref_right);
        double angle_err_right = fabs(angle_recorded_right - angle_ref_right)/fabs(angle_ref_right);
        double curvature_err_right = fabs(curvature_right - curvature_ref_right)/fabs(curvature_ref_right);

        fprintf(fp, "%lf %lf %lf %lf\n", t, position_left, position_right, position_ref_right);
        // We log 'angle_recorded_...' which is the value we calculated at start of step
        fprintf(fp2, "%lf %lf %lf %lf\n", t, 0., angle_recorded_left, 0.); 
        fprintf(fp3, "%lf %lf %lf %lf %lf\n", t, 0., angle_recorded_right, 0., angle_ref_right);
        fprintf(fp4, "%lf %lf %lf %lf %lf %lf\n", t, curvature_left, curvature_right, curvature_ref_right, kappa_num, my_dt);
        fprintf(fp5, "%lf %lf %lf %lf\n", t, position_err_right, angle_err_right, curvature_err_right);

        fflush(fp); fflush(fp2); fflush(fp3); fflush(fp4); fflush(fp5);
    }
}