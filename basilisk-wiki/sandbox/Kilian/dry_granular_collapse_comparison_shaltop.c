/**
# Dry granular flow with a SHALTOP-type topographic correction

## Introduction

<p align="justify">
This code was developed as part of a first-year Master's internship on the
numerical modelling of granular landslides. It forms one part of a broader
study devoted to shallow-layer models, friction laws, model validation and
parametric investigations of deposit mobility.
</p>

<p align="justify">
This program simulates in one dimension a **pile of glass beads** initially
held back by a gate, then released onto an inclined slope followed by a
horizontal section. The geometry and initial dimensions correspond to the
dry 35° case of Poulain *et al.* (2023), for glass beads of
diameter $d=4$ mm. This code is designed to compare a Basilisk implementation
of a **SHALTOP-type granular model** with the reference profiles reported in
the article.
</p>

<p align="justify">
The conservative transport step is still handled by Basilisk's
`saint-venant.h` solver, but the variables and source terms are rewritten so
that they mimic the one-dimensional SHALTOP formulation: the avalanche
thickness is measured **normal to the topography**, the velocity is tangent to
the bed, and the curvature-dependent correction of the normal pressure is kept
in the basal friction term.
</p>

<p align="justify">
This is therefore not the official SHALTOP code. It is a Basilisk
implementation kept as close as possible to the previously documented
[HySEA-comparison version](https://basilisk.fr/sandbox/Kilian/dry_landslide_comparison_hysea.c),
so that the numerical framework and comparison workflow remain almost
unchanged.
</p>

## Geometry sketch

<p align="justify">
Figure 1 introduces the main geometric notations used in the following, in particular those associated with the initial pile, the topographic transition and the slope.
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/kilian-gthub/basilisk-ressources/main/ressources/pente_droite.png"
       width="70%">
</p>

<p align="center">
*Figure 1 - Sketch of the initial geometry used in the simulations: horizontal bottom, smoothed transition, slope of angle $\theta$, initial granular pile and main geometric quantities of the domain.*
</p>

## SHALTOP local geometry

<p align="justify">
The physical topography is denoted by $z_b(x)$. The SHALTOP model is written
using the avalanche thickness $h_n$ measured in the direction normal to this
bed. The local geometric factor and its derivative are
</p>

$$
\displaystyle
c
=
\cos\theta_{\rm loc}
=
\frac{1}{\sqrt{1+\left(\partial_x z_b\right)^2}},
\qquad
c_x
=
-\left(\partial_x z_b\right)
 \left(\partial_{xx}z_b\right)c^3.
$$

<p align="justify">
In order to keep Basilisk's conservative Saint-Venant transport step, the code
uses the variables
</p>

$$
\displaystyle
q=\frac{h_n}{c},\qquad v=c u_t,
$$

<p align="justify">
where $u_t$ is the tangential velocity. The physical vertical height displayed
in the output is therefore
</p>

$$
\displaystyle
h_{\rm vert}=h_n c=c^2 q,
\qquad
z_s=z_b+h_{\rm vert}=z_b+c^2q.
$$

<p align="center">
  <img src="https://raw.githubusercontent.com/kilian-gthub/basilisk-ressources/main/ressources/schema_topographie_shaltop.png"
       width="40%">
</p>

<p align="center">
*Figure 2 - Local SHALTOP geometry. The avalanche thickness $h_n$ is measured
normal to the bed $z_b(x)$, while $u_t$ is the tangential velocity. The local
angle determines the geometric factor $c=\cos\theta_{\rm loc}$.*
</p>

<p align="justify">
The local cosine is not specific to SHALTOP. The granular Froude number used
in the Pouliquen--Forterre law already contains the factor $\cos\theta$ for a
uniform inclined plane. In the
[HySEA-comparison code](https://basilisk.fr/sandbox/Kilian/dry_landslide_comparison_hysea.c),
this angle is evaluated locally from the bed slope. In the present
SHALTOP-type formulation, however, the same geometric factor $c$ also enters
the definition of the normal thickness, the velocity variables, the pressure
terms and the curvature correction.
</p>

## Solved equations

<p align="justify">
The conservative mass equation becomes
</p>

$$
\displaystyle
\frac{\partial q}{\partial t}+
\frac{\partial(qv)}{\partial x}=0.
$$

<p align="justify">
The momentum equation implemented here is the one-dimensional SHALTOP-type
correction written on the horizontal velocity $v$:
</p>

$$
\displaystyle
\frac{\partial v}{\partial t}+v\frac{\partial v}{\partial x}
= c\,c_x\,u_t^2
-Gc^2\frac{\partial}{\partial x}\left(z_b+c^2q\right)
-G\mu(h_n,Fr)c^2\,\mathrm{sign}(v)
\left[
1+
\frac{(\partial_{xx}z_b)u_t^2}{Gc}
\right]_+.
$$

<p align="justify">
The first two terms on the right-hand side are applied in the `shaltop_geometry`
event. The last term is the basal Pouliquen--Forterre friction law, including
the SHALTOP curvature correction of the normal pressure.
</p>

# Code

Import of Basilisk and C libraries.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grid/cartesian1D.h"
#include "saint-venant.h"

/**
## Physical, numerical and geometric parameters
*/

// Physical reference scales and gravitational acceleration.

const double physical_length_scale = 1.;
double physical_time_scale = 1.;
const double physical_gravity_acceleration = 9.81;

// Physical experimental geometry: dry 35-degree case, 4 mm beads.

const double slope_angle = 35.*pi/180.;
const double physical_gate_position = 0.;
const double physical_gate_bottom_height = 0.15;
const double physical_initial_pile_length = 0.131;
const double physical_initial_pile_height = 0.103;
const double physical_smoothing_width = 0.10;

// Same geometry, then nondimensionalized in main().

double gate_position;
double gate_bottom_height;
double initial_pile_length;
double initial_pile_height;
double smoothing_width;

// Parameters of the Pouliquen--Forterre law for 4 mm glass beads.

const double delta1 = 22.1*pi/180.;
const double delta2 = 31.8*pi/180.;
const double delta3 = 23.3*pi/180.;
const double beta_pf = 0.136;
const double xi_pf = 1e-3;
const double physical_grain_diameter = 4e-3;
double epsilon_pf = 1e-2;

// Numerical and output parameters.

const int default_grid_size = 1024;
const double physical_end_time = 3.;
const double physical_output_dt = 0.02;
const double physical_front_height_threshold = 1e-4;
const double physical_domain_start = -0.65;
const double physical_domain_length = 0.85;

double simulation_end_time;
double front_height_threshold;

// Internal nondimensional variables, computed in main().

double slope_gradient;
double break_position, smoothing_start, smoothing_end;
double grain_scale_length;
double friction_1, friction_2, friction_3;

// The physical SHALTOP topography is recorded here for clarity and output.
// The Saint-Venant bottom zb[] is set to zero so that the Cartesian source
// term of saint-venant.h is not added on top of the SHALTOP source terms.

scalar physical_bottom[];

// Auxiliary functions defined later in the resources section.

static void fetch_reference_files (void);
static void write_stderr_profile (double physical_time);
static void write_comparison_profile (const char * filename);
static void write_animation_profile (int frame_index);

/**
## Smoothed topography and SHALTOP geometric quantities

<p align="justify">
The physical bed $z_b(x)$ is horizontal downstream and then connected to a
$35^\circ$ slope by a quadratic transition. For SHALTOP, the code also needs
its first and second derivatives:
</p>

$$
\displaystyle
\partial_x z_b,
\qquad
\partial_{xx} z_b.
$$

<p align="justify">
The local geometric factors are
</p>

$$
\displaystyle
c
=
\frac{1}{\sqrt{1+\left(\partial_x z_b\right)^2}},
\qquad
c_x
=
-\left(\partial_x z_b\right)
 \left(\partial_{xx}z_b\right)c^3.
$$

<p align="justify">
Here $z_b(x)$ always denotes the physical topography. In the C code,
`physical_bottom[]` stores this physical bed, and
`shaltop_bottom_geometry()` evaluates the same geometry for the explicit
SHALTOP source terms. By contrast, `zb[]` is the internal bottom field of
`saint-venant.h` and is deliberately kept equal to zero. Assigning the
physical topography to `zb[]` would add the Cartesian topographic source a
second time.
</p>
*/

// This function defines z_b(x), its slope and its curvature using a smoothed transition.

static void shaltop_bottom_geometry (double xp,
                                     double * bottom,
                                     double * bottom_slope,
                                     double * bottom_curvature)
{
  if (xp <= smoothing_start) {                                         // Horizontal part of the bed.
    *bottom = 0.;
    *bottom_slope = 0.;
    *bottom_curvature = 0.;
    return;
  }

  if (xp >= smoothing_end) {                                           // Inclined part of the bed.
    *bottom = slope_gradient*(xp - break_position);
    *bottom_slope = slope_gradient;
    *bottom_curvature = 0.;
    return;
  }

  double s = (xp - smoothing_start)/smoothing_width;                   // Relative position inside the smoothed transition.

  *bottom = slope_gradient*smoothing_width*sq(s)/2.;                   // Parabolic C1 transition.
  *bottom_slope = slope_gradient*s;                                    // Linear slope inside the transition.
  *bottom_curvature = slope_gradient/smoothing_width;                  // Constant curvature inside the transition.
}

// This function returns c = cos(theta_local).

static double local_cosine (double bottom_slope)
{
  return 1./sqrt(1. + sq(bottom_slope));
}

// This function returns c_x from the first and second derivatives of z_b.

static double local_cosine_slope (double bottom_slope,
                                  double bottom_curvature,
                                  double bottom_cosine)
{
  return -bottom_slope*bottom_curvature*bottom_cosine*bottom_cosine*bottom_cosine;
}

/**
## Pouliquen--Forterre friction law

<p align="justify">
The same granular Froude-number definition derived from the Pouliquen
experiments is used here, evaluated with the tangential velocity $u_t$ and the
normal thickness $h_n$. The local Froude number is
</p>

$$
\displaystyle
Fr
=
\frac{|u_t|}{\sqrt{G h_n c}},
\qquad
c=\cos\theta_{\rm loc}.
$$

The law contains three regimes:

$$
\mu(h_n,Fr)=
\begin{cases}
\displaystyle
\mu_1+
\dfrac{\mu_2-\mu_1}{1+\beta h_n/(LFr)},
& Fr\geq\beta,\\[1.2em]
\displaystyle
\mu_3+
\dfrac{\mu_2-\mu_1}{1+h_n/L}
+
\left[
\mu_1+
\dfrac{\mu_2-\mu_1}{1+h_n/L}
-\mu_3-
\dfrac{\mu_2-\mu_1}{1+h_n/L}
\right]
\left(\dfrac{Fr}{\beta}\right)^\xi,
& 0<Fr<\beta,\\[1.2em]
\displaystyle
\min\left(
\mu_3+
\dfrac{\mu_2-\mu_1}{1+h_n/L},
\left|\partial_x\left(z_b+c h_n\right)\right|
\right),
& Fr=0.
\end{cases}
$$

<p align="justify">
Near $Fr=0$, a small numerical regularization is added so that the transition
between the static state and the intermediate regime is continuous.
</p>
*/

// This function regularizes the transition between the static and intermediate regimes.

static double smoothstep (double s)
{
  s = max(0., min(1., s));                                             // Clamp s inside the transition interval.

  return sq(s)*(3. - 2.*s);                                            // Smooth interpolation between the two regimes.
}

// This function defines the Pouliquen--Forterre law using the normal thickness.

static double mu_pouliquen_forterre (double Fr, double normal_height,
                                     double free_surface_slope)
{
  double height_over_grain = normal_height/grain_scale_length;          // Ratio between local normal thickness and grain length scale.

  double mu_stop =                                                     // Stopping friction threshold.
    friction_1 + (friction_2 - friction_1)/(1. + height_over_grain);

  double mu_start =                                                    // Starting friction threshold.
    friction_3 + (friction_2 - friction_1)/(1. + height_over_grain);

  double mu_static = min(mu_start, free_surface_slope);                 // Static branch of the law: Fr = 0.

  if (Fr <= 0.)
    return mu_static;

  double mu_intermediate =                                             // Intermediate branch of the law: 0 < Fr < beta.
    mu_start + (mu_stop - mu_start)*pow(Fr/beta_pf, xi_pf);

  if (Fr < epsilon_pf)                                                 // Transition near Fr = 0 to avoid a numerical discontinuity.
    return mu_static + smoothstep(Fr/epsilon_pf)*(mu_intermediate - mu_static);

  if (Fr < beta_pf)
    return mu_intermediate;

  // Dynamic branch: Fr >= beta.
  return friction_1 + (friction_2 - friction_1)/(1. + beta_pf*normal_height/(grain_scale_length*Fr));
}

/**
## Main program and nondimensionalization

<p align="justify">
The main program initializes the nondimensional quantities used by the solver.
Choosing the length scale $L_{\mathrm{ref}}=1~\mathrm{m}$ keeps the same
numerical values for lengths expressed in meters. The time scale is chosen so
that the nondimensional gravity is $G=1$.
</p>

We define

$$
\displaystyle
x=L_{\mathrm{ref}}x^*,\qquad
h=L_{\mathrm{ref}}h^*,\qquad
t=T_{\mathrm{ref}}t^*,\qquad
u=\frac{L_{\mathrm{ref}}}{T_{\mathrm{ref}}}\nu^*,
\qquad
T_{\mathrm{ref}}=\sqrt{\frac{L_{\mathrm{ref}}}{g}}.
$$

<p align="justify">
The optional first command-line argument can be used to change the grid size,
and the optional second one can be used to change the very small
regularization interval $\varepsilon_{\rm PF}$.
</p>
*/

int main (int argc, char * argv[])
{
  // Reference scales.

  physical_time_scale = sqrt(physical_length_scale/physical_gravity_acceleration);
  G = physical_gravity_acceleration*sq(physical_time_scale)/physical_length_scale;

  slope_gradient = tan(slope_angle);                                   // Bed slope written as a slope coefficient.

  // Conversion of physical geometric parameters into nondimensional variables.

  gate_position = physical_gate_position/physical_length_scale;
  gate_bottom_height = physical_gate_bottom_height/physical_length_scale;
  initial_pile_length = physical_initial_pile_length/physical_length_scale;
  initial_pile_height = physical_initial_pile_height/physical_length_scale;
  smoothing_width = physical_smoothing_width/physical_length_scale;

  // Compute the break position and the bounds of the transition.

  break_position = gate_position - gate_bottom_height/slope_gradient;
  smoothing_start = break_position - smoothing_width/2.;
  smoothing_end = break_position + smoothing_width/2.;

  // Nondimensional parameters of the basal friction law.

  grain_scale_length = 1.3*physical_grain_diameter/physical_length_scale;
  friction_1 = tan(delta1);
  friction_2 = tan(delta2);
  friction_3 = tan(delta3);

  // Convert the thresholds and times used by the output events.

  front_height_threshold = physical_front_height_threshold/physical_length_scale;
  simulation_end_time = physical_end_time/physical_time_scale;

  // Physical domain, grid size and maximum time step.

  X0 = physical_domain_start/physical_length_scale;
  L0 = physical_domain_length/physical_length_scale;
  N = argc > 1 ? atoi(argv[1]) : default_grid_size;

  if (argc > 2)
    epsilon_pf = atof(argv[2]);

  if (epsilon_pf <= 0. || epsilon_pf >= beta_pf) {
    fprintf(stderr,
            "ERROR: epsilon_pf must satisfy 0 < epsilon_pf < beta_pf = %g.\n",
            beta_pf);
    return 1;
  }

  DT = 1e-4;

  fetch_reference_files();
  run();
}

/**
## Geometry and initial condition

<p align="justify">
In this SHALTOP-type formulation, `h[]` stores
$q=h_n/c$, not the vertical height. The initial physical free surface is still
horizontal and equal to
</p>

$$
\displaystyle
z_s(0)=z_{\rm gate}+H_{\rm pile}.
$$

<p align="justify">
Since $z_s=z_b+c^2q$, the initial value imposed in the pile is
</p>

$$
\displaystyle
q=\frac{z_s-z_b}{c^2}.
$$
*/

event init (i = 0)
{
  foreach() {
    double bottom, bottom_slope, bottom_curvature;
    shaltop_bottom_geometry(x, &bottom, &bottom_slope, &bottom_curvature);

    double bottom_cosine = local_cosine(bottom_slope);

    physical_bottom[] = bottom;                                        // Physical topography z_b stored separately from Basilisk's zb[].
    zb[] = 0.;                                                         // Disable the Cartesian topographic source term of saint-venant.h.

    double inside_initial_pile =                                       // Cell located inside the initial pile footprint.
      (x >= gate_position && x <= gate_position + initial_pile_length);
    double initial_surface_level =                                     // Elevation of the initial horizontal free surface.
      gate_bottom_height + initial_pile_height;

    h[] = inside_initial_pile ?                                        // Conserved variable q = h_n/c.
      max((initial_surface_level - bottom)/sq(bottom_cosine), 0.) : 0.;

    u.x[] = 0.;                                                        // v = c*u_t: pile initially at rest.
  }

  boundary ({zb, physical_bottom, h, u});                              // Update boundary conditions.
}

/**
## SHALTOP geometric correction

<p align="justify">
With `zb[] = 0`, the internal Cartesian topographic source is disabled and the
Saint-Venant solver advances
</p>

$$
\displaystyle
\partial_t v+v\partial_xv=-G\partial_x q.
$$

<p align="justify">
The SHALTOP-type pressure and geometric terms require instead
</p>

$$
\displaystyle
\partial_t v+v\partial_xv=
-Gc^2\partial_x(z_b+c^2q)+c c_x u_t^2.
$$

<p align="justify">
The event below therefore adds only the missing correction
</p>

$$
\displaystyle
G\partial_xq-Gc^2\partial_x(z_b+c^2q)+c c_x u_t^2.
$$
*/

event shaltop_geometry (i++)
{
  foreach() {
    if (h[] <= dry) {
      u.x[] = 0.;
      continue;
    }

    double bottom, bottom_slope, bottom_curvature;
    shaltop_bottom_geometry(x, &bottom, &bottom_slope, &bottom_curvature);

    double bottom_cosine = local_cosine(bottom_slope);
    double bottom_cosine_slope =
      local_cosine_slope(bottom_slope, bottom_curvature, bottom_cosine);

    // Physical surface elevation E = z_b + c^2 q, evaluated with neighboring local cosines.

    double bottom_minus, slope_minus, curvature_minus;
    double bottom_plus, slope_plus, curvature_plus;
    shaltop_bottom_geometry(x - Delta, &bottom_minus, &slope_minus, &curvature_minus);
    shaltop_bottom_geometry(x + Delta, &bottom_plus, &slope_plus, &curvature_plus);

    double cosine_minus = local_cosine(slope_minus);
    double cosine_plus = local_cosine(slope_plus);

    double surface_minus = bottom_minus + sq(cosine_minus)*h[-1];
    double surface_plus = bottom_plus + sq(cosine_plus)*h[1];

    double conserved_height_slope = (h[1] - h[-1])/(2.*Delta);         // d_x q.
    double surface_slope = (surface_plus - surface_minus)/(2.*Delta);  // d_x(z_b + c^2 q).

    double tangential_velocity = u.x[]/(bottom_cosine + 1e-30);         // u_t = v/c.

    double correction =
        G*conserved_height_slope
      - G*sq(bottom_cosine)*surface_slope
      + bottom_cosine*bottom_cosine_slope*sq(tangential_velocity);

    u.x[] += dt*correction;
  }

  boundary ({u});
}

/**
## Basal friction

<p align="justify">
For SHALTOP, the basal friction is applied on the tangential motion and
contains a curvature correction of the normal pressure:
</p>

$$
\displaystyle
\partial_t v
=-G\mu(h_n,Fr)c^2\,\mathrm{sign}(v)
\left[
1+
\frac{(\partial_{xx}z_b)u_t^2}{Gc}
\right]_+.
$$

<p align="justify">
The factor in brackets is the positive part of the curvature-modified normal
pressure. On a straight slope it is equal to one.
</p>
*/

event friction (i++)
{
  foreach() {
    if (h[] <= dry) {                                                  // Dry-cell treatment.
      u.x[] = 0.;
      continue;
    }

    double bottom, bottom_slope, bottom_curvature;
    shaltop_bottom_geometry(x, &bottom, &bottom_slope, &bottom_curvature);

    double bottom_cosine = local_cosine(bottom_slope);

    double normal_height = bottom_cosine*h[];                          // h_n = c q.
    double tangential_velocity = u.x[]/(bottom_cosine + 1e-30);         // u_t = v/c.
    double tangential_speed = fabs(tangential_velocity);

    double froude_number =
      tangential_speed/(sqrt(G*normal_height*bottom_cosine) + 1e-30);   // Fr = |u_t|/sqrt(G h_n c).

    // Static branch: slope of the physical free surface z_b + c^2 q.

    double bottom_minus, slope_minus, curvature_minus;
    double bottom_plus, slope_plus, curvature_plus;
    shaltop_bottom_geometry(x - Delta, &bottom_minus, &slope_minus, &curvature_minus);
    shaltop_bottom_geometry(x + Delta, &bottom_plus, &slope_plus, &curvature_plus);

    double cosine_minus = local_cosine(slope_minus);
    double cosine_plus = local_cosine(slope_plus);

    double surface_minus = bottom_minus + sq(cosine_minus)*h[-1];
    double surface_plus = bottom_plus + sq(cosine_plus)*h[1];
    double free_surface_slope = fabs((surface_plus - surface_minus)/(2.*Delta));

    double effective_friction =
      mu_pouliquen_forterre(froude_number, normal_height, free_surface_slope);

    double normal_pressure_factor = max(0.,
      1. + bottom_curvature*sq(tangential_velocity)/(G*bottom_cosine + 1e-30));

    double horizontal_speed = fabs(u.x[]);
    double horizontal_deceleration =
      G*effective_friction*sq(bottom_cosine)*normal_pressure_factor;

    double reduced_speed = max(horizontal_speed - horizontal_deceleration*dt, 0.);

    u.x[] = horizontal_speed > 1e-30 ?
      reduced_speed*u.x[]/horizontal_speed : 0.;                       // Reduce the magnitude without changing the sign.
  }

  boundary ({u});
}

/**
## Plots: profiles, front position and animation

<p align="justify">
Profiles are saved at $t=0$, $0.2$, $0.4$ and $3$ s for the final comparison.
The front position is tracked in parallel. The profiles used for the animation
are saved every `physical_output_dt` seconds and then read back by Gnuplot in
the Main results section.
</p>

<p align="justify">
For compatibility with the previous Python scripts, the profiles printed to
`stderr` keep the same five-column format:
</p>

```
x  h_vertical  physical_time  physical_bottom  tangential_velocity
```
*/

// Initial profile used for comparison with the figure from Poulain et al.
event profile_t0 (t = 0.)
{
  write_comparison_profile ("profile_dry35_t0.dat");
  write_stderr_profile (t*physical_time_scale);
}

// Profile at t = 0.2 s.
event profile_t02 (t = 0.2/physical_time_scale)
{
  write_comparison_profile ("profile_dry35_t02.dat");
  write_stderr_profile (t*physical_time_scale);
}

// Profile at t = 0.4 s.
event profile_t04 (t = 0.4/physical_time_scale)
{
  write_comparison_profile ("profile_dry35_t04.dat");
  write_stderr_profile (t*physical_time_scale);
}

// Final profile at t = 3 s.
event profile_t30 (t = 3./physical_time_scale)
{
  write_comparison_profile ("profile_dry35_t30.dat");
  write_stderr_profile (t*physical_time_scale);
}

// Track the front position as the distance travelled from the gate.
event front (t = 0.; t <= simulation_end_time; t += physical_output_dt/physical_time_scale)
{
  double front_position = gate_position;

  foreach (reduction(min:front_position)) {
    double bottom, bottom_slope, bottom_curvature;
    shaltop_bottom_geometry(x, &bottom, &bottom_slope, &bottom_curvature);
    double bottom_cosine = local_cosine(bottom_slope);
    double vertical_height = sq(bottom_cosine)*h[];

    if (vertical_height > front_height_threshold)
      front_position = min(front_position, x);                         // Farthest front toward decreasing x.
  }

  fprintf (stdout, "%.12g %.12g\n",
           t*physical_time_scale, (gate_position - front_position)*physical_length_scale);
}

// Save the profiles that will be read back by the GIF Gnuplot block.
event save_animation_profiles (t = 0.; t <= simulation_end_time; t += physical_output_dt/physical_time_scale)
{
  int frame_index = (int) lround((t*physical_time_scale)/physical_output_dt);
  write_animation_profile (frame_index);
}

event stop (t = simulation_end_time)
{
  return 1;
}

/**
## Comparison resources and saved profiles

The four PNG backgrounds containing the reference SHALTOP profiles from the
article are downloaded at startup from the GitHub repository. The CSV file
containing the experimental profiles is also retrieved automatically. The
first three experimental profiles are extracted into separate `.dat` files so
that Gnuplot can overlay them on the final comparison.
*/

// Files used for the results.

const char * github_resources_url =
  "https://raw.githubusercontent.com/kilian-gthub/basilisk-ressources/main/ressources";

const char * comparison_background_files[] = {
  "t0_red_plus_black.png",
  "t0_2_red_plus_black.png",
  "t0_4_red_plus_black.png",
  "t3_red_plus_black.png"
};

const char * experimental_csv_url =
  "https://basilisk.fr/sandbox/Kilian/ressource/data_article_poulain.csv";

const char * experimental_csv_file = "data_article_poulain.csv";

const char * experimental_profile_files[] = {
  "exp_profile_t0.dat",
  "exp_profile_t02.dat",
  "exp_profile_t04.dat"
};

static void extract_experimental_profiles (const char * csv_filename);

// Download the comparison backgrounds and the experimental CSV file.
static void fetch_reference_files (void)
{
  for (int k = 0; k < 4; k++) {
    char command[2048];

    // The PNG backgrounds are fetched from GitHub. curl is preferred; wget is used as a fallback.
    snprintf (command, sizeof(command),
              "rm -f '%s'; "
              "(curl -fsSL --retry 2 --connect-timeout 15 -o '%s' '%s/%s' "
              "|| wget -q -O '%s' '%s/%s') "
              ">> reference_files.log 2>&1",
              comparison_background_files[k],
              comparison_background_files[k], github_resources_url, comparison_background_files[k],
              comparison_background_files[k], github_resources_url, comparison_background_files[k]);

    if (system(command) != 0)
      fprintf (stderr, "WARNING: unable to download %s. See reference_files.log.\n",
               comparison_background_files[k]);
  }

  char csv_command[2048];
  // The experimental CSV is downloaded and then converted into three Gnuplot profiles.
  snprintf (csv_command, sizeof(csv_command),
            "rm -f '%s'; "
            "(curl -fsSL --retry 2 --connect-timeout 15 -o '%s' '%s' "
            "|| wget -q -O '%s' '%s') "
            ">> reference_files.log 2>&1",
            experimental_csv_file,
            experimental_csv_file, experimental_csv_url,
            experimental_csv_file, experimental_csv_url);

  if (system(csv_command) != 0)
    fprintf (stderr, "WARNING: unable to download %s. See reference_files.log.\n",
             experimental_csv_file);
  else
    extract_experimental_profiles(experimental_csv_file);
}

// Remove newline characters at the end of a string.
static void trim_newline (char * s)
{
  if (!s) return;
  size_t n = strlen(s);
  while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r'))
    s[--n] = '\0';
}

// Read a CSV line of the form label,value1,value2,...

static int parse_csv_numbers (char * line, double * values, int max_values)
{
  int count = 0;
  char * token = strtok(line, ",");
  int column = 0;

  while (token) {
    if (column > 0 && count < max_values) {
      while (*token == ' ' || *token == '\t')
        token++;
      if (*token != '\0')
        values[count++] = atof(token);
    }
    token = strtok(NULL, ",");
    column++;
  }
  return count;
}

// Write an experimental profile in the figure coordinates, with the sign of x reversed.

static void write_experimental_profile (const char * filename,
                                        double * physical_x_values, int nx,
                                        double * physical_z_values, int nz)
{
  FILE * fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "WARNING: unable to write %s.\n", filename);
    return;
  }

  int n = nx < nz ? nx : nz;
  for (int i = 0; i < n; i++)
    fprintf(fp, "%.12g %.12g\n", -physical_x_values[i], physical_z_values[i]);

  fclose(fp);
}

// Extract the first three experimental times from the Poulain et al. CSV file.

static void extract_experimental_profiles (const char * csv_filename)
{
  FILE * fp = fopen(csv_filename, "r");
  if (!fp) {
    fprintf(stderr, "WARNING: unable to read %s.\n", csv_filename);
    return;
  }

  char line[65536];
  double physical_x_t1[4096], physical_y_t1[4096];
  double physical_x_t2[4096], physical_y_t2[4096];
  double physical_x_t3[4096], physical_y_t3[4096];
  int nx1 = 0, ny1 = 0, nx2 = 0, ny2 = 0, nx3 = 0, ny3 = 0;

  if (!fgets(line, sizeof(line), fp)) {
    fclose(fp);
    return;
  }

  while (fgets(line, sizeof(line), fp)) {
    trim_newline(line);
    char raw_line[65536];
    strcpy(raw_line, line);

    char * label = strtok(raw_line, ",");
    if (!label)
      continue;

    if (!strcmp(label, "x_t1"))
      nx1 = parse_csv_numbers(line, physical_x_t1, 4096);
    else if (!strcmp(label, "y_t1"))
      ny1 = parse_csv_numbers(line, physical_y_t1, 4096);
    else if (!strcmp(label, "x_t2"))
      nx2 = parse_csv_numbers(line, physical_x_t2, 4096);
    else if (!strcmp(label, "y_t2"))
      ny2 = parse_csv_numbers(line, physical_y_t2, 4096);
    else if (!strcmp(label, "x_t3"))
      nx3 = parse_csv_numbers(line, physical_x_t3, 4096);
    else if (!strcmp(label, "y_t3"))
      ny3 = parse_csv_numbers(line, physical_y_t3, 4096);
  }

  fclose(fp);

  write_experimental_profile(experimental_profile_files[0], physical_x_t1, nx1, physical_y_t1, ny1);
  write_experimental_profile(experimental_profile_files[1], physical_x_t2, nx2, physical_y_t2, ny2);
  write_experimental_profile(experimental_profile_files[2], physical_x_t3, nx3, physical_y_t3, ny3);
}

// Write the standard five-column profile to stderr for external scripts.
static void write_stderr_profile (double physical_time)
{
  foreach() {
    double bottom, bottom_slope, bottom_curvature;
    shaltop_bottom_geometry(x, &bottom, &bottom_slope, &bottom_curvature);
    double bottom_cosine = local_cosine(bottom_slope);
    double vertical_height = sq(bottom_cosine)*h[];
    double tangential_velocity = u.x[]/(bottom_cosine + 1e-30);

    fprintf (stderr, "%.12g %.12g %.12g %.12g %.12g\n",
             x*physical_length_scale, vertical_height*physical_length_scale, physical_time,
             bottom*physical_length_scale, tangential_velocity*physical_length_scale/physical_time_scale);
  }
}

// Columns: x as displayed in the article, z_b, z_b+h_vertical, h_vertical and tangential velocity.
static void write_comparison_profile (const char * filename)
{
  FILE * fp = fopen (filename, "w");

  if (!fp) {
    fprintf (stderr, "WARNING: unable to write %s.\n", filename);
    return;
  }

  foreach() {
    double bottom, bottom_slope, bottom_curvature;
    shaltop_bottom_geometry(x, &bottom, &bottom_slope, &bottom_curvature);
    double bottom_cosine = local_cosine(bottom_slope);
    double vertical_height = sq(bottom_cosine)*h[];
    double tangential_velocity = u.x[]/(bottom_cosine + 1e-30);

    fprintf (fp, "%.12g %.12g %.12g %.12g %.12g\n",
             -x*physical_length_scale, bottom*physical_length_scale,
             (bottom + vertical_height)*physical_length_scale,
             vertical_height*physical_length_scale,
             tangential_velocity*physical_length_scale/physical_time_scale);
  }

  fclose (fp);
}

// One file is written per animation frame. It is kept so that
// the ~~~gnuplot block in the Basilisk documentation can then generate
// and display the GIF directly on the sandbox page.
static void write_animation_profile (int frame_index)
{
  char filename[128];
  snprintf (filename, sizeof(filename), "frame_dry35_shaltop_%04d.dat", frame_index);

  FILE * fp = fopen (filename, "w");
  if (!fp) {
    fprintf (stderr, "WARNING: unable to write %s.\n", filename);
    return;
  }

  foreach() {
    double bottom, bottom_slope, bottom_curvature;
    shaltop_bottom_geometry(x, &bottom, &bottom_slope, &bottom_curvature);
    double bottom_cosine = local_cosine(bottom_slope);
    double vertical_height = sq(bottom_cosine)*h[];

    fprintf (fp, "%.12g %.12g %.12g %.12g\n",
             -x*physical_length_scale, bottom*physical_length_scale,
             (bottom + vertical_height)*physical_length_scale,
             vertical_height*physical_length_scale);
  }

  fclose (fp);
}


/**
# Main results

## Collapse animation

This animation shows the granular pile starting to **move on the slope**.
The topography is plotted in black, while the upper contour of the pile is
plotted in sand yellow.

~~~gnuplot Figure 3 - Time animation of the dry granular pile with the SHALTOP-type correction
set terminal gif animate delay 5 loop 0 size 850,672 enhanced font "Liberation Sans,12"
set output "animate_dry35_muv_shaltop.gif"

set xlabel "x (m)"
set ylabel "z (m)"
set xrange [0.65:-0.20]
set yrange [-0.05:0.38]
set grid
set key top left

slope_fill = "#e2d1c4"
pile_fill  = "#fff1b8"
sand_line  = "#d8a72e"

# t = 0, 0.02, ..., 3.00 s: 151 profiles saved by the C code.
do for [k=0:150] {
    time = 0.02*k
    datafile = sprintf("frame_dry35_shaltop_%04d.dat", k)

    set title sprintf("SHALTOP-type Basilisk model - t = %.2f s", time)
    plot datafile using 1:(-0.05):2 with filledcurves lc rgb slope_fill title "slope", \
         datafile using 1:2:3 with filledcurves lc rgb pile_fill title "granular pile", \
         datafile using 1:2 with lines lw 2.2 lc rgb "#000000" title "topography", \
         datafile using 1:($4 > 1.e-4 ? $3 : 1/0) with lines lw 2.4 lc rgb sand_line title "pile outline"
}

unset output
~~~

## Comparison with the experiment and reference SHALTOP profiles

Figure 4 compares the Basilisk SHALTOP-type computation with the **reference
SHALTOP profiles extracted** from Poulain *et al.*. The Basilisk topography is
plotted in black, the SHALTOP-type Basilisk pile outline in sand yellow, the
experimental profile in blue and the reference SHALTOP result from the article
in dashed red.

~~~gnuplot Figure 4 - Superposition of Basilisk SHALTOP-type, experimental and reference SHALTOP profiles
reset session
set terminal pngcairo size 850,672 enhanced font "Liberation Sans,10"
set output "comparison_dry35_shaltop.png"

unset key
unset border
unset tics
unset title
unset xlabel
unset ylabel
unset grid
unset colorbox
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
set style line 1 lc rgb "#000000" lw 2.2
set style line 2 lc rgb "#d8a72e" lw 2.6 dt 2
set style line 3 lc rgb "#40549b" lw 2.4
set style line 4 lc rgb "#d62728" lw 2.4 dt 2

# Global margins: the four panels are slightly reduced to place
# the axis labels and the legend outside the article images.
set object 100 rect from screen 0.075,0.905 to screen 0.900,0.985 fc rgb "white" fillstyle solid 0.85 noborder front
set label 100 "x (m)" at screen 0.525,0.035 center front font "Liberation Sans,12"
set label 101 "z (m)" at screen 0.030,0.490 center rotate by 90 front font "Liberation Sans,12"

set arrow 201 from screen 0.105,0.965 to screen 0.155,0.965 nohead front lw 2.2 lc rgb "#000000"
set label 201 "Basilisk topography" at screen 0.160,0.955 left front font "Liberation Sans,9"
set arrow 202 from screen 0.420,0.965 to screen 0.470,0.965 nohead front lw 2.6 dt 2 lc rgb "#d8a72e"
set label 202 "Basilisk SHALTOP-type" at screen 0.475,0.955 left front font "Liberation Sans,9"
set arrow 203 from screen 0.105,0.935 to screen 0.155,0.935 nohead front lw 2.4 lc rgb "#40549b"
set label 203 "Experimental profile" at screen 0.160,0.925 left front font "Liberation Sans,9"
set arrow 204 from screen 0.420,0.935 to screen 0.470,0.935 nohead front lw 2.4 dt 2 lc rgb "#d62728"
set label 204 "SHALTOP (Poulain et al.)" at screen 0.475,0.925 left front font "Liberation Sans,9"

set multiplot
set size 0.455,0.410

# t = 0 s: top-left panel
set origin 0.070,0.490
set xrange [0:1545]
set yrange [0:1221]
plot "t0_red_plus_black.png" binary filetype=png with rgbimage notitle, \
     "profile_dry35_t0.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(2.5 + ($2 - 0.0)*(1053.5 - 2.5)/(0.3 - 0.0)) with lines ls 1 notitle, \
     "profile_dry35_t0.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):($4 > 1.e-4 ? (2.5 + ($3 - 0.0)*(1053.5 - 2.5)/(0.3 - 0.0)) : 1/0) with lines ls 2 notitle, \
     "exp_profile_t0.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(2.5 + ($2 - 0.0)*(1053.5 - 2.5)/(0.3 - 0.0)) with lines ls 3 notitle

# t = 0.2 s: top-right panel
set origin 0.525,0.490
set xrange [0:1545]
set yrange [0:1219]
plot "t0_2_red_plus_black.png" binary filetype=png with rgbimage notitle, \
     "profile_dry35_t02.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(0.5 + ($2 - 0.0)*(1051.5 - 0.5)/(0.3 - 0.0)) with lines ls 1 notitle, \
     "profile_dry35_t02.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):($4 > 1.e-4 ? (0.5 + ($3 - 0.0)*(1051.5 - 0.5)/(0.3 - 0.0)) : 1/0) with lines ls 2 notitle, \
     "exp_profile_t02.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(0.5 + ($2 - 0.0)*(1051.5 - 0.5)/(0.3 - 0.0)) with lines ls 3 notitle

# t = 0.4 s: bottom-left panel
set origin 0.070,0.080
set xrange [0:1543]
set yrange [0:1221]
plot "t0_4_red_plus_black.png" binary filetype=png with rgbimage notitle, \
     "profile_dry35_t04.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(0.5 + ($2 - 0.0)*(1053.36270491803 - 0.5)/(0.3 - 0.0)) with lines ls 1 notitle, \
     "profile_dry35_t04.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):($4 > 1.e-4 ? (0.5 + ($3 - 0.0)*(1053.36270491803 - 0.5)/(0.3 - 0.0)) : 1/0) with lines ls 2 notitle, \
     "exp_profile_t04.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(0.5 + ($2 - 0.0)*(1053.36270491803 - 0.5)/(0.3 - 0.0)) with lines ls 3 notitle

# t = 3 s: bottom-right panel
set origin 0.525,0.080
set xrange [0:1543]
set yrange [0:1217]
plot "t3_red_plus_black.png" binary filetype=png with rgbimage notitle, \
     "profile_dry35_t30.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(0.5 + ($2 - 0.0)*(1049.90983606557 - 0.5)/(0.3 - 0.0)) with lines ls 1 notitle, \
     "profile_dry35_t30.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):($4 > 1.e-4 ? (0.5 + ($3 - 0.0)*(1049.90983606557 - 0.5)/(0.3 - 0.0)) : 1/0) with lines ls 2 notitle

unset multiplot
unset output
~~~

*/

/**
# References

* POLGE POULICHET, K. (2026). *Numerical modelling of landslides under different mechanical conditions*. First-year Master's internship report, Sorbonne Université, Institut Jean le Rond d'Alembert. [Full report (PDF)](https://raw.githubusercontent.com/kilian-gthub/basilisk-ressources/main/ressources/Rapport_POLGE_POULICHET_Kilian.pdf)
* [P. Poulain *et al.*, *Performance and limits of a shallow-water model for landslide-generated tsunamis: from laboratory experiments to simulations of flank collapses at Montagne Pelée (Martinique)*, Geophysical Journal International, 233, 796--825, 2023.](https://www.ipgp.fr/~mangeney/poulain-etal_gji-2023.pdf)
* [O. Pouliquen & Y. Forterre, *Friction law for dense granular flows: application to the motion of a mass down a rough inclined plane*, Journal of Fluid Mechanics, 453, 133--151, 2002.](https://yoelforterre.wordpress.com/wp-content/uploads/2016/09/jfmmoire02.pdf)
* [M. Peruzzetto *et al.*, *Topography curvature effects in thin-layer models for gravity-driven flows without bed erosion*, Journal of Geophysical Research: Earth Surface, 126(4), e2020JF005657, 2021.](https://doi.org/10.1029/2020JF005657)
*/
