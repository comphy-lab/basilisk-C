/**
# Dry granular flow on a slope

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
diameter $d=4$ mm.This code is also designed to compare the Basilisk results with the
**non-hydrostatic HySEA code solution**, using the
same initial configuration and the same output times.
</p>

<p align="justify">
The solver is Basilisk's Saint-Venant equations, to which a **basal friction
law** from Pouliquen & Forterre (2002) is added.
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

## Solved equations

<p align="justify">
The unknowns are the vertical thickness $h(x,t)$ and the depth-averaged
horizontal velocity $u(x,t)$. The **Saint-Venant** system with topography
$z_b(x)$ reads
</p>

$$
\displaystyle
\frac{\partial h}{\partial t}
+ \frac{\partial (hu)}{\partial x} = 0,
$$

$$
\displaystyle
\frac{\partial (hu)}{\partial t}
+ \frac{\partial}{\partial x}
(hu^2 + \frac{g h^2}{2})
= -gh\frac{\partial z_b}{\partial x}
  -gh\mu(h,Fr)\frac{u}{|u|}.
$$

# Code

Import of Basilisk and C libraries
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

// Parameters of the Pouliquen--Forterre law.

const double delta1 = 22.1*pi/180.;
const double delta2 = 31.8*pi/180.;
const double delta3 = 23.3*pi/180.;
const double beta_pf = 0.136;
const double xi_pf = 1e-3;
const double physical_grain_diameter = 4e-3;
double epsilon_pf = 1e-2;

// Numerical and output parameters.

const int grid_size = 1024;
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

// Auxiliary functions defined later in the resources section.
static void fetch_reference_files (void);
static void write_comparison_profile (const char * filename);
static void write_animation_profile (int frame_index);

/**
## Smoothed topography

If $x_{break}$ is the theoretical break position and $\ell_t$ the width
of the transition, the topography is:

$$
z_b(x)=
\begin{cases}
\displaystyle
0, & x\leq x_b-\ell_t/2,\\
\tan\theta\,\ell_t\,s^2/2,
& x_b-\ell_t/2 < x < x_b+\ell_t/2,\\
\tan\theta\,(x-x_b), & x\geq x_b+\ell_t/2,
\end{cases}
\qquad
\displaystyle
s=\frac{x-(x_b-\ell_t/2)}{\ell_t}.
$$

The slope is therefore zero at the beginning of the transition and equal to
$\tan\theta$ at its end.
*/

// This function defines the bed topography using a smoothed transition.

static double smooth_bottom (double xp)
{
  if (xp <= smoothing_start)                                            // Horizontal part of the bed: z_b = 0.
    return 0.;

  if (xp >= smoothing_end)                                              // Inclined part of the bed, after the transition.
    return slope_gradient*(xp - break_position);

  double s = (xp - smoothing_start)/smoothing_width;                    // Relative position inside the smoothed transition, between 0 and 1.

  return slope_gradient*smoothing_width*sq(s)/2.;                       // Parabolic transition between the horizontal bed and the slope.
}

// This function defines the local bed slope.

static double smooth_bottom_slope (double xp)
{
  if (xp <= smoothing_start)
    return 0.;

  if (xp >= smoothing_end)
    return slope_gradient;                                              // Constant slope on the inclined part.

  double s = (xp - smoothing_start)/smoothing_width;

  return slope_gradient*s;                                              // Linear slope inside the transition.
}

/**
## Pouliquen--Forterre friction law

<p align="justify">
The quantities $\mu_i=\tan\delta_i$ are defined as friction coefficients depending on the slope angle $\theta$, $L=1.3d$ is a characteristic length of the order of the glass-bead diameter, and $\beta$ is a constant that also depends on the granular material. The empirical parameter $\xi$ is used to **regularize** the law between the static and dynamic regimes.
</p>

The law contains three regimes:

$$
\mu(h,Fr)=
\begin{cases}
\displaystyle
\mu_1 + \frac{\mu_2-\mu_1}{1+\beta h/(LFr)}, & Fr\geq\beta,\\[0.8em]
\displaystyle
\mu_3 + \frac{\mu_2-\mu_1}{1+h/L}
+ \left(\mu_1 + \frac{\mu_2-\mu_1}{1+h/L}-\mu_3 - \frac{\mu_2-\mu_1}{1+h/L}\right)
\left(\frac{Fr}{\beta}\right)^\xi,
& 0<Fr<\beta,\\[0.8em]
\displaystyle
\min\left(\mu_3 + \frac{\mu_2-\mu_1}{1+h/L},
\left|\partial_x(z_b+h)\right|\right), & Fr=0.
\end{cases}
$$

<p align="justify">
Near $Fr=0$, a small numerical regularization is added so that the transition between the static state and the intermediate regime is continuous. A polynomial is used that is zero at $s=0$, one at $s=1$, and has zero derivative at both ends. It provides a **smooth connection** from the static state to the intermediate branch in the interval $0<Fr<\varepsilon_{\rm PF}$.
</p>
*/

// This function regularizes the transition between the static and intermediate regimes.

static double smoothstep (double s)
{
  s = max(0., min(1., s));                                              // Clamp s inside the transition interval.

  return sq(s)*(3. - 2.*s);                                             // Smooth interpolation between the two regimes.
}

// This function defines the Pouliquen--Forterre law.

static double mu_pouliquen_forterre (double Fr, double local_height,
                                     double bottom_slope, double height_slope)
{
  double height_over_grain = local_height/grain_scale_length;           // Ratio between local height and grain length scale.

  double mu_stop =                                                      // Stopping friction threshold.
    friction_1 + (friction_2 - friction_1)/(1. + height_over_grain);

  double mu_start =                                                     // Starting friction threshold.
    friction_3 + (friction_2 - friction_1)/(1. + height_over_grain);

  double free_surface_slope = fabs(bottom_slope + height_slope);
  double mu_static = min(mu_start, free_surface_slope);                 // Static branch of the law: Fr = 0.

  if (Fr <= 0.)
    return mu_static;

  double mu_intermediate =                                              // Intermediate branch of the law: 0 < Fr < beta.
    mu_start + (mu_stop - mu_start)*pow(Fr/beta_pf, xi_pf);

  if (Fr < epsilon_pf)                                                  // Transition near Fr = 0 to avoid a numerical discontinuity.
    return mu_static + smoothstep(Fr/epsilon_pf)*(mu_intermediate - mu_static);

  if (Fr < beta_pf)
    return mu_intermediate;

  // Dynamic branch: Fr >= beta.
  return friction_1 + (friction_2 - friction_1)/(1. + beta_pf*local_height/(grain_scale_length*Fr));
}

/**
## Main program and nondimensionalization

<p align="justify">
The main program initializes the nondimensional quantities used by the solver. Choosing the length scale $L_{\mathrm{ref}}=1~\mathrm{m}$ makes it possible to keep the same numerical values for lengths expressed in meters, while passing them to Basilisk in **nondimensional form**. The time scale is chosen so that the nondimensional gravity is $G=1$, which simplifies the form of the solved equations.
</p>

We define

$$
\displaystyle
x=L_{\mathrm{ref}}x^*,\qquad
h=L_{\mathrm{ref}}h^*,\qquad
t=T_{\mathrm{ref}}t^*,\qquad
u=\frac{L_{\mathrm{ref}}}{T_{\mathrm{ref}}}u^*,
\qquad
T_{\mathrm{ref}}=\sqrt{\frac{L_{\mathrm{ref}}}{g}}.
$$

With this choice,

$$
\displaystyle
G=g\frac{T_{\mathrm{ref}}^2}{L_{\mathrm{ref}}}=1.
$$

<p align="justify">
The geometric parameters, the front-detection threshold and the final time are therefore converted before calling `run()`. The outputs are then converted back to physical units for positions, heights, times and velocities.
</p>
*/

int main (void)
{
  // Reference scales.

  physical_time_scale = sqrt(physical_length_scale/physical_gravity_acceleration);
  G = physical_gravity_acceleration*sq(physical_time_scale)/physical_length_scale;

  slope_gradient = tan(slope_angle);                                    // Bed slope written as a slope coefficient.

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

  // Taille du domaine physique, nombre de mailles et pas de temps.

  X0 = physical_domain_start/physical_length_scale;
  L0 = physical_domain_length/physical_length_scale;
  N = grid_size;
  DT = 1e-4;

  fetch_reference_files();
  run();
}

/**
## Geometry and initial condition

<p align="justify">
The bed is horizontal downstream and then connected to a $35^\circ$ slope by a **quadratic transition** of width `physical_smoothing_width`. This transition is continuous in both elevation and slope, which avoids an abrupt topographic break.
</p>

<p align="justify">
The pile is initially at rest behind the gate, between $x=x_{\rm gate}$ and $x=x_{\rm gate}+L_{\rm pile}$. Its free surface is horizontal and set at the elevation `physical_gate_bottom_height + physical_initial_pile_height`. The chosen values match the dimensions reported by Poulain *et al.* (2023). Here `h[]` denotes a vertical height above the bed.
</p>
*/

event init (i = 0)
{
  foreach() {
    zb[] = smooth_bottom(x);                                            // Nondimensional topography at the cell center.

    double inside_initial_pile =                                        // Cell located inside the initial pile footprint.
      (x >= gate_position && x <= gate_position + initial_pile_length);
    double initial_surface_level =                                      // Elevation of the initial horizontal free surface.
      gate_bottom_height + initial_pile_height;

    // Initial vertical height above the bed.
    h[] = inside_initial_pile ? max(initial_surface_level - zb[], 0.) : 0.;

    u.x[] = 0.;                                                         // Pile initially at rest.
  }

  boundary ({zb, h, u});                                                // Update boundary conditions.
}

/**
## Basal friction

<p align="justify">
At each iteration, the bed slope $\partial_x z_b$ is computed with `smooth_bottom_slope` and the height slope $\partial_x h$ is computed with a **centered difference**. These two terms give the free-surface slope $\partial_x(z_b+h)$ used in the static branch of the Pouliquen--Forterre law.
</p>

<p align="justify">
The local Froude number is then evaluated as:
</p>

$$
\displaystyle
Fr=\frac{|u|}{\sqrt{G\,h\,\cos\alpha}},
\qquad
\cos\alpha=\frac{1}{\sqrt{1+(\partial_x z_b)^2}}.
$$
*/

event friction (i++)
{
  foreach() {
    if (h[] <= dry) {                                                   // Condition sur le sol sec.
      u.x[] = 0.;
      continue;
    }

    double bottom_slope = smooth_bottom_slope(x);                       // Analytical bed slope.
    double height_slope = (h[1] - h[-1])/(2.*Delta);                    // Discrete height slope.

    double bottom_cosine = 1./sqrt(1. + sq(bottom_slope));              // Cosine of the local bed angle.
    double flow_speed = fabs(u.x[]);                                    // Magnitude of the longitudinal velocity.
    double froude_number = flow_speed/sqrt(G*h[]*bottom_cosine);        // Local Froude number.

    double effective_friction =                                         // Effective friction coefficient.
      mu_pouliquen_forterre(froude_number, h[], bottom_slope,
                            height_slope);

    double reduced_speed = max(flow_speed - G*effective_friction*dt,    // Deceleration due to basal friction.
                               0.);

    u.x[] = flow_speed > 1e-30 ? reduced_speed*u.x[]/flow_speed : 0.;   // Reduce the magnitude without changing the sign.
  }

  boundary ({u});
}

/**
## Plots: profiles, front position and animation

<p align="justify">
Profiles are saved at $t=0$, $0.2$, $0.4$ and $3$ s for the final comparison.
The front position is tracked in parallel. The profiles used for the animation
are saved every `physical_output_dt` seconds and then read back by Gnuplot
in the Main results section.
</p>
*/

// Initial profile used for comparison with the figure from Poulain et al.
event profile_t0 (t = 0.)
{
  write_comparison_profile ("profile_dry35_t0.dat");

  foreach()
    fprintf (stderr, "%.12g %.12g %.12g %.12g %.12g\n",
             x*physical_length_scale, h[]*physical_length_scale, t*physical_time_scale,
             zb[]*physical_length_scale, u.x[]*physical_length_scale/physical_time_scale);
}

// Profile at t = 0.2 s.
event profile_t02 (t = 0.2/physical_time_scale)
{
  write_comparison_profile ("profile_dry35_t02.dat");

  foreach()
    fprintf (stderr, "%.12g %.12g %.12g %.12g %.12g\n",
             x*physical_length_scale, h[]*physical_length_scale, t*physical_time_scale,
             zb[]*physical_length_scale, u.x[]*physical_length_scale/physical_time_scale);
}

// Profile at t = 0.4 s.
event profile_t04 (t = 0.4/physical_time_scale)
{
  write_comparison_profile ("profile_dry35_t04.dat");

  foreach()
    fprintf (stderr, "%.12g %.12g %.12g %.12g %.12g\n",
             x*physical_length_scale, h[]*physical_length_scale, t*physical_time_scale,
             zb[]*physical_length_scale, u.x[]*physical_length_scale/physical_time_scale);
}

// Final profile at t = 3 s.
event profile_t30 (t = 3./physical_time_scale)
{
  write_comparison_profile ("profile_dry35_t30.dat");

  foreach()
    fprintf (stderr, "%.12g %.12g %.12g %.12g %.12g\n",
             x*physical_length_scale, h[]*physical_length_scale, t*physical_time_scale,
             zb[]*physical_length_scale, u.x[]*physical_length_scale/physical_time_scale);
}

// Suivi du front : distance parcourue depuis la porte.
event front (t = 0.; t <= simulation_end_time; t += physical_output_dt/physical_time_scale)
{
  double front_position = gate_position;

  foreach (reduction(min:front_position))
    if (h[] > front_height_threshold)
      front_position = min(front_position, x);                          // Farthest front toward decreasing x.

  fprintf (stdout, "%.12g %.12g\n", t*physical_time_scale, (gate_position - front_position)*physical_length_scale);
}

// Save the profiles that will be read back by the GIF Gnuplot block.
event save_animation_profiles (t = 0.; t <= simulation_end_time; t += physical_output_dt/physical_time_scale)
{
  // The frame index is obtained directly from the physical time.
  int frame_index = (int) lround((t*physical_time_scale)/physical_output_dt);
  write_animation_profile (frame_index);
}

event stop (t = simulation_end_time)
{
  return 1;
}

/**
## Comparison resources and saved profiles

<p align="justify">
The four PNG backgrounds are downloaded at startup from the GitHub repository.
The CSV file containing the experimental profiles is also retrieved automatically.
The first three experimental profiles are extracted into separate `.dat` files
so that Gnuplot can overlay them on the final comparison.
</p>
*/

// Files used for the results.

const char * github_resources_url =
  "https://raw.githubusercontent.com/kilian-gthub/basilisk-ressources/main/ressources";

const char * hysea_background_files[] = {
  "t0_green_dashed_plus_black.png",
  "t0_2_green_dashed_plus_black.png",
  "t0_4_green_dashed_plus_black.png",
  "t3_green_dashed_plus_black.png"
};

const char * experimental_csv_url =
  "https://basilisk.fr/sandbox/Kilian/ressource/data_article_poulain.csv";

const char * experimental_csv_file = "data_article_poulain.csv";

const char * comparison_profile_files[] = {
  "profile_dry35_t0.dat",
  "profile_dry35_t02.dat",
  "profile_dry35_t04.dat",
  "profile_dry35_t30.dat"
};

const char * experimental_profile_files[] = {
  "exp_profile_t0.dat",
  "exp_profile_t02.dat",
  "exp_profile_t04.dat"
};

static void extract_experimental_profiles (const char * csv_filename);

// Download the HySEA-NF backgrounds and the experimental CSV file.
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
              hysea_background_files[k],
              hysea_background_files[k], github_resources_url, hysea_background_files[k],
              hysea_background_files[k], github_resources_url, hysea_background_files[k]);

    if (system(command) != 0)
      fprintf (stderr, "WARNING: unable to download %s. See reference_files.log.\n",
               hysea_background_files[k]);
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

// Supprime les retours chariot en fin de ligne.
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

// Columns: x as displayed in the article, z_b, z_b+h, h and horizontal velocity.
static void write_comparison_profile (const char * filename)
{
  FILE * fp = fopen (filename, "w");

  if (!fp) {
    fprintf (stderr, "WARNING: unable to write %s.\n", filename);
    return;
  }

  foreach()
    fprintf (fp, "%.12g %.12g %.12g %.12g %.12g\n",
             -x*physical_length_scale, zb[]*physical_length_scale, (zb[] + h[])*physical_length_scale,
             h[]*physical_length_scale, u.x[]*physical_length_scale/physical_time_scale);

  fclose (fp);
}

// One file is written per animation frame. It is kept so that
// the ~~~gnuplot block in the Basilisk documentation can then generate
// and display the GIF directly on the sandbox page.
static void write_animation_profile (int frame_index)
{
  char filename[128];
  snprintf (filename, sizeof(filename), "frame_dry35_%04d.dat", frame_index);

  FILE * fp = fopen (filename, "w");
  if (!fp) {
    fprintf (stderr, "WARNING: unable to write %s.\n", filename);
    return;
  }

  foreach()
    fprintf (fp, "%.12g %.12g %.12g %.12g\n",
             -x*physical_length_scale, zb[]*physical_length_scale, (zb[] + h[])*physical_length_scale, h[]*physical_length_scale);

  fclose (fp);
}


/**
# Main results

## Collapse animation

<p align="justify">
Figure 2 shows the time evolution of the dry 35 degrees granular collapse.
The black line is the bed, the
pale yellow area is the granular mass and the sand-yellow line is the free
surface of the pile. At the beginning, the material is **compact near the gate**.
Then it spreads along the slope, the front moves downslope and the maximum
thickness decreases. This is the expected behaviour: **gravity stretches** the
pile while basal friction progressively slows it down.
</p>

~~~gnuplot Figure 2 - Dry 35 degrees granular collapse computed with Basilisk: time evolution of the pile over the smoothed topography (black bed, pale-yellow pile, sand-yellow free surface)
set terminal gif animate delay 5 loop 0 size 850,672 enhanced font "Liberation Sans,12"
set output "animate_dry35_muv_hysea.gif"

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
    datafile = sprintf("frame_dry35_%04d.dat", k)

    set title sprintf("Basilisk dry 35 degrees collapse - t = %.2f s", time)
    plot datafile using 1:(-0.05):2 with filledcurves lc rgb slope_fill title "slope", \
         datafile using 1:2:3 with filledcurves lc rgb pile_fill title "granular pile", \
         datafile using 1:2 with lines lw 2.2 lc rgb "#000000" title "topography", \
         datafile using 1:($4 > 1.e-4 ? $3 : 1/0) with lines lw 2.4 lc rgb sand_line title "pile free surface"
}

unset output
~~~

## Comparison with the experiment and HySEA

<p align="justify">
Figure 3 compares the result obtained with Basilisk against experimental profiles and the non-hydrostatic HySEA result reported by Poulain et al. The objective is not only to show the final position but also to verify whether the model reproduces the same **overall motion** at various time steps. The fluid mass simulated with Basilisk remains virtually identical to that simulated with HySEA, thereby **validating the code**. However, both codes still diverge from the experimental data—particularly in **the transitional phase**—which is to be expected for a simplified non-hydrostatic shallow-water model.
</p>

~~~gnuplot Figure 3 - Dry 35 degrees granular collapse: Basilisk shallow-water result compared with experimental profiles and non-hydrostatic HySEA backgrounds at t = 0, 0.2, 0.4 and 3 s
reset session
set terminal pngcairo size 850,672 enhanced font "Liberation Sans,10"
set output "comparison_dry35_hysea_nf.png"

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
set style line 4 lc rgb "#32a852" lw 2.4 dt 2

# Global margins: the four panels are slightly reduced to place
# the axis labels and the legend outside the article images.
set object 100 rect from screen 0.075,0.905 to screen 0.900,0.985 fc rgb "white" fillstyle solid 0.85 noborder front
set label 100 "x (m)" at screen 0.525,0.035 center front font "Liberation Sans,12"
set label 101 "z (m)" at screen 0.030,0.490 center rotate by 90 front font "Liberation Sans,12"

set arrow 201 from screen 0.105,0.965 to screen 0.155,0.965 nohead front lw 2.2 lc rgb "#000000"
set label 201 "Basilisk topography" at screen 0.160,0.955 left front font "Liberation Sans,9"
set arrow 202 from screen 0.420,0.965 to screen 0.470,0.965 nohead front lw 2.6 dt 2 lc rgb "#d8a72e"
set label 202 "Basilisk shallow-water" at screen 0.475,0.955 left front font "Liberation Sans,9"
set arrow 203 from screen 0.105,0.935 to screen 0.155,0.935 nohead front lw 2.4 lc rgb "#40549b"
set label 203 "Experimental profile" at screen 0.160,0.925 left front font "Liberation Sans,9"
set arrow 204 from screen 0.420,0.935 to screen 0.470,0.935 nohead front lw 2.4 dt 2 lc rgb "#32a852"
set label 204 "Non-hydrostatic HySEA" at screen 0.475,0.925 left front font "Liberation Sans,9"

set multiplot
set size 0.455,0.410

# t = 0 s: top-left panel
set origin 0.070,0.490
set xrange [0:1545]
set yrange [0:1221]
plot "t0_green_dashed_plus_black.png" binary filetype=png with rgbimage notitle, \
     "profile_dry35_t0.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(2.5 + ($2 - 0.0)*(1053.5 - 2.5)/(0.3 - 0.0)) with lines ls 1 notitle, \
     "profile_dry35_t0.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):($4 > 1.e-4 ? (2.5 + ($3 - 0.0)*(1053.5 - 2.5)/(0.3 - 0.0)) : 1/0) with lines ls 2 notitle, \
     "exp_profile_t0.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(2.5 + ($2 - 0.0)*(1053.5 - 2.5)/(0.3 - 0.0)) with lines ls 3 notitle

# t = 0.2 s: top-right panel
set origin 0.525,0.490
set xrange [0:1545]
set yrange [0:1219]
plot "t0_2_green_dashed_plus_black.png" binary filetype=png with rgbimage notitle, \
     "profile_dry35_t02.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(0.5 + ($2 - 0.0)*(1051.5 - 0.5)/(0.3 - 0.0)) with lines ls 1 notitle, \
     "profile_dry35_t02.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):($4 > 1.e-4 ? (0.5 + ($3 - 0.0)*(1051.5 - 0.5)/(0.3 - 0.0)) : 1/0) with lines ls 2 notitle, \
     "exp_profile_t02.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(0.5 + ($2 - 0.0)*(1051.5 - 0.5)/(0.3 - 0.0)) with lines ls 3 notitle

# t = 0.4 s: bottom-left panel
set origin 0.070,0.080
set xrange [0:1543]
set yrange [0:1221]
plot "t0_4_green_dashed_plus_black.png" binary filetype=png with rgbimage notitle, \
     "profile_dry35_t04.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(0.5 + ($2 - 0.0)*(1053.36270491803 - 0.5)/(0.3 - 0.0)) with lines ls 1 notitle, \
     "profile_dry35_t04.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):($4 > 1.e-4 ? (0.5 + ($3 - 0.0)*(1053.36270491803 - 0.5)/(0.3 - 0.0)) : 1/0) with lines ls 2 notitle, \
     "exp_profile_t04.dat" using (96 + ($1 - 0.6)*(1255 - 96)/(0.0 - 0.6)):(0.5 + ($2 - 0.0)*(1053.36270491803 - 0.5)/(0.3 - 0.0)) with lines ls 3 notitle

# t = 3 s: bottom-right panel
set origin 0.525,0.080
set xrange [0:1543]
set yrange [0:1217]
plot "t3_green_dashed_plus_black.png" binary filetype=png with rgbimage notitle, \
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
*/
