/**

# Introduction

# The Code

## How to execute?

If we want the stratified case, add `-DDR=1`

~~~
rm -r kolmogorov/
mkdir kolmogorov
qcc -autolink -O -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -DMTRACE=3 -DDR=1 -o kolmogorov/kolmogorov kolmogorov.c -lm
cd kolmogorov
srun kolmogorov
~~~

## Set Up

We use the multilayer solver with the non-hydrostatic extension.

*/
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"

#include "layered/perfs.h"

  /**
We also add an option to make the simulation with a stratified flow.
*/
#if DR
#define drho(T) (-1.*(T))
#include "layered/dr.h"
#endif
  /**

We define the four physical parameters and give them default values. Even if both Froude and Prantl numbers are used only in the stratified case, we define them even if we don't use it for practicity in the main code.
*/
double Re = 1.;
double ratio = 0.5;
double Fr = 1.;
double Pr = 1.;
  /**

We also define as a global variable for the spatial resolution `N_res`. We will use it later to define the resolution in both directions in the following way : $N_z = N_{res}$ and $N_x =  N_{res}/\alpha$. Such a definition for the resolutions ensure that $\Delta x = \Delta z$.

*/
int N_res = 8;
  /**

Now, we define some numerical variable that will control the simulation, which runs as follow:

* the simulation is initialized with the *theoretical* laminar solution, which is not the *numerical* laminar solution (they are the same up to order $\Delta x^2$). This profile will diffuse with the characteristic timescale of diffusion $\tau = Re$ (in advective time) toward the *numerical* laminar solution.
* once we reached `T_perturb` $=$ `wait_X_tau_laminar` $\tau$, we add a (divergence-free) perturbation to the velocity field
* again, we wait `T_fit` $=$ `wait_X_tau_perturb` $\tau$ to let the most increasing perturbation dominate the other ones, and we measure its exponential variation during `T_growth`

*/
double Tmax = 0.;
double T_growth = 20.;

double wait_X_tau_laminar = 3.;
double wait_X_tau_perturb = 2.;

double T_fit;
double T_perturb;
  /**
  
The following booleans state if we want to get:

* an animation of the velocity fields as output i.e. u.x(x,z)(t) and u.z(x,z)(t) (`do_movie = true`)
* an animation of the velocity profile along the streamwise direction i.e. u.x(z)(t) (`do_profile = true`)
* a plot of the exponential behavior of the perturbation (`do_plot_growth = true`)
* snapshots (under the form of datafile) of $u_x$, $u_z$ and $T$ in the stratified case (`do_snapshot = true`)
* an animation of the velocity field as a vector field (`do_arrowplot = true`)

*/
bool do_movie = false;
bool do_profile = false;
bool do_plot_growth = false;
bool do_snapshot = false;
bool do_arrowplot = false;
  /**

Here, we define some universal properties of the flow that will stay the same throughout all simulations such as the vertical size of the domain `H0`, the number of wavelength for the forcing in the vertical direction `N_osc` and the amplitude of the perturbation `epsilon`.

*/
int N_osc = 1;
double H0;

double epsilon = 0.001;

/**

Finally, we define some global variables needed in various functions.

*/
char simulation_type[16];
char extension_phy[64];
char extension_num[64];
double A_cos;
bool perturb_added;

/**

### General before-main functions

At this point, and just before the `main()` function, we need to define functions that will be executed after some runs, such as the plot of the stability diagram.

We define a general function to set a title with gnuplot depending on both physical and numerical parameters. Currently, these parameters are the Reynolds number `Re`, the aspect ratio of the domain `ratio`, the amplitude of the perturbation `epsilon`, the spatial resolution `n_resolution` and the maximal timestep `DT`.

*/
void gnu_title (FILE * fp, char * title) {
  char gnu_command[512];
  snprintf(gnu_command, sizeof(gnu_command),
           "set title \" %s - %s \\nRe = %g, {/Symbol a} = %g, Fr = %g , {/Symbol e} = %g \\n n_res = %d, DT = %g \"\n",
           title, simulation_type, Re, ratio, Fr, epsilon, N_res, DT);
  fprintf(fp, gnu_command); //fputs(gnu_command, fp); ?
}
  /**

We define a function to draw the stability diagram (which will be used at the very end of the code execution).

*/
void plot_stability() {
  FILE * fp = popen("gnuplot 2> /dev/null", "w");
  fprintf(fp,
          "set ylabel '{/Symbol a}'\n"
          "set xlabel 'Re'\n"
          "set grid\n"
          "set title \" Stability Diagram - %s\\nn_res = %d, DT = %g, Fr = %g\"\n",
          simulation_type, N_res, DT, Fr
  );

  fprintf(fp,
          "files = system(\"ls sigma_%s\")\n"
          "N = words(files)\n"
          "set palette defined (0 \"blue\", 1 \"red\")\n"
          "set cbrange [0:1]\n"
          "set pointsize 2\n"
          "unset colorbox\n"
          "p for [i=1:N] word(files,i) u 2:1:(($3 > 0) ? 1 : 0) w p linetype 7 lc palette t '', '../Rec_freeslip' u 1:2 w l t 'Theory - G.M.'\n"
          ,extension_num
  );
  fprintf(fp,
          "set term pngcairo size 1024,720\n"
          "set output 'ReCRatio_%s.png'\n"
          "replot\n"
          , extension_num);
  fflush(fp);

}
  /**

The last before-main function we define is a plot of the growth rate $\sigma$ as a function of the Reynolds number for various values of the aspect ratio.

TODO : change this by a gnuplot-basilisk code at the end

*/
void plot_sigma_all (FILE * fp) {
  gnu_title(fp, "");
  fprintf(fp,
          "set xlabel '%s'\n"
          "set ylabel '{/Symbol s}'\n"
          "set grid\n",
          "Re"/*, CFL_H*/); ///////////////////////// var_3
  fprintf(fp,
          "files = system(\"ls sigma_%s\")\n"
          "N = words(files)\n"
          , extension_num);

  if (false) {
  fprintf(fp,
          "p for [i=1:N] word(files,i) u 2:3:4 w errorbars t substr(word(files,i),7,32) \n"
          "set term pngcairo size 1024,720\n"
          "set output 'sigmaAll_out_%s.png'\n"
          "replot\n"
          , extension_num);
  }
  // same plot but filtering the fit with the RMS residual //
  fprintf(fp,
          "p for [i=1:N] word(files,i) u 2:(($5<0.15) ? $3 : NaN):4 w errorbars t substr(word(files,i),7,32) \n"
          "set term pngcairo size 1024,720\n"
          "set output 'sigma_out_%s.png'\n"
          "replot\n"
          , extension_num);

  fflush(fp);
}
  /**

If one wants to study the flow near the transition, we also define a rough approximation of the critical Reynolds number of the linear instability of the laminar flow.

*/
double Re_approx (double alpha) {
  return (sqrt(2)/(1 - sq(alpha)));
}
  /**

## Run

*/
int main () {
  cell_lim = mono_limit;
  
  /**
We define the vertical size of the domain such that the wavelength of the forcing is $2\pi$.
*/
  H0 = 2.*pi*N_osc;
  /**

We change the slip length default value of the Navier boundary condition at the bottom in order to get a free-slip boundary condition (see [multilayer.h](/src/multilayer.h))
*/
  lambda_b[] = {HUGE};
  /**

As mentionned before, we use the `rigid lid` option in order to get a free-slip boundary condition (and then symetrize the flow). So we set `rigid = true`, and in order to avoid numerical issues, we also set $\theta_H$ to some value which is not 1/2 (as done in [gulf-stream.c](/src/examples/gulf-stream.c)).

*/  
  theta_H = 0.51;
  rigid = true;
  /**

The left-right boundary conditions are set to be periodic:

*/
  periodic(right);
  /**

*/
  DT = 1.; // need to be non-dimensional in order to work with non-dimensional system
  TOLERANCE = 1e-3;

  /**
We set two loops to iterate over the numerical resolution (both temporal and spatial) of the simulation, which can be used to check numerical convergence.
  
*/  
  for (DT = 0.01; DT >= 0.01; DT /= 10.) {
    for (N_res = 32; N_res <=32; N_res *=2) {
      snprintf(extension_num, sizeof(extension_num), "DT%g_N%d", DT, N_res);
  /**


*/
      for (ratio = 0.1; ratio <= 0.9; ratio += 0.1) {
        for (Re = 5.; Re <= 5.; Re += 1.) {
  /**
We set an option if we want to do simulation close to the line $Re = Re_c(\alpha)$ in the parameter space. To do so, we compute an approximation of the critical Reynolds number (at given aspect ratio $\alpha$), and we do the simulation only if the Reynolds number of the simulation is not too far.
*/
          double Re_c_th = Re_approx(ratio);
          if (true || ((Re_c_th - 1 < Re) && (Re < 2*(Re_c_th + 1) + 1))) {

#if DR
            snprintf(extension_phy, sizeof(extension_phy), "Re%g_alpha%g_Fr%g", Re, ratio, Fr);
            snprintf(simulation_type, sizeof(simulation_type), "stratified");
#else
            snprintf(extension_phy, sizeof(extension_phy), "Re%g_alpha%g", Re, ratio);
#endif
  /**
  We use the physical parameters to define the horizontal size of the domain, and the numerical resolution (which are thus kind of intensive quantities).
*/
            L0 = 2.*pi/ratio;
            nl = N_res; //number of layer
            N = floor(N_res/ratio); // number of points in the non layered direction.s

  /**            
  And we define the non-dimensionnalized viscosity and gravity.
*/
            nu = 1./Re;
            G = 1.;//  /sq(Fr); //useful only in stratified cases

  /**
  Finally, we set the times of some events as discussed above, and run the simulation.
*/
            
            T_perturb = wait_X_tau_laminar*Re;
            T_fit = T_perturb + wait_X_tau_perturb*Re;
            //if (Tmax == 0.) { //if we want such a thing, we should define bool adaptative_time = true, because here, T_max won't change after the first simulation since it won't be null anymore
            Tmax = T_fit + T_growth;
            //}
            fprintf(stderr, "\nRunning %s %s\n", extension_phy, extension_num);
#if DR
            fprintf(stderr, "Boussinesq\n");
#endif
            fprintf(stderr, "T_perturb = %g, T_fit = %g, Tmax = %g\n", T_perturb, T_fit, Tmax);
            A_cos = 0.;
            run();
            //sleep(10.);
          } //endif
        } //end of Re loop
      } //end of ratio loop
      plot_stability();
    } //end of N_res loop
  } //end of DT loop
  fprintf(stderr, "\nSimulations ended succesfully\n");
}
  /**
We initialize the flow either with an input file (wip), or with the laminar flow if such file does not exist.
Note that we also have to initialize the height of each layer.
  
*/
event init (i = 0) {
  /**
We check for the existence of an init file.
  
*/
  char filename[64];
  snprintf(filename, sizeof(filename),"../mean_profile_%d",N_res);
  FILE * fp = fopen(filename, "r");                         
  /**
If such file does not exist, we initialize with the theoretical laminar flow: $\mathbf{u} = cos(z) \mathbf{e_x}$. And with no deviation of the stratification with respect to the background one.

*/
  if (!fp) {                                                
    fprintf(stderr, "No file to initialize with. Initializing with the theoretical laminar profile\n");
    foreach() {
      zb[] = 0.;
      double z = 0.;
      foreach_layer() {                                     
        h[] = H0/nl;                                        
        z += h[]/2.;                                         
        u.x[] = 1.*cos(z); // the laminar profile
        w[] = 0.;
        phi[] = 0.;
#if DR
        T[] = 0; // no deviation
#endif
        z += h[]/2.;                                        
      }                                            
    }                         
  }  
  /**
If an init file exist, we use it.

*/
  else {                                                    
    double u_lam[nl];                                       
    for (int l = 0; l < nl; l++) {                          
      fscanf(fp, "%lf", &u_lam[l]);                         
    }                                                       
    foreach() {                                             
      zb[] = 0.;                                            
      for (int l = 0; l < nl; l++) {                        
        h[0,0,l] = H0/nl;                                   
        u.x[0,0,l] += u_lam[l];
        // TO DO : ADD T
      }                                                     
    }                                                       
  }
  fprintf(stderr, "Initialized\n");  
  perturb_added = false;                                                      
}
  /**
 
 
*/
void fit_mean_profile () {
  double norm = 0;
  foreach(serial) {
    double z = zb[];
    foreach_layer() {
    z += h[]/2.;
    norm += fabs(u.x[]-cos(z)); //maybe patch with some volume average
    z += h[]/2.;
    }
  }
  norm /= 1.*N*nl;
  // fprintf(stderr, "degug - norm = %g\n", norm);
  A_cos = 1. + norm/(1. - exp(-wait_X_tau_laminar));
  fprintf(stderr, "mean profile has an approximate amplitude of %g\n", A_cos);
}
  /**
Then, we add a divergence-free perturbation.
Since we wanted to add this perturbation at a different time depending on the Reynold number, I've encountered some issue with the compliation, hence the presence of the function comparator.

*/
bool comparator (double _time, double _max) {
  //fprintf(stderr, "%g\n", _max);
  return (_time >= _max);
}

double degug_timeperturbation (double _time) {
  return _time;
}
                                                            
event perturbation (i++) {
  bool is_time = comparator(t, T_perturb);
  if (is_time && (perturb_added == false)) {
    fprintf(stderr, "Perturbation added at t = %g\n", t);
    perturb_added = true;
/**
Just before we add the perturbation, we fit the mean profile
*/
    fit_mean_profile();
/**
We want to add a small divergence-free perturbation, exciting the mode $\alpha$ in the $x$ direction, that respect the boundary conditions.

Mathematically, we will use a stream function $\psi$ (such that $\mathbf{u} = \mathbf{rot}(\psi \mathbf{e_y})$) defined by : $\psi = \epsilon \cos (\alpha x) \textrm{regul}(z)$ where $\textrm{regul}$ is a smooth function that goes to 0 up to its third derivative at $z = 0$ or $z = H0$ in order to respect the free-slip boundary condition.
We choose to use $\textrm{regul}(z) = \exp\left( \dfrac{1}{\left(\dfrac{z - \pi}{\pi}\right)^2 ~ - 1}\right) e$
*/
    double regul;
    scalar psi = new scalar[nl];
    foreach() {
      double z = zb[];
      for (int l = 0; l<nl; l++) {
        z += h[0,0,l]/2.;
        regul = exp(1/(sq((z-pi)/(1.*pi)) -1))*exp(1.);
        psi[0,0,l] = epsilon*cos(ratio*x)*regul;
#if DR
        T[0,0,l] += -1.*epsilon*cos(ratio*x)/sq(Fr);
#endif
        z+= h[0,0,l]/2.;
    } //end of the layer loop
  } //end of foreach
/**
Then, we compute the 2D rotational of this stream function:
$\mathbf{u} = \boldsymbol{\nabla}\land(\psi \mathbf{e_y}) = \begin{pmatrix}
-\partial_z \psi \\
0 \\
+\partial_x \psi \\
\end{pmatrix}$

The derivatives are here approximated by a centered scheme. Despite looking for a free-slip perturbation, our definition of the `regul` function implies that the perturbation is no-slip, hence the 0 velocity added to the top and bottom layers.

*/
    foreach() {
      for (int l = 1; l<nl-1; l++) {
          u.x[0,0,l] += (psi[0,0,l+1] - psi[0,0,l-1])/(2.*h[0,0,l]) ;
	  w[0,0,l] += -(psi[1,0,l] - psi[-1,0,l])/(2.*Delta);
      } //end of the layer loop
    } //end of foreach
  } //endif
}
  /**
We define exterior force field $\mathbf{f} = \dfrac{1}{Re} \cos(z) \mathbf{e_x}$ which is implemented in this solver with $\mathbf{f} = \textrm{ha}/\textrm{h}$.
  
*/
event acceleration (i++) {
  foreach_face(x) {
    double z = zb[];
    foreach_layer() {
      z += h[]/2.;
      ha.x[] += h[]*1./Re*cos(z);
      z += h[]/2.;
    }
  }
}

  /**
We also add the viscosity for the horizontal component of the velocity, and particle diffusion in the stratified case.

See [diffusion.h](https://basilisk.fr/src/layered/diffusion.h#vertical_diffusion).

A way to rewrite the Robin boundary condition to transform it into a Neumann BC is the following:

$$\partial_z T_b = \dfrac{T_b - T_0}{\lambda_b}$$

So we can make $\lambda_b \longrightarrow + \infty$ keeping the following ratio : $T_0/\lambda_b = -1$
  
*/
event viscous_term (i++) {
  horizontal_diffusion({u.x, w}, 1./Re, dt);
#if DR
  horizontal_diffusion({T}, 1./(Re*Pr), dt) ;
  foreach()
    vertical_diffusion (point, h, T, dt, 1./(Re*Pr), 1, -1*HUGE, HUGE);
//    vertical_diffusion (point, h, T, dt, 1./(Re*Pr), 0, 0, HUGE);
#endif
}

  /**
As we defined the flotability as a deviation of some background stratification, we should add the non-linear term caused by this background stratification.

Not that by adding manualy this term, the stratified version of the solver can be correct only up to the order $1$ in $\Delta t$.
*/

#if DR
event density_advection (i++) {
  foreach() {
    foreach_layer() {
      T[] -= dt*w[]/sq(Fr); 
    }
  }
}
#endif

  /**

## Utilities

First, as done in [utils.h](/src/utils.h), we define a fuction computing the vorticity field in the 2D multilayer case.

As we use free-slip BC, the vorticity is null at the boundaries, and we could use this condition to implement a general centered scheme using finite difference for the computation in the 'bulk' (i.e. not at the boundaries).

*/
void vorticity_layered (vector u, scalar omega) {
  scalar dxuy[] = 0.;
  scalar dyux[] = 0.;
  foreach() {
    for (int l = 0; l < nl; l++) {
      if (l != 0 && l != nl-1 ) {
	dxuy[0,0,l] = (w[1,0,l] - w[-1,0,l])/(2.*Delta);
	dyux[0,0,l] = (u.x[0,0,l+1] - u.x[0,0,l-1])/(2.*h[]);
	omega[0,0,l] = dxuy[0,0,l] - dyux[0,0,l];
      }
      else {
	omega[0,0,l] = 0.;
      }
    }
  }
}
  /**
  
We also define a function to get the projection of a field at some altitude `z0` on the Fourier mode `kx`. Mathematically, we compute the following quantities:

$f_{c}(k_x,z) = \int_{L_x} f(x,z) \cos(k_x x) dx$

$f_{s}(k_x,z) = \int_{L_x} f(x,z) \sin(k_x x) dx$

$f(k_x,z) = \sqrt{f_{c}^2 + f_{s}^2}$
*/
double FT_1D_mode (scalar field, double kx, double z0) {
  double field_kc = 0.;
  double field_ks = 0.;
  foreach(reduction(+:field_ks) reduction(+:field_ks)) {
    double z = zb[];
    for (int l = 1; l < nl; l++) {
      if (l==1) {
        z += h[0,0,0];
      }
      if ((z-h[0,0,l-1]/2 <= z0) && (z0 < z+h[0,0,l]/2)) {
	// double approx_mode = field[0,0,l];
        double approx_mode = (field[0,0,l]*(z0 - (z-h[0,0,l-1]/2.)) + field[0,0,l-1]*((z + h[0,0,l]/2.) - z0 ) )/(h[0,0,l-1]/2. + h[0,0,l]/2.);
        field_kc += approx_mode*cos(kx*x)*Delta;
        field_ks += approx_mode*sin(kx*x)*Delta;
      }
      z += h[0,0,l];
    }
  }
  double field_k = sqrt(sq(field_kc) + sq(field_ks));
  return field_k;
}

/**

We need an indicator reflecting how non-laminar is the flow. We choose the L1 norm of uy

 */

double fourier_1D (scalar field, double k) {
  double field_cosk[nl];
  double field_sink[nl];
  double field_k[nl];
  foreach () {
    for (int l = 0; l<nl; l++){
      field_cosk[l] = field[0,0,l]*Delta*cos(x);
      field_sink[l] = field[0,0,l]*Delta*sin(x);
    }
  }
  for (int l = 0; l<nl; l++){
    field_k[l] = sqrt(sq(field_cosk[l]) + sq(field_sink[l]));
  }

  double mu = 0.;
  double mu_2 = 0.;
  double sigma = 0.;

  for (int l = 0; l<nl; l++){
    mu += field_k[l];
    mu_2 += sq(field_k[l]);
  }

  mu /= nl;
  mu_2 /= nl;
  sigma = sqrt(mu_2 - sq(mu));

  return mu;//, sigma} ;
}

/**

We will need an indicator reflecting how non-laminar is the flow. We choose the L1 norm of uy. 

 */

double indicator (scalar field) {
  double abssum = 0.;
  foreach(reduction(+:abssum)) {
    foreach_layer(){
      abssum += fabs(field[])*h[]*Delta; //maybe *ratio / N_osc
    }
  }
  return abssum;
}

/**
  And also an indicator reflecting how stratified the flow is. For this, we also choose the L1 norm of the total stratification, that is :
  $\mathcal{S} = b' + (z - H0/2)/Fr^2$
*/

double indicator_lin (scalar field, double slope) {
  double abssum = 0.;
  foreach(reduction(+:abssum)) {
    double z = zb[];
    foreach_layer(){
      z += h[]/2.;
      abssum += fabs(field[] - slope*(z - H0/2.))*h[]*Delta; //maybe *ratio / N_osc
      z += h[]/2.;
    }
  }
  return abssum;
}

/**
## Outputs

### Field Animation
See [lee.h](/src/examples/lee.h) for another example

Starting from the end, we make the movie by putting all the snapshots in an `.mp4` file.
*/
event moviemaker (t = end) {
  if (do_movie) {
    sleep(10.);
    char cmd[512];
    snprintf(cmd, sizeof(cmd),
             "for f in plot-*.png; do convert $f ppm:- ; rm -f $f; done | "
             "ppm2mp4 out_movie_%s_%s.mp4", extension_phy, extension_num);
    system(cmd);
  }
}
  /**

We need a function to draw a 2D field with gnuplot. Knowing that we want to see the perturbation, on $u_x$, we will need to substract the laminar flow. Hence, the function not only draws the field `f`, but can plot the field $f(x,z) - \textrm{substr\_cos} * \cos(z)$.

*/

void plot_field (FILE * fp, scalar f, char * f_name, double substr_cos, double substr_cte, double substr_slp) {
  gnu_title(fp,f_name);
  fprintf (fp,
           "set pm3d map\n"// interpolate 2,2\n"
           "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
           //"set cbrange [-5e-6:5e-6]\n"
           );
  fprintf(fp,
          "set xlabel 'x'\n"
          "set ylabel 'z'\n"
          "set xrange [0:%g]\n"
          "set yrange [0:%g]\n"
          ,L0, H0);
  fprintf(fp, "sp '-' u 1:2:3\n");
  foreach (serial) {
    double z = zb[];
    foreach_layer() {
      z += h[]/2.;
      fprintf (fp, "%g %g %g\n", x, z, f[] - substr_cos*cos(z) - substr_cte - substr_slp*z);
      z += h[]/2.;
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n"); 
}

/**

Thus, we define the event which will call the previous function every timestep (this can be changed, but having a snapshot every timestep has proved to be very useful for debugging).

*/
event gnuplot (i++) {
  double start_snap = T_perturb - DT;
  if (do_movie && (t>=start_snap)) {
    static FILE * fp = popen("gnuplot 2> /dev/null", "w");  

  /**
We can set which field we want to output by setting its index to `1`.
*/
    
    int index_ux = 1;
    int index_uy = 1;
    int index_omega = 0;


    int total_index = 2;

#if DR
    int index_T = 1;
    total_index +=1;
#endif

  /**

We set the multiplot, and plot each field, with an option to substract the laminar flow.

*/
      // this setup works with gnuplot v4.6 but may bug with v5.2
      fprintf(fp,
              "set term pngcairo size 1000,800\n"
              "set output 'plot-%06d.png'\n"
              "set multiplot layout %d,1\n",
              i, total_index);
      if (index_ux!=0) {
        char field_name[32];
        snprintf(field_name, sizeof(field_name), "u_x'(t = %g)", t);
        plot_field(fp, u.x, field_name, 0., 0., 0.);
//        plot_field(fp, u.x, field_name, A_cos, 0., 0.); // to substract the laminar flow
      }
      if (index_uy!=0) {
        char field_name[32];
        snprintf(field_name, sizeof(field_name), "u_y'(t = %g)", t);
        plot_field(fp, w, field_name, 0., 0., 0.);
      }
      if (index_omega!=0) {
        scalar omega = new scalar[nl];
        vorticity_layered(u, omega);
        plot_field(fp, omega, "omega", 0., 0. ,0.);
      }
#if DR
      if (index_T!=0) {
        char field_name[32];
        snprintf(field_name, sizeof(field_name), "b_tot(t = %g)", t);
        plot_field(fp, T, field_name, 0., H0/2./sq(Fr), -1./sq(Fr)); // adding the background flotability
      }
#endif
      fprintf (fp,
               "unset multiplot\n"
               "set output\n"
               );

    // fflush(fp);
  }
}
/**

### Profile u(z)

Similarly to the movie, we want to animate the profile at the end.

*/
event animate_profile (t = end) {
  sleep(10); // wait for the buffers of profile to be written
  if (do_profile) {
    char cmd[512];
    snprintf(cmd, sizeof(cmd),
             "for f in profile_%s_%s*.png; do convert $f ppm:- ; rm -f $f; done | "
             "ppm2mp4 profile_movie_%s_%s.mp4", extension_phy, extension_num, extension_phy, extension_num);
    system(cmd);
  }
}

/**

WIP : to move
We compute the Lp norm of the numerical residual when converging to the laminar flow.

*/

event residual (i++) { // t < T_perturb
  int norm_Lp = 1;
  double norm = 0.;
  foreach(serial) {
    double z = zb[];
    foreach_layer() {
    z += h[]/2.;
    norm += fabs(pow(u.x[]-cos(z),norm_Lp));
    z += h[]/2.;
    }  
  }

  char file[256];
  snprintf(file, sizeof(file), "residual_norm_%s_%s", extension_phy, extension_num);
  
  static FILE * fp = fopen(file, "w");
  fprintf(fp, "%g %g\n", t, pow(norm, 1./norm_Lp)/(N*nl));
  fflush(fp);
}

/**

We want a function that can do a snapshot of the profile. In order to do so, we will need two things :
- Write down the profile in a logfile
- Use the logfile in gnuplot

*/

event profile (i++) {
  if (do_profile) {
    FILE * fp = fopen("profile_log", "w");
    foreach(serial) {
      double z = zb[];
      foreach_layer() {
        z += h[]/2.;
        fprintf(fp, "%g %g %g\n", x, z, u.x[]);
        z += h[]/2.;
      }
      // fprintf(fp, "\n"); //separate each x
    } 
    fflush(fp);
    
    char cmd[256];
    snprintf(cmd, sizeof(cmd), "mv profile_log profile_log_%s_%s", extension_phy, extension_num);
    system(cmd);

    // plot //
    static FILE * pp = popen("gnuplot 2> /dev/null", "w");
    fprintf(pp,
            "reset\n"
            "set xrange [0.:%g]\n" //x lim
            "set grid\n"
            "set palette defined(0 \"red\", 1 \"blue\")\n"
            "unset colorbox\n"
            ,H0);

    // fit //
    if (false) {
    fprintf(pp,
           "set fit errorvariables\n"
           "f(x) = a*cos(x)\n"
           "a = 1\n"
           "fit f(x) 'profile_log_%s_%s' u 2:3 via a\n"
           "set print 'profile_fit_%s_%s'\n"
           "print a\n"
           // "print a_err\n"
           "unset print\n",
           extension_phy, extension_num, extension_phy, extension_num);
    }
    // plot
    char field_name[32];
    snprintf(field_name, sizeof(field_name), "ux(z, t = %g)", t);
    gnu_title(pp, field_name);
    fprintf(pp,
    "set yrange [-1.2:1.2]\n"
    "plot 'profile_log_%s_%s' u 2:3:1 w points lc palette, cos(x)\n"
    "set term pngcairo size 1024,720\n"
    "set output 'profile_%s_%s_%06d.png'\n"
    "replot\n",
    extension_phy, extension_num, extension_phy, extension_num, i);
    fflush(pp);
  }
}

event animate_profile (t = end) {
  //sleep(10); // wait for the buffers of profile to be written
  if (do_arrowplot) {
    char cmd[512];
    snprintf(cmd, sizeof(cmd),
             "for f in plotarrow_%s_%s*.png; do convert $f ppm:- ; rm -f $f; done | "
             "ppm2mp4 arrow_movie_%s_%s.mp4", extension_phy, extension_num, extension_phy, extension_num);
    system(cmd);
  }
}
/**

Output the final fields

*/

event snapshot (t+=1.) {
  if ((do_arrowplot || do_snapshot) && perturb_added) {
    char file[256];
    if (do_snapshot) {
      snprintf(file, sizeof(file), "snapshot_%s_%s_%06d", extension_phy, extension_num, i);
    }
    else {
      snprintf(file, sizeof(file), "snapshot");
    }
    //sleep(1);

    FILE * fp = fopen(file, "w");
    foreach(serial) {
      double z = zb[];
      foreach_layer() {
        z += h[]/2.;
#if DR
        fprintf(fp, "%g %g %g %g %g\n", x, z, u.x[], w[], T[]);
#else
        fprintf(fp, "%g %g %g %g\n", x, z, u.x[], w[]);
#endif
        z += h[]/2.;
      }
    }

    if (do_arrowplot) {

      fflush(fp);
      FILE * pp = popen("gnuplot 2> /dev/null", "w");
      char field_name[32];
      snprintf(field_name, sizeof(field_name), "u(t = %g)", t);
      gnu_title(pp, field_name);
      fprintf(pp,
              "set grid\n"
              "set xlabel 'x'\n"
              "set ylabel 'z'\n"
              "set xrange [0:%g]\n"
              "set yrange [0:%g]\n"
              , L0, H0
      );
      int N_plot = 31;
      int skip_plot = 1 + N_res/N_plot;
      double scale = 5.;
      //fprintf(stderr, "degug - %s, %d, %d, %g\n", file, N_res, skip_plot, scale);
      //char gnu_cmd[512];
      //snprintf(gnu_cmd, sizeof(gnu_cmd), "p '%s' u (((int($0)%% %d)%% %d == 0)*((int($0)/%d)%% %d == 0) ? $1 : NaN):2:($3/%g):($4/%g) w vectors notitle\n"
//              , file, N_res, skip_plot, N_res, skip_plot, scale, scale);

//      fprintf(stderr, "debug - %s\n", gnu_cmd);

      fprintf(pp,
              "p '%s' u (((int($0)%% %d)%% %d == 0)*((int($0)/%d)%% %d == 0) ? $1 : NaN):2:($3/%g):($4/%g) w vectors notitle\n"
              , file, N_res, skip_plot, N_res, skip_plot, scale, scale
      );
      fprintf(pp,
              "set term pngcairo size 1024,720\n"
              "set output 'plotarrow_%s_%s_%06d.png'\n"
              "replot\n",
               extension_phy, extension_num, i);
      fflush(pp);

    }
  }
}

  /**
### Growth Rate of the Perturbation

We generate a datafile containing the indicator of how non-laminar the flow is, and the indicator of how stratified it is if needed.

*/
event datafile (i++) {
  static FILE * fp = fopen("data", "w");
  double E = indicator(w);
#if DR
  double b_tot = indicator_lin(T, -1./sq(Fr));
  fprintf(fp, "%d %g %g %g\n", i, t, E, b_tot);
#else
  fprintf(fp, "%d %g %g\n", i, t, E);
#endif
  fflush(fp); //write the buffer immediately, so that growth will always have the complete data set

  /**
  We migh want to end the simulation if a non-(linear state is reached)
*/
  if (false && (E > 1.)) { // stop condition
    return 1.;
  }
  /**
  Or if the flow is completely homogeneized.
*/
#if DR
  if (false) {
    fprintf(stderr, "Mixed");
  }
#endif
}
/**

  Using the previous datafile, we generate figures of our indicators, and fit some informations like the growth rate, (WIP and the turbulent mixing timescale).

*/
event growth (t = end) {
  char cmd[256];
  snprintf(cmd, sizeof(cmd), "mv data data_%s_%s", extension_phy, extension_num);
  system(cmd);

/**
  We fit the growth rate
*/
  
  static FILE * fp = popen("gnuplot 2> /dev/null", "w");
  gnu_title(fp, "");
  // fit
  fprintf(fp,
          "set fit errorvariables\n"
          "FIT_CONVERGED = 0\n" // is initialized only if the fit is done, and be set to 1 if it converges 
          "a = 0\n"
          "b = 0\n"
          "f(x) = a + b*x\n"
          "fit [%g:][-8*log(10):-0*log(10)] f(x) 'data_%s_%s' u 2:(log($3)) via a,b\n" // we keep only the points after the projection t~>20, and the points that are not noise amp >e-8, and still linear (a not too big)
          "if (!FIT_CONVERGED) {a = NaN; b = NaN; b_err = NaN; FIT_STDFIT = NaN}\n"
          "set print 'sigma_%s' append\n" // give some error values if not converged
          "print %g, %g, b, b_err, FIT_STDFIT, %d, %g\n"
          "unset print\n",
           T_fit, extension_phy, extension_num, extension_num, ratio, Re, N, DT);
/**
  If set to do so, we generate the plots  
*/
  if (do_plot_growth) {
    fprintf(fp,
          "set logscale y\n"
          "set xlabel 't'\n"
          "set ylabel 'abs(uy_k)'\n"
          "set yrange [:1000]\n"
          "set grid\n"
          "plot 'data_%s_%s' u 2:3 w p t 'uy_L1', exp(a+b*x) t 'linear fit' \n"
  	  "set term pngcairo size 1024,720\n"
          "set output 'amp_growth_%s_%s.png'\n"
          "replot\n"
          , extension_phy, extension_num, extension_phy, extension_num);
#if DR
    fprintf(fp,
          "reset\n"
//          "set logscale y\n"
          "set xlabel 't'\n"
          "set ylabel 'abs(b_tot)'\n"
          "set yrange [0:]\n"
          "set grid\n"
          "plot 'data_%s_%s' u 2:4 w p t 'b_tot'\n"
          "set term pngcairo size 1024,720\n"
          "set output 'flo_growth_%s_%s.png'\n"
          "replot\n"
          , extension_phy, extension_num, extension_phy, extension_num);
#endif
  }
  fflush(fp);
/**
  In the stratified case, we also want to get the typical time taken by the turbulent flow to break the stratification. WIP
*/
#if DR
  fprintf(fp, "reset\n");
  double strat_0 = pi*sq(H0)/(2*ratio);
  fprintf(fp,
	  "set fit errorvariables\n"
	  "ksi(x) = ksiinf\n"
	  "ksiinf = 1\n"
	  "FIT_CONVERGED = 0\n"
	  "fit [%g:] ksi(x) 'data_%s_%s' u 2:4 via ksiinf\n"
	  , Tmax - 10, extension_phy, extension_num
         );

  fprintf(fp,
	  "ksi_r(x) = (x - ksiinf)/(%g - ksiinf)\n"
          "T_1 = 0\n"
	  "stats 'data_%s_%s' using (abs(ksi_r($4) - 1 ) > 0.05 ? $2 : 1/0) nooutput\n"
	  "T_1 = STATS_min\n"
          "T_2 = 0\n"
	  "stats 'data_%s_%s' using (abs(ksi_r($4) - 0 ) > 0.05 ? $2 : 1/0) nooutput\n"
	  "T_2 = STATS_max\n"
          "if (!FIT_CONVERGED) {T1 = NaN; T2 = NaN; ksiinf = NaN; ksiinf_err = NaN; FIT_STDFIT = NaN}\n"
	  "set print 'strat_%s' append\n"
	  "print %g, %g, %g, T_2 - T_1, T_1, T_2, ksiinf, ksiinf_err, FIT_STDFIT, %d, %g\n"
	  , strat_0, extension_phy, extension_num, extension_phy, extension_num, extension_num, Re, ratio, Fr, N, DT
         );
#endif
}

/**

### Others

*/
event logfile (i+=1) {
  fprintf(stdout, "%d %g\n", i, t);
}

/**

*/
event stop(t = Tmax) {
  return 1.;
}

/**

# Outputs

if needed, we can add a filtering on the second column to select a specific Reynolds number

 //we can choose a resolution by adding extension_num

~~~gnuplot Spectre

set xlabel 'k'
set ylabel '{/Symbol s}'
set grid
files = system("ls sigma_*")
N = words(files)

p for [i=1:N] word(files,i) u 1:(($5<0.15) ? $3 : NaN):4 w errorbars t substr(word(files,i),7,32)
~~~

*/