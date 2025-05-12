/**
# A brief overview of the Lamb-Oseen vortex

## Main characteristics

+ Non-steady planar flow $\longrightarrow \partial_t \neq 0$ and $(r, \theta)$ coordinates ;
+ Rotational symmetry for velocity field $\longrightarrow \boldsymbol{U} = u_\theta(r,t) \, \boldsymbol{e_\theta}$ ;
+ Newtonian fluid: the shear stress is proportional to the velocity gradient;
+ Homogeneous and incompressible flow $\longrightarrow \rho = C^{te}$ ;
+ Volumic forces are negligeable $\longrightarrow \boldsymbol{f_v} = \boldsymbol{0}$.

Therefore, the *Navier-Stokes* equations for the *Lamb-Oseen* vortex read as:

$$
\left\{\begin{matrix}
\boldsymbol{\nabla} \cdot \boldsymbol{U} &=& \boldsymbol{0} \\ 
\partial_t \boldsymbol{U} + \left(\boldsymbol{\boldsymbol{U} \cdot \boldsymbol{\nabla}} \right ) \boldsymbol{U} &=& - \dfrac{1}{\rho}\boldsymbol{\nabla}p + \nu \Delta \boldsymbol{U}  
\end{matrix}\right.
$$

## Local equations for pressure and azimuth velocity

Using the properties listed above for the considered flow, we do have the following simplifications:

$$
\left\{\begin{array}{lcl}
\dfrac{\partial p}{\partial r} & = & \rho \dfrac{{u_\theta}^2}{r}\\ 
& & \\
\dfrac{\partial u_\theta}{\partial t} & = & \nu \left(
\dfrac{\partial^2 u_\theta}{\partial r^2} 
+ \dfrac{1}{r} \dfrac{\partial u_\theta}{\partial r}
- \dfrac{{u_\theta}}{r^2}
\right ) 
\end{array}\right.
$$

As we will solve these equations in *cartesian coordinates* (and not in cylindrical ones), the azimuth velocity has to be projected:

$$
u_x = - u_\theta \sin(\theta) = - u_\theta \frac{y}{\sqrt{x^2 + y^2}} \quad ; \quad u_y = u_\theta \cos(\theta) = u_\theta \frac{x}{\sqrt{x^2 + y^2}}
$$




# Simulation in the physical space


## Setting general parameters

First, we import the `navier-stokes/centered` library. The full domain is of size $L_0 = 10$, with a grid resolution of $2^7 = 128$. Since kinematic viscosity is introduced in our problem, we enable it with the value $10^{-3}$ (by default, density is equal to $1.0$, according to its definition in the used solver). Parameters dedicated to the initial velocity profile are defined there, along with the *maximum level of refinement* for the mesh.

*/


#include "navier-stokes/centered.h"

#define LEVEL 7 // grid size
#define MAX_LEVEL 9 // refinement level
#define EPSILON 1.e-2
#define CORE_SIZE 1.e-1


// MAIN 
int main(){
    L0 = 10;
    origin(-L0/2., -L0/2.);
    init_grid(1 << LEVEL);

    // viscosity
    const face vector muc[] = {1e-3, 1e-3};
    mu = muc;

    run();
}

// ADAPTIVE REFINEMENT
event adapt (i += 50){
    adapt_wavelet ((scalar *){u}, (double []){5e-5, 5e-5}, maxlevel = MAX_LEVEL);
}


/**
The choice of parameters was constrained by the numerical rule of thumb of having *at least* 10 points in the smallest scale of the simulation. 
Here, with a vortex core radius $d=0.1$, the grid resolution should be less than $\Delta \leqslant 2d/10 = 2.10^{-2}$.

Since $L_0 = 2^N \times \Delta$, then:

$$
N_{max} = \left \lfloor \log_2 \left(\frac{L_0}{\Delta} \right ) \right \rfloor = 9
$$



## Initial Guess for the velocity profile

Given a *core radius* for the vortex of size $d = 0.1$ and a *desingularization term* $\varepsilon = 0.01$, an initial guess for the azimuth velocity profile is chosen for its generality (simulating a range of initial velocity profiles from *Rankine* to *Lamb-Oseen*), such as:

$$
u_\theta(r, t=0) 
=
\frac{r \left(1-\tanh \left(\frac{r-d}{\varepsilon }\right)\right)}{2 d}
+ \frac{d \left(\tanh \left(\frac{r-d}{\varepsilon }\right)+1\right)}{2 (r+\varepsilon )}
- \frac{d \left(1-\tanh \left(\frac{d}{\varepsilon }\right)\right)}{2 \varepsilon } 
$$

where the first term is the *linear part* corresponding to the solid revolving core of the vortex; the second one refers to the *inverse part* outside of the vortex core (**irrotational** flow is considered outside). Finally, the last term is nothing more than a *correction part*, whose purpose is having a *zero-value* in $(0,0)$, whatever was the choice of $(\varepsilon,d)$ couple.



*/


// INITIALIZATION OF VELOCITY FIELD
double velocity_profile(double epsilon, double d, double z){
    // epsilon: desingularization term; d: half-core size ; z: coordinate
    double linear_part = (z / (2. * d)) * (1. - tanh((z - d) / epsilon) ) ;
    double inverse_part = (d / (z + epsilon)) * (1./2.) * (1. + tanh((z - d) / epsilon) ) ;
    double correction_part = d * (1. - tanh(d / epsilon)) / (2. * epsilon) ;

    return (linear_part + inverse_part - correction_part) ;
}

    // boundary conditions are set up to 'symmetry' by default --> no action needed there
    // GAMMA = 1. and nu = 1.e-3 --> bc. rho = 1. by default so nu = mu / rho = 1.e-3
event init (t = 0){
    foreach() {
        double radius = sqrt(x*x + y*y) ;                                   
        u.x[] =  -1. * velocity_profile(EPSILON, CORE_SIZE, radius) * y/radius ;
        u.y[] =  velocity_profile(EPSILON, CORE_SIZE, radius) * x/radius ;
    }
}



/**
## Outputs

To observe a sufficient diffusion of the vorticity (calculated below), we need to go as further as $t = 3000$.

*/
  
// END OF SIMULATION
event end (t = 3000) {
    printf("\ni = %d t = %g\n", i, t);
}

/**
We then store in files `yprof_[t_i]`, for different time values $t_i$, the interpolated vertical velocity profiles, on the whole width of the domain, when $y = 0$.

*/

// OUTPUTS
event profiles (t = {0, 20, 40, 80, 100, 200, 300, 400, 500, 750, 1000, 3000}){
    char filename[20];
    sprintf(filename, "yprof_%g", t);
    FILE * fp = fopen(filename, "w");
    for (double x = -L0/2. ; x <= L0/2. ; x += 0.01)
        fprintf (fp, "%g %g\n", x, interpolate (u.y, x, 0));
    fclose (fp);
}


/**
~~~gnuplot Vertical Velocity Profiles in the Physical Space  
set xlabel 'x'
set ylabel 'u_y'
set xrange [-5:5]

plot "yprof_0" using 1:2 with lines title "t = 0", \
"yprof_40" using 1:2 with lines title "t = 40", \
"yprof_80" using 1:2 with lines title "t = 80", \
"yprof_100" using 1:2 with lines title "t = 100", \
"yprof_200" using 1:2 with lines title "t = 200", \
"yprof_300" using 1:2 with lines title "t = 300", \
"yprof_400" using 1:2 with lines title "t = 400", \
"yprof_500" using 1:2 with lines title "t = 500", \
"yprof_600" using 1:2 with lines title "t = 600", \
"yprof_700" using 1:2 with lines title "t = 700", \
"yprof_800" using 1:2 with lines title "t = 800", \
"yprof_900" using 1:2 with lines title "t = 900", \
"yprof_1000" using 1:2 with lines title "t = 1000", \
"yprof_3000" using 1:2 with lines title "t = 3000"
~~~
*/


/**
After computing the *vorticity* (`Basilisk C` function), it is then possible to output its scalar field and the different refinement levels of the grid, whose animations will allow us to see the evolution of the vortex. 

*/

event images (t += 10){
    // vorticity is computed there:
    scalar omega[];
    vorticity (u, omega);

    static FILE * fpa = fopen ("out_adapt.ppm", "w");
    output_ppm (omega, fpa, linear = true);

    foreach()
        omega[] = level;
    static FILE * fpg = fopen ("grid.ppm", "w");
    output_ppm (omega, fpg, min = 0, max = MAX_LEVEL);
}





/**
# Towards self-similarity

## The analytical solution (in cylindrical space)

The exact solution of the azimuth velocity for the *Lamb-Oseen* vortex is well-known^[[Wikipedia page (in german) of the *Lamb-Oseen* vortex: analytical solution](https://de.wikipedia.org/wiki/Hamel-Oseenscher-Wirbel#Umfangsgeschwindigkeit)] and is given by:

$$
\boxed{
u_\theta(r,t) = \frac{\Gamma}{2 \pi r} \left(1 - \mathrm{e}^{-r^2 / 4\nu t} \right )  
}
$$

where $\Gamma$ is the circulation of the flow.

Indeed, one can actually observe that:

$$
u_\theta^{sim}(r,t) = \frac{\Gamma}{r} f \left( \frac{r}{\sqrt{\nu t }}  \right) = \frac{\Gamma}{r} f \left( \eta  \right) 
\quad , \quad \eta = \frac{r}{\sqrt{\nu t }} 
\quad\quad\quad\quad (1)
$$ 

if injected in the simplified *Navier-Stokes* equations given at the beginning of this page, leads to the far more simple ODE:

$$
f'' + \left(\frac{\eta}{2} - \frac{1}{\eta}\right)f' = 0
$$

whose solution is the boxed one. 


In the following section will be discussed how was obtained the formula $(1)$, but for the sake of brevity, the method used will be directly applied in cartesian coordinates, as `Basilisk` works in such space.





## The Scale-Invariance Method directly applied in cartesian coordinates 


## Outputs in the self-similar space


~~~gnuplot Vertical Velocity Profiles in the Self-Similar Space 
reset

array times[16]
times[1] = 0
times[2] = 20
times[3] = 40
times[4] = 60
times[5] = 80
times[6] = 100
times[7] = 200
times[8] = 300
times[9] = 400
times[10] = 500
times[11] = 600
times[12] = 700
times[13] = 800
times[14] = 900
times[15] = 1000
times[16] = 3000
t0 = -2.5 # time of singularity

set xlabel 'xi = x / sqrt(nu * (t - t_0))'
set ylabel 'u_y sqrt(nu * (t - t_0)) / Gamma'
set xrange [-20:20]
set key left top # legend location

f_y(x) = (1. / (2. * pi)) * (x/(x * x + 0 * 0)) * (1. - exp(-(x * x + 0 * 0)/4.))

plot for [i=3:16] sprintf('yprof_%i',times[i]) using ($1)/sqrt(1e-3*(times[i]-t0)):($2)*sqrt(1e-3*(times[i]-t0)) with lines title sprintf('t = %i', times[i]), \
f_y(x)*0.627064 dt 2 lc rgb "red" lw 2 title "f_y(xi, 0)
# 0.627064 = GAMMA circulation value linked to the initial guess for 
# EPSILON = 1.e-2 & CORE_SIZE 1.e-1
~~~


## But what is this $0.627064$ thing in the `gnuplot` script ?? An explanation.




*/


/**
# APPENDIX: Convergence Studies

*/

















