/**
# Reinitialization of the level-set function

Redistancing function with subcell correction, see the work of [Russo & Smereka,
1999](#russo_remark_2000) with corrections by [Min & Gibou](#Min2007) and by 
[Min](#Min2010).

Let $\phi$ be a function close to a signed function that has been perturbed by
numerical diffusion (more precisely, a non-zero tangential velocity). By
iterating on this eikonal equation :
$$
\left\{\begin{array}{ll}
\phi_t + sign(\phi^{0}) \left(\left| \nabla \phi\right| - 1 \right) = 0\\ 
\phi(x,0) = \phi^0(x)
\end{array}
\right.
$$
we can correct or redistance $\phi$ to make it a signed function.

We use a [Godunov Hamiltonian](alex_functions.h#godunov-hamiltonian) approximation for
$\left| \nabla \phi\right|$:
$$
\left| \nabla \phi \right|_{ij} = H_G(D_x^+\phi_{ij}, D_x^-\phi_{ij}, D_y^+\phi_{ij}, D_y^-\phi_{ij})
$$
where $D^\pm\phi_{ij}$ denotes the one-sided ENO difference finite difference in
the x- direction:
$$
D_x^+ = \dfrac{\phi_{i+1,j}-\phi_{i,j}}{\Delta} - \dfrac{\Delta}{2}minmod(D_{xx}\phi_{ij}, D_{xx}\phi_{i+1,j})
$$
$$
D_x^- = \dfrac{\phi_{i,j}-\phi_{i-1,j}}{\Delta} + \dfrac{\Delta}{2}minmod(D_{xx}\phi_{ij}, D_{xx}\phi_{i+1,j})
$$
here $D_{xx}\phi_{ij}= (\phi_{i-1,j} - 2\phi{ij} + \phi_{i+1,j})/\Delta^2$.

The minmod function is zero when the two arguments have different signs, and
takes the argument with smaller absolute value when the two have the same
sign.

The Godunov Hamiltonian $H_G$ is given as:
$$
H_G(a,b,c,d) = \left\{ \begin{array}{ll}
\sqrt{max((a^-)^2,(b^+)^2 + (c^-)^2,(d^+)^2)} \text { when } sgn(\phi^0_{ij})
\geq 0\\
\sqrt{max((a^+)^2,(b^-)^2 + (c^+)^2,(d^-)^2)} \text { when } sgn(\phi^0_{ij}) < 0
\end{array}
\right.
$$
with:
$$
x^+ = max(0, x)\\
x^- = min(0, x)\\
$$
 */

/**
We use a minmod and the second order approximation for the derivatives that I redefined in my sandbox [weno2.h]().
*/
#include "weno2.h"
//#include "alex_functions.h"
/**

## Precaculation of the inputs of the hamiltonian:
$$
D^+\phi_{ij} = \dfrac{\phi_{i+1,j} - \phi_{ij}}{\Delta x} - \dfrac{\Delta}
{2}minmod(D_{xx}\phi_{ij},D_{xx}\phi_{i+1,j})
$$
$$
D^-\phi_{ij} = \dfrac{\phi_{i,j} - \phi_{i-1,j}}{\Delta x} - \dfrac{\Delta}
{2}minmod(D_{xx}\phi_{ij},D_{xx}\phi{i-1,j})
$$
where:
$$
D_{xx}\phi_{ij} = \dfrac{\phi_{i-1,j} - 2\phi_{ij} + \phi_{i+1,j}}{\Delta x^2}
$$
*/


void prehamil(Point point, coord  * grapl, coord * gramin, scalar s){
  foreach_dimension(){
    grapl->x  = WENOdiff_x(point,s,1);
    gramin->x = WENOdiff_x(point,s,-1);
  }
}

/**
## Godunov Hamiltonian
$$
H_G(a,b,c,d) = \left\{ \begin{array}{ll}
\sqrt{max((a^-)^2,(b^+)^2 + (c^-)^2,(d^+)^2)} \text { when } sgn(\phi^0_{ij})
\geq 0\\
\sqrt{max((a^+)^2,(b^-)^2 + (c^+)^2,(d^-)^2)} \text { when } sgn(\phi^0_{ij}) < 0
\end{array}
\right.
$$
*/

double hamiltonian (Point point, scalar s0, coord grapl, coord  gramin)
{
  double hamil = 0;
  if(s0[] > 0){
    foreach_dimension(){
      double a = min(0.,grapl.x); 
      double b = max(0.,gramin.x);
      hamil += max(sq(a),sq(b));
    }
    return sqrt(hamil);
  }
  else{
    foreach_dimension(){      
      double a = max(0.,grapl.x);
      double b = min(0.,gramin.x);
      hamil += max(sq(a),sq(b));
    }
  }
  return sqrt(hamil);
}

/**
$$
My function
$$
*/

double sign2 (double x) {
  return(x > 0. ? 1. : x<0 ? -1. : 0.);
}

/**
## Root extraction for the subcell fix near the interface:
$$
\Delta x^+ = \left\{ \begin{array}{ll}
\Delta x \cdot \left( \dfrac{\phi^0_{i,j}-\phi^0_{i+1,j}-sgn(\phi^0_{i,j}-\phi^0_{i+1,j})\sqrt{D}}{}\right) \text{ if } \left| \phi^0_{xx}\right| >\epsilon \\
\Delta x \cdot \dfrac{\phi^0_{ij}}{\phi^0_{i,j}-\phi^0_{i+1,j}} \text{ else.}\\
\end{array}
\right.
$$
with:
$$
\phi_{xx}^0 = minmod(\phi^0_{i-1,j}-2\phi^0_{ij}+\phi^0_{i+1,j}, \phi^0_{i,j}-2\phi^0_{i+1j}+\phi^0_{i+2,j}) \\
D = \left( \phi^0_{xx}/2  - \phi_{ij}^0 - \phi_{i+1,j} \right)^2  - 4\phi_{ij}^0\phi_{i+1,j}^0
$$

For the $\Delta x^-$ calculation, replace all the $+$ subscript by $-$, this
is dealt with properly with the `dir` variable in the following function.
*/

foreach_dimension()
static inline double root_x(Point point, scalar s, double eps, int dir)
{
  // dir == 1 or -1 offsets the position of the interface
  double phixx = my_minmod(s[2*dir] + s[] - 2*s[dir], 
                           s[1] + s[-1] - 2*s[]);
  if(fabs(phixx) > eps){
    double D = sq(phixx/2.-s[] - s[dir])-4*s[]*s[dir];
    // fprintf(stderr, "%g %g %g\n", D, phixx, eps);
    return 1/2.+( s[] - s[dir] - sign2(s[] - s[dir])*sqrt(D))/phixx;
  }
  else{
    return s[]/(s[]- s[dir]);
  }
}

/**
## Forward Euler Integration

Simple Euler integration for the LS_reinit() function.

*/

double ForwardEuler(scalar dist, scalar temp, scalar dist0, double dt){
  double res=0.;

  foreach(reduction(max:res)){
    double delt =0.;
/**
Near the interface, *i.e.* for cells where:
$$
\phi^0_i\phi^0_{i+1} \leq 0 \text{ or } \phi^0_i\phi^0_{i-1} \leq 0
$$

The scheme must stay truly upwind, meaning that the movement of the 0
level-set of the function must be as small as possible.

The cells which contain the interface are tagged with `flag` ($<0$ for
interfacial cells). 
*/

    double flag = 1.;
    foreach_dimension(){
      flag = min (flag, dist0[-1]*dist0[]);
      flag = min (flag, dist0[ 1]*dist0[]);
    }

    coord grapl, gramin;
    prehamil(point, &grapl, &gramin, temp);

    if(flag < 0.){ // the cell contains the interface

/**
Near the interface, *i.e.* for cells where:
$$
\phi^0_i\phi^0_{i+1} \leq 0 \text{ or } \phi^0_i\phi^0_{i-1} \leq 0
$$

The scheme must stay truly upwind, meaning that the movement of the 0
level-set of the function must be as small as possible. Therefore the upwind
numerical scheme is modified to:

$$
D_x^+ = \dfrac{0-\phi_{ij}}{\Delta x^+} - \dfrac{\Delta x^+}{2} minmod(D_
{xx}\phi_{ij},D_{xx}\phi_{i+1,j}) \text{ if } \phi_{ij}\phi_{i+1,j} < 0
$$
$$
D_x^- = \dfrac{\phi_{ij}-0}{\Delta x^-} + \dfrac{\Delta x^-}{2} minmod(D_
{xx}\phi_{ij},D_{xx}\phi_{i-1,j}) \text{ if } \phi_{ij}\phi_{i+1,j} < 0
$$
correction by Min & Gibou 2006.
*/
      double size = 1.e10;
      foreach_dimension(){
        if(dist0[]*dist0[1]<0){
          double dx = Delta*root_x(point, dist0, 1.e-10, 1);
          double sxx1 = (temp[2] + temp[] - 2*temp[1])/sq(Delta);
          double sxx2 = (temp[1] + temp[-1] - 2*temp[])/sq(Delta);
          if(dx !=0.)
            grapl.x = -temp[]/dx - dx* my_minmod(sxx1,sxx2)/2.;
          else 
            grapl.x = 0.;
          size = min(size, dx);
        }
        if(dist0[]*dist0[-1]<0){
          double dx = Delta*root_x(point, dist0, 1.e-10, -1);
          // if(dx>10.){
            // fprintf(stderr, "%g %g %g\n", dist0[],dist0[-1,0],dx);
          // }
          double sxx2 = (temp[1] + temp[-1] - 2*temp[])/sq(Delta);
          double sxx3 = (temp[-2] + temp[0] - 2*temp[-1])/sq(Delta);
          if(dx!=0.)
            gramin.x = temp[]/dx + dx* my_minmod(sxx3,sxx2)/2.;
          else 
            gramin.x = 0.;
          size = min(size, dx);
        }
      }
      delt = sign2(dist0[]) * min(dt,fabs(size)/2.) *
      (hamiltonian(point, dist0, grapl,gramin) - 1);
      dist[] -= delt;
    }
    else{ 

/**
Far from the interface, we use simply the Hamiltonian defined earlier.
*/
      delt = sign2(dist0[]) * 
      (hamiltonian(point, dist0, grapl, gramin)- 1);
      dist[] -= dt*delt;
    }
    res = max (res,fabs(delt));
  }


  boundary({dist});
  restriction({dist});

  return res;
}




struct LS_reinit {
  scalar dist;
  double dt;
  int it_max;
};

/**
## LS_reinit() function
*/
int LS_reinit(struct LS_reinit p){
  scalar dist = p.dist; 

  double dt   = p.dt;     // standard timestep (0.5*Delta)
  int it_max  = p.it_max;// maximum number of iteration (100)

/**
In 2D, if no specific timestep is set up by the user, we take the most
restrictive one
with regards to the CFL condition :
$$
\Delta t = 0.5 * \Delta x
$$
*/
  
  //if(dt == 0) dt = 0.5 * L0/(1 << grid->maxdepth);
  if(dt == 0) dt = 0.1 * L0/(1 << grid->maxdepth);


/**
Default number of iterations is 20 times, which is sufficient to have the first
10 neighbor cells to the 0-level-set properly redistanced.
*/

  if(it_max == 0)it_max = 100;

  vector gr_LS[];
  int i ;

/**
Convergence is attained is residual is below $dt\times 10^{-6}$
*/  
  double eps = dt*1.e-6;

/**
We create `dist0[]` which will be a copy of the initial level-set function
before the iterations and `temp[]` which will be $\phi^{n}$ used for the
iterations.
*/
  scalar dist0[];
  foreach(){
    dist0[] = dist[] ;
  }
  boundary({dist0});

/**
 Time integration iteration loop.

One can choose between Runge Kutta 2 and Forward Euler temporal integration.
*/
  for (i = 1; i<=it_max ; i++){
    double res = 0;

/**

## RK3
We use a Runge Kutta 3 compact version taken from [Shu and Osher](#Shu1988)
made of 3 backward Euler steps:

* Step1-2
$$
\frac{\widetilde{\phi}^{n+1}  - \phi^n}{\Delta t}  = \text{RHS}^n\\
\dfrac{\widetilde{\phi}^{n+2}  - \widetilde{\phi}^{n+1}}{\Delta t}  = \widetilde{RHS}^{n+1} 
$$
with :
$$
RHS =  sgn (\phi_{ij}^0)\cdot \left[ H_G\left( D_x^+\phi_{ij}^n, D_x^-\phi_{ij}^n, D_y^+\phi_{ij}^n,
D_y^-\phi_{ij}^n \right)\right]
$$
*/
    scalar temp[],temp1[], temp2[];
    foreach(){
      temp[] = dist[] ;
      temp1[] = dist[] ;
    }
    boundary({temp,temp1});
    ForwardEuler(temp1,temp,dist0,dt);
    foreach(){
      temp2[] = temp1[] ;
    }
    boundary({temp2});
    ForwardEuler(temp2,temp1,dist0,dt);
/**
* Intermediate value
$$
\widetilde{\phi}^{n+1/2}  = \dfrac{3}{4}\widetilde{\phi}^{n} + \dfrac{1}{4}\widetilde{\phi}^{n+2}
$$
*/
    foreach(){
      temp1[] = 3./4*dist[] + temp2[]/4.;
      temp2[] = temp1[];
    }
    boundary({temp1,temp2});

/**
* Step 3
$$
\widetilde{\phi}^{n+3/2} - \widetilde{\phi}^{n+1/2} = \widetilde{RHS}^{n+1/2}
$$
*/
    ForwardEuler(temp2,temp1,dist0,dt);
/**
* Final Value
$$
\widetilde{\phi}^{n+1} = \widetilde{\phi}^{n} + \dfrac{2}{3}\widetilde{\phi}^{n+3/2}
$$
*/
    foreach(reduction(max:res)){
      res = max(res, 2./3.*fabs(dist[] - temp2[]));
      dist[] = dist[]/3. + temp2[]*2./3.;
    }
    boundary({dist});
    restriction({dist});
/**
Iterations are stopped when $L_1 = max(|\phi_i^{n+1}-\phi_i^n|) < eps$
*/
    if(res<eps){
      return i;
    }
  }
  return it_max;
}


/**
## References

~~~bib

@Article{Shu1988,
  author        = {Chi-Wang Shu and Stanley Osher},
  title         = {Efficient implementation of essentially non-oscillatory shock-capturing schemes},
  year          = {1988},
  volume        = {77},
  pages         = {439-471},
  issn          = {0021-9991},
  __markedentry = {[limare:6]},
  doi           = {10.1016/0021-9991(88)90177-5},
}

@article{russo_remark_2000,
  title = {A remark on computing distance functions},
  volume = {163},
  number = {1},
  journal = {Journal of Computational Physics},
  author = {Russo, Giovanni and Smereka, Peter},
  year = {2000},
  pages = {51--67}
}

@article{Min2007,
  author        = {Chohong Min and Frédéric Gibou},
  title         = {A second order accurate level set method on non-graded adaptive cartesian grids},
  year          = {2007},
  volume        = {225},
  pages         = {300-321},
  issn          = {0021-9991},
  doi           = {10.1016/j.jcp.2006.11.034},
}

@article{Min2010,
  author        = {Chohong Min},
  title         = {On reinitializing level set functions},
  year          = {2010},
  volume        = {229},
  pages         = {2764-2772},
  issn          = {0021-9991},
  doi           = {10.1016/j.jcp.2009.12.032},
}

~~~
*/
