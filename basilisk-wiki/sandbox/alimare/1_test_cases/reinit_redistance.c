/**
#LS_reinit() test case

This case is extracted from [Russo et al.,1999](#russo_remark_2000) we
initialize a perturbed distance field, where the zero level-set is an ellipse of
the form:
$$
\phi (x,y,0) = f(x,y) \times g(x,y)
$$
where the perturbation is :
$$
f(x,y) = \epsilon  + (x - x_0)^2 +(y - y_0)^2
$$
and the ellipse is :
$$
g(x,y) = \left( \sqrt{\frac{x^2}{A^2}+\frac{y^2}{B^2}} -R \right)
$$
with $A=2$, $B=1$, $R = 1$ , $x_0 = 3.5$, $y_0 = 2.$.

We want to recover a perfect distance field, \textit{i.e.} remove the initial
perturbation.

*/

#include "../../popinet/distance_point_ellipse.h"
#include "../../popinet/redistance.h"
#include "../alex_functions.h"
#include "../basic_geom.h"
#include "view.h"

double perturb (double x, double y, double eps, coord center){
  return eps + sq(x - center.x) + sq(y - center.y);
}

void draw_isolines(scalar s, double smin, double smax, int niso, int w){
  vertex scalar vdist[];
  foreach_vertex()
  vdist[] = interpolate (s, x , y);
  boundary ({vdist});
  for (double sval = smin ; sval <= smax; sval += (smax-smin)/niso){
    isoline ("vdist", sval, lw = w);
  }
}


#define Pi 3.141592653589793
double eval_err_LS(double A, double B, int nb_point, scalar s){

/**
We will evaluate the flux on the theoretical 0-level-set of our
function with $N_\Sigma$ point $\mathbf{x_p}=(x_p, y_p)$ exactly on the
interface and defined such that:
$$
x_p = A \cos(\frac{2\pi p}{N_\Sigma}) \text{\hspace{1cm}} y_p = B \sin(\frac{2\pi
p}{N_\Sigma}) 
$$
Therefore, we want this function to calculate :
$$
E_\Sigma = \int_\Sigma |\phi(\mathbf(x),t)|ds
$$
that we discretize into:
$$
E_\Sigma = \frac{1}{2}\sum_{p=1}^{N_\Sigma}\left( |\tilde\phi(\mathbf
{x_p},t)| + |\tilde\phi(\mathbf{x_{p+1}},t)|\right) \left| \mathbf{x_{p+1}} - \mathbf{x_{p}} \right|
$$
Where $\tilde\cdot$ is the bilinear interpolation of $\phi$ at the desired
position, ideally it should be at least one order above the desired order of
accuracy of our method. In our case, this would require a biquadratic
interpolation or a hermitian one.

*/
  double sum = 0.;
  for(int i = 1; i<=nb_point; i++){
// 1rst point
    coord pLS = {A*cos(2.*Pi*i/nb_point), B*sin(2.*Pi*i/nb_point)};
    Point p = locate (pLS.x, pLS.y);
    int n = 1 << p.level;
    coord ccenter;
    ccenter.x = X0 + (p.i - 1) * L0/n;
    ccenter.y = X0 + (p.i - 1) * L0/n;

    int Stencil[2];
    double coeff[4];
    coord p_interp = {pLS.x - ccenter.x, pLS.y - ccenter.y};
    InterpStencil(p_interp, Stencil);

    double temp = fabs(mybilin(p, s, Stencil, p_interp, coeff));

    coord pLS2 = {A*cos(2.*Pi*(i+1)/nb_point), B*sin(2.*Pi*(i+1)/nb_point)};
    p = locate (pLS2.x, pLS2.y);
    n = 1 << p.level;
    ccenter.x = X0 + (p.i - 1) * L0/n;
    ccenter.y = X0 + (p.i - 1) * L0/n;
    p_interp.x = pLS2.x - ccenter.x;
    p_interp.y = pLS2.y - ccenter.y;
    InterpStencil(p_interp, Stencil);
    double temp2 = fabs(mybilin(p, s, Stencil, p_interp, coeff));

    sum += (temp + temp2)*sqrt( sq(pLS.x-pLS2.x) + sq(pLS2.y - pLS2.y))/2.;
  }
  return sum;
}


scalar dist[];
scalar * level_set = {dist};

int main() {

  origin (-5., -5.);
  L0 = 10;

  int MAXLEVEL = 6;
  for(MAXLEVEL = 6; MAXLEVEL < 9; MAXLEVEL++){
    init_grid (1 << MAXLEVEL);  
  
    double A = 4., B=2.;
    coord  center_perturb = {-3.5,-2.};
    foreach(){
      double a,b;
      dist[] = DistancePointEllipse(A,B,x,y,&a, &b)*
      perturb(x,y, 0.1, center_perturb);
    }
    boundary({dist});
  
    squares ("dist", map = cool_warm, min = -2, max = 2);
    draw_isolines(dist, -2., 2., 20, 1);
    save("dist_init.png");
  
    redistance(dist);
    
    scalar err[],Logerr[];
    foreach(){
      double a,b;
      err[] = fabs(DistancePointEllipse(A,B,x,y,&a, &b)-dist[]);
      Logerr[] = log(err[]/fabs(dist[]+1.e-10));
    }
    boundary({err});
    norm n = normf(err);
  
    double err2 = eval_err_LS(A, B, 200, dist);
    fprintf(stderr, "%d %g %g %g %g\n",1<<MAXLEVEL, n.avg, n.rms,n.max,
      err2*L0/(1 << MAXLEVEL));



    if(MAXLEVEL == 8){
      squares ("Logerr", map = cool_warm, min = -10, max = -4);
      save("err.png");
    }
  }

  squares("dist", map = cool_warm, min = -2, max = 2);
  draw_isolines(dist, -2., 2., 20, 1);
  save("dist_final.png");
  exit(1);
}


/**
We show here the initial and final level-set for the same isovalues.

![Initial level-set](reinit_redistance/dist_init.png) 

![Final level-set](reinit_redistance/dist_final.png)

We plot now:
$$\frac{|dist_{theo} - dist|}{|dist|}$$

![Relative error. Logscale between $10^{-10}$ and $10^{-4}$](reinit_redistance/err.png)

~~~gnuplot error analysis
f(x) = a + b*x
f1(x) = a1 + b1*x
f2(x) = a2 + b2*x
fit f(x) 'log' u (log($1)):(log($2)) via a,b
fit f1(x) 'log' u (log($1)):(log($3)) via a1,b1
fit f2(x) 'log' u (log($1)):(log($4)) via a2,b2
ftitle(a,b) = sprintf("%.3f-%4.3f x", a, -b)
set logscale xy
set xrange [32:512]
set xtics 32,2,512
plot 'log' u 1:2 t 'avg', exp(f(log(x))) t ftitle(a,b), \
    'log' u 1:3 t 'rms', exp(f1(log(x))) t ftitle(a1,b1), \
    'log' u 1:4 t 'max', exp(f2(log(x))) t ftitle(a2,b2)
~~~
We have first order convergence globally on the distance function.

~~~gnuplot Error on the 0-level-set.
unset logscale
unset xrange
f3(x) = a3 + b3*x
fit f3(x) 'log' u (log($1)):(log($5)) via a3,b3
ftitle(a,b) = sprintf("%.3f-%4.3f x", a, -b)
set logscale xy
set xrange [32:512]
set xtics 32,2,512
plot 'log' u 1:5 t '0-level-set error', exp(f3(log(x))) t ftitle(a,b)
~~~
Here we study the value of the level-set function on a set of points where it is
theoretically 0, we show that we have also a 1rst-order convergence.

## References

~~~bib

@article{russo_remark_2000,
  title = {A remark on computing distance functions},
  volume = {163},
  number = {1},
  journal = {Journal of Computational Physics},
  author = {Russo, Giovanni and Smereka, Peter},
  year = {2000},
  pages = {51--67}
}

~~~
*/