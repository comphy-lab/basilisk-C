//Pbl 7.1.2-4


/**We are interested by the (inner) solution of the boundary value problem $7.1.2-4$ of Polyanin's book page 474. It reads for $R\leq r < \infty$

$$
    w\left(r,\theta\right) = \frac{a_0}{2k} + \sum\limits_{n=1}^{\infty} \frac{R}{kR +n }\left(\frac{R}{r}\right)^n\bigg(a_n \cos(n\theta) + b_n\sin(n\theta)\bigg),  
$$
with 
$$
a_n = \frac{1}{\pi}\int\limits_0^{2\pi}f(\psi)\cos(n\psi)d \psi~~~n = 0,1,2,\cdots\\
b_n = \frac{1}{\pi}\int\limits_0^{2\pi}f(\psi)\sin(n\psi)d \psi~~~n = 1,2,3\cdots
$$

where the function $f(\psi)$ appearing above is chosen such that at the boundary at $r=R$ the following Robin boundary condition is imposed :
$$
k w + \partial_n w = k w +\nabla w \cdot n = f\left(\theta\right)
$$
Beware to the sign convention in front of the derivative term ! 
We solve the problem for this set of parameters : $k = 3$, $R=1.63568$ and $f(\theta) = cos(\theta)$. 
For this specific case, it comes $a_n = 0$, for $n = 0,.....$ except $n=1$ where $a_1=1$ and $b_n = 0$ for each $n$. Then the solution is $w(r,\theta)=\frac{R}{Rk+1}\left(\frac{R}{r}\right)cos(\theta)$
*/

#include "embed.h"
#include "poisson.h"
#include "run.h"
#include "view.h"

#define R0 1.63568
#define k 3.
#define r sqrt(sq(x)+sq(y))
#define f1(x,y) atan2(y,x)
#define circle(x, y, R) -(sq(x) + sq(y) - sq(R))


static
double dirichlet_homogeneous_bc (Point point, Point neighbor, scalar s, int * data) {
  return 0.;
}
double exact(double x, double y){
  return x/(k*R0+1);
}

scalar u_exact[],rhs[],u_cs[]=0.;

int NN=16;
int main(){
  for (N = 16; N <= 256; N *= 2) {
    NITERMAX=1000;
    L0=6.;
    origin(-L0/2.,-L0/2.);
    init_grid(N);
    solid(cs,fs,circle(x,y,R0));

#if TREE
    cs.refine = cs.prolongation = fraction_refine;
#endif
    restriction ({cs,fs});
    cm = cs;
    fm = fs;
    scalar u[];


    u[embed]=robin(k,1.,x/r);
    u.third = true;

    foreach() {
      u[]=0.;
      rhs[]=0.;
    }
    struct Poisson p;
    p.alpha = fs;
    p.lambda = zeroc;
    p.embed_flux = embed_flux;
    scalar res[],err[];
    double maxp = residual ({u}, {rhs}, {res}, &p), maxf = 0.;

    mgstats s = poisson (u, rhs,alpha=fs, tolerance = 1e-12, minlevel = 4);


    scalar e[], ep[], ef[];
    foreach() {
      if (cs[] == 0.)
	ep[] = ef[] = e[] = u_exact[] = nodata;
      else {
	u_exact[] = exact(x,y);
	e[] = u[] - u_exact[];

	ep[] = cs[] < 1. ? e[] : nodata;
	ef[] = cs[] >= 1. ? e[] : nodata;
      }
    }
    norm n = normf (e), np = normf (ep), nf = normf (ef);
    fprintf (stdout, "%d %.3g %.3g %.3g %.3g %.3g %.3g %d %d\n",
	     N, n.avg, n.max, np.avg, np.max, nf.avg, nf.max,s.i, s.nrelax);

    if (N==32){
	clear();
	/* view (tx = -0., ty = -0., width = 2000, height = 2000); */
	squares("u_exact", linear=false, map=jet);
	draw_vof ("cs",filled=0, fc = {0.,0.8,0.8}, lw = 2);
	save ("analyt.png");

	clear();
	squares("e", linear=false, map=jet);
	save ("error.png");

	clear();
	/* view (tx = -0., ty = -0., width = 2000, height = 2000); */
	squares("u", linear=false, map=jet);
	draw_vof ("cs",filled=0, fc = {0.,0.8,0.8}, lw = 2);
	save ("u.png");
	cells();
      }
      
  }
}

/** 
## Result

###
![Solution on the finest mesh](valid_robin_pbl7/u.png)

![Error on the finest mesh](valid_robin_pbl7/error.png)

Maximum residual convergence.

~~~gnuplot Maximum residual convergence
set xrange [*:*]
ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)
f(x) = a + b*x
fit f(x) '< grep maxres log' u (log($2)):(log($3)) via a,b
f2(x) = a2 + b2*x
fit f2(x) '' u (log($2)):(log($4)) via a2,b2
set xlabel 'Resolution'
set logscale
set cbrange [1:2]
set xtics 16,2,1024
set grid ytics
set ytics format "% .0e"
set xrange [16:1024]
set yrange [1e-6:]
set key bottom left
~~~

Error convergence.

~~~gnuplot Error convergence
set xrange [*:*]
ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)
f(x) = a + b*x
f2(x) = a2 + b2*x
fit f(x) 'out' u (log($1)):(log($3)) via a,b
fit f2(x) '' u (log($1)):(log($4)) via a2,b2
f3(x) = a3 + b3*x
fit f3(x) '' u (log($1)):(log($6)) via a3,b3
set xrange [16:1024]
set ylabel 'Error'
set yrange [*:*]
plot '' u 1:3 pt 6 t 'max all cells', exp(f(log(x))) t ftitle(a,b), \
     '' u 1:4 t 'avg partial cells', exp(f2(log(x))) t ftitle(a2,b2), \
     '' u 1:6 t 'avg full cells', exp(f3(log(x))) t ftitle(a3,b3)
~~~

*/
