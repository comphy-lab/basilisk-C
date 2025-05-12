/**
# Vector induction by line

This is a solver for the equation:

$$\mathbf F(\mathbf{x}) = \int_C \mathrm{d}\mathbf{F} = \int_C
f(r)\mathrm{d}\mathbf{l} .$$

Where, $\mathbf{F}$ is a vector field, $C$ is a line or loop, and $f$
is a function that that takes as argument the distance of the point
$\mathbf{x}$ and the line segment with length and direction
$\mathrm{d}\mathbf{l}$.

## An example

Consider an infinite straight line from $\mathbf{x_1} = \{ 0,0,- \infty
\}$ to $\mathbf{x_2} = \{ 0,0,\infty \}$ and $f(r) =
\frac{\Gamma}{\pi^{3/2}a^3} e^{-\frac{r^2}{a^2}}$. Then, using the
relevant coordinates $\rho = \sqrt{x^2 + y^2} \text{\ and\ } z$

$$\mathbf{F}(\mathbf{x}) = \mathbf{F}(\rho) = \int_{-\infty}^{\infty}
\Gamma e^{-\frac{\rho^2 - z^2}{a^2}} d\mathrm{z} = $$ 
$$ \Gamma e^{-\frac{\rho^2}{a^2}}\int_{-\infty}^{\infty}e^{-\frac{z^2}{a^2}} \mathrm{d}z
\mathbf{e_z} = \frac{\sqrt{\pi}a}{\pi^{3/2}a^3}\Gamma e^{-\frac{\rho^2}{a^2}}\mathbf{e_z} = \\ 
\frac{\Gamma}{\pi a^2}e^{-\frac{\rho^2}{a^2}} \mathbf{e_z}.$$

Which is closely related to the Lamb-Oseen vortex or Gaussian vortex
of size $a$ and strength $\Gamma$.

## Implementation

We follow A. Castillo's filaments.h/c functions, and compute $F$ from
a discrete parameteric curve $\mathbf{C}(t)$ and its parametric
direction vector $\mathrm{d}\mathbf{l}(t)$. Noting that the magnitude
of the direction vector can be used to have a varying strength along
the line.
 */
trace
coord get_vorticity (coord X, int n_seg, coord * c, coord * dl, double (* fun)(double r)) {
  coord F = {0, 0, 0};
  for (int i = 0; i < n_seg; i++) {
    double r = 0;
    foreach_dimension()
      r += sq(X.x - c[i].x);
    double f = fun(r);
    foreach_dimension()
      F.x += dl[i].x*f;
  }
  return F;
}
trace
void get_vor_vector (vector omg, coord (* C)(double t), double t0,
		     double te, int n_seg, double (* fun)(double r)) {
  coord c[n_seg], dl[n_seg];
  double Tl = te - t0, Dt = Tl/(double)n_seg;
  int i = 0;
  double sqrt3 = sqrt(3);
  for (double ti = t0 + Dt/2; ti < te; ti += Dt) {
    c[i] = C(ti);
    coord d1 = C(ti - Dt/2), d2 = C(ti + Dt/2);
    foreach_dimension()
      dl[i].x = d2.x - d1.x;
    i++;     
  }
  foreach() {
    foreach_dimension()
      omg.x[] = 0.;
    coord cc = {x, y, z}, vort = {0, 0, 0};
    foreach_child() {
      coord ccc;
      foreach_dimension()
	ccc.x = cc.x + child.x*Delta/sqrt3;
      coord vor = get_vorticity (ccc, n_seg, c, dl, fun);
      foreach_dimension()
	vort.x += vor.x/(double)(1 << dimension);
    }
    foreach_dimension()
      omg.x[] = vort.x;
  }
  boundary ((scalar*){omg});
}

double Gamma_LO = 1, a_LO = 1;
double Lamb_Oseen (double sqr) {
  return Gamma_LO/(pow(pi, 3./2.)*cube(a_LO))*exp(-(sqr/sq(a_LO)));
}

/**
## Test

[Lamb Oseen vortex in 3D](test_filaments.c)
 */
