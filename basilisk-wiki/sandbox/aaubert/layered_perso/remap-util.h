/**
This file implement the Piecewise Polynomial Reconstruction (PPR) used
by the remap function for the multilayer. It is mostly based on the
[Fortran code](/src/ppr) that was previously used by Basilisk to do
this remap. As the name says, it is based on a piecewise polynomial
reconstruction of the function on each interval to be able to compute
the mean value on the new mesh. We implement here a parabolic method
that consists of fitting order 2 polynomial on each interval.

The limiter used by the remap function */

double minmodremap (double a, double b) {
  return a*b <= 0. ? 0. :
    fabs(a) < fabs(b) ? a : b;
}

static void remap_central (double s_left, double s_right,
			   int npos, const double xpos[npos], const double fdat[npos],
			   int k,
			   double matrix_coeff[npos-1][3])
{
  double xk = xpos[k], dx = xpos[k+1] - xk;
  if ((fdat[k+1] - fdat[k])*(fdat[k] - fdat[k-1]) < 0.) {
    s_left = fdat[k];
    s_right = fdat[k];
  }
  else {
    // edge values limiter
    double sigma_left = 2.*(fdat[k] - fdat[k-1])/dx;
    double sigma_center = 2.*(fdat[k+1] - fdat[k-1])/(xpos[k+2] - xpos[k-1] + dx);
    double sigma_right = 2.*(fdat[k+1] - fdat[k])/dx;
    double sigma = minmodremap(sigma_center, minmodremap(sigma_left, sigma_right));
    if ((s_left - fdat[k-1])*(fdat[k] - s_left) < 0.)
      s_left = fdat[k] - 1./2.*dx*sigma;
    if ((s_right - fdat[k])*(fdat[k+1] - s_right) < 0.)
      s_right = fdat[k] + 1./2.*dx*sigma;
  }

  double coeff[3] = {
    3.*(s_right + s_left - 2.*fdat[k]),
    6.*fdat[k] - 2.*s_right - 4.*s_left,
    s_left
  };
  
  if ((s_left + s_right != 0.) && (fabs(s_left - s_right)/fabs(s_left + s_right) < 1.e-3)) {
    coeff[0] = 0.;
    coeff[1] = 0.;
    coeff[2] = fdat[k];
  }
  else if (coeff[0] != 0. && - coeff[1]/(2.*coeff[0]) > 1.e-10 && -coeff[1]/(2.*coeff[0]) < 1. - 1.e-10) {
    if (- coeff[1]/(2.*coeff[0]) < 0.5)
      s_right = 3.*fdat[k] - 2.*s_left;
    else
      s_left = 3.*fdat[k] - 2.*s_right;

    coeff[0] = 3.*(s_right + s_left - 2.*fdat[k]);
    coeff[1] = 6.*fdat[k] - 2.*s_right - 4.*s_left;
    coeff[2] = s_left;

    if (s_left + s_right != 0. && fabs(s_left - s_right)/fabs(s_left + s_right) < 1e-3) {
      coeff[0] = 0.;
      coeff[1] = 0.;
      coeff[2] = fdat[k];
    }
  }

  matrix_coeff[k][0] = coeff[0]/sq(dx);
  matrix_coeff[k][1] = (coeff[1] - 2.*coeff[0]*xk/dx)/dx;
  matrix_coeff[k][2] = coeff[2] - xk*(coeff[1] - coeff[0]*xk/dx)/dx;
}

static void remap3 (int npos, const double xpos[npos], const double fdat[npos], int k,
		    double matrix[3][3], double s_left, double s_right,
		    double matrix_coeff[npos-1][3])
{
  for (int j = 0; j < 3; j++)
    matrix[0][j] = 1./(3. - j);  //dzeta1=1 and dzeta0=0

  matrix_inverse3 (matrix);

  double b[3];
  b[0] = fdat[k];
  b[1] = s_left;
  b[2] = s_right;

  double coeff[3] = {0.,0.,0.};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      coeff[i] += matrix[i][j]*b[j];

  double xk = xpos[k], dx = xpos[k+1] - xk;
  matrix_coeff[k][0] = coeff[0]/sq(dx);
  matrix_coeff[k][1] = (coeff[1] - 2.*coeff[0]*xk/dx)/dx;
  matrix_coeff[k][2] = coeff[2] - xk*(coeff[1] - coeff[0]*xk/dx)/dx;
}

static void remap_bottom (double s_right,
			  int npos, const double xpos[npos], const double fdat[npos],
			  double f_b, double lambda_b,
			  double matrix_coeff[npos-1][3])
{
  const int k = 0;

  double matrix[3][3];
  matrix[1][0] = 0.;       //left condition
  matrix[1][1] = lambda_b;
  matrix[1][2] = 1.;

  matrix[2][0] = 1.;       //right condition
  matrix[2][1] = 1.;
  matrix[2][2] = 1.;

  remap3 (npos, xpos, fdat, k, matrix, f_b, s_right, matrix_coeff);
}

static void remap_top (double s_left,
		       int npos, const double xpos[npos], const double fdat[npos], int k,
		       double f_t, double lambda_t,
		       double matrix_coeff[npos-1][3])
{
  double matrix[3][3];
  matrix[1][0] = 0.;               //left condition
  matrix[1][1] = 0.;
  matrix[1][2] = 1.;

  matrix[2][0] = 1. + 2.*lambda_t; //right condition
  matrix[2][1] = 1. + lambda_t;
  matrix[2][2] = 1.;

  remap3 (npos, xpos, fdat, k, matrix, s_left, f_t, matrix_coeff);
}

void my_remap_c (int npos, int nnew,
		 const double xpos[npos], const double xnew[nnew],
		 const double fdat[npos], double fnew[nnew],
		 double f_b, double lambda_b, double f_t, double lambda_t)
{
  
  //boundary condition : bottom Navier with coefficient f_b, lambda_b
  //     top Navier with coefficient f_t, lambda_t

  // this matrix contains all the coefficients of the polynomial in the different intervals
  double matrix_coeff[npos - 1][3];

  /**  
  We first start by reconstructing the polynom on each interval. For
  this, we need the mean value of the function on this interval, and a
  left and a right value at the boundary. We thus need to solve the
  system (denoting $y_L$ and $y_R$ the left and right boundary $q$ the
  function that we need to remap and $Q$ the quadratic polynom that we
  want to fit)
  $$
  \begin{cases}
  Q(y_L)=q_L \\
  Q(y_R)=q_R \\
  \frac{1}{y_R-y_L}\int_{y_L}^{y_R}Q(y)dy=\bar{q}
  \end{cases}
  $$
  
  We create a local variable $\xi$ to normalize for every interval width. We choose 
  $$
  y=y_L+\xi\left(y_R-y_L\right)
  $$
  such that $\xi_L=0$ and $\xi_R=1$ (nota : in the Fortran code, they choose $\xi_L=-1$).
  
  This modify the polynom $Q(y)$ into a polynom $\tilde{Q}(\xi)$ with
  the new coefficient linked to those of $Q$
  $$
  Q(y)=ay^2+y+c ~,~~\tilde{Q}(\xi)=\tilde{a}y^2+\tilde{b}y+\tilde{c} \\
  a=\frac{\tilde{a}}{\left(y_R-y_L\right)^2} \\
  b=\frac{\tilde{b}}{y_R-y_L}-\frac{2\tilde{a}y_L}{\left(y_R-y_L\right)^2} \\
  c=\tilde{c}-\frac{\tilde{b}y_L}{y_R-y_L}+\frac{\tilde{a}y_L^2}{\left(y_R-y_L\right)^2}
  $$
  
  We thus need to solve the system
   $$
  \begin{cases}
  \tilde{Q}(0)=q_L \\
  \tilde{Q}(1)=q_R \\
  \int_{0}^{1}\tilde{Q}(\xi)dy=\bar{q}
 \end{cases}
  $$
  
  We obtain the linear system
  $$
  \begin{cases}
  \tilde{c}=q_L \\
  \tilde{a}+\tilde{b}+\tilde{c}=q_R \\
  \frac{\tilde{a}}{3}+\frac{\tilde{b}}{2}+\tilde{c}=\bar{q}
  \end{cases}
  $$
  
  The solution of this system is simply given by
  $$
  \tilde{c}=q_L \\
  \tilde{b}=6\bar{q}-4q_L-2q_R \\
  \tilde{a}=3\left(q_L+q_R-2\bar{q}\right)
  $$
  
  To conclude the reconstruction part, we need to construct the edge
  estimate $q_L$ and $q_R$. Close to the boundary, we can directly
  used the value of the boundary condition (this might change the
  linear system if the condition is not a dirichlet). To preserve the
  order 2 approximation, we need to use at least an order 3 method to
  construc this edge value and we will use a polynom of order 3 to
  interpolate the edge value from the neighbour mean values.
  
  We are thus looking for a polynom $P(y)$ such that
  $$
  \begin{cases}
  \frac{1}{y_{-1}-y_{-2}}\int_{y_{-2}}^{y_{-1}}P(y)dy=\bar{q}_{-2} \\
  \frac{1}{y_{0}-y_{-1}}\int_{y_{-1}}^{y_{0}}P(y)dy=\bar{q}_{-1} \\
  \frac{1}{y_{1}-y_{0}}\int_{y_{0}}^{y_{1}}P(y)dy=\bar{q}_{0}  \\
  \frac{1}{y_{2}-y_{1}}\int_{y_{1}}^{y_{2}}P(y)dy=\bar{q}_{1}
  \end{cases}
  $$
  for the edge value at $y=y_0$.
 
  As previously, we normalize the coordinate by introducing $\xi$ such
  that $y=y_0+\xi\left(y_1-y_0\right)$ and looking for a polynom
  $\tilde{P}(\xi)$ where the new coefficients are linked with the
  previous ones
  $$
  Q(y)=ay^3+by^2+cy+d ~,~~\tilde{Q}(\xi)=\tilde{a}y^3+\tilde{b}y^2+\tilde{c}y+\tilde{d} \\
  a=\frac{\tilde{a}}{\left(y_1-y_0\right)^3} \\
  b=\frac{\tilde{b}}{\left(y_1-y_0\right)^2}-\frac{3\tilde{a}y_0}{\left(y_1-y_0\right)^3} \\
  c=\frac{\tilde{c}}{y_1-y_0}-\frac{2\tilde{b}y_0}{\left(y_R-y_L\right)^2}+\frac{3\tilde{a}y_0^2}{\left(y_R-y_L\right)^3} \\
  d=\tilde{d}-\frac{\tilde{c}y_0}{y_1-y_0}+\frac{\tilde{b}y_0^2}{\left(y_1-y_0\right)^2}-\frac{\tilde{a}y_0^3}{\left(y_1-y_0\right)^3}
  $$
  
  The system that we need to solve is
  $$
  \begin{cases}
  \frac{y_1-y_0}{y_{-1}-y_{-2}}\int_{\xi_{-2}}^{\xi_{-1}}\tilde{P}(\xi)dy=\bar{q}_{-2} \\
  \frac{y_1-y_0}{y_{0}-y_{-1}}\int_{\xi_{-1}}^{\xi_{0}}\tilde{P}(\xi)dy=\bar{q}_{-1} \\
  \frac{y_1-y_0}{y_{1}-y_{0}}\int_{\xi_{0}}^{\xi_{1}}\tilde{P}(\xi)dy=\bar{q}_{0}  \\
  \frac{y_1-y_0}{y_{2}-y_{1}}\int_{\xi_{1}}^{\xi_{2}}\tilde{P}(\xi)dy=\bar{q}_{1}
  \end{cases}
  $$
 
  Which give the linear system
  $$
  \begin{cases}
  \frac{y_1-y_0}{y_{-1}-y_{-2}}\left[\tilde{a}\frac{\xi_{-1}^4-\xi_{-2}^4}{4}+\tilde{b}\frac{\xi_{-1}^3-\xi_{-2}^3}{3}+\tilde{c}\frac{\xi_{-1}^2-\xi_{-2}^2}{2}+\tilde{d}\left(\xi_{-1}-\xi_{-2}\right)\right]=\bar{q}_{-2}\\
  \frac{y_1-y_0}{y_{0}-y_{-1}}\left[\tilde{a}\frac{\xi_{0}^4-\xi_{-1}^4}{4}+\tilde{b}\frac{\xi_{0}^3-\xi_{-1}^3}{3}+\tilde{c}\frac{\xi_{0}^2-\xi_{-1}^2}{2}+\tilde{d}\left(\xi_{0}-\xi_{-1}\right)\right]=\bar{q}_{-1}\\
  \frac{y_1-y_0}{y_{1}-y_{0}}\left[\tilde{a}\frac{\xi_{1}^4-\xi_{0}^4}{4}+\tilde{b}\frac{\xi_{1}^3-\xi_{0}^3}{3}+\tilde{c}\frac{\xi_{1}^2-\xi_{0}^2}{2}+\tilde{d}\left(\xi_{1}-\xi_{0}\right)\right]=\bar{q}_{0} \\
  \frac{y_1-y_0}{y_{2}-y_{1}}\left[\tilde{a}\frac{\xi_{2}^4-\xi_{1}^4}{4}+\tilde{b}\frac{\xi_{2}^3-\xi_{1}^3}{3}+\tilde{c}\frac{\xi_{2}^2-\xi_{1}^2}{2}+\tilde{d}\left(\xi_{2}-\xi_{1}\right)\right]=\bar{q}_{1}
  \end{cases}
  $$
 
  The boundary conditions for the second layer (and the second to
  last) will change this system a bit.
  
  Finally, we apply a limiter for this edge values.
  
  First, if
  $\left(\bar{q}_i-\bar{q}_{i-1}\right)\left(\bar{q}_i-\bar{q}_{i+1}\right)<0$,
  meaning that we are at an extrema, we take a constant polynom.
  
  Second, if the edge estimate create an extrema, we interpolate this
  edge value linearly using a slope $\sigma$ computed from
  $$
  \sigma=minmod(\sigma_C,minmod(\sigma_L,\sigma_R))
  $$
  where $\sigma_C$, $\sigma_L$ and $\sigma_R$ are respectively the
  centered, left and right estimated slope
  $$
  \sigma_C=2\frac{\bar{q}_{i+1}-\bar{q}_{i-1}}{h_{i_1}+2h_i+h_{i+1}} ~~ 
  \sigma_L=2\frac{\bar{q}_i-\bar{q}_{i_1}}{h_i} ~~ \sigma_R=2\frac{\bar{q}_{i+1}-\bar{q}_i}{h_i}
  $$
  and the limiter is
  $$
  minmod(a,b)=\begin{cases}
  a~\text{if}~ab>0~\text{and}~|a|<|b| \\
  b~\text{if}~ab>0~\text{and}~|b|<|a| \\
  0~\text{otherwise}
  \end{cases}
  $$
  
  Lastly, we do not want to add a new extrema in the interval we are
  reconstructing. For a second order polynom, the extrema (if
  $\tilde{a}\neq0$) is at
  $\xi_{ext}=-\frac{\tilde{b}}{2\tilde{a}}$. With the value of the
  coefficient of the polynom $\tilde{Q}$, we have
  $$
  \xi_{ext}=-\frac{6\bar{q}-2q_R-4q_L}{6q_R+6q_L-12\bar{q}}
  $$
  
  If $\xi_{ext}\in[0,0.5]$, we change $q_R$ such that
  $\xi_{ext}=0$. The new value is then $q_R=3\bar{q}-2q_L$.
  
  If $\xi_{ext}\in[0.5,1]$, we change $q_L$ such that
  $\xi_{ext}=1$. The new value is then $q_L=3\bar{q}-2q_R$.
 
  Note that the different limiters will break the continuity of the
  reconstruction. */

  double s_left = 0;
  for (int k = 0; k < npos - 2; k++) {
    double s_right, matrix[4][4], b[4];

    double xk = xpos[k], dx = xpos[k+1] - xk;
    for (int i = 0; i < 4; i++) {
      double dzeta1 = (xpos[i+k-1] - xk)/dx, dz1n = dzeta1;
      double dzeta0 = (xpos[i+k] - xk)/dx, dz0n = dzeta0;
      for (int j = 3; j >= 0; j--, dz1n *= dzeta1, dz0n *= dzeta0)
	matrix[i][j] = 1./(4. - j)*(dz0n - dz1n);
      b[i] = (xpos[k+i] - xpos[k-1+i])*fdat[k-1+i]/dx;
    }
    
    if (k == 0) { // bottom boundary condition
      matrix[0][0] = 0.;  // take into account the boundary condition
      matrix[0][1] = 0.;  // instead of the integral at the layer -1
      matrix[0][2] = lambda_b;  // at this position, dzeta=0
      matrix[0][3] = 1.;
      b[0] = f_b;
    }
    else if (k == npos - 3) { // top boundary condition
      double dzetab = (xpos[k+2] - xk)/dx;
      matrix[3][0] = cube(dzetab) + 3.*sq(dzetab)*lambda_t;  // take into account the boundary condition
      matrix[3][1] = sq(dzetab) + 2.*dzetab*lambda_t;  // instead of the integral at the layer -1
      matrix[3][2] = dzetab + lambda_t;  // at this position, dzeta=(x_0-x_1)/(x_2-x_1)
      matrix[3][3] = 1.;
      b[3] = f_t;
    }
      
    assert (smatrix_inverse (4, matrix, 1.e-30) > 1.e-20);
    double coeff[4] = {0.,0.,0.,0.};
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
	coeff[i] += matrix[i][j]*b[j];
    s_right = coeff[3];
    for (int i = 0; i < 3; i++)
      s_right += coeff[i];   //right dzeta=1

    if (k == 0)
      remap_bottom (s_right, npos, xpos, fdat, f_b, lambda_b, matrix_coeff);
    else
      remap_central (s_left, s_right, npos, xpos, fdat, k, matrix_coeff);
    
    s_left = s_right;    
  }
  remap_top (s_left, npos, xpos, fdat, npos - 2, f_t, lambda_t, matrix_coeff);

  /**  
  We now need to compute the mean value on the new interval using the
  reconstruction previously obtained. We start from the bottom and
  work all the way to the top. We first fill the integral on each
  interval and we divide by the height at the end. */

  int inew = 0, ipos=0;
  double xcur = xpos[0], xmax = xnew[nnew - 1];

  for (int k = 0; k < nnew - 1; k++)
    fnew[k] = 0.;
  
  while ((inew < nnew - 1) && (xcur < xmax)) {
    if (ipos < npos - 1 && xpos[ipos + 1] <= xcur)
      ipos++;
    int j = ipos, jn = inew;
    double x1 = xpos[ipos + 1] < xnew[inew + 1] ? xpos[++ipos] : xnew[++inew];
    double xcurn = xcur, xpn = x1;
    for (int i = 2; i >= 0; i--, xcurn *= xcur, xpn *= x1)
      fnew[jn] += matrix_coeff[j][i]/(3. - i)*(xpn - xcurn);
    xcur = x1;
    if (ipos >= npos-1 && inew < nnew-1) {
      //if xnew[nnew-1]>xpos[npos-1], we need to finish with the last information that we have
      double xcurn = xcur, x1 = xnew[inew+1], xpn = x1;
      for (int i = 2; i >= 0; i--, xcurn *= xcur, xpn *= x1)
        fnew[inew] += matrix_coeff[ipos-1][i]/(3. - i)*(xpn - xcurn);
      xcur = xnew[++inew];
    }
  }

  for (int k = 0; k < nnew - 1; k++)
    fnew[k] /= xnew[k+1] - xnew[k];
}
