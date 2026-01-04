/**
# Combined compact finite difference

We can solve two coulped implicit schemes for the first and second derivative of a signal. This can yield 12th-order accurate first and second derivatives on a 5-point stencil. Indeed, in the convergence zone, doubling the resoltution results in a $2^{12} = 4096$ fold reduction of the error. 

~~~gnuplot Error convergence for the 12th-order scheme and a 2nd order scheme. 
set term svg size 700,500
set size square
set xr[16:1024]
set yr[1e-15:10]
set grid
set logscale x 2
set logscale y
set key outside
set ylabel 'L2 error'
set xlabel 'N'
plot 'out' t '2nd-order 1st der', '' u 1:3 t '2nd-order 2nd der' , '' u 1:4 t '12th order 1st der', '' u 1:5 t '12th-order 2nd der', 1e-15*x t 'Scaling for binary limit', 1e15*x**(-12) w l lw 2 t '12th order convergence'
~~~

By adding more coupled equations, for even higher derivatives, an arbitrary order can be achieved.

Say, you are not interested in the logarithm of the error, one would plot their values on a linear axis.

~~~gnuplot Even on a linear axis the 12th-order scheme seems to compare well. It also better represents the fact that the schemes actually converge for high $N$
reset
set term svg size 700,500
set size square
set xr[30:600]
set yr[-0.1:1.1]
set grid
set logscale x 2
set key outside
set ylabel 'L2 error'
set xlabel 'N'
plot 'out' t '2nd-order 1st der' w l lw 2, '' u 1:4 w l lw 2 t'12th order 1st der' 
~~~
*/

int main() {
  periodic (left);
  L0 = 20;
  X0 = -L0/2.;
  for (N = 32; N <= 512; N *= 2) {
    init_grid (N);
    scalar s[], ds[], dds[];
    foreach() 
      s[] = exp(-sq(x));
    foreach() {
      ds[] = -(s[1] - s[-1])/(2*Delta);
      dds[] = (s[1] - 2*s[] + s[-1])/sq(Delta);
    }
    double err = 0, err2 = 0;
    foreach(reduction(+:err) reduction(+:err2)) {
      err += Delta*sq(-ds[] + 2*x*exp(-sq(x)));
      err2 += Delta*sq(dds[] - (4*sq(x) - 2)*exp(-sq(x)));
    }
    printf ("%d %g %g", N, sqrt(err), sqrt(err2));
  
    // See Octave script below for these coefficients
    for (int j = 0; j < 50; j++) {
      foreach() {
	ds[] = (+445./2592.*(s[2] - s[-2])/Delta
		- 100./81.*(s[1] - s[-1])/Delta
		+ 23./432.*(ds[2] + ds[-2])
		- 4./9.*(ds[1] + ds[-1])
		+ 1./216.*(dds[2] - dds[-2])*Delta
		- 4./27.*(dds[1] - dds[-1])*Delta); 
      }
      foreach() {
	dds[] = (-65/324.*(s[2] + s[-2])/sq(Delta)
		 + 320./81.*(s[1] + s[-1])/sq(Delta)
		 - 15./2.*(s[])/sq(Delta)
		 - 25/432.*(ds[2] - ds[-2])/Delta
		 + 40./27.*(ds[1] - ds[-1])/Delta
		 - 1./216.*(dds[2] + dds[-2])
		 + 8./27.*(dds[1] + dds[-1])); 
      }
    }
    err = 0, err2 = 0;
    foreach(reduction(+:err) reduction(+:err2)) {
      err += Delta*sq(-ds[] + 2*x*exp(-sq(x)));
      err2 += Delta*sq(dds[] - (4*sq(x) - 2)*exp(-sq(x)));
    }
    printf (" %g %g\n", sqrt(err), sqrt(err2));
  }
}

/**
Coefficient finder script (using Octave)

~~~octave
format rat

N = 13
n = [0 : N - 1];
% Taylor expansion arrays
f = 1./factorial(n);
d1 = (-1).^n;
d2 = 2.^n;
% f
B(13,:) = d2.*f;
B(12,:) = f;
B(11,:) = (n == 0)
B(10,:)  = d1.*f;
B(9,:)  = d1.*d2.*f;
% f'
B(8,2:end) = B(13,1:N-1);
B(7,2:end) = B(12,1:N-1);
B(6,2:end) = B(10,1:N-1);
B(5,2:end) = B(9,1:N-1);
$ f''
B(4,3:end) = B(13,1:N-2);
B(3,3:end) = B(12,1:N-2);
B(2,3:end) = B(10,1:N-2);
B(1,3:end) = B(9,1:N-2);

B = transpose(B);
Bi = inv(B);
coefs_eq1 = transpose(Bi(:,2))
coefs_eq2 = transpose(Bi(:,3))
~~~


*/