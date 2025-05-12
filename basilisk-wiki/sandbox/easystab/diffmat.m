%{

*This document belongs to the [easystab](http://basilisk.fr/sandbox/easystab/README) project, check the main page of the project to understand the general philosophy of the project.*

*The present program was specifically designed as a support for the course ["Introductions to hydrodynamical instabilities"](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) for the "DET" Master's cursus in Toulouse. For the underlying theory please see [Lecture notes for chapter 3](http://basilisk.fr/sandbox/easystab/NumericalMethodsForEigenvalueProblems.md)*

*To use this program: click "raw page source" in the left column, copy-paste in the Matlab/Octave editor, and run. You should have previously downloaded the package [easypack.zip](http://basilisk.fr/sandbox/easystab/easypack.zip).*

# Presentation of the differentiation matrix

The differentiation matrix is at the heart of all the codes that we will be using here. 
This script is a way to show how it works and give a means to test it. It is discussed from a general point of view in [pedagogy#differentiation-matrices](). 

In this code, I give an example and a test based on finite differences. For most of the codes, we embed the creation of the differentiation matrices in a special function [dif1D.m]() in 1D (and [dif2D.m]() in 2D). Please see [diffmat_dif1D.m]() for the same code as here but using [dif1D.m]() to produce the differentiation matrices.

We suppose that we have a function
`f` that is discretized on a grid `x`, let's say first a simple equidistant grid with `N` gridpoints. 
Then in the computer,
`f` is described as a column vector of the `N` values of the function `f` at each point of the grid.
To compute the derivative `f'` of this function, knowing its values at the gridpoints, the simplest way is to use 
finite differences:
$$
f'_i=\frac{f_{i+1}-f_{i-1}}{2h},
$$
where the index `i` refers to the number of the cell and h is the distance betwen grid points. This formula is an approximation of the exact derivative:
$$
f'(x)=\lim_{h \to 0} \frac{f(x+h)-f(x)}{h}
$$

The differentiation is a linear operation, so it can be described by a matrix-vector product. 
$$
f'=Df
$$
So how can we build the matrix?
We will do this simply by stacking in a matrix the operations of the finite difference.

$$
\left(
\begin{array}{l}
f'_{1}\\ f'_{2}\\ \vdots \\  \\ f'_{N}
\end{array}\right)
=
\frac{1}{2\Delta x}
\left(\begin{array}{ccccc}
-3 & 4 & -1 &\dots & 0\\
-1 & 0 & 1 &  &  \\
&\ddots&\ddots&\ddots& \\
&&-1&0&1\\
0&\dots&1&-4&3
\end{array}\right)
\left(\begin{array}{l}
f_1\\ f_2\\ \vdots \\  \\ f_N
\end{array}\right)
$$

You see that we where in the need to change the formula for the derivative at the first and last
grid points because we have used centered stencils that require the knowledge of the values of f 
before and after the grid point i. Thus at the first and last gridpoints we have used uncentered
stencils.

%}

clear all; clf

% parameters
L=2*pi; % domain length
N=15; % number of points

% the grid
x=linspace(0,L,N)';
h=x(2)-x(1); % the grid size

% first derivative
D=zeros(N,N);
D(1,1:3)=[-3/2, 2, -1/2]/h;
for ind=2:N-1
    D(ind,ind-1:ind+1)=[-1/2, 0, 1/2]/h;
end
D(end,end-2:end)=[1/2, -2, 3/2]/h;

% second derivative
DD=zeros(N,N);
DD(1,1:3)=[1, -2, 1]/h^2;
for ind=2:N-1
    DD(ind,ind-1:ind+1)=[1, -2, 1]/h^2;
end
DD(end,end-2:end)=[1, -2, 1]/h^2;


%{

Now that the differentiation matrices are built and stored in D and DD, we can test how they work
by comparing to the derivative of a function for which we know the exact derivative, 
for instance the cosinus(x). In the graph below we plot a cosinus on the grid, the its derivative: -sin in red
    and compare it to D*f in red with a dashed line. The same comparison is also done for the second derivative, in magenta.

%}

% test of the derivatives
f=cos(x);
fp=-sin(x)';

plot(x,cos(x),'b.-',x,-sin(x),'r.-',x,D*cos(x),'r+--',x,-cos(x),'m+-',x,DD*cos(x),'mx--');
legend('cos','-sin','D*cos(x)','-cos','DD*cos')
print('-dsvg','diffmat.svg'); % save the figure

%{

# test of the computation
And here is the figure that is produced by the code:

![Comparison of the numerical and exact derivative for the first and second derivative of a cosinus](diffmat/diffmat.svg)

First look at the two red curves. This is the comparison of the exact and numerical derivatives. They are close but they are not the same. The discrepancy is the numerical approximation. To improve the approximation, you can increase the number of gridpoints or use a higher order formula for the differentiation matrix. in Easystab we use often pseudospectral differentiation matrices, based for instance on Chebychev polynomials, see for instance [???]().

# Links

Here we have shown the differentiation matrix for a 1D function. We can as well build differentiation matrices for 2D and 3D functions (the function $f$ depends for instance of three spatial direction $f(x,y,z)$). See a pedagogical introduction to higher dimensionality in [pedagogy#1D,-2D-and-3D]() and the implementation in the codes [diffmat_2D.m]() and [diffmat_3D.m](). 

We use differentiation matrices to solve differential equations, for instance see [pedagogy#boundary-conditions]() for a general introduction to the way to enforce boundary conditions, and [differential_equations.m]() for the implementation. In [poisson1D.m](), [poisson2D.m]() and [poisson3D.m]() we solve Poisson equations in 1D, 2D and 3D.   

# Exercices/contributions

* Please do a convergence study in the number of gridpoints, showing how the error tends to zero when the number of points increase
* Please build the differentiation matrices using 5-point stencils for the finite differences (for this you can use the function [ufdwt.m]() which gives you the finite difference weights)
* Please test the accuracy of the differentiation on a different function than a cosinus
* Please code the differentiation matrices for a periodic domain (for which there is no need to use uncentered stencils at the boundaries, but instead use the value of the function "on the other side" of the domain -> [periodicals boudaries.m]()
* Please test the differentiation matrix based on Chebychev spectral discretisation ([chebdif.m]())like done in most of the other codes, for instance [blasius.m]().
* Please test the differentiation matrix based on Fourier spectral discretisation for periodic domain [fourdif.m]().
* Please build the differentiation matrix for higher order derivatives. [diffmat_thirdorder.m](stab2014/diffmat_thirdorder.m).
* Please do the comparison between the second derivative computed by multiplication twice by the first derivative matrix or once by the second derivative matrix. [comparison](stab2014/diffmat_comparison.m)


%}
