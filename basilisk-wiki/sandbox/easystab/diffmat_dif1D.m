%{

*This document belongs to the [easystab](http://basilisk.fr/sandbox/easystab/README) project, check the main page of the project to understand the general philosophy of the project.*

*The present program was specifically designed as a support for the course ["Introductions to hydrodynamical instabilities"](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) for the "DET" Master's cursus in Toulouse. For the underlying theory please see [Lecture notes for chapter 3](http://basilisk.fr/sandbox/easystab/NumericalMethodsForEigenvalueProblems.md)*

*To use this program: click "raw page source" in the left column, copy-paste in the Matlab/Octave editor, and run. You should have previously downloaded the package [easypack.zip](http://basilisk.fr/sandbox/easystab/easypack.zip).*

# Differentiation matrices using dif1D.m 

This code is just like [diffmat.m]() but instead of building explicitely the differentiation matrices, we use the function [dif1D.m](). This is the way we do in most of the codes, since producing the differentiation matrices is one of the fondamental elements of Easystab.

%}

clear all; clf

% parameters
L=2*pi; % domain length
N=15; % number of points

%{

# dif1D.m

The output arguments are:

* D: the first derivative differentiation matrix
* DD: the second derivative differentiation matrix
* wx: the integration weights (they are used to integrate a function on the grid, please see [integration_2D.m]() to learn about how to build and understand these weights).
* x: the location of grid cells (This location is increasing with cell number, that is, $x(i+1)>x(i)$).

The intput arguments are:

* 'fd': the method of interpolant. Here 'fd' means finite differences. The other possible choices are: 'fp' for periodic finite differences (see [periodicals%20boudaries.m]() for explanations), 'cheb' for Chebychev polynomial, 'fou' for Fourier periodic. A number of other choices are also implemented, see direcly in the program [dif1D.m]().
* 0: the start of the grid.
* L: the end of the grid.
* N: the number of grid cells.
* 3: the number of element in the finite difference stencil. Here, 3 means centered stencils that use the value on the right and the value on the left of the point at which we want to compute the derivative. If you chose 5, you will use two points at the left and two points at the right (this means an approximation of higher order). For other choices than finite difference, this number is not used, since spectral differentiation matrices are dense (there are no zeros...)

To get a better idea of how this is done, please see [dif1D.m]().

%}

% building the differentiation matrices
[D,DD,wx,x]=dif1D('fd',0,L,N,3);

%{

# Validation of the computation

%}

% test of the derivatives
f=cos(x);
fp=-sin(x)';

plot(x,cos(x),'b.-',x,-sin(x),'r.-',x,D*cos(x),'r.--',x,-cos(x),'m.-',x,DD*cos(x),'m.--');
legend('cos','-sin','D*cos(x)','-cos','DD*cos')
%print('-dsvg','diffmat.svg'); % save the figure
saveas(gcf,'diffmat','svg');
%{

And here is the figure that is produced by the code:

![Comparison of the numerical and exact derivative for the first and second derivative of a cosinus](diffmat_dif1D/diffmat.svg)

# Links

You can look at [diffmat.m]() where we explicitely build the differentiation matrice for finite differences. To see how we use these 1D differentiation matrices to compute the derivatives on a 2D grid please see [diffmat_2D.m](). That code uses [dif2D.m]() which takes as input the output of [dif1D.m]() to combine them into 2D differentiation matrices. You can as well differentiate in 3D, see [diffmat_3D.m]().

If you want to learn about Chebychev and Fourier differentiation matrices, please see the paper by the guys who wrote the functions [matlabdifmatsuite.pdf]().

For the finite difference, I wrote the function [fddif.m]() which uses the stencils from [ufdwt.m](), writen by Greg von Winckel.

# Exercices/contributions

* Please
* Please
* ...

%}
