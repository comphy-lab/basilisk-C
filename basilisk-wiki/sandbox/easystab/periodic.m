%{

This is a code that I use to test the periodic boundaries in 1D with different kind of discretization. The Problem I solve is the same as was solved in 2D in [poisson2D.m]() and in 3D in [poisson3D.m]().

To learn more about the differenciation matrices made specially for periodical boundaries, please see [periodicals boudaries.m]().

This code is also an exemple of using the function [dif1D.m]() to build simply the differentiation matrices. This is useful to quickly go from one type of discretization to the other. The possible choices by now are

* *fou* for fourier
* *cheb* for chebychev
* *fd* for finte difference
* *fp* for periodic finite difference

%}

clear all; clf

% parameters
N=50; % number of gridpoints
L=1; % domain length

Z=zeros(N,N); I=eye(N); II=eye(2*N);

%{
# Fourier
%}
% differentiation matrices
[d.x,d.xx,d.wx,x_fou]=dif1D('fou',0,L,N);

% System matrices
A=d.xx;

% Forcing
k=2;
b=-pi^2*(k^2)*sin(pi*k*x_fou);

%{
Since periodicity is already enforced, we only need to impose the value of the function somewhere, here at the first gridpoint. 
%}
% boundary conditions
loc=[1];
A(loc,:)=I(loc,:);
b(loc)=0;

% solving the linear system
f_fou=A\b;

% compute error
fe_fou=sin(pi*k*x_fou);
err_fou=norm(f_fou-fe_fou,2)

%{
# Chebychev
%}

% differentiation matrices
[d.x,d.xx,d.wx,x_cheb]=dif1D('cheb',0,L,N);

% System matrices
A=d.xx;

% Forcing
k=2;
b=-pi^2*(k^2)*sin(pi*k*x_cheb);

%{
Here we need of course to impose two boundary conditions. The first one is the eriodicity, that the value of the function at the first and last gridpoints should be the same. This is written

    f(1)=f(N)
    
but since we do not manipulate $f$ but instead write in the form of a matrix/vector product the constraints that $f$ should satisfy, this is written
 
    I(1,:)*f-I(N,:)*f=0

or again

    (I(1,:)-I(N,:))*f=0

The second condition, just like for above, is the fact that the value of the solutino at the first gridpoint should be zero.
%}
% boundary conditions
loc=[1,N];
A(loc,:)=[I(1,:)-I(N,:); I(1,:) ];
b(loc)=0;

% solving the linear system
f_cheb=A\b;

% compute error
fe_cheb=sin(pi*k*x_cheb);
err_cheb=norm(f_cheb-fe_cheb,2)


%{
# Periodic finite difference

For finite difference, there is one additional input argument for the [dif1D.m]() function, the number of points in the finite-difference stencil. Here we chose 5. Choose more if you want a higher order scheme (7, 9 ...)

%}

% differentiation matrices
[d.x,d.xx,d.wx,x_fp]=dif1D('fp',0,L,N,5);

% System matrices
A=d.xx;

% Forcing
k=2;
b=-pi^2*(k^2)*sin(pi*k*x_fp);

%{
Just for Fourier, here the periodicity is enforced by default in the differentiation matrix, so there is only the need for one boundary condition. (and you can check in fact that the last point in the grid is not at 1, just like for Fourier).
%}
% boundary conditions
loc=[1];
A(loc,:)=[I(1,:)];
b(loc)=0;

% solving the linear system
f_fp=A\b;

% compute error
fe_fp=sin(pi*k*x_fp);
err_fp=norm(f_fp-fe_fp,2)


%{
# Finite difference
%}

% differentiation matrices
[d.x,d.xx,d.wx,x_fd]=dif1D('fd',0,L,N,5);

% System matrices
A=d.xx;

% Forcing
k=2;
b=-pi^2*(k^2)*sin(pi*k*x_fd);

%{
The boundary conditions are just like for Chebychev above.
%}
% boundary conditions
loc=[1,N];
A(loc,:)=[I(1,:); I(1,:)-I(N,:)];
b(loc)=0;

% solving the linear system
f_fd=A\b;

% compute error
fe_fd=sin(pi*k*x_fd);
err_fd=norm(f_fd-fe_fd,2)


% plotting the results
plot(x_fou,f_fou,'b',x_cheb,f_cheb,'m',x_fp,f_fp,'k');
xlabel('x'); ylabel('f'); legend('Fourier','Chebychev','finite difference periodic')
title('Poisson problem');


%{

The code gives the screen output:

    err_fou =
       4.0309e-13
    err_cheb =
       3.5058e-13
    err_fp =
       1.3834e-05
    err_fd =
       8.8751e-06

![And the figure](periodic.png)

# Exercices/contributions

* Please do a convergence study with the grid to compare the different discretizations


%}