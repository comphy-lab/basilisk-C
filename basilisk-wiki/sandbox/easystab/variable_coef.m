%{
In [differential_equation.m]() we use the differentiation matrices to solve a 1D differential equation with constant coefficients. Many physical systems have variable coefficients, this is for instance the case for the stability of fluid flows when there is a base flow as for instance in the classical stability of the Poiseuille profile in the channel in [poiseuille_uvp.m]().

Here we show an example of solving such equation and the way we use a diagonal matrix to account for this variation.
%}

clear all; clf

% parameters
L=1; % domain length
N=25; % number of points

% 1D differentiation matrices
[d.x,d.xx,d.wx,x]=dif1D('cheb',1,L,N,3);
I=eye(N); Z=zeros(N,N);
%{
# The differential equation

We seek to solve 
$$
x^2 f_{xx}-f_x=x^2+x
$$
on $x\in[1,2]$ with boundary conditions
$$
f(1)=1, f_x(1)=5
$$
this is exercice B -1.1 7 of [Zachary S. Tseng, math course on differential equations](http://www.math.psu.edu/tseng/class/Math251/Notes-2nd%20order%20ODE%20pt1.pdf). The solution exist on $(0,\infty)$ and is
$$
f(x)=\frac{x^3}{3}+\frac{7x^2}{4}+\frac{x^2}{2}log(x)-\frac{13}{12}.
$$
As compared to [differential_equation.m]() we have two new things, the fact that $f_{xx}$ is multiplied by a function of $x$, this is the main difficulty, and the fact that the nonhomogeneous term is also a function of $x$.

Since the coefficient is different for each grid cell, we may have to use a vector to input this informatino in the equation. In fact since this coefficient multiplication is not a variable but an operator, it has to be a matrix. Thus this will be a diagonal matrix whose diagonal elements are the values of the coefficient for each grid cell as illustrated on the figure. Here we have called $c(x)$ this variable coefficient.  

![Matrix representation of a variable coefficient](variable_coeffig.png)

%}

% I build the linear differential equation
A=diag(x)*d.xx-d.x;

% boundary conditions
I=eye(N);
loc=[1,N];
A(loc,:)=[I(1,:); d.x(1,:)];
b=x.^2+x; b(loc)=[1,5];

% solve the system
f=A\b;
solexact=x.^3/3+7*x.^2/4+x.^2/2.*log(x)-13/12;

%{
# Validation
%}

% plotting
plot(x,f,'b.-',x,solexact,'r-');
xlabel('x');ylabel('f');
legend('numerical','theory')
xlim([1,1+L]); grid on

set(gcf,'paperpositionmode','auto')
print('-dpng','-r75','variable_coef.png')

% the approximation error
err=norm(f-solexact)

%{
The screen output gives:

    err =
       1.8393e-12

And here is the figure:

![Theory and computation](variable_coef.png)

# Links

In [poisson2D.m]() and [poisson3D.m]() we solve differential equations with constant coefficients. 

In [poiseuille_uvp.m](), [orr_sommerfeld_squire.m]() for instance, we compute eigenmodes of linear systems with variable coefficients. The variable corefficients are the velocity profile of the base flow and its derivatives. 

In [blasius.m]() the nonlinear equation that we want to solve leads to variable coefficients once linearized about the actual guess for the solution in the steps of the Newton iterations. Another nonlinear equation which gives variable coefficients once linearized is in [meniscus.m]().

In [thermal_entrance.m]() the variable coefficients in the Jacobian (the matrix representation of the linearization of the system) come from the fact that we are in a cylindrical coordinate system.

In [venturi.m]() this is the Navier-Stokes equations in 2D, with as well variable coefficients comming from the base flow arround which the equations are linearized.

# Exercices/Contributions
    
* Please write the codes that solve exercices 1 to 9 on page 14 of [Zachary S. Tseng, math course on differential equations](http://www.math.psu.edu/tseng/class/Math251/Notes-2nd%20order%20ODE%20pt1.pdf). Doing that is good practice for you because you don't need to change much in the present code (be aware of the domain of validity of the solution!).
* Please find in the litterature a differential equation in 2D and 3D with variable coefficients which has an analytical solution, and solve it with easystab, basing your code on [poisson2D.m]() and [poisson3D.m]().

%}
