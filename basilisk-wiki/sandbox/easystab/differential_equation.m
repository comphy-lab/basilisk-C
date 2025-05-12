%{

# Solving a linear differential equation

This code is built upon [diffmat_dif1D.m]() where the differentiation matrices
are built and validated. Here we solve a linear non-homogeneous differential 
equation in 1D with non-homogeneous boundary conditions.
    
We show here how to use the differentiation matrix to code the differential equation as a linear system and impose Dirichlet and Neuman boundary conditions.

To learn more about differentiation matrices, see [pedagogy#differentiation-matrices]() and for boundary conditions, please see [pedagogy#boundary-conditions]().
%}

    clear all; clf

    % parameters
    L=2*pi; % domain length
    N=25; % number of points

    % 1D differentiation matrices
    [D,DD,wx,x]=dif1D('fd',0,L,N,3);
    
%{
# The differential equation

We solve this
$$
f_{xx}=1
$$
with the boundary conditions 
$$
f(1)=0, f_x(L)=1
$$

So you see that this is a non-homogeneous equation, since there is a forcing term. And we have both Dirichlet condition (at $x=0$) and Neumann boundary condition (at $x=L$). You see that we have used the first order differentiation matrix also when imposing the boundary condition. Off course this is very natural since these matrices are the way we have transformed the operation of differentiation for our discrete systems.

What we have done is the following. We have replaced the first equation (the first line of the system matrix *A* by a new equation that tells that the value of *f* at the first gridpoint should be equal to 0. We did this using the first line of the identity matrix *I*. Indeed, this first line multiplied by the vector of *f* is the first value in *f*. We then did the same thing for the Neumann boundary condition at the last gridpoint. We have replaced the line in *A* by the last line of the differentiation matrix *D*. Indeed, this last line multiplied by the vector of *f* will give you the derivative at $L$.
    %}

    % I build the linear differential equation
    A=DD;
    
    % boundary conditions
    I=eye(N);
    A([1,N],:)=[I(1,:); D(N,:)];
    b=1+zeros(N,1); b([1,N])=[0,1];

    % solve the system
    f=A\b;

    % plotting
    plot(x,f,'b.-',x,x.^2/2+(1-L)*x,'r-');
    xlabel('x');ylabel('f');
    legend('numerical','theory')
    xlim([0,L]); grid on

    set(gcf,'paperpositionmode','auto')
    print('-dsvg','-r75','plot.svg')

    %{
    
![The results](differential_equation/plot.svg)

# Links

We do more or less the same things but in 2D in [poisson2D.m]() and in 3D in [poisson3D.m](). 

In [variable_coef.m]() we solve in 1D a differential equation with variable coefficients.

    # Exercices/Contributions
    
    * Please put the two boundary conditions on the same side of the domain
    (for instance, impose the value of $f$ and its derivative at $x=0$) --> [differential_equation_sameside.m]()
    * Please solve $f_{xx}=\sin(x)$ --> [differential_equation_sin.m]()
    * Please check the solution and convergence of the system a*Df+b*DDf=c, where a,b stands for scalars and c a vector --> [differential_equation_secondorder.m]()
    * Please check the convergence with the number of gridpoints --> [convergence_of_the_solution.m]()
    * Please use Chebychev differentiation matrix [chebdif.m]() and compare with the finite
    differences -> [differential_equation_chebychev.m](http://basilisk.fr/sandbox/easystab/stab2014/differential_equation_chebychev.m)
    * Please use a non-constant coefficient i your differential equation
    $a(x)f_{xx}=0$ (first find the way to code that using differentiation
    matrices)
    *Please check convergence of the system aDDf+bDf+cf+d=0 [classical_differential_equation.m]()
    %}
