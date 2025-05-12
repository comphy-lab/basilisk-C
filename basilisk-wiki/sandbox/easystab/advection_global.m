%{
# Advection equation without time marching

Here is a nice formulation, where we remove the idea of the time marching, and we replace with the global solution of the system with treating the time derivatve just like we do with the space derivative. This becomes in fact a 2D problem where time is one of the directions. This is not very useful for the case shown here, but it will be very usefull in nonlinear to get periodic orbits, and even in linear to get steady response to periodic forcing, without computing the transient.

Conceptually this is difficult to understand, but realize that there is no formal difference between a spatial derivative and a temporal derivative in a differential equation. So a 1D system with in addition a time evolution becomes a 2D system. The same way, and 2D system with a time evolution is in fact a 3D system. 

Since I want to compute periodic orbits for nonlinear systems in 2D, they will become to searching for a steady solution of the nonlinear equations, with some special time boundary conditions: I say that the response is periodic in time, so to do this, I impose periodic boundary conditions in time. This is easily done using Fourier differentiation matrices for the time differentiation.

It is also useful for "final value problems" instead of "initial value problem". usually we tell what the initial condition is, and then we march forward in time. but maybe in some cases we want to impose the final solution and see what the initial condition will be. The global formulation that I test here will make this very easy (at the cost of an additional dimension of the system, which means just more computational time because the matrices are bigger). This is done for the advection equation in [advection_global_reversed.m]() and for the diffusion in [diffusion_global_reversed.m]().


This code is related to [burgers_global.m]() where we do the same thing but with a nonlinear equation. In that case, we can as well solve globally in time, but we need to do Newton iterations using the Jacobian.
%}

clear all; clf

% parameters
Nx=41; % number of gridpoints
Nt=40; % gridpoints in time
Lx=10; % domain length in x
Lt=4; % duration in time
U=1; % advection velocity
xpos=2; % position of the initial condition

%{
We build the differentiation matrices, just like we did in [venturi.m](). And we rename the $y$ direction by a $t$. This is just a renaming, there is no changing.
%}

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,Lx,Nx);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Lt,Nt,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% rename y to time derivative
D.t=D.y;
D.tt=D.yy;
t=y;
T=Y;

%{
here we change the definition of the location vectors for the boundary conditions, for ease of the initial conditions. 
%}
l.left=[l.cbl;l.left;l.ctl];
l.right=[l.cbr;l.right;l.ctr];
l.start=l.bot;
l.end=l.top;


%{
# System matrices

This is the advection equation
$$
f_t+Uf_x=0
$$
so the system matrix is
$$
A=\partial_t+U\partial_x
$$
and the system to be solved is
$$
Af=b
$$
where $b$ is just zeros, and will be uefull to impose the nonhomogeneous boundary conditions.
%}

% system matrices
A=D.t+U*D.x;
b=zeros(NN,1);

%{
As we often do, we build an initial guess $f_0$ which satisfies all the boundary conditions, so that these boundary conditions can just be expressed
$$
Cf=Cf_0
$$
with $C$ the constraint matrix.
%}
% initial guess
f0=exp(-((X-xpos)/0.5).^2);

%{
# Boundary conditions

They are, $f=0$ at the left boundary, and a bump at the "start boundary", the initial condition.
%}

% boundary conditions
loc=[l.start; l.left];
C=I(loc,:);
 
b(loc)=C*f0(:); 
A(loc,:)=C;

%{
# Solve system
%}

f=A\b;
f=reshape(f,Nt,Nx);

%{
# Validation

For the validation, we compare the final time (the "top boundary") withthe initial condition just translated of $UL_t$, that is the motion of moving to the right at velocity $U$ during time $L_t$.
%}

% show evolution of f
subplot(2,1,1);
surf(X,T,f); shading interp; view(2)
xlabel('x'); ylabel('t');
title('time evolution of the string');

% compare final time with theory
subplot(2,1,2);
xx=linspace(0,Lx,100);
ftheo=exp(-((xx-(xpos+U*Lt))/0.5).^2);

plot(X(l.end),f(l.end),'b.',xx,ftheo,'r-');
legend('numerical','theory'); title('solution at final time')
xlabel('x'); ylabel('f(x,t=Lt)');


%{
Here is the figure at time *tmax*:
![Validation](advection_global.png)

# Exercises/contributions

* Please do this also for other equations: diffusion and advection diffusion, and vibrating string. -->[diffusion](stab2014/diffusion_global.m), [advection diffusion](stab2014/advection_diffusion_global.m)
* Please do that for a system forced periodically, using a Fourier differentiation matrix in time. This allows to skip the computation of the initial transient and solve directly for the steady state periodic response.
* Please do this for a nonlinear system (an initial condition problem from the gobal oint of view)
* Please do this for a nonlinear system forced periodically (with Fourier in time).
%}