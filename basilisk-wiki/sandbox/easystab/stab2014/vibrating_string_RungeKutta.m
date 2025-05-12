Comment by Jerome: You put a figure that is not produced by the code, I don't know how you made this figure, so I cannot help you. This is not reproductible.



%{
# March in time using the second-order Runge Kutta scheme 
In this code, I will use the original code [/sandbox/easystab/vibrating_string.m]() but using the second-order Runge Kutta scheme instead of the Crank Nicolson scheme for the time marching.
For the second-order Runge Kutta scheme, I use the approximation:
$$
y_{n+1}-y_n=dt.f(t_n+dt/2,y_{n+1/2})
$$
  with
$$
y_{n+1/2}=y_n+dt/2.f(t_n,y_n)
$$
Applying this formula to the model of vibrating string (with the same notations), we obtaint the march in time matrix:
$$
q(t+dt/2)=q(t)+dt/2.A.q(t)=(II+A.dt/2)q(t)
$$

$$
E(q(t+dt)-q(t))=dt.A.q(t+dt/2)=dt.A(II+A.dt/2)q(t)
$$
we rearrange this equation to build the march-in-time matrix
$$
q(t+dt)=E^{-1}(E+dt.A(II+A.dt/2))q(t)
$$
So we have the march-in-time matrix by using the second order Runge Kutta scheme
$$
M=E^{-1}(E+dt.A(II+A.dt/2))
$$
%}

clear all; clf

% parameters
N=50; % number of gridpoints
L=2*pi; % domain length
c=1; % wave velocity
dt=0.001; % time step
tmax=10; % final time

% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 
II=eye(2*N);

%{
# System matrices

Here we build the matrices that give the discretization of the equations. The equations is
$$
f_{tt}+c^2f_{xx}=0
$$ 
This is the classical vibrating string equation, also known as D'Alembert's equation. This is also known as the *wave equation*, here just in 1D. The vibration of the string will depend upon how heavy it is (the mass per unit length), and how tensed it is (the forced applied onto the string by the two attachment points). The heavier it is, the slower will th wave be, and the more tensed, the faster, so these two parameters compete for the wave speed *c*. Here we have written the  equations already with the parameterization that let *c* appear.

This is an equation with two time derivatives and to apply always the same methods, we will transform it into an equation with a single derivative, simply by augmenting the state vector: instead of describing the state of the system by only the position of the string, we will store the position of the string and also the velocity of the string (for every gridpoint). The state vector is thus
$$
q=\left(
\begin{array}{c}
v \\
f \\
\end{array}\right)
$$
once this done, we need to tell the system that $v$ is indeed the velocity, that is
$$
f_t=v
$$ 
and then the systems dynamics in terms of both $v$ and $f$
$$
v_t+c^2f_{xx}=0
$$
we thus have the array representation for the system
$$
\left(
\begin{array}{cc}
I & 0 \\
0 & I \\
\end{array}\right)
\left(\begin{array}{c}
v_t \\
f_t \\
\end{array}\right)=
\left(
\begin{array}{cc}
0 & -c^2D_{xx} \\
I & 0 \\
\end{array}\right)
\left(\begin{array}{c}
v \\
f \\
\end{array}\right)
$$
Here we have put explicitely a large identity matrix on the left because this will be useful to impose the boundary conditions. We thus have a linear system of the form
$$
Eq_t=Aq
$$
%}


% system matrices
E=[I,Z; Z,I];
A=[Z,I; c^2*dxx,Z]; 

% locations in the state
v=1:N;
f=N+1:2*N;

%{
# Boundary conditions
The original system has one variable and two time derivatives, so we should apply two boundary conditions. We will simply say that the position of the string at its both ends should be zero: two attachment points. This is a homogeneous Dirichlet condition applied on the first and last point of the state vector for the state position.

But there is one problem with the matrix E. As you can see in the next part, we construct the march-in-time matrix by using the inverse of matrix E. But if we impose the boundary condition as we do in the original code, matrix E will be non invertible (cause it has one row of 0).

To solve this problem, we will use the zero derivative boundary condition instead of the homogenous Dirichlet condition. This will give us the same result thanks to the characteristic of the Runge Kutta scheme. This scheme is a forward scheme, it means that it will calculate the value of the function at time $t$ by the previous value at time $t-dt$. So if we impose the derivative in time of the function is zero and the initial condition is zero, we will always have the value of the function at that points is equal to zero independent of time $t$ which means we have the homogenous Dirichlet condition.

As a result, we can impose the boundary condition like below:

%}

% boundary conditions
loc=[f(1),f(N)];
A(loc,:)=0;
E(loc,:)=II(loc,:);

%{
# March in time
Like showing in the first part of this page, for the march-in-time matrix here, I use the formular that I built before for the march in time
%}
% march in time matrix by using second order Runge-Kutta scheme
M=E\(E+A*(II+dt/2*A)*dt);
         

%{
# Initial condition
We need to describe the state of the system at initial time. We say here that the string is initially deformed as a sinus (this satisfies the boundary conditions), and that the velocity is zero. From that state we then perform a loop to repeatedly advance the state of a short forward step in time. We use *drawnow* to show the evolution of the simulation as a movie when running the code. We store for validation the string position at the midle of the domain, to do this without worrying about the the grid points are, we interpolate $f$ with the function *interp1*. 
%}

% initial condition
q=[zeros(N,1); sin(pi*x/L)]; 

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

for ind=1:Nt    
    q=M*q; % one step forward
    e(ind)=interp1(x,q(f),L/2); % store center point
    
    % plotting
    subplot(1,2,1);
    plot(x,q(f),'b',x,q(v),'r--');    
    axis([0,L,-2,2])    
%     drawnow
end
legend('position','velocity'); title('Vibrating string')
xlabel('x'); ylabel('f');

%{
# Validation
Now we show that the code does what we think it does, and also give ourselves the means to tell how precise and robust it is. The period of oscillation of this wave should be the time that it takes for a wave to travel at speed $c$ along its wavelength. Here we have chosen an initial condition whose wave length is twice the domain length $2L$ (this is "mode $1/2$"), thus the period of oscillations is 
$$
T=2L/c
$$
thus the central point of the string should evolve like
$$
cos(2\pi t/T)
$$
Here we build a second time vector, to plot the theoretical solution with large resolution in time to have a smooth plot even when the time step of the computation is large.
%}

% time evolution of central point
subplot(1,2,2);
Ttheo=2*L/c;
tt=linspace(0,tmax,500);
etheo=cos(2*pi*tt/Ttheo);
plot(tvec,e,'b.-',tt,etheo,'r-');
legend('numerical','theory');
xlabel('time'); ylabel('f(L/2)');


%{
Here is the figure at time *tmax*:

![Validation](vibrating_string_RungeKutta.jpg)

From the figure above, we can confirm that our numeric results approach the theorical solution.
%}
