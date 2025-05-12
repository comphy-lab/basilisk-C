%{
# March in time of the vibrating string

In this code, we simulate the evolution in time of an initial condition of a model of a guitar string, which has one point fixed and the other attached to a spring. The spring will have an effect of dissipation for the system. What is interesting with this system from the technical point of view is that the dynamic equation has a second derivative in time and we show here how to transform this into a larger system with just a single time derivative.

We use the code of [the vibrating string](../vibrating_string.m) and by changing the boundary conditions to add one spring at one end of the string. To run this code you need the code Chebychev differentiation, [here](../chebdif.m)

%}

clear all; clf

% parameters
N=50; % number of gridpoints
L=2*pi; % domain length
c=1; % wave velocity
dt=0.05; % time step
tmax=50; % final time
mu=1; % unit length density 
T=c^2*mu; % tension 
k=0.01; % elasticity coefficient
a=0.1; % damping coefficient

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
Here we build the matrices that give the discretization of the equations. And the last point is attached to a spring. So the differential equations for this condition is
$$
dx*ftt=-c^2*fx-k/mu*f-a/mu*ft
$$
k/mu and c/mu are two constants, we denote them c2 and c3.

we thus have the array representation for the system
$$
\left(
\begin{array}{cc}
I & 0 & . & .\\
0 & I & . & .\\
. & . & dx & c3\\
. & . & 0 & 1\\
\end{array}\right)
\left(\begin{array}{c}
v_t \\
f_t \\
vs_t \\
fs_t \\
\end{array}\right)=
\left(
\begin{array}{cc}
0 & c^2*dxx & . & . \\
I & 0 & . & .\\
. & -c^2*dx & 0 & -c2\\
. & . & 1 & 0\\
\end{array}\right)
\left(\begin{array}{c}
v \\
f \\
vs \\
fs \\
\end{array}\right)
$$
Here we have put explicitely a large identity matrix on the left because this will be useful to impose the boundary conditions. We thus have a linear system of the form
$$
Eq_t=Aq
$$
%}

% locations in the state
v=1:N;
f=N+1:2*N;
vs=2*N+1;
fs=2*N+2;

for j=1:2 % A loop is used to operate the validation   
    if j==2
       c=1;k=1000;a=0;%specific parameters for infinite elasticity coefficient and null damping coefficient
    end
    
c2=k/mu; % k/mu
c3=a/mu; % a/mu

% system matrices
E=[I,Z,Z(:,[1,2]); Z,I,Z(:,[1,2]);Z([1,2],:),Z([1,2],:),[L/N,c3;0,1]];
A=[Z,c^2*dxx,Z(:,[1,2]); I, Z,Z(:,[1,2]);Z([1,2],:),[-c^2*dx(N,:);Z(1,:)],[0,-c2;1,0]];

%{
# Boundary conditions
The original system has one variable and two time derivatives, so we should apply two boundary conditions. We will simply say that the first point of the string is fixed. This is a homogeneous Dirichlet condition applied on the first point of the state vector for the state position.
%}
% boundary conditions
loc=[f(1),f(N)];
E(loc,:)=0; 
A(loc,1:2*N)=II(loc,1:2*N);
A(f(N),fs)=-1;
%{
# March in time

To write a relation between the state at a given time and the state a little later, we will do a numerical approximation of the time derivative. We use the Crank-Nicolson scheme, which is quite simple and robust. We just integrate in time the evolution equation from time $t$ to time $t+dt$
$$
\int_t^{t+dt} E q_t dt=\int_t^{t+dt} Aq dt
$$ 
the integral on the left hand side is just $E(q(t+dt)-q(t))$ since $E$ does not change in time in this example, and the right hand side we approximate with the trapezoidal rule for integration
$$
E(q(t+dt)-q(t))\approx A dt (q(t+dt)+q(t))/2
$$ 
since now we want to express the new state $q(t+dt)$ as a function of the old state $q(t)$, we rearrange
$$
(E-Adt/2)q(t+dt)\approx (E+A dt/2) q(t)
$$ 
and we build the march-in-time matrix by inverting the left hand side matrix to transform this implicit system into an explicit system
$$
q(t+dt)\approx (E-Adt/2)^{-1}(E+A dt/2) q(t)
$$ 
This matrix inverse is well-defined since we have already imposed the propoer boundary conditions. 
%}

% march in time matrix 
M=(E-A*dt/2)\(E+A*dt/2);

%{
# Initial condition
We need to describe the state of the system at initial time. We say here that the string is initially deformed as a sinus. The velocity is zero. From that state we then perform a loop to repeatedly advance the state of a short forward step in time. We use drawnow to show the evolution of the simulation as a movie when running the code.
%}

% initial condition
q=[zeros(N,1); sin(pi*x/L); 0 ; 0]; 

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

for ind=1:Nt    
    q=M*q; % one step forward
    
    if j==1
    plot(x,q(f),'b',x,q(v),'r--',L,q(fs),'blacko');    
    axis([0,L,-2,2])    
    title('Vibrating string')
    xlabel('x'); ylabel('f');
    drawnow
    elseif j==2
    figure(2);
    plot(x,q(f),'.b',x,cos(c*pi/L*ind/Nt*tmax)*sin(x*pi/L),'r');
    title('infinite elasticity coefficient and null damping coefficient');
    axis([0,L,-2,2]);
    drawnow
    end
 end

% Legends are displayed after calculation as they strongly slow down dynamical graphics
    if j==1
        figure(1);legend('position','velocity');
    elseif j==2
        figure(2);legend('numerical solution','theorical solution');
    end
end

%{
# Validation
We test it with infinite elasticity coefficient and null damping coefficient.
It is supposed to behave as a fixed point, thus we should see a strings vibrating fixed at both ends.
In mode 1 with adapted boundary and initial conditions, each one should follow this equation
$$
f(x,t)=cos(c\pi t/L)sin(\pi x/L)
$$

%}

%{
Here is the figure at time *tmax*:
<center>
![](/sandbox/easystab/stab2014/vibrating_string_spring_dissipation.png)
</center>
And the figure we test with infinite elasticity coefficient
<center>
![](/sandbox/easystab/stab2014/vibrating_string_spring_dissipation_valid.png)
</center>
%}