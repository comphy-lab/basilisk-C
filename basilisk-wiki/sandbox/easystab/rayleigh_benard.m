%{
# Rayleigh-Benard instability

This is the instability of a fluid layer heated from below. The fluid is sandwiched between two horizontal plates, and the bottom plate is hot and the top plate is cold. Since the hot fluid is lighter, it wants to raise and the cold fluid want to fall. This is hindered by two diffusive effect: the fluid viscosity slows down the motion, and the thermal diffusivity will smear out the temperature.

The important parameters are the fluid viscosity, the distance betwen the two plates, the temperature difference between the two plates, the fluid density, the gravity, the thermal diffusivity and the thermal dilatation (how much the fluid becomes lighter when it is heated).

%}

clear all; clf;

% parameters
N=50; % number of gridpoints
rho=1; % fluid density
mu=0.001; % fluid viscosity
g=1; % gravity
alpha=2.5; % wavenumber in x
d=1; % thermal dilatation
k=1; % thermal diffusivity
Ty=-1; % vertical gradient of temperature
L=1; % domain height

% differentiation matrices
scale=-2/L;
[y,DM] = chebdif(N,2); 
D=DM(:,:,1)*scale;    
DD=DM(:,:,2)*scale^2;    
y=y/scale; 
Z=zeros(N,N); I=eye(N); 

% renaming the differentiation matrices
dy=D; dyy=DD;
dx=i*alpha*I; dxx=-alpha^2*I;
Delta=dxx+dyy;

%{
# System matrices

We start with the Navier-Stokes equations for the fluid. Here we consider only 2D waves so there is no $w$ component
$$
\begin{array}{l}
\rho(u_t+uu_x+vu_y)=-p_x+\mu\Delta u\\
\rho(v_t+uv_x+vv_y)=-p_y+\mu\Delta v-\rho g\\
u_x+v_y=0\\
\end{array}
$$
We see that there is the volume force in the $v$ equation, dependent on the density, this will be the term through which the dilatation will affect the flow (the "floattability term", telling how much the fluid "floats" when it is heated). We linearize these equation about the base flow, but since the base flow is zero (initially the fluid is at rest), the equation become very simple, simply removing nonlinear terms
$$
\begin{array}{l}
\rho u_t=-p_x+\mu\Delta u\\
\rho v_t=-p_y+\mu\Delta v-\rho g\\
u_x+v_y=0\\
\end{array}
$$
now we consider the equation for the temperature, it is simply an advection-diffusion equation, where the temperature $T$ is advected by the flow and diffused by the thermal diffusivity
$$
T_t+uT_x+vT_y=k\Delta T
$$
we linearize this equation about the steady base flow and the constant vertical temperature gradiant. For the base state: $T_x=0, T_y=cste$ thus for the fluctuation $\theta$ of temperature about the base temperature we have the dynamics
$$
\theta_t+vT_y=k\Delta\theta
$$
where we see how the flow affects the evolution of the temperature: when the is an upward flow ($v$ positive), the base temperature is advected upward thus giving a source term in the equation for $\theta$. 

We need now to describe how the temperature affects the density of the fluid. For this we say
$$
\rho=\rho_r[1-d(T-T_r)]
$$
where $\rho_r$ is the reference density of the fluid at temperature $T_r$ about which we have some small variation. We could replace $\rho$ everywhere in the flow equation by this formula, but it does in fact not play a big role for inertia: the variation of density is so small that the only place where it will have a big impact is in the floatability term, saying how the gravity will affect the flow on the right hand side of the $v$ equation. When we linearize this term, it becomes $-d\rho_rg\theta$ and we drop the subscript $r$. 

We thus have the system of four equations
$$
\begin{array}{l}
\rho u_t=-p_x+ \mu\Delta u,\\
\rho v_t=-p_y+\mu\Delta v-d\rho g \theta,\\
u_x+v_y=0\\
\theta_t+vT_y=\Delta \theta
\end{array}
$$
We can rewrite this system of equations in a matrix form
$$
Eq_t=Aq
$$
with the matrices
$$
q=\left(\begin{array}{c}
u \\ v\\ p\\ \theta
\end{array}\right)
, \quad
A=\left(\begin{array}{cccc}
\mu\Delta&0&-\partial_x&0\\
0&\mu\Delta&-\partial_y&\rho g d\\
\partial_x&\partial_y&0&0\\
0&-T_y&0&k\Delta
\end{array}\right)
, \quad
E=\left(\begin{array}{cccc}\rho&0&0&0\\0&\rho&0&0\\0&0&0&0\\0&0&0&1\end{array}\right)
$$
%}

% system matrices
A=[mu*Delta, Z, -dx, Z; ...
   Z, mu*Delta, -dy, rho*g*d*I;  ...
   dx, dy, Z, Z;  ...
   Z, -Ty*I, Z, k*Delta];
E=blkdiag(rho*I,rho*I,Z,I);

%{
# Boundary conditions

The natural things would be that $u$ and $v$ are zero at the walls, and the temperature perturbation $\theta$ is also zero at the walls (the temperature is imposed at the wall, without perturbations). But in fact for here we have an analytical solution for the stability when the flow is allowed to slip along the walls, so we will instead here impose that
$$
\begin{array}{l}
u_y|_0=0\\
u_y|_L=0\\
v|_0=0\\
v|_L=0\\
\theta|_0=0\\
\theta|_L=0\\
\end{array}
$$
thus the boundary conditions are expressed $Cq=0$, with the constraint matrix
$$
C=\left(\begin{array}{cccc}
\partial_y|_0&0&0&0\\
\partial_y|_L&0&0&0\\
0&I|_L&0&0\\
0&I|_0&0&0\\
0&0&0&I|_L\\
0&0&0&I|_0\\
\end{array}\right)
$$

%}


% boundary conditions
II=eye(4*N); ddy=blkdiag(dy,dy,dy,dy);
u0=1; uL=N; v0=N+1; vL=2*N; T0=3*N+1; TL=4*N;
loc=[u0,uL,v0,vL,T0,TL]; 
C=[ddy([u0,uL],:); II([v0,vL,T0,TL],:)];
A(loc,:)=C;
E(loc,:)=0; 

%{
# Computing eigenmodes

For the linear system $Eq_t=Aq$ we have imposed the boundary conditions in a way that $E$ is not invertible, so we solve a generalized eigenvalue problem which will have infinite eigenvalues corresponding to the boundary conditions constraints. We remove them form the result and we plot only the five eigenvalues that have largest growth rate. These modes will be relevant for the stability analysis. here the wave speed of the mode is zeor, these are stationnary mode, they have no propagation in the $x$ direction, so we plot only the growth rate.
%}

% computing eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

plot(alpha,real(s(1:5)),'b.')
xlabel('wavenumber \alpha'); ylabel('exponential growth rate')
grid on; hold on

%{
# Validation

We can get an analytical solution for this problem if we have slip at the walls. The growth rate of the eigenmodes should be
$$
s=-\frac{1}{2}(1+Pr)(j^2\pi^2+\alpha^2)\pm \sqrt{\frac{1}{4}(Pr-1)^2(j^2\pi^2+\alpha^2)^2+\alpha^2\frac{RaPr}{j^2\pi^2+\alpha^2}}
$$
where $j$ is the number of the mode, and the Prandtl and Rayleigh (nondimensional) numbers are
$$
Pr=\frac{\mu}{\rho k}, \quad Ra=\frac{\rho g \Delta T L^3}{\mu k}
$$
In the code, we plot these theoretical growth rates for the first five modes ($j$ goes from 1 to 5), and for varying $\alpha$.
%}

% validation against the theory
Pr=mu/(rho*k);
Ra=rho*g*(-Ty*L)*L^3/(mu*k);

j=1:5;
alphavec=linspace(0,8,100);
Stheo=zeros(length(j),length(alphavec));
for ind=1:length(alphavec);
    alpha=alphavec(ind);
    Stheo(:,ind)=-0.5*(1+Pr)*(j.^2*pi^2+alpha^2)+sqrt(0.25*(Pr-1)^2*(j.^2*pi^2+alpha^2).^2+alpha^2*Ra*Pr./(j.^2*pi^2+alpha^2));
end
plot(alphavec,Stheo,'k-')
print('-dsvg','plot.svg');


%{

![The results](rayleigh_benard/plot.svg)


# Exercices/Contributions

* Please show that the computation fits the theory for many values of the parameters
* Please draw the neutral curve, showing the boundary between the flow situations that are stable and those who are unstable
* Please change the boundary conditions to have a wall with no-slip and find a some results in the litterature to validate your computations
* Please show that the flow is always stable when the hot plate is at the top
* Please draw the velocity field and the temperature field corresponding to the five most unstable modes.
* Please write a code that does the advection of tracer particles by the velocity field of the first eigenmode in a stable case and in an unstable case
* Extension considering surface tension and a free surface, also called [Rayleigh Benard Marangoni convection](stab2014/rayleigh_benard_marangoni.m)
%}
