%{
# Gravity free surface waves

Here is a code to get the waves at the surface of the sea. The interesting technicality here is do do a stability case with a moving free surface, so this will induce some new type of boundary conditions, and we have an additional variable, which is the position of the free surface.

The physical parameters are the gravity, the fluid density, the depth of liquid and the wavenumer of the wave we are looking at ($2\pi$ divided by the wavelength of the wave). We also add the fluid viscosity.
%}

clear all; clf;

n=100;      % number of gridpoints
alpha=1;    % wavenumber in x
L=1;        % Fluid height in y
rho=1;      % fluid density
mu=0.0001;    % fuid viscosity
g=1;        % gravity

% differentiation matrices
scale=-2/L;
[y,DM] = chebdif(n,2); 
D=DM(:,:,1)*scale;    
DD=DM(:,:,2)*scale^2;    
y=(y-1)/scale; 
I=eye(n); Z=zeros(n,n);

% renaming the matrices
dy=D; dyy=DD;
dx=i*alpha*I; dxx=-alpha^2*I;
Delta=dxx+dyy;

%{
# System matrices

The fluid motion is described by the Navier--Stokes equations and the continuity equation, linearized about a zero base flow, we suppose that appart from the waves, the fluid is at rest (contrary to the free surface when the fluid is flowing on an incline). Thus we have
$$
\begin{array}{l}
\rho u_t=-p_x+\mu\Delta u\\
\rho v_t=-p_y+\mu\Delta v-\rho g\\
u_x+v_y=0\\
\end{array}
$$ 
We describe the position of the free surface by $L+\eta(x,t)$. We need an equation for the evolution of $\eta$ when the fluid is in motion, this is simply an advection equation of $\eta$ by the velocity field $u,v$ at the position of the free surface
$$
\eta_t+u(x,y=L+\eta)\eta_x=v(x,y=L+\eta)
$$
you see that this is a nonlinear equation because $\eta$ appears through the arguments of $u$ and $v$. We will linearize this by doing an taylor expansion: express for instance $u(x,y=L+\eta)$ as a function of $u(x,y=L)$
$$
u(x,y=L+\eta)=u(x,y=L)+\eta u_y(x,y=L)
$$
thus the nonlinearity is now explicit. We linearize this about a zero base flow, which means that $u,v$ are small and also $\eta$ is small, thus if we take away the quadratic terms, we finally simply have a vertical advection of $\eta$ by the vertical velocity component.
$$
\eta_t=v(x,y=L)
$$
this is the evolution equation for $\eta$.

if we now put the system in matrix form $Eq_t=Aq$ we have
$$
q=\left(\begin{array}{c}
u \\ v\\ p\\ \eta
\end{array}\right)
, \quad
A=\left(\begin{array}{cccc}
\mu\Delta&0&-\partial_x&0\\
0&\mu\Delta&-\partial_y&0\\
\partial_x&\partial_y&0&0\\
0&I_L&0&0
\end{array}\right)
, \quad
E=\left(\begin{array}{cccc}\rho&0&0&0\\0&\rho&0&0\\0&0&0&0\\0&0&0&1\end{array}\right)
$$
but we should be careful when building the matrices, that contrary to [rayleigh_benard.m]() where the fourth additional variable was the perturbation of temperature, was defined on the grid $y$, here $\eta$ is just a scalar. The last line of the system matrices shows how the time derivative of $\eta$ is equal to the value of $v$ at the top grid point, here the gridpoint number $n$.
%}

% System matrices
A=[mu*Delta, Z, -dx, Z(:,1); ...
   Z, mu*Delta, -dy, Z(:,1); ...
   dx, dy, Z, Z(:,1); ...
   Z(1,:),I(n,:),Z(1,:),0];

E=blkdiag(rho*I,rho*I,Z,1);

%{
# Boundary conditions

At the wall we have the classical no-slip boundary conditions that $u$ and $v$ and the first gridpoint must be 0. Then we have to additional cnstraints comming from the continuity of the stress. 

The first one is the most tricky. The pressure at the free surface should be equal to the atmospheric pressure (0)
$$
p(L+\eta)=0
$$
and again here to express this as a boundary condition at $y=L$ we need to do a Taylor expansion
$$
p(L+\eta)\approx p(L)+\eta p_y(L)=0
$$
thus $p(L)=-\eta p_y(L)$, and here we will linearize this. There is no base flow, but because of the gravity there is an hydrostatic distribution of pressure, so there is a base pressure gradient to be accounted for in the linearizeation.
$$
p(L)=-\eta (-\rho g)
$$
which is now our boundary condition on the perturbation pressure.

**Note that this way to write the boundary condition misses one term coming from the normal stress. To be implemented and validated (DF, 01/2020)**

Now for the continuity of the viscous stress we have
$$
u_y+v_x=0
$$
at $y=L$. Here there is no base flow so the linearization is easy, there is no base flow term introduced by the Taylor expansino from $L+\eta$ to $L$.

To impose these boundary conditions, we replace the bottom and top gridpoint equations by these equations.

The constraint matrix $C$ is thus
$$
C=\left(\begin{array}{cccc}
I|_0&0&0&0\\
0&I|_0&0&0\\
0&0&-I|_L& \rho g\\
\partial_y|_L&\partial_x|_L&0&0\\
\end{array}\right)
$$
%}

% boundary conditions
loc=[1,n,n+1,2*n];
C=[I(1,:),Z(1,:),Z(1,:),0; ... 
   Z(1,:),I(1,:),Z(1,:),0; ...
   Z(1,:),Z(1,:),-I(n,:),rho*g; ...
   dy(n,:),dx(n,:),Z(1,:),0]; 

E(loc,:)=0;  
A(loc,:)=C;

% compute eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

%{
# Validation

We can get a theory for the inviscid equivalent of this situation, using a streamfunction to describe the velocity field
$$
u=-\phi_y, v=\phi_x
$$
so the continuity equation becomes $\Delta \phi=0$. Assuming a wave in $x$ this equation becomes
$$
-\alpha^2\hat{\phi}+\hat{\phi}_yy=0
$$
which has as solution
$$
\hat{\phi}=A \cosh(\alpha y)+B \sinh(\alpha y)
$$
The first boundary condition is no penetration at the bottom wall which gives $A=0$. At the top we have the advection of the interface by $v$
$$
s\eta=\hat{\phi}_x(L)=i\alpha B \sinh(\alpha y)
$$
Now we need a relation betwen the velocity and the pressure at the interface, which we get from the $u$ component of the Euler equation
$$
\rho u_t=-p_x
$$
thus 
$$
-s\rho \hat{\phi}_y(L+\eta)=-i\alpha \hat{p}(L+\eta)
$$
which after our Taylor expansion in $y$ gives
$$
-s\rho \hat{\phi}_y(L) \approx -i\alpha \hat{\eta} P_y(L)
$$
where the base pressure gradient $P_y=-\rho g$ comes into play.

Combining this and the advection of the interface, we get the solution for the eigenvalue $s$
$$
s^2=-\alpha g \tanh(\alpha L)
$$
which shows that we have two purely imaginary solutions for $s$, one positive and one negative. The wave velocity is equal to minus the imaginary part of $s$ divided by the wavenumber, so we see that we have two waves, one going to the right and one going to the left, and both with celerity
$$
c=\sqrt{\frac{g \tanh(\alpha L)}{\alpha}}
$$
this is the expression that we use to validate our viscous code for small values of $\mu$.
%}

% validation
alphavec=linspace(0,5,100);
ctheo=sqrt(g*tanh(alphavec*L)./alphavec);

cnum=abs(imag(s(1)))/alpha;
plot(alphavec,ctheo,'r-',alpha,cnum,'b.');
xlabel('alpha');ylabel('wave velocity'); title('validation')
grid on

set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','free_surface_gravity.png');


%{

![wave velocity as a function of the wavenumber $\alpha$](/sandbox/easystab/free_surface_gravity.png)

# Exercices/Contributions

* Please check that in the limit of short wavelength compared to the fluid depth, (deep water limit), we have a dispersive system of waves. Both theoretially and numerically ---> [free_surface_gravity-short_wavelength](http://basilisk.fr/sandbox/easystab/free_surface_gravity-short_wavelength.m)
* Please check that in the limit of long wavelength compared to the fluid depth, (shallow water limit), we have a non-dispersive system of waves. Both theoretially and numerically
* Please write a theory for a very thin and viscous layer (the lubrication limit) and show that our codes works well also in this limit.
* Please invert the direction of the gravity and show that we have a system of unstable waves without propagation (the Rayleigh-Taylor instability), show that the theory works still well.
* Please add tension surface and remove gravity and make a theory for these capillary waves just like we did, but with a pressure jump through the interface. Code this and show that the theory and the code fit well
* Please do the Rayleigh-taylor instability with gravity and surface tension and show how the surface tension is stabilizing.
* Please change the code to have two fluids on top of each other, and do as well the theory for validation ---> [free_surface_gravity_twofluids.m]()


%}