%{
# Free surface oscillations in 2D (with surface tension) 

We have a perfect fluid on a rectangular box with the free surface at the top. Here we chose to not have the Navier-Stokes equations for the fluid but just describe its flow with a velocity potential (and the Bernoulli equation for the relation between the flow and the pressure). Thus, all the dynamics is just the dynamics of the free surface, and the flow is just a constraint for the evolution of this free surface. On the contrary to [free_surface_gravity.m](), here the oscillations of the free surface are not gravity waves, they are just capillary waves, so the theory is a little different from (and complementary to) [free_surface_gravity.m#validation]().

See [free_surface_2D_functions.m]() for the same thing but writen with the functions [dif1D.m]() and [dif2D.m]().

Dependency:

* [chebdif.m]() for the Chebychev differentiation matrices

%}


clear all; clf

% parameters
Nx=20; % gridpoints in x
Ny=20; % gridpoints in y 
Lx=1; % domain size in x
Ly=1; % domain size in y
pts=5; % number of points in finite difference stencils
sigma=1; % surface tension
rho=1; % fluid density

% 1D and 2D differentiation matrices
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,Lx,Nx,3);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,Ly,Ny,3);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

%{
# Location vectors

As it is for [venturi.m](), and maybe here a little bit more, one of the central aspects of building the system matrices and imposing the boundary conditions has a lot to knowing where the different variables are in the state vector and in the matrices. This is done using the "location vectors", which are nothing else than a memory of the indices that correspond to the different parts of the state variable $q$. Here we have the velocity potential $\phi(x,y,t)$ defined on a 2D grid, and we have $\eta$ the position of the free surface, defined just along the top of the rectagular box (it is a 1D variable). So one nice aspect of the present example is that we have to combine variables of different dimensionalities.

Here the *top* of the box are the grid points at the top and excluding the corners. The corners are instead included to the *right* and *left* location vectors. This is good because we impose boundary conditions for the free surface at the corners, so we can use *top* just to access the inside points of the top of the mesh. 

Implicitely, the *top*, *bottom*, *left* and *right* location vectors pertain to the velocity potential $\phi$ which thus occupies the first Nx*Ny elements of the state vector. We then have the values of $\eta$ following.

%}

% vectors for coordinate selections 
NN=Nx*Ny;
l.left=[l.cbl; l.left; l.ctl];
l.right=[l.cbr; l.right; l.ctr];
l.phi=(1:NN)';
l.eta=((NN+1):(NN+Nx))';

% useful matrices
Z=zeros(NN,NN); ZZ=zeros(NN+Nx,NN+Nx);    zx=zeros(Nx,Nx);
I=eye(NN);      II=eye(NN+Nx);            ix=eye(Nx);

%{
# System matrices

The continuity equation for $\phi$ is nothing else than its Laplacian equal zero
$$
\Delta \phi=0
$$
the velocity field is recovered from $\phi$ by derivations
$$
u=\phi_x, v=\phi_y
$$

This equation for $\phi$ is not an evolution equation, there is no time derivative. The ingredient of evolution comes for the dynamics of the free surface. Just as for [free_surface_gravity.m#boundary-conditions](), the interface is locally advected vertically by $v$ thus
$$
\eta_t=\phi_y(y=L)
$$

To build the dynamic matrix *A* we stack the Laplacian of the 2D grid and a matrix of zeros of the size of the $\eta$ grid-that is-with *Nx* points. And then we add the block that connects the vertical velocity at the top to $\eta$. For this we use the location vectors. The state vector is thus
$$
q=\left(\begin{array}{c}
\phi\\ \eta 
\end{array}\right)
$$

The left hand-side matrix $E$ is full of zeros for the part of $\phi$ (its equation has no time-derivative: it is just that $\Delta \phi=0$). The $\eta$ part of $E$ is the identity matrix.
%}

% system matrices
A=blkdiag(D.xx+D.yy,zx);
A(l.eta(2:end-1),l.phi)=D.y(l.top,:);
E=blkdiag(Z,ix);

%{
# Boundary conditions

We use no-penetration (slip) of the fluid (since it is not viscous, we cannot impose no-slip). This means that on the lateral boundaries we have $u=\phi_x=0$, this is the first equation in *c*

> Dx([left;right],:)*II(phi,:)

At the right, there is the "big" identity matrix *II* for which we extract only the rows corresponding to $\phi$ such that 

> II(phi,:)*q=q(phi)

and then we multiply this with the $x$ differentiation matrix *Dx* but with just the rows that correspond to the *left* and *right* gridpoints. Thus this line of the constraint matrix tells that all this must be zero whatever happens.

At the bottom boundary we have $v=\phi_y=0$, this is the second equation in *c*

> Dy(bot,:)*II(phi,:)

We are now done with the boundary conditions for $\phi$ and we can continue with the boundary conditions for $\eta$. We tell that it should have a zeros derivative at $x=0$ and $x=Lx$. This is good here because we can write by hand in [#validation]() a theory for the validation of this case. This is writen in the third row of *c* where we use the 1D differentiation matrix in $x$ *dx*, and selecting only the first and last rows of it

> dx([1,Nx],:)*II(eta,:)

We then need to impose that the interface induces no variation of volume compared to the square domain-that is-its integral should be zero. This is just a scalar constraint. For this we use the integration vector *ix* built at the time of building the 1D matrices. This is the mast row of *c*

> d.wx*II(eta,:)

Now we need yet one last constraint. We have writen in the dynamic matrix how the velocity field affects the interface (vertical advection), but we did not write how the interface affects the fluid. Across a free surface tensed by surface tension, there is the "Laplace pressure jump" proportional to the interface curvature times the value $\sigma$ of the surface tension
$$
p_{in}-p_{out}=\sigma/R
$$
where $R$ is the radius of curvature of the interface. We assume that the outside pressure is zero, and if the interface is close to flat, the expression of the curvature is simple (linear)
$$
p_{in}=-\sigma\eta_{xx}.
$$
If the second derivative of $\eta$ is positive, the concave side of the interface is below, so the pressure should be less in the fluid than outside, this is why there is the minus sign. 

Now we need to relate $\phi$ to the pressure inside the domain, for this we use the instationnary Bernoulli equation linearized about a zero base flow
$$
\rho \phi_t+p=0 
$$
thus we have
$$
\rho \phi_t=\sigma \eta_{xx}
$$
which is our last boundary condition. Note that this is not just a simple boundary conditions as we had until now, it is a little more subtle because there is a time derivative of $\phi$ in this boundary condition. This means that in the numerical expression of this constraint, we will need to have a term in the left-hand-side of 
$$
Eq_t=Aq
$$
so there will not only be zeros in $E$. This is coded in *Ca* and *Ce*. *Ca* is the part of the constraint that goes into *A* and *Ce* is the part of it that goes into *E*. 

%}

% boundary conditions
loc=[l.top;l.left;l.bot;l.right;l.eta([1,2,Nx])];

c=[D.x([l.left;l.right],:)*II(l.phi,:); ...
   D.y(l.bot,:)*II(l.phi,:); ...
   d.x([1,Nx],:)*II(l.eta,:); ...
   d.wx*II(l.eta,:)];

Ca=[c; ...
    sigma*d.xx(2:end-1,:)*II(l.eta,:)];

Ce=[0*c; ...
    rho*II(l.top,:)]; 

A(loc,:)=Ca;
E(loc,:)=Ce;

%{
The structure of these matrices is quite sophisticated so it is interesting to see what they look like
%}

% show the matrices
figure(1); 
subplot(1,2,1); spy(E); title('E');
subplot(1,2,2); spy(A); title('A');
set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','free_surface_2D_spy.png');

%{
![spy of the dynamic matrices](free_surface_2D_spy.png)
%}

%{
# Eigenmodes

There is no dissipation and no instability in this flow, so all the modes should be neutral (the real part of the eigenvalues $s$ should be 0), but there are oscillations so the imaginary parts should not be zero). The spectrum should be right-left symetric because the system itself is, so we remove the eigenmodes with negative imaginary parts.
%}

% compute eigenmodes
disp('computing eigenmodes');
[U,S]=eig(full(A),full(E));
s=diag(S);  [t,o]=sort(abs(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
rem=imag(s)<0; s(rem)=[]; U(:,rem)=[];

%{
# Validation

Just as for the example in [wave_like.m#validation](), we have chosen boundary conditions on $\eta$ so that by chance our numerical solution should be the same as a wavelike solution periodic in $x$. We do here a theory similar to that of [free_surface_gravity.m#validation]() but with surface tension (pressure jump through the interface) instead of gravity (pressure gradient in the fluid). And also as a difference we use velocity potential instead of stream-function (for a change...).

We have the wavelike assumption
$$
\phi(x,y,t)=\hat{\phi}(y)\exp(i\alpha x+st)+C.C.
$$
so the Laplacian of $\phi$ becomes
$$
-\alpha^2\hat{\phi}+\hat{\phi}_{yy}=0
$$
and the solution should be like
$$
\hat{\phi}=A\cosh(\alpha y)+B\sinh(\alpha y)
$$
The bottom boundary condition is $\hat{\phi}_y(y=0)=0$ thus $A=0$.

The advection of the interface gives
$$
s\hat{\eta}=\hat{\phi}(y=L)=\alpha B \sinh(\alpha L)
$$
and the presure jump through the interface gives
$$
s\hat{\phi}=-\alpha^2\sigma \hat{\eta}
$$
thus
$$
\hat{\eta}=-\frac{s}{\alpha^2\sigma}B\cosh(\alpha L)
$$
combining the two equations we get the *dispertion relation* for this system
$$
s^2=-\alpha^3\tanh(\alpha L)\sigma
$$
thus
$$
s=\pm i\sqrt{\sigma \alpha^3\tanh(\alpha L)}
$$
So we see that the eigenvalues are purely imaginary (neutral modes) and oscillate faster and faster as the wavenumber $\alpha$ increases (the wavelength $\lambda$ decreases). When the domain is very deep, $L$ is large compared to $2\pi/\alpha$ so $\alpha L$ becomes large, we have the limit of deep water waves 
$$
s=\pm i \sqrt{\sigma \alpha ^3}
$$

To validate our computations, we plot the spectrum (the eigenvalues plotted in the complex plane) and we add this theory. For this we have to chose the values of $\alpha$ of the wavelike theory that are relevant for our rectangular box. We say that thanks to the boundary conditions on $\eta$, there can be half a wavelength in *Lx*, one full wavelength and so on, that is
$$
\lambda=L_x/(k/2), k=1,2,3 ...
$$
We show this in the left plot of the figure below, showing a good agreement between the theory and the computation. To get even beter agreement, please increase the grid resolution (the computation is a little longer and there will be too much vectors on the velocity field plot). 

%}

% validation
subplot(1,4,1)
lambda=Lx./([1:15]/2); 
alpha=2*pi./lambda;

stheo=i*sqrt(sigma*alpha.^3.*tanh(alpha*Ly))';
plot(real(s),imag(s),'b.',real(stheo),imag(stheo),'ro');
axis([-1,1,-10,100]);
xlabel('real part of eigenvalue'); ylabel('imaginary part of eigenvalue');
title('spectra'); legend('numeric','theory','location','north')
grid on

%{
# Velocity fields

It is nice to see what the flow looks like for the different modes of oscillation of our fluid/interface system. We need to get the velocity field and also the connection with the interface. We plot both the real part (in blue) and the imaginary part (in red). 

You can recognize directly by eye that these configurations corresponds to stable oscillations because at the top of the domain, the velocity field is in the direction opposite to the deflection of the free surface: the flow pushes back the surface toward rest (but it will overshoot and go to the other side, so the fluid will push back again and so on). If this was unstable, the flow would instead push the interface away from its rest position.

%}

% show velocity field and free surface
sel=[2,3,4];
for ind=1:length(sel);
    subplot(1,4,ind+1)
    
    % select eigenvector
    q=U(:,sel(ind)); 
    
    % extract u and v and reshape
    u=reshape(D.x*q(l.phi),Ny,Nx);
    v=reshape(D.y*q(l.phi),Ny,Nx);
    quiver(x,y,real(u),real(v),'b'); hold on 
    quiver(x,y,imag(u),imag(v),'r');
    
    % free surface
    e=q(l.eta);
    plot(x,real(e)+Ly,'b-',x,imag(e)+Ly,'r-')
    xlabel('x'); ylabel('y'); title(['mode ' num2str(sel(ind))]);
    axis equal; axis([0,Lx,0,2*Ly]);
    grid on


end
hold off

% print figure
set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','free_surface_2D.png');

%{
![The spectrum and the velocity field of the first three eigenmodes](free_surface_2D.png)

You can see on this figure that we have a numerical eigenvalue $s=0$, why so? it is just because $\phi$ comes into the equations only through its derivatives, so it is immaterial to the system if we add a constant to $\phi$. This means that if we build a state $q$ where $\phi$ is constant and $\eta$ is 0, then
$$
Aq=0
$$
so $s=0$ is the associated "eigenvalue". To remove this uninteresting (so-called "spurious") mode, we sould add a constraint on the state, just like we do in [peristalsis.m]() or [venturi.m]() for the pressure, saying that the value of $\phi$ somewhere in the domain should be zero.

# Exercices/Contribution

* Please change the boundary conditions on $\eta$ to have fixed ends of the interface $\eta(x=0)=0, \eta(x=Lx)=0$ and compare how it changes the eigenmodes of oscillations.
* Please do like in [free_surface_gravity_particles.m](), show the time evolution of the eigenmodes using the exponential in time. An then if you have the courage, advect also particles using this time dependent velocity field to give a realistic impression of the fluid/interface motion.
* Please code all this using Navier-Stokes instead of velocity potential (for Masters! ... master students?)
* Please code gravity waves instead of capillary waves
* Please using the gravity waves as above, look at the stable/unstable eigenmodes and re-evaluate what was said here about the flow pushing the interface back to rest (stable) or moving it away from rest (unstable).
* Please add a constraint to remove the spurious mode with $s=0$

%}
