%{
# Free-surface: displaying the eigenmodes using particles

In this code, we take the code that computes the eigenmodes for free surface waves [free_surface_gravity.m]() and we visualize the flow motion of the first eigenmode by expanding the solution to physical space and advecting tracer particles, like we did in a simple case in [particles.m]().
%}

clear all; clf;

n=100;      % number of gridpoints
alpha=1;    % wavenumber in x
L=2;        % Fluid height in y
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

% System matrices
A=[mu*Delta, Z, -dx, Z(:,1); ...
   Z, mu*Delta, -dy, Z(:,1); ...
   dx, dy, Z, Z(:,1); ...
   Z(1,:),I(n,:),Z(1,:),0];

E=blkdiag(rho*I,rho*I,Z,1);

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
# The animation

We show the position of the free surface stored in *q(eta)*, and the particles are initialy set as a mesh, and stretched in *y* to fit with the initial position of the free surface for the top prticles.
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% particle animation of the eigenmodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
nper=1      % number of periods of oscillation
nt=30       % number of time steps per period
nx=40       % number of points in x
modesel=2;  % which more do animate

% select the eigenmode
u=1:n; v=u+n; p=v+n; eta=3*n+1;
q=U(:,modesel); 
lambda=s(modesel);

% The time and x extent
tvec=linspace(0,nper*2*pi/abs(lambda),nper*nt); 
dt=tvec(2)-tvec(1);
Lx=2*pi/alpha; 
x=linspace(0,Lx,nx);

% scale mode amplitude 
q=0.05*q/abs(q(eta)); 

% initialize tracer particles
[px,py]=meshgrid(linspace(0,Lx,60),linspace(0,L,30)); 
py=py.*(1+2*real(exp(lambda*tvec(1))*q(eta)*exp(i*alpha*px))/L);
px=px(:);py=py(:);  

% time loop
for ind=1:nper*nt
    
    
%{
# Expand the solution to physical space

We have assumed a wavelike behaviour of the flow solution in the $x$ direction, thus we have writen for all the flow variables, for instance for the $u$ velocity
$$
u(x,y,t)=\hat{u(y)}\exp(\alpha x+s t)+CC
$$
where $CC$ is the complex conjugate of the other term. This way we have transformed a variable $u$ depending on three coordinates, to a variable $\hat{u}$ depending on one single coordinate, thus transforming everything into a 1D problem. What we do now is to transform back the 1D problem into a full original problem by using the exponential. 

A complex quantity plus its complex conjugate gives twice the real part of it, this is what we code here.

In the code, we do the expansion to hysical space for all the variables $(u,v,p,\eta)$ at once in a single command.

%}
    % expand mode to physical space 
    qq=2*real(exp(lambda*tvec(ind))*q*exp(i*alpha*x));
        
    % plot pressure
    surf(x,y,qq(p,:)-10,'facealpha',0.3); view(2); shading interp; hold on
    
    % plot free surface
    plot(x,L+qq(eta,:),'k-',x,0*x+L,'k--'); hold on
    
    %%%% plot the particles
    plot(mod(px,Lx),py,'k.');

%{
# Interpolation of the particle positions

The particles are not on a cartesian grid, they can be wherever the field has brought them, the only thing that matter here for the interpolation is where they are on the $y$ grid, and then we can do the expansion in $x$ as a second step, this is less computations.
%}
    % compute velocity at position of particles
    pu=interp1(y,q(u),py);
    pv=interp1(y,q(v),py);

%{
# The Taylor expansion

When we wrote the boundary conditionat the free surface for setting up the linear model, we have used Taylor expansions to express the velocity and pressure at the position of the moving interface as a function of the values at the fixed position $L$. here we do again the same thing for the particles that are above $L$.

%}

    % For particles above L, use Taylor expansion for velocity
    pu(py>L)=q(u(n))+q(eta)*dy(n,:)*q(u);
    pv(py>L)=q(v(n))+q(eta)*dy(n,:)*q(v);
    
    % expand to physical space
    puu=2*real(exp(lambda*tvec(ind))*pu.*exp(i*alpha*px));
    pvv=2*real(exp(lambda*tvec(ind))*pv.*exp(i*alpha*px));
 
    % advect particles
    px=px+puu*dt;
    py=py+pvv*dt;
    
    xlabel('x');    
    ylabel('y'); 
    axis equal; axis([0,Lx,0,1.3*L]); 
    grid off
    hold off
    drawnow
end

set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','free_surface_gravity_particles.png');

%{
![The figure at final iteration of the time loop, shwing the free surface at the top, the particles in black, and the disturbance pressure field in color ](/sandbox/easystab/free_surface_gravity_particles.png)

# Exercices/Contributions

* Please make a movie instead of just an image to display on the page --> [free_surface_gravity_particles_gif.m]()
* Please play with the parameters of the particles, and the amplitude of the eigenmode
* Please check that the first and second eigenmodes correspond to waves propagating in opposite directions -->[free_surface_gravity_first_modes.m](/sandbox/easystab/stab2014/free_surface_gravity_first_modes)
* Please use this animation to check the theory for the wave speed of gravity waves
* Please plot also the vector field using quiver in a way that confirms the motion of the particles that we observe -->[free_surface_gravity_vector_field.m](/sandbox/easystab/stab2014/free_surface_gravity_vector_field)
* Please investigate the motion of the other eigenmodes (which correspond to damped stationnary waves, for now, there is a bug that makes things strange for the non propagative modes...)
* Please consider the unstable case with *g* negative and find the nice scaling of the mode amplitude and the time of the animation so that we nicely see the flow in this instability
* Please investigate how the Taylor expansino is not accurate when the surface deformation is large (you can see this with the particles that are close to the interface)
* Please do the animation showing the total pressure instead of only the pressure disturbance (that is, add the hydrostatic pressure field as well to see what this changes)[free_surface_gravity_total_pressure](/sandbox/easystab/stab2014/free_surface_gravity_total_pressure)
* Please do the same visualisation with particle for the other flow instabilities
%}

