%{
# Potential flow in 2D

This is a code to solve the equation for potential flow in a square. 


%}


clear all; clf

% parameters
Nx=10; % gridpoints in x
Ny=10; % gridpoints in y 
Lx=2; % domain size in x
Ly=2; % domain size in y
pts=5; % number of points in finite difference stencils
method='fd'; % discretization

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D(method,-1,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D(method,-1,Ly,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

%{
The equation for the velocity potential $\phi(x,y)$ is
$$
\phi_{xx}+\phi_{yy}=0.
$$
%}

% system matrices
A=D.xx+D.yy;
b=zeros(NN,1);

%{
We start with the known analystical solution of a potential flow
$$
\phi=xy
$$
thus $u(x,y)=y$ and $v(x,y)=x$.
%}

% initial guess
phi0=X.*Y; 

%{
The boundary conditions are the velocity normal to the walls. The components of the velocity are
$$u=\phi_x$$
and $$v=\phi_y.$$

We use the analytical expression of the flow to impose the value of the nonhomogeneous boundary conditions.
%}

% boundary conditions
loc=[l.top;l.left;l.bot;l.right];

C=[D.x([l.left;l.right],:); ...
   D.y([l.bot;l.top],:)];

A(loc,:)=C;
b(loc,:)=C*phi0(:);

% solve system
phi=A\b;

% validation
u=reshape(D.x*phi,Nx,Ny);
v=reshape(D.y*phi,Nx,Ny);

utheo=reshape(D.x*phi0(:),Nx,Ny);
vtheo=reshape(D.y*phi0(:),Nx,Ny);

% show velocity field
subplot(1,2,1);
quiver(X,Y,u,v,'b');
hold on
quiver(X,Y,utheo,vtheo,'r');
axis equal; axis([-1 -1+Lx -1 -1+Ly])
grid on; 
xlabel('x'); ylabel('y'); title('velocity field')

% show error
subplot(1,2,2);
phi=phi-phi(1)+phi0(1);
mesh(X,Y,abs(reshape(phi,Ny,Nx)-phi0));
xlabel('x'); ylabel('y'); title('phi error')

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','potential.png');

%{

The error is close to machine accuracy even thought the resolution is low since the solution $\phi$ is linear.
The two velocity fields (numerical and computed) are on top of each other.
![The figure](potential.png)
%}
