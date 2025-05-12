%{

# Kelvin-Helmholtz instability of a shear layer

This is probably the first instability of fluid mechanics. Helmholtz wrote "the surfaces of discontinuity behave to some extent like bodies in unstable equilibrium". These "surfaces of discontinuity" are his discontinuous idealization of the progressive zone through which the velocity goes from one constant to an other one (the "shear layer" or also the "mixing layer".

This code resembles very much [poiseuille_uvp.m](), except for two things: 

* the base flow is now *tanh(y)* 
* we use a differentiation based on Hermite polynomials instead of Chebychev polynomial.

Chebychev polynomials have dense collocation nodes (gridpoints) close to the boundaries so they are good for boundary value problems. On the other hand, Hermite polynomial collocation nodes are dense in the center, precisely where we have our shear layer. Hermite polynomials naturally assume that the variable tends to zero at infinity, so there is no boundary condition to be enforced. Hermite polynomials go from $-\infty$ to $+\infty$, so there are no "boundaries" indeed.

Dependency:

* [herdif.m]() That builds the differentiation matrices based on Hermite polynomials
* [herroots.m]() That is used by herdif.m
* [poldif.m]() That is used by herdif.m

%}


clear all; clf;

alpha=0.5;    % the wave number
L=30;        % the height of the domain, from -1 to 1
Re=10000;    % the Reynolds number
N=100;      % the number of grid points

%{

In herdif, the first input is the number of nodes, the second one is the desired order of derivation and the third one is a scaling parameter for the distance between the farthest nodes. Make it smll for a "large domain".

%}

% Hermite differentiation matrices
[y, DM] = herdif(N,2,1);
D=DM(:,:,1);
DD=DM(:,:,2);
Z=zeros(N,N); I=eye(N); 

% renaming the differentiation matrices
dy=D; dyy=DD;
dx=i*alpha*I; dxx=-alpha^2*I;
Delta=dxx+dyy;

% base flow
U=tanh(y); 
Uy=(1-tanh(y).^2);

%{
# System matrices

The matrices $E$ and $A$ are the same as for [poiseuille_uvp.m]().

%}

% the matrices
S=-diag(U)*dx+Delta/Re;
A=[ ...
    S, -diag(Uy), -dx; ...
    Z, S, -dy; ...
    dx, dy, Z];
E=blkdiag(I,I,Z);

% computing eigenmodes 
[UU,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); UU=UU(:,o);
rem=abs(s)>1000; s(rem)=[]; UU(:,rem)=[];

%{
# Validation

We validate our computation based on figure 4.27 page 238 in the book "Drazin and Reid, Hydrodynamic stability, Cambridge university press". Here we use the wavelike formulation
$$
u(x,y,t)=\hat{u}(y)\exp(i\alpha x+st)
$$
but the classical one (used in the book) is
$$
u(x,y,t)=\hat{u}(y)\exp(i\alpha(x-ct))
$$
such that the real part of $c$ is the wave speed, and the exponential growth rate is the imaginary part of $c$. This is what is plotted on the graph.

%}

% loading the scanned figure for comparison and plotting
subplot(1,2,1);
a=imread('kelvin_helmholtz_ref.png'); ss=size(a);
xx=linspace(0,5,ss(2));
yy=linspace(01,-4,ss(1));
image(xx,yy,a); axis xy; hold on

% validation
plot(alpha,real(s(1))/alpha,'r.','markersize',20)
grid on
xlabel('wavenumber alpha'); 
ylabel('exponential growth rate/alpha')
title('Validation');

%{
# Velocity field

We show the velocity field of the unstable eigenmode and as well the associated pressure field as a colormap.
%}

% showing the velocity field
Nx=20;
u=1:N; v=u+N; p=v+N; 
q=UU(:,1); 
Lx=2*pi/alpha;  x=linspace(-Lx/2,Lx/2,Nx);

% expand to physical space
qphys=2*real(q*exp(i*alpha*x));

% add the base flow to the perturbations
uu=qphys(u,:);
vv=qphys(v,:);
pp=qphys(p,:);

% show the velocity field
sely=1:2:N;
subplot(1,2,2);
quiver(x,y(sely),uu(sely,:),vv(sely,:),'k'); hold on
surf(x,y,pp-10,'facealpha',0.5); shading interp;
axis([x(1),x(end),y(1),y(end)]);
xlabel('x'); ylabel('y'); title('Velocity field of the initial condition');

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','kelvin_helmholtz_hermite.png');

%{

![The figure](kelvin_helmholtz_hermite.png)


# Exercices/Contributions

* Please find a way to validate the figure of the book as well for $\alpha>1$ where the mode becomes stable
* Please find the book and try to understand what is the second curve below and compare it to our code
* Please do the comparison of our code to the following figure in the book where the is the neutral curve for viscous shear layers, showing the dividing line between unstable wavelengthes and stable ones as the Reynolds number is changed
* Please use the eigenmode from this code as an initial condition in Gerrris and validate the coparison linear/nonlinear ----------->[kelvin_helmholtz_gerris.m]()


%}