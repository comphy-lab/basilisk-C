%{
# Free surface with no-slip wall boundary

We use the code of free surface (with surface tension) [free_surface_2D.m](/sandbox/easystab/free_surface_2D.m) and modify the boundary condition to see what happens to the eigenmodes of ocsillations if we impose no-slip wall boundary at the ends of the interface $\eta(x=0)=0, \eta(x=Lx)=0$ instead of zero derivative boundary condition in the original code.

The theory and the matrix systems are the same as the original code, here we just modify the boundary condition of the problem. 

Dependency:

* [chebdif.m](sandbox/easystab/chebdif.m) for the Chebychev differentiation matrices
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

%1D differentiation matrices
scale=-2/Lx;
[x,DM] = chebdif(Nx,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=x/scale; x=x-x(1); 
intx=([diff(x)',0]+[0,diff(x)'])/2;

scale=2/Ly;
[y,DM] = chebdif(Ny,2); 
dy=DM(:,:,1)*scale;    
dyy=DM(:,:,2)*scale^2;    
y=y/scale; y=y-y(end); 
inty=([diff(y)',0]+[0,diff(y)'])/2;

% 2D differentiation matrices
Dx=kron(dx,eye(Ny));
Dxx=kron(dxx,eye(Ny));
Dy=kron(eye(Nx),dy);
Dyy=kron(eye(Nx),dyy);
[X,Y]=meshgrid(x,y);


% vectors for coordinate selections 
NN=Nx*Ny;
dom=reshape(1:NN,Ny,Nx);
top=dom(1,2:end-1); top=top(:);
bot=dom(end,2:end-1); bot=bot(:);
left=dom(:,1); left=left(:);
right=dom(:,end); right=right(:);
phi=(1:NN)';
eta=((NN+1):(NN+Nx))';

% useful matrices
Z=zeros(NN,NN); ZZ=zeros(NN+Nx,NN+Nx);    zx=zeros(Nx,Nx);
I=eye(NN);      II=eye(NN+Nx);            ix=eye(Nx);

% system matrices
A1=blkdiag(Dxx+Dyy,zx);
A1(eta(2:end-1),phi)=Dy(top,:);
A2=blkdiag(Dxx+Dyy,zx);
A2(eta(2:end-1),phi)=Dy(top,:);
E=blkdiag(Z,ix);

%{
# Boundary conditions

In the first case, we tell that it should have a zeros derivative at $x=0$ and $x=Lx$. This is good here because we can write by hand in [#validation]() a theory for the validation of this case. This is writen in the third row of *c* where we use the 1D differentiation matrix in $x$ *dx*, and selecting only the first and last rows of it

> dx([1,Nx],:)*II(eta,:)

Secondly, to impose the fixed ends for the interface, we will select the first and last elements of vector $eta$ 

%}

% boundary conditions
loc=[top;left;bot;right;eta([1,2,Nx])];

c1=[Dx([left;right],:)*II(phi,:); ...
   Dy(bot,:)*II(phi,:); ...
   dx([1,Nx],:)*II(eta,:); ...      %zero derivative boundary condition
   intx*II(eta,:)];

Ca1=[c1; ...
    sigma*dxx(2:end-1,:)*II(eta,:)];

c2=[Dx([left;right],:)*II(phi,:); ...
   Dy(bot,:)*II(phi,:); ...
   II(eta([1,Nx],:),:);...     %No-slip wall condition
   intx*II(eta,:)];

Ca2=[c2; ...
    sigma*dxx(2:end-1,:)*II(eta,:)];

Ce=[0*c1; ...
    rho*II(top,:)]; 

A1(loc,:)=Ca1;
A2(loc,:)=Ca2;
E(loc,:)=Ce;

%{
# Eigenmodes

There is no dissipation and no instability in this flow, so all the modes should be neutral (the real part of the eigenvalues $s$ should be 0), but there are oscillations so the imaginary parts should not be zero). The spectrum should be right-left symetric because the system itself is, so we remove the eigenmodes with negative imaginary parts.
%}

% compute eigenmodes for zeroo derivative boundary condition (A1,E)
disp('computing eigenmodes');
[U,S]=eig(A1,E);
s=diag(S);  [t,o]=sort(abs(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
rem=imag(s)<0; s(rem)=[]; U(:,rem)=[];

% validation
figure(1)
subplot(1,2,1)
lambda=Lx./([1:15]/2); 
alpha=2*pi./lambda;

stheo=i*sqrt(sigma*alpha.^3.*tanh(alpha*Ly))';
plot(real(s),imag(s),'b.',real(stheo),imag(stheo),'ro');
axis([-1,1,-10,100]);
xlabel('real part of eigenvalue'); ylabel('imaginary part of eigenvalue');
title('spectra'); legend('numeric','theory','location','north')
grid on

% compute eigenmodes for no-slip wall boundary condition (A2,E)
disp('computing eigenmodes');
[U,S]=eig(A2,E);
s=diag(S);  [t,o]=sort(abs(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
rem=imag(s)<0; s(rem)=[]; U(:,rem)=[];

% validation
subplot(1,2,2)
lambda=Lx./([1:15]/2); 
alpha=2*pi./lambda;

plot(real(s),imag(s),'b.');
axis([-1,1,-10,100]);
xlabel('real part of eigenvalue'); ylabel('imaginary part of eigenvalue');
title('spectra'); legend('numeric')
grid on


% print figure
set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','free_surface_noslip_boundary_condition.png');

%{
![Eigenmodes of the zero derivative boundary condition (left) and the fixed ends boundary condition (right)](/sandbox/easystab/stab2014/free_surface_noslip_bc.png)

You can see on this figure that we have a numerical eigenvalue $s=0$, why so? it is just because $\phi$ comes into the equations only through its derivatives, so it is immaterial to the system if we add a constant to $\phi$. This means that if we build a state $q$ where $\phi$ is constant and $\eta$ is 0, then
$$
Aq=0
$$
so $s=0$ is the associated "eigenvalue". To remove this uninteresting (so-called "spurious") mode, we sould add a constraint on the state, just like we do in [peristalsis.m]() or [venturi.m]() for the pressure, saying that the value of $\phi$ somewhere in the domain should be zero.
%}

