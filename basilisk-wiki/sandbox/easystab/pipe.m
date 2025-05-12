%{
# Stability of the pipe flow

Here is the same code as [pipe_sym.m]() but with double the number of points because we do not explicitely account for symmetry. So the code is simpler but has double as large matrices.

For explanations on the technicalities see, [pipe_sym.m]().

Please keep an even number of gridpoints *N* so that you won't have a gridpoint at $r=0$ where $1/R$ is ill-behaved.
%}

clear all; clf; format compact

% parameters
N=100;      % number of grid points  
Re=5000;     % Reynolds number
n=1;        % azymuthal wavenumber
k=1;        % axial wavenumber

% differentiation matrices
scale=1; 
[r,DM] = chebdif(N,2);
dr=DM(:,:,1)/scale;
drr=DM(:,:,2)/scale^2;
r=r*scale;
Z=zeros(N,N); I=eye(N);

dt=i*n*I; dtt=-n^2*I;
dz=i*k*I; dzz=-k^2*I;

% base flow
W=1-r.^2; Wr=-2*r;
  
% usefull
rm1=diag(1./r);
rm2=diag(1./(r.^2));

% Laplacian
lap=drr+rm1*dr+rm2*dtt+dzz;

% System matrices
A=[-diag(W)*dz+(lap-rm2)/Re, -2*rm2*dt/Re, Z, -dr; ...
  2*rm2*dt/Re, -diag(W)*dz+(lap-rm2)/Re, Z, -rm1*dt; ...
  -diag(Wr), Z, -diag(W)*dz+lap/Re, -dz; ...
  rm1+dr, rm1*dt, dz, Z];

E=blkdiag(I,I,I,Z);

% Boundary conditions
II=eye(4*N);
loc=[1,N,N+1,2*N,2*N+1,3*N]; 
A(loc,:)=II(loc,:); 
E(loc,:)=0; 

% computing eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

% loading the scanned figure for comparison
a=imread('pipe_spectra.png'); ss=size(a);
x=linspace(0,1,ss(2)); y=linspace(0,-1,ss(1));
image(x,y,a); axis xy; 
hold on

plot(-imag(s)/k,real(s),'ro')
xlabel('wave speed'); ylabel('exponential growth rate')
grid on; hold on

%{
# Exercices/Contributions

* Please look at the suggestions in [pipe_sym.m]()

%}