%{
# Pipe flow stability using symmetry

This is to compute the eigenmodes of the Hagen-Poiseuille flow (cylindrical pipe with flow due to an axial pressure gradient). In this version of the code, we build differentiation matrices that account for the symmetry of the different variables and we treat only half of the pipe diameter (for positive $r$). This is good to avoid troubles with factors like $1/r$ in the Navier-Stokes equations in axisymmetric form.

%}

clear all; clf; format compact

% parameters
N=50;      % number of grid points  
Re=5000;     % Reynolds number
n=1;        % azymuthal wavenumber
k=1;        % axial wavenumber

% differentiation matrices
scale=1; 
[r,DM] = chebdif(2*N,2);
dr=DM(:,:,1)/scale;
drr=DM(:,:,2)/scale^2;
r=r*scale;
Z=zeros(N,N); I=eye(N);

%{
The radial coordinate is $r$, the azymuthal one is $\theta$ (thus the postfix $t$ for the differentiation matrix in $\theta$), and the axial coordinate is $z$. The velocity profile has only components in $z$ and is called $W$. here we have the easy differentiation matrices in $\theta$ and in $z$, which we get from the wavelike assumption in these directions (as opposed to $r$ where the quantities are not wavelike)
%}

dt=i*n*I; dtt=-n^2*I;
dz=i*k*I; dzz=-k^2*I;

%{
# Enforcing symmetry

Here is the core of the technical part of this code where we tell in fact that we will only deal with the unknown for positive $r$, and build matrices differently for the variables that are even (symmetric) or odd (antisymmetric) with respect the the axis. Here we allow to compute for different azymuthal wavenumbers *n*, and unfortunately the symmetry of the variables changes depending on the parity of $n$, but this is all accounted here using the variable *s*.

For instance for $n=0$ (the perturbation does not depend on $\theta$), the radial velocity is even, so its values for the gridpoints with $r>0$ are equal to its values for the gridpoints with $r<0$. So to get the radial derivative we do (here the ABCD matrix is the radial differentiation matrix) 
$$
\left(\begin{array}{c}
u_r(r>0) \\ u_r(r<0) 
\end{array}\right)
=
\left(\begin{array}{cc}
A &B \\ C & D  
\end{array}\right)
\left(\begin{array}{c}
u(r>0) \\ u(r<0) 
\end{array}\right)
$$
and since we are only interested in the values for $r>0$ w can write
$$
u_r(r>0)=Au(r>0)+Bu(r<0)
$$
and since on the grid, $u(r<0)$ is just the same as the vector $u(r>0)$ but fliped up-down, we just need to flip up-down the sub-matrix $B$ and multiply it to $u(r>0)$, thus
$$
u_r(r>0)=(A+flipud(B))u(r>0)
$$
where $flipud$ is a function to do the up-down flip of the matrices. In fact in the code, instead of using this function, I flip directly by flipping the vector of indices *sn* which then serves both to select the $B$ matrix and to flip it.

For odd functions the only difference is that the values for $r<0$ must be fliped but also negated. This is done using the *s* variable in the code. 
%}

% enforce symmetry about axis
sp=1:N; sn=2*N+1-sp; % selection and flipping vectors
s=(-1)^mod(n,2); % symmetry changes depending on n

dre=dr(sp,sp)+s*dr(sp,sn);
drre=drr(sp,sp)+s*drr(sp,sn);
dro=dr(sp,sp)-s*dr(sp,sn);
drro=drr(sp,sp)-s*drr(sp,sn);
r=r(sp); 

% base flow
W=1-r.^2; Wr=-2*r;
 
%{
Here we build the $1/r$ and $1/r^2$ factors that come often in the equations.
%}

% usefull
rm1=diag(1./r);
rm2=diag(1./(r.^2));

% Laplacian
lape=drre+rm1*dre+rm2*dtt+dzz;
lapo=drro+rm1*dro+rm2*dtt+dzz;

%{
Please here describe the Jacobian for the axisymmetric Navier--Stokes equations.
%}

% System matrices
A=[-diag(W)*dz+(lapo-rm2)/Re, -2*rm2*dt/Re, Z, -dre; ...
  2*rm2*dt/Re, -diag(W)*dz+(lapo-rm2)/Re, Z, -rm1*dt; ...
  -diag(Wr), Z, -diag(W)*dz+lape/Re, -dz; ...
  rm1+dro, rm1*dt, dz, Z];

E=blkdiag(I,I,I,Z);

% Boundary conditions
II=eye(4*N);
loc=[1,N+1,2*N+1]; 
A(loc,:)=II(loc,:); 
E(loc,:)=0; 

% computing eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

%{
# Validation

We compute the eigenmodes and we compare this to figure 3.3 in "Stability and transition in shear flows, Springer, by P Schmid and D henningson", page 67 of my version. The fit is good.
%}

% loading the scanned figure for comparison
a=imread('pipe_spectra.png'); ss=size(a);
x=linspace(0,1,ss(2)); y=linspace(0,-1,ss(1));
image(x,y,a); axis xy; 
hold on

plot(-imag(s)/k,real(s),'ro')
xlabel('wave speed'); ylabel('exponential growth rate')
grid on; hold on

set(gcf,'paperpositionmode','auto')
print('-dpng','-r100','pipe_sym.png')

%{
![](pipe_sym.png)

# Exercices/Contributions

* Please find some litterature to compare for other *n* than 1
* Please write a version of the code without accounting for symmetry in the differentiation matrices, having a computational domain from $r=0$ to $1$. You will have to be carefull when imposing the boundary conditions at $r=0$, to account correctly the symmetry of the different variables.
* Please do the same as above with $r$ from $r=-1$ to $1$ (but don't put a grid point on the axis where $1/r$ is infinite, use an even number of grid points!) ---> [pipe.m]()
%}


