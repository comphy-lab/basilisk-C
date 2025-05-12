%{

# The Rayleigh-Plateau instability

**This is an ongoing work! Not yet validated!**

This is the instability of an infinite cylinder of liquid due to surface tension. The system is unstable for wavelengthes longer than its perimeter.  

The code is based on [pipe_sym.m](), but I have changed the names of the different velocity components. Here there is no base flow. I also use a way like in that code to enforce symmetry about the axis of the cylinder, so that I only need to enforce the boundary conditions at the free-surface.

%}

clear all; clf; format compact

% parameters
N=80;      % number of grid points  
Re=1000;     % Reynolds number
alpha=0.5;        % axial wavenumber
sig=1; % surface tension

% differentiation matrices
[d.y,d.yy,d.wy,y]=dif1D('cheb',-1,2,2*N);

% enforce symmetry about axis (even or odd)
sp=1:N; sn=2*N+1-sp; % selection and flipping vectors
d.ye=d.y(sp,sp)+d.y(sp,sn);
d.yye=d.yy(sp,sp)+d.yy(sp,sn);
d.yo=d.y(sp,sp)-d.y(sp,sn);
d.yyo=d.yy(sp,sp)-d.yy(sp,sn);
y=y(sp); 

Z=zeros(N,N); I=eye(N); 
d.x=i*alpha*I; d.xx=-alpha^2*I;

%{
# The state vector $q$

The state is the axial and radial velocity components $U$ and $v$. And there is as well pressure $p$ and interface displacement $\eta$ (which is a scalar wheres $u$ and $v$ are defined as a function of $y$).
$$
q=\begin{pmatrix}
u\\v\\p\\\eta\end{pmatrix}
$$

%}
% useful
l.u=(1:N)'; l.v=l.u+N; l.p=l.v+N; l.e=l.p(end)+1; % location vectors
II=eye(3*N+1);
Iu=II(l.u,:); Iv=II(l.v,:);Ip=II(l.p,:);Ie=II(l.e,:); % selection matrices

% Laplacian
lape=d.xx+d.yye+diag(1./y)*d.ye;
lapo=d.xx+d.yyo+diag(1./y)*d.yo;

%{
# System matrices


The system equations are the Navier-Stokes equations in cylindrical coordinates linearized about a zero base flow, and the advection of the interface $\eta$ by the flow (which once inearized about a zero ase flow reduces to the vertical advection of $\eta$ as $\eta_t=v$) $Eq_t=Aq$
$$
\left(
\begin{array}{rcl}
u_t&=&-p_x+(u_{xx}+u_{yy}+u_y/y)/Re\\
v_t&=&-p_y+(v_{xx}+v_{yy}+v_y/y-v/y^2)/Re\\
0&=&u_x+v_y+v/y\\
\eta_t&=&v|_1
\end{array}
\right.
$$
%}

% System matrices
A=[lape/Re*Iu-d.x*Ip; ...
  (lapo-diag(1./y.^2))/Re*Iv-d.ye*Ip; ...
  d.x*Iu+(d.yo+diag(1./y))*Iv; ...
  Iv(N,:); ...
  ];
E=[Iu; Iv; 0*Ip; Ie];

%{
# Boundary conditions

They are: zero tangential stress at the free surface, and also at the free surface: normal stress equal to the pressure jump due to the curvature of the interface
$$
\left(\begin{array}
0&=&u_y|_1+v_x|_1 \\
0&=&p-\sigma \left[\frac{1}{\eta(1+\eta_x^2)^{1/2}}-\frac{\eta_xx}{1+\eta_x^2}\right]
\end{array}\right.
$$
The first one is already linear, and we need to linearize the second one about an interface at $y=1+\eta$, which gives for the perturbations
$$
\left(\begin{array}
0&=&u_y|_1+v_x|_1 \\
0&=&p+\sigma (1-\alpha^2)\eta
\end{array}\right.
$$

%}

% Boundary conditions
loc=[l.u(N);l.v(N)];
C=[d.ye(N,:)*Iu+d.x(N,:)*Iv; ... % uy+vx=0 at free surface
    Ip(N,:)+sig*(1-alpha^2)*Ie]; % jump in normal stress 

A(loc,:)=C; 
E(loc,:)=0; 

% computing eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

s(1:10)

% % loading the scanned figure for comparison
% a=imread('pipe_spectra.png'); ss=size(a);
% x=linspace(0,1,ss(2)); y=linspace(0,-1,ss(1));
% image(x,y,a); axis xy; 
% hold on
% 
% plot(-imag(s)/k,real(s),'ro')
% xlabel('wave speed'); ylabel('exponential growth rate')
% grid on; hold on

set(gcf,'paperpositionmode','auto')
print('-dpng','-r100','rayleigh_plateau.png')

%{
![](pipe_sym.png)

# Exercices/Contributions

* Please
* Please
* Please

%}
