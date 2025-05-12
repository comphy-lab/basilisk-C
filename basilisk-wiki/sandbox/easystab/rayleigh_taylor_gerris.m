%{
# Free-surface waves with two fluids: gravity and surface tension

This is basically the same thing than [free_surface_gravity.m]() but here there are two fluids on top of each other, so this is a little more complicated, but we do not have to neglect the density of the top fluid so this is more general. Technically it works like that, we build two dynamic systems for the two fluid layers with different physical properties, then we need to connect this two fluid layers through the moving interface.

We take a very low viscosity because we will compare with a theory for a perfect fluid.
%}

clear all; 

Lx=4;      % wavelength in x
alpha=2*pi/Lx;    % wavenumber in x
g=1;        % gravity
sigma=0.1;    % surface tension

% for domain 1 (top)
N1=50;      % number of gridpoints
L1=2;        % Fluid height in y
rho1=1;      % fluid density
mu1=0.00001;    % fuid viscosity

% for domain 2 (bottom)
N2=50;      % number of gridpoints
L2=2;        % Fluid height in y
rho2=0.5;      % fluid density
mu2=0.00001;    % fuid viscosity


% domain 1
% differentiation matrices
scale=-2/L1;
[y,DM] = chebdif(N1,2); 
D=DM(:,:,1)*scale;    
DD=DM(:,:,2)*scale^2;    
y1=(y-1)/scale+L2; 
I1=eye(N1); Z1=zeros(N1,N1);
% renaming the matrices
dy1=D; dyy1=DD;
dx1=i*alpha*I1; dxx1=-alpha^2*I1;
Delta1=dxx1+dyy1;

% domain 2
% differentiation matrices
scale=-2/L2;
[y,DM] = chebdif(N2,2); 
D=DM(:,:,1)*scale;    
DD=DM(:,:,2)*scale^2;    
y2=(y-1)/scale; 
I2=eye(N2); Z2=zeros(N2,N2);
% renaming the matrices
dy2=D; dyy2=DD;
dx2=i*alpha*I2; dxx2=-alpha^2*I2;
Delta2=dxx2+dyy2;



%{
# System matrices

The motion of each fluid is described by the Navier-Stokes equations and continuity, linearized about a zero base flow, here written in a matrix form 
$$
e_1q_{1t}=a_1q_1
$$ 
for the top fluid and  
$$
e_2q_{2t}=e_2q_2
$$
for the bottom fluid. We also have 
$$
\eta_t=v_1(0)
$$ 
for the advection of the interface by the velocity next to the interface. Note that we could have as well written $\eta_t=v_2(0)$ since anyway we will see for the boundary conditions that $v$ is continuous across the interface so $v_1(0)=v_2(0)$.
The total state of the system is now called $q$
$$
q=\left(\begin{array}{c}
q_1 \\ q_2 \\ \eta
\end{array}\right)
$$
with all the information: both of the top and bottom fluids and the interface position $\eta$. Thus the system is $Eq_t=Aq$ with
$$
E=\left(\begin{array}{ccc}
\rho&0&0\\
0&\rho&0\\
0&0&1
\end{array}\right)
, \quad
A=\left(\begin{array}{ccc}
e1 & 0 & 0\\
0 & e2 & 0\\
I_{v_1(0)}  & 0 & 0\\
\end{array}\right)
$$
where $I_{v_1(0)}$ is the vector that multiplied by $q_1$ gives the value of $v_1$ at the interface (at $y=0$). This way we have in the dynamics the motion fo fluids 1 and 2, and the motion of the interface.

Here a little discussion on the different ways to build a matrix. What we have used here  is a way to build the matrices by stacking blocks together. The function *blkdiag* means "block diagonal", it build a larger matrices by stacking the input matrix arguments together to build a block diagonal matrix as output. An other way is like is done below for *a1* and *a2* is by concatenation of the differentiation matrices, the matrix full of zeros *Z* and the identity matrix *I*. 

We can also access the elements of the big matrix using the indices of submatrices. In most of the codes we use also this, especially for the boundary conditions. So we build here the vectors of indices that correspond in the big matrix to the different components of the state: *u1=1:N1* means that the *u* velocity of the fluid 1 is stored on the elements from 1 to *N1* of the state $q$, and so on. Then we see also that the position $\eta$ of the interface is the last element of the state $q$. 

So this is why when we want to insert the connection between $\eta$ and the vertical velocity at the bottom of fluid 1, we write 

> *A(eta,v1(1))=1;*

We will make use even more of these vectors of indices when implementing the boundary conditions below.

%}

% location vectors
u1=1:N1; v1=u1+N1; p1=v1+N1;
u2=p1(end)+(1:N2); v2=u2+N2; p2=v2+N2;
eta=p2(end)+1;

% System matrices
a1=[mu1*Delta1, Z1, -dx1; ...
   Z1, mu1*Delta1, -dy1; ...
   dx1, dy1, Z1];

a2=[mu2*Delta2, Z2, -dx2; ...
   Z2, mu2*Delta2, -dy2; ...
   dx2, dy2, Z2];

e1=blkdiag(rho1*I1,rho1*I1,Z1);
e2=blkdiag(rho2*I2,rho2*I2,Z2);

A=blkdiag(a1,a2,0);
E=blkdiag(e1,e2,1);

% we add the equation for eta
A(eta,v1(1))=1;

%{
# Boundary conditions

At the wall we have the classical no-slip boundary conditions. This means that $u1$ and $v1$ must be zero at the top of domain 1. here the top correspond the the first index, so $u_1(1)=0, v_1(1)=0$. For the fluid 2, the wall is at the bottom, so $u_2(N2)=0, v_2(N2)=0$. These are our 4 homogeneous Dirichlet boundary conditions. 

We then have 4 additional constraints comming from the connection of the two fluids through the interface. We have the continuity of vertical velocity $v1(N1)=v2(1)$, the continuity of $u$: $u1(N1)=u2(1)$, and then we have the continuity of the stress, that is the continuity of the pressure  
$$
p_1(0)-\rho_1 g\eta=p_2(0)-\rho_2 g\eta 
 $$
where we have done Taylor expantions to express the pressure at $y=\eta$ using the pressure at $y=0$ like we did in the simpler case of [free_surface_gravity.m#boundary-conditions](). Now, instead of having the left hand side of this equation be zero (the atmospheric pressure), we have that the left hand side equals the right-hand side from the pressure of the other fluid at $y=0$.

The last constraint is the continuity of tangential stress
$$
\mu_1(u_{1y}(0)+v_{1x}(0))=\mu_2(u_{2y}(0)+v_{2x}(0).
$$

Since we have many variables, we should find a nice and simple way to write all these constraints with the least chance of making mistakes. For this we make an extensive use of identity matrices and location indices. 

The difficult thing with matrix representation of the boundary conditions is that we do not manipulate the state itself $q$, but we manipulate matrices that multiply this state vector. A way to work arround this difficulty is to use an identity matrix $II$ of the same size than the state $q$, and say that for instance here for the $u$ velocity of fluid 1

> $u_1$=q(u1)=II(u1,:)*q

and an other example for $\eta$

> $\eta$ = q(eta)= II(eta,:)*q

where $u1$ and $eta$ are the index vectors forthe location of $U_1$ and $\eta$ in the state vector $q$. here a little trickier example for the contiuity of the viscous stress where we need to compute the $y$ derivative of the velocity components

> $u_{1y}(0)$=dy(N,:)*II(u1,:)*q

where from right to left of the multiplications we first extract the $u_1$ velocity from the state $q$ using the identity matrix, then we multiply by the line number $N$ of the $y$ differentiation matrix to extract the derivative at position $N$, that is the $y$ derivative at $y=0$. Now all these can be combined to get in compact form the constraint for viscous stress which is coded below

> mu1*(dy1(1,:)*II(u1,:)+dx1(1,:)*II(v1,:))-mu2*(dy2(end,:)*II(u2,:)+dx2(end,:)*II(v2,:)) 

this line produces an array that multiplying $q$ from the left must give the result 0.
%}

% boundary conditions
dir=[u1(end),v1(end),u2(1),v2(1)]; % no-slip at both walls
loc=[dir,u1(1),v1(1),u2(end),v2(end)];
II=eye(3*N1+3*N2+1);

C=[II(dir,:); ... % no-slip at walls
   II(v1(1),:)-II(v2(end),:); ... % continuity of v across interface
   II(u1(1),:)-II(u2(end),:); ... % continuity of u across interface
  (II(p1(1),:)-rho1*g*II(eta,:))-(II(p2(end),:)-rho2*g*II(eta,:))+sigma*alpha^2*II(eta,:); ... % Continuity of pressure
   mu1*(dy1(1,:)*II(u1,:)+dx1(1,:)*II(v1,:))-mu2*(dy2(end,:)*II(u2,:)+dx2(end,:)*II(v2,:))]; % continuity of shear stress

E(loc,:)=0;  
A(loc,:)=C;

% compute eigenmodes
disp('compute eigenmodes');
[UU,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); UU=UU(:,o);
rem=abs(s)>1000; s(rem)=[]; UU(:,rem)=[];

% validation
alphavec=linspace(0,3,100);
stheo=sqrt((alphavec*(rho1-rho2)*g-alphavec.^3*sigma)./(rho1./tanh(alphavec*L1)+rho2./tanh(alphavec*L2)));

snum=real(s(1));
plot(alphavec,real(stheo),'r-',alpha,snum,'b.');
xlabel('alpha');ylabel('exponential growth rate'); 
title('validation')
grid on

%
%{
# Initial condition for Gerris

From now on, we prepare the 2D field for the initial condition of Gerris. *Nx* is the nimber o grid points in *x*. And here we select the first eigenmode, which is the least stable one. 
%}

% showing the velocity field
Nx=100;
q=UU(:,1); 
x=linspace(-Lx/2,Lx/2,Nx);

%{
Here we chose the relative amplitude of the eigenmode and the base flow (the base flow has max u velocity equal to 1). Chosing an amplitude too large means that the evolution will be nonlinear right from the beginning, which we don't want because we would like to compare the evolution of the energy in Gerris and see that as long as the perturbatino is small, it is the same as in this code. If we start on the other hand with an amplitude too small, grid resolution in Gerris and round-off error may disturb the evolution of our eigenmode.

To expand to physical space, we do like in [kelvin_helmholtz_hermite.m#velocity-field]().
%}
% scale mode amplitude 
q=0.05*q/abs(q(eta)); 

% expand each flow domain into the other domain using Taylor expantion
y=[y2(1:end-1);y1]; N=length(y);
qu1=[q(u1(1))+0*(y2(1:end-1)-L2)*dy1(1,:)*q(u1); q(u1) ];
qv1=[q(v1(1))+0*(y2(1:end-1)-L2)*dy1(1,:)*q(v1); q(v1) ];
qp1=[q(p1(1))+0*(y2(1:end-1)-L2)*dy1(1,:)*q(p1); q(p1) ];

qu2=[q(u2); q(u2(end))+0*(y1(2:end)-L2)*dy2(end,:)*q(u2); ];
qv2=[q(v2); q(v2(end))+0*(y1(2:end)-L2)*dy2(end,:)*q(v2); ];
qp2=[q(p2); q(p2(end))+0*(y1(2:end)-L2)*dy2(end,:)*q(p2); ];

%plot(real(qp2),y,'b.-',imag(qp2),y,'r.-')

% expand to physical space
qphysu1=2*real(qu1*exp(i*alpha*x));
qphysv1=2*real(qv1*exp(i*alpha*x));
qphysp1=2*real(qp1*exp(i*alpha*x));

qphysu2=2*real(qu2*exp(i*alpha*x));
qphysv2=2*real(qv2*exp(i*alpha*x));
qphysp2=2*real(qp2*exp(i*alpha*x));

e=2*real(q(eta)*exp(i*alpha*x))+L2;

% combine domain 1 and 2 according to the interface position
Y=repmat(y,1,Nx); E=repmat(e,N,1);
sel1=Y>E; % select gridpoints where y is larger than eta
uu=sel1.*qphysu1+(~sel1).*qphysu2;
vv=sel1.*qphysv1+(~sel1).*qphysv2;
pp=sel1.*qphysp1+(~sel1).*qphysp2;

% show the velocity and pressure field
figure(1); clf
selx=1:5:Nx; sely=1:2:N; 
quiver(x(selx),y(sely),uu(sely,selx),vv(sely,selx),'k'); hold on
surf(x,y,pp-10,'facealpha',0.5); shading interp;
plot(x,real(e),'b');
axis equal;axis([x(1),x(end),1.5,2.5]);
xlabel('x'); ylabel('y'); title('Velocity field and interface of the initial condition');

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','rayleigh_taylor_gerris_field.png');

%{
# Saving to disc
There are several ways to save fields for gerris, here we chose the simplest one, the [cgd format (cartesian grid data)](http://gerris.dalembert.upmc.fr/gfsfunction.html#Cartesian_Grid_Data_.28CGD.29_files)
%}

% translate to have interface position at y=0 for gerris
y=y-L2;
e=e-L2; 

% save to cartesian grid data format for gerris (.cgd)
disp('saving initial condition')
% for u
fid=fopen('u.cgd','w');
fprintf(fid,'%s\n','2 x y');
fprintf(fid,'%u %u\n',Nx,N);
fprintf(fid,'%f ',x);fprintf(fid,'\n');
fprintf(fid,'%f ',y);fprintf(fid,'\n');
fprintf(fid,'%f\n',uu);  
fclose(fid);
% for v
fid=fopen('v.cgd','w');
fprintf(fid,'%s\n','2 x y');
fprintf(fid,'%u %u\n',Nx,N);
fprintf(fid,'%f ',x);fprintf(fid,'\n');
fprintf(fid,'%f ',y);fprintf(fid,'\n');
fprintf(fid,'%f\n',vv);  
fclose(fid);

% for eta
% create a closed contour
xeta=[-100 x 100 100 -100];
yeta=[e(1) e e(end) 100 100];
% write to file
fid=fopen('rayleigh_taylor_shape','w');
for ind=1:length(xeta)
   fprintf(fid,'%f %f \n',[xeta(ind),yeta(ind)]);
end
fclose(fid);
% create gts file for the interface
!shapes rayleigh_taylor_shape > rayleigh_taylor_shape.gts


%{
# Starting Gerris

The parameter file for the Gerris simulation is in [rayleigh_taylor.gfs]().

We can start Gerris directly from here. The command to be executed in the shell is stored in the string *command*. The character "!" at the start of a command tells that the command should be executed in the shell. This command, stored in a string, is then evaluated using the function *eval*.

In Gerris it is possible to define some of the parameter values directly from the command line, like we do here to set the value of the *boxsize*. To set it automatically to the value of *Lx* in the code, we use *sprintf* which will insert in the string the value of *Lx* with two digits after the comma.

The Gerris simulation saves on the disc a file *intevo* which memorizes the time in the first column and the value of the energy of the perturbation as a second column.
%}

% start gerris simulation
disp('starting gerris')
cparam=['-Dboxsize=' sprintf('%3.2f',Lx) ' ' ...
	    '-Drho1=' sprintf('%3.2f',rho1) ' ' ...
        '-Drho2=' sprintf('%3.2f',rho2) ' ' ...
        '-Dsigma=' sprintf('%3.2f',sigma) ' '];
command=['!gerris2D -m ' cparam ' rayleigh_taylor.gfs | gfsview2D rt.gfv'];
disp(command);
eval(command)

%{

![A snapshot from the vorticity of the Gerris simulation](rayleigh_taylor_gerris_snapshot.png)

# Validation

To validate these computations, we compare the evolution of the energy here and in Gerris. The eigenvalue of the mode we use here is stored in *s(1)*. Since the energy is the square of the mode amplitude, and the mode amplitudes grows in time like 
$$
\exp(s_rt)
$$
where $s_r$ is the real part of the eigenvalue $s$, then the energy grows like
$$
\exp(2s_rt)
$$
We scale this energy so that the linear one and the nonlinear one are initialy equal. 
%}

% validation
figure(2);clf;
a=load('intevo');
t=a(:,1); e=a(:,2);
etheo=exp(real(s(1))*t); 
etheo=etheo/etheo(1)*e(1); 

semilogy(t,e,'b.-',t,etheo,'r');
xlabel('time'); ylabel('amplitude of perturbation')
legend('gerris','easystab');
title('Comparison Gerris/easystab')
grid on

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','rayleigh_taylor_validation.png')

%{
![Comparison of easystab and Gerris](rayleigh_taylor_validation.png)

# Exercices/Contributions

* Please validate the code for several values of the viscosity and of the wavelength
* Please do these validation with Gerris for other instabilities, like [poiseuille_uvp.m]() or [couette_uvwp.m](), or [pipe_sym.m]() or even [rayleigh_benard.m]() for which you will need to initiate as well the temperature field. 


%}
