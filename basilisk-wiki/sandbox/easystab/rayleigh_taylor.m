%{
# Free-surface waves with two fluids: gravity and surface tension

This is basically the same thing than [free_surface_gravity.m]() but here there are two fluids on top of each other, so this is a little more complicated, but we do not have to neglect the density of the top fluid so this is more general. Technically it works like that, we build two dynamic systems for the two fluid layers with different physical properties, then we need to connect this two fluid layers through the moving interface.

We take a very low viscosity because we will compare with a theory for a perfect fluid.
%}

clear all; clf;

alpha=1.4;    % wavenumber in x
g=1;        % gravity
sigma=0.1;    % surface tension

% for domain 1 (top)
N1=50;      % number of gridpoints
L1=1;        % Fluid height in y
rho1=1;      % fluid density
mu1=0.00001;    % fuid viscosity

% for domain 2 (bottom)
N2=50;      % number of gridpoints
L2=1;        % Fluid height in y
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
#Boundary conditions

The boundary conditions are just as like in [free_surface_gravity_twofluids.m#boundary-conditions](), with just one single difference: the pressure is nolonger continuous across the interface, it has a pressure jump due to the combination of interface curvature and surface tension (as described in [free_surface_2D.m]()). Thus the pressure constraint across the interface becomes
$$
[p_1(0)-\rho_1 g\eta]-[p_2(0)-\rho_2 g\eta]=-\sigma\alpha^2\eta 
$$
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
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

%{
# validation

Now we are in an unstable configuration with the top fluid heavier than the bottom fluid, so the least stable eigenmode should be purely real (no oscillations), and positive (amplitude growing in time). We can do a model for the growth rate just like we did for [free_surface_gravity.m#validation]() for just gravity and [free_surface_2D.m#validation]() for just surface tension (but only one fluid). The derivation of the model is very similar and we get the exponential growth rate:
$$
s=\sqrt{\frac{\alpha(\rho_1-\rho_2)g-\alpha^3\sigma}{\frac{\rho_1}{\tanh(\alpha L_1)} + \frac{\rho_2}{\tanh(\alpha L_2)} }}
$$

%}

% validation
alphavec=linspace(0,3,100);
stheo=sqrt((alphavec.*(rho1-rho2)*g-alphavec.^3*sigma)./(rho1./tanh(alphavec*L1)+rho2./tanh(alphavec*L2)));

snum=real(s(1));
plot(alphavec,real(stheo),'r-',alpha,snum,'b.');
xlabel('alpha');ylabel('exponential growth rate'); 
title('validation')
grid on
break
set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','rayleigh_taylor.png');


%{
![Growth rate of the Rayleigh-Taylor instability](rayleigh_taylor.png)

# Exercices/Contributions

* Please investigate what happens for $\alpha$ large enough that the growth rate becomes zero
* Please compare this result with high viscosity and thin layer for the top fluid with the lubrication limit
* Please use the unstable eigenmode as an initial condition in Gerris and compare the nonlinear evolution and the linear one ------------> [rayleigh_taylor_gerris.m]()


%}