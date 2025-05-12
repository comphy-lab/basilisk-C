%{
# Gravity free-surface waves with two fluids

This is basically the same thing than [free_surface_gravity.m]() but here there are two fluids on top of each other, so this is a little more complicated, but we do not have to neglect the density of the top fluid so this is more general. Technically it works like that, we build two dynamic systems for the two fluid layers with different physical properties, then we need to connect this two fluid layers through the moving interface.

We take a very low viscosity because we will compare with a theory for a perfect fluid.
%}

clear all; clf;

% for domain 1 (top)
N1=50;      % number of gridpoints
L1=1;        % Fluid height in y
rho1=0.5;      % fluid density
mu1=0.00001;    % fuid viscosity

% for domain 2 (bottom)
N2=50;      % number of gridpoints
L2=1;        % Fluid height in y
rho2=1;      % fluid density
mu2=0.00001;    % fuid viscosity

alpha=2;    % wavenumber in x
g=1;        % gravity

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
  (II(p1(1),:)-rho1*g*II(eta,:))-(II(p2(end),:)-rho2*g*II(eta,:)); ... % Continuity of pressure
   mu1*(dy1(1,:)*II(u1,:)+dx1(1,:)*II(v1,:))-mu2*(dy2(end,:)*II(u2,:)+dx2(end,:)*II(v2,:))]; % continuity of shear stress

E(loc,:)=0;  
A(loc,:)=C;

% compute eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

%{
# validation

For a theory we will de the same derivation as for [free_surface_gravity.m#validation](), with the difference that we have an other fluid at the top. We call $\phi_1$ and $\phi_2$ the stream functions in domain 1and 2. They satisfy the Lalace equations
$$
\Delta \phi_1=0, \quad \Delta \phi_2=0
$$
thus we have
$$
\hat{\phi}_1=A_1 \cosh (\alpha (y-L_1))+B_1 \sinh (\alpha (y-L_1))
$$
and the no-penetration at the top wall at height $y=L_1$ gives $A_1=0$ thus
$$
\hat{\phi}_1=B_1 \sinh (\alpha (y-L_1))
$$
we have chosen to parameterize the hyperbolic functions with $y-L_1$ just to simplify the imposition of the boundary conditions. For domain 2 we have
$$
\hat{\phi}_2=A_2 \cosh (\alpha (y+L_2))+B_2 \sinh (\alpha (y+L_2))
$$
and the no-penetration at the bottom wall at height $y=-L_2$ gives $A_2=0$ thus
$$
\hat{\phi}_2=B_2 \sinh (\alpha (y+L_2))
$$
The advection of the interface gives 
$$
s\hat{\eta}=v_1(0)=\phi_{1x}(0)=i\alpha \phi_1(0)=I\alpha A_1\sinh(-\alpha L_1)
$$

We now use the continuity of $v$ across the interface
$$
\begin{array}{l}
v_1(\eta)=v_2(\eta)\\
v_1(0)=v_2(0)\\
i\alpha \hat{\phi}_1(0)=i\alpha \hat{\phi}_2(0)\\
A_1\sinh(-\alpha L_1)=A_2\sinh(\alpha L_2)\\
A_2=A_1\sinh(-\alpha L_1)/\sinh(\alpha L_2)\\
\end{array}
$$


The continuity of the pressure at the interface gives 
$$
\begin{array}{l}
p_1(\eta)=p_2(\eta) \\
p_1(0)+P_1(0)+\eta P_{1y}(0)=p_2(0)+P_2(0)+\eta P_{2y}(0) \\
p_1(0)-p_2(0)=(\rho_1-\rho_2)g\eta\\
\end{array}
$$
where we have used the fact that the base pressure $P$ is continuous accross the interface, the base pressure gradient is  $P_y=-\rho g$ with $\rho$ different in the two fluids, and that $\eta p_y$ is too small to be kept.

we need now to express the pressure as a function of $\phi$, for this we use the $u$ component of the Euler equation at $y=0$
$$
\begin{array}{l}
\rho u_t=-p_x\\
\rho s\hat{u}(0)=-i\alpha \hat{p}(0)\\
\end{array} 
$$
thus $\hat{p}(0)=\rho s\hat{\phi}_y/i\alpha$. We can thus express the continuity of the pressure using the streamfunctions
$$
i(\hat{p}_1-\hat{p}_2)=\rho_1A_1\cosh(-\alpha L_1)-\rho_2 A_2 \cosh(\alpha L_2)
$$
we can now use the relation between $A_1$ and $A_2$ to get
$$
sA_1 \sinh(\alpha L_1)\left[ \frac{\rho_1}{\tanh(\alpha L_1)}+\frac{\rho_2}{\tanh(\alpha L_2)}\right]=i(\rho_1-\rho_2)g \eta
$$
if we now use the relation obatined from the advection of $\eta$ we obtain
$$
s^2=\frac{\alpha(\rho_1-\rho_2)g}{\frac{\rho_1}{\tanh(\alpha L_1)} + \frac{\rho_2}{\tanh(\alpha L_2)}}
$$
and we have used the fact that $\cosh$ is even and $\sinh$ is uneven. This is the final relation telling the behaviour of the eigenmodes of the flow.

We see that the flow is unstable if the top fluid is heavier than the bottom fluid (Rayleigh-Taylor instability), since then there are two solutions for $s$, one positive and one negative. The positive one correspond to the unstable eigenmode. 

If now $\rho_1$ tends to zero, we recover the solution from [free_surface_gravity.m#validation]()
$$
s^2=\alpha \tanh(\alpha L_2) g
$$

For the validation we will compare the wave velocity from teh code to this. The wave speed is minus the imaginary part of the eigenvalue divided by the wavenumber $\alpha$. The theory for this is
$$
c=\sqrt{\frac{-(\rho_1-\rho_2)g}{\alpha \left[\frac{\rho_1}{\tanh(\alpha L_1)} + \frac{\rho_2}{\tanh(\alpha L_2)}\right] }}
$$
%}

% validation
alphavec=linspace(0,5,100);
ctheo=sqrt(-(rho1-rho2)*g./(alphavec.*(rho1./tanh(alphavec*L1)+rho2./tanh(alphavec*L2))));

cnum=max(imag(s(1:50)))/alpha;
plot(alphavec,ctheo,'r-',alpha,cnum,'b.');
xlabel('alpha');ylabel('wave velocity'); 
title('validation')
grid on

set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','free_surface_gravity_twofluids.png');


%{

![The wave velocity as a function of the wavenumber $\alpha$](free_surface_gravity_twofluids.png)

# Exercices/Contributions

* Please check that the codes behaves well in the limit of low density of the top fluid (thus its results tend to that of [free_surface_gravity.m]())
* Please check what happens when the top fluid is heavier than the bottom fluid (Rayleigh_taylor instability)-------------> [rayleigh_taylor.m]()
* Please add the surface tension at the interface (the pressure is nolonger continuous at the interface, it is equal to the Laplace pressure jump $\Delta p=\sigma/R$ where $\sigma$ is the intensty of the surface tension and $R$ is the radius of curvature of the interface)-------------> [rayleigh_taylor.m]()
* Please compare the results of this code to a simulation with Gerris by using the first eigenmode as an initial condition in Gerris.



%}
