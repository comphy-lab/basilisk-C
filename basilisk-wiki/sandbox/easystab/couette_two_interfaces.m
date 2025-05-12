%{
# Capillary/advection instability in 2D

Here is an idea of Luca Biancofiore, Eyal Heifetz and François Gallaire: two capillary interfaces inside of a plane Couette flow. Based on the paper of Heifetz, there can be some instability when the left traveliing capillary wave of the top interface has the same velocity as the right-traveling capillary wave of the bottom interface. They interact through the vertical velocity, phase lock, then grow or decay depending on the wave speed (wavenumber) 

We take a very low viscosity because we will compare with a theory for a perfect fluid.
%}

clear all; clf;

% parameters
N1=50; N2=30; N3=50;     % number of gridpoints
L=8;        % bottom is at -L, top at L
rho=1;      % density
mu=0.001;    % fuid viscosity
sigma=4.5;    % surface tension
alpha=0.4;    % wavenumber in x

mu1=mu; mu2=mu; mu3=mu;

%{
We build differentiation matrices for each of the fluid domain. Domain number increase from bottom to top, just like the indices in the coordinate. From $-L$ to $-1$ is domain 1, from $-1$ to $1$ is domain 2, and from $1$ to $L$ is domain 3. The interface displacement of the bottom interface is called $\eta_1$ and the one for the top interface is $\eta_2$.
%}
% differentiation matrices
[d1.y,d1.yy,d1.wy,y1]=dif1D('cheb',-L,L-1,N1); I1=eye(N1); Z1=zeros(N1,N1);
[d2.y,d2.yy,d2.wy,y2]=dif1D('cheb',-1,2,N2); I2=eye(N2); Z2=zeros(N2,N2);
[d3.y,d3.yy,d3.wy,y3]=dif1D('cheb',1,L-1,N3); I3=eye(N3); Z3=zeros(N3,N3);

d1.x=i*alpha*I1; d1.xx=-alpha^2*I1; d1.lap=d1.xx+d1.yy;
d2.x=i*alpha*I2; d2.xx=-alpha^2*I2; d2.lap=d2.xx+d2.yy;
d3.x=i*alpha*I3; d3.xx=-alpha^2*I3; d3.lap=d3.xx+d3.yy;

% location vectors
u1=1:N1; v1=u1+N1; p1=v1+N1;
u2=p1(end)+(1:N2); v2=u2+N2; p2=v2+N2;
u3=p2(end)+(1:N3); v3=u3+N3; p3=v3+N3;
eta1=p3(end)+1;
eta2=eta1+1;

%{
The base flow is a simple Couette flow with $u=y$.
%}

% base flow
U1=y1; U1y=1+0*y1; S1=-diag(U1)*d1.x+mu1*d1.lap;
U2=y2; U2y=1+0*y2; S2=-diag(U2)*d2.x+mu2*d2.lap;
U3=y3; U3y=1+0*y3; S3=-diag(U3)*d3.x+mu3*d3.lap;

%{
The system matrices are just like for any planar parrallel flow, like for instance [poiseuille_uvp.m]().
%}

% System matrices
a1=[S1,-diag(U1y),-d1.x;Z1,S1,-d1.y;d1.x,d1.y,Z1]; e1=blkdiag(rho*I1,rho*I1,Z1);
a2=[S2,-diag(U2y),-d2.x;Z2,S2,-d2.y;d2.x,d2.y,Z2]; e2=blkdiag(rho*I2,rho*I2,Z2);
a3=[S3,-diag(U3y),-d3.x;Z3,S3,-d3.y;d3.x,d3.y,Z3]; e3=blkdiag(rho*I3,rho*I3,Z3);

%{
To understand the connection between the different domains and the connection between the vertical velocity and the interface displacement, please see [rayleigh_taylor.m]().
%}

A=blkdiag(a1,a2,a3,-i*alpha*U1(N1),-i*alpha*U2(N2));
E=blkdiag(e1,e2,e3,1,1);

% we add the connection eta1-v1 and eta2-v2
A(eta1,v1(N1))=1;
A(eta2,v2(N2))=1;

%{
# Boundary conditions

The boundary conditions are like the one for [rayleigh_taylor.m#boundary-conditions](), except that instead of having one capillary interface, we have two of them. Also, here we do not put the gravity. The way it is formulated here in primitive variables, there is no base flow component in the boundary conditions. 
%}

% boundary conditions
dir=[u1(1),v1(1),u3(N3),v3(N3)]; % no-slip at both walls
loc=[u1([1,N1]),v1([1,N1]),u2([1,N2]),v2([1,N2]),u3([1,N3]),v3([1,N3])];
II=eye(3*(N1+N2+N3)+2);

C=[II(dir,:); ... % no-slip at bottom and top walls
   II(v1(N1),:)-II(v2(1),:); ... % continuity of v across eta1 
   II(v2(N2),:)-II(v3(1),:); ... % continuity of v across eta2
   II(u1(N1),:)-II(u2(1),:); ... % continuity of u across eta1
   II(u2(N2),:)-II(u3(1),:); ... % continuity of u across eta2
   II(p1(N1),:)-II(p2(1),:)-sigma*alpha^2*II(eta1,:); ... % Continuity of pressure across eta 1
   II(p2(N2),:)-II(p3(1),:)-sigma*alpha^2*II(eta2,:); ... % Continuity of pressure across eta 2
   mu1*(d1.y(N1,:)*II(u1,:)+d1.x(N1,:)*II(v1,:))-mu2*(d2.y(1,:)*II(u2,:)+d2.x(1,:)*II(v2,:)); ... % continuity of shear stress across eta 1
   mu2*(d2.y(N2,:)*II(u2,:)+d2.x(N2,:)*II(v2,:))-mu3*(d3.y(1,:)*II(u3,:)+d3.x(1,:)*II(v3,:)); ... % continuity of shear stress across eta 2
   ]; 

E(loc,:)=0; A(loc,:)=C;

% compute eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

%{
# Validation

For the theory, we assume an inviscid fluid described by its streamfunction $\phi$. If there was just one interface, there would be two waves, one traveling to the right and one traveling to the left. Just like fro [rayleigh_taylor.m]() we get the result for the eigenvalue $s$
$$
s=\sqrt{\frac{\alpha(\rho_1-\rho_2)g-\alpha^3\sigma}{\frac{\rho_1}{\tanh(\alpha L_1)} + \frac{\rho_2}{\tanh(\alpha L_2)} }}
$$
but here the two lfuids have the sale density, there is no gravity, and we assume that the walls are very far $L=\infty$, so we have
$$
s=\pm i\sqrt{\frac{\alpha^3\sigma}{2\rho}}
$$
An other thing, if the interface is at a height $y_0$ where the base flow velocity is non zero, we have the advection effect coming in the advection equation for the interface
$$
\eta_t+U(y_0)\eta_x=v(y_0)
$$
so the eigenvalue becomes
$$
s+iU(y_0)\alpha=\pm i\sqrt{\frac{\alpha^3\sigma}{2\rho}}
$$
with just an added imaginary part describing the temporal phase shift due to the advection.

So now we have two of these interfaces at the two heights $y=-1$ for $\eta_1$ and $y=1$ for $\eta_2$. The wave can coexist nicely such that their individual structure is unchanged, except for a coupling coming throught the fact that each interface is advected vertically by its own wave, but also by the wave traveling on the other interface. This coupling is two-way. This is the coupling that may induce the instability. here are the advections of each interfaces
$$
\begin{array}{l}
(s+i\alpha)\eta_2=v_1(1)+v_2(1)\\
(s-i\alpha)\eta_1=v_1(-1)+v_2(-1)\\
\end{array}
$$
where we see the vertical velocity from wave 1 acting on interface 1 and 2 and vertical velocity from wave 2 acting on interfaces 1 and 2.

From the continuity of $v$ through the interface and the pressure jump, individually for each wave we have 
$$
\begin{array}{l}
v_1(1)+v_2(1)=i\alpha (A_1+A_2\exp(-2\alpha))\\
v_1(-1)+v_2(-1)=i\alpha (\exp(-2\alpha)A_1+A_2)
\end{array}
$$
where we see that the coupling comes through an amplitude term $\exp(-2\alpha)$ which becomes very small when the two interfaces are far apart compared to the wavelength of the waves. The $A_1$ and $A_2$ are the amplitudes of the waves and are
$$
\begin{array}{l}
A_1=\frac{i\alpha^2\eta_1}{2\rho}\frac{1}{s-i\alpha}\\
A_2=\frac{i\alpha^2\eta_2}{2\rho}\frac{1}{s+i\alpha}
\end{array}
$$
Combining all this gives
$$
\begin{array}{l}
(s+i\alpha)^2\eta_2=-K(E\frac{s+i\alpha}{s-i\alpha}\eta_1+\eta_2)\\
(s-i\alpha)^2\eta_1=-K(\eta_1+E\frac{s-i\alpha}{s+i\alpha}\eta_2)
\end{array}
$$
with $K=\alpha^3\sigma/\rho$ and $E=\exp(-2\alpha)$. 

To write this as an eigenvalue problem, we do a change of variable
$$
\begin{array}{l}
\tilde{\eta_2}=(s-i\alpha)\eta_2\\
\tilde{\eta_1}=(s+i\alpha)\eta_1\\
\end{array}
$$
to get
$$
\begin{array}{l}
(s+i\alpha)^2\tilde{\eta_2}=-K(E\tilde{\eta_1}+\tilde{\eta_2})\\
(s-i\alpha)^2\tilde{\eta_1}=-K(\tilde{\eta_1}+E\tilde{\eta_2})
\end{array}
$$
which is a quadratic eigenvalue problem that we transform into a linear one by the transformation
$$
\begin{array}{l}
e_2=s\tilde{\eta_2}\\
e_1=s\tilde{\eta_1}\\
\end{array}
$$
to get the system
$$
s
\left(
\begin{array}{cccc}
2i\alpha&1&0&0\\
1&0&0&0\\
0&0&-2i\alpha &1\\
0&0&1&0
\end{array}\right)
\left(\begin{array}{c}
\tilde{\eta_2}\\
e_2\\
\tilde{\eta_1}\\
e_1\\
\end{array}\right)
=
\left(\begin{array}{cccc}
-K+\alpha^2&0&-KE&0\\
0&1&0&0\\
-KE&0&-K+\alpha^2&0\\
0&0&0&1\\
\end{array}\right)
\left(\begin{array}{c}
\tilde{\eta_2}\\
e_2\\
\tilde{\eta_1}\\
e_1\\
\end{array}\right)
$$
And this is what we compare to the numerical calculations. We do just an other thing on the plots: we show the eigenvalues just for the decoupled waves.

%}

% theory with 4 coupled waves
subplot(1,2,1);
for a=[linspace(0.1,0.8,300),alpha];
   K=a^3*sigma/2; E=exp(-2*a);
   Ec=[2*i*a,1,0,0; ...
       1,0,0,0; ...
       0,0,-2*i*a,1; ...
       0,0,1,0; ...
       ];
   Ac=[-K+a^2,0,-K*E,0; ...
       0,1,0,0; ...
       -K*E,0,-K+a^2,0; ...
       0,0,0,1; ...
       ];
   [Uc,Sc]=eig(Ac,Ec); sc=diag(Sc);  
   plot(a,max(real(sc)),'r.'); hold on
end
plot(alpha,max(real(s)),'b.');
xlabel('alpha'); ylabel('growth rate');
title('max growth rate as function of alpha');
grid on

% theory with decoupled interface
sd1=i*(alpha+[1,-1]*sqrt(sigma*alpha^3/(2*rho)));
sd2=i*(-alpha+[1,-1]*sqrt(sigma*alpha^3/(2*rho)));
sd=[sd1,sd2];

% compare the spectrum
subplot(1,2,2);
plot(real(s),imag(s),'b.',real(sc),imag(sc),'ro',real(sd),imag(sd),'co');
xlabel('real part'); ylabel('imaginary part');
title('spectrum'); legend('numerical','coupled theory','decoupled theory');
grid on; axis([-0.1,0.1,-2,2]);


set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','couette_two_interfaces.png');


%{

![Results and validation](couette_two_interfaces.png)

# Exercices/Contributions

* 



%}
