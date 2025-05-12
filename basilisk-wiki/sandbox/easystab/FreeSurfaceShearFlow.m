%{
# Stability of a plane shear flow with a free surface

This programs investigates the stability of a plane shear flow with profile
$\bar{U}(y)$ in a domain $y \in [O,H]$. The bottom ($y=0$) is a wall while the upper surface ($y=H$) 
is a *free surface* characterized by surface
tension $\gamma$. The domain above is a gaz with
negligible density and viscosity and uniform pressure $P_0$.
The gravity field is $g$, in the downward (y<0) direction.
 
This code is adapted from two programs in the easystab project :  
[sandbox/easystab/TS_PlanePoiseuille.m]() for a shear flow without free surface
[sandbox/easystab/free_surface_gravity.m]() for a free surface without shear flow

** NOTE : this program remains to be validated !**


## Equations

### Flow decomposition 

Noting $y = \eta(x,t)$ the equation of the free surface, the flow is thus expanded as follows :

$$ 
\underbrace{\left[
\begin{array}{c} u \\ v \\  p \\ \eta
\end{array} 
\right]}_{q} 
\,
= 
\, 
\underbrace{\left[
\begin{array}{c} \bar{U}(y) \\ 0 \\ \bar{P}(y) \\ H  
\end{array} 
\right]}_{q_0}
\quad + \quad  \epsilon 
\underbrace{{\left[
\begin{array}{c} \hat{u}(y) \\ \hat{v}(y)\\ \hat{p}(y) \\ \hat{\eta} 
\end{array} 
\right]}}_{\hat{q}} e^{i k x - i \omega t}
$$


### Base flow

If neglecting the effect of viscosity on the baseflow, the law $\bar{U}(y)$
can be arbitrary, and the associated pressure distribution is hydrostatic :
$\bar{P}(y) = P_0 + \rho g (H-y)$.

(NB if viscosity is retained and we want the base flow to be a *steady* solution of NS, 
the only solution is the parabolic profile
$U(y) = y/H( 2- y/H)$ and an axial component of $g$ must be present, i.e.
the plate must be tilted by some angle $\alpha$).

### Linear equations : bulk equations

Within the bulk ($0<y<H)$, linearisation of the Navier-stokes equations
leads to:

$$
\begin{array}{rcl}
-i \omega \hat{u} &=&  -i k \bar{U} \hat{u} - \bar{U}' \hat{v} - i k \hat{p}
+\mu (-k^2 \hat{u}+\partial_y^2 \hat{u})
\\
- i \omega  \hat{v}&=&-i k \bar{U} \hat{v} -\partial_y \hat{p}
+\mu(-k^2 \hat{v}+\partial_y^2 \hat{v}) \\
0 &=& i k \hat{u}+\partial_y \hat{v} \\
\end{array}
$$

### Boundary conditions.

#### Bottom 
The boundary conditions at the bottom are $\hat{u}(0) = \hat{v}(0) = 0$.

#### Kinematic
The kinematic boundary condition at free surface comes from linearization
of $[d /dt ]\eta(x,t) = v(x,\eta(x,t),t)$. After linearization:

$$
i (k  \bar{U}(H) - \omega)  \hat{\eta} = \hat{v}(H)
$$

#### Dynamic

The dynamic condition is more complicate to express. In starting variables
we have to write 

$$ (p  - P_0 ) \vec{n} - 2 \mu {\bf D} \cdot \vec{n} = \gamma K \vec{n}
\quad \mathrm{ at } \quad y = \eta(x,t)
$$ 

where $\bf D$ is the rate-of-strain tensor,  $\vec n$
the normal vector at the interface, and  $K$ the curvature.

All these terms have to be developed in powers of $\epsilon$ :
$$
\vec{n} = \frac{1}{\sqrt{1+(\partial_x \eta)^2}} {\left[
\begin{array}{c} -\partial_x \eta \\ 1
\end{array} 
\right]} 
$$



$$ \equiv
\underbrace{\left[
\begin{array}{c} 0 \\ 1
\end{array}
\right]}_{\vec{n}_0} + \epsilon 
\underbrace{\left[
\begin{array}{c} i k \hat{\eta} \\ 0 
\end{array}
\right]}_{\vec{n}_1}
e^{i k x - i \omega t}
+ {\mathcal O}(\epsilon ^2)
$$



$$
{\bf D} = 
\underbrace{\left[
\begin{array}{cc} 
0 & \bar{U}_y \\
\bar{U}_y & 0 
\end{array}
\right]}
_{{\bf D}_0} + \epsilon 
\underbrace{\left[
\begin{array}{cc} 
2 i k \hat{u} & \partial_y \hat{u} + i k \hat{v} \\
 \partial_y \hat{u} + i k \hat{v} & 2 \partial_y \hat{v}
\end{array}
\right]}
_{{\bf D}_1}
e^{i k x - i \omega t}
+ {\mathcal O}(\epsilon ^2)
$$

$$
{\left[{\bf D} \cdot \vec{n}\right]}_{y=\eta}  
= 
{\left[{\bf D}_0 \cdot \vec{n}_0\right]}_{y=H}
+ \epsilon \left\{ {\bf D}_0 \cdot \vec{n}_1 + {\bf D}_1 \cdot \vec{n}_0 +
\partial_y  \left[ {\bf D}_0 \cdot \vec{n}_0 \right] \hat{\eta}   \right\}
e^{i k x - i \omega t}
+ {\mathcal O}(\epsilon ^2)
$$


$$ 
{[(p  - P_0 ) \vec{n}]}_{y=\eta}  = \epsilon \left[ \hat{p}(H) - \rho g \hat{\eta}  \right] \vec{n}_0 
e^{i k x - i \omega t}
+ {\mathcal O}(\epsilon ^2)
$$

$$ 
K = - \epsilon k^2 \hat{\eta} e^{i k x - i \omega t} +  {\mathcal O}(\epsilon ^2)
$$

At the end, introducing everything in the dynamic BC and taking the $y$ and
$x$ components leads to



$$
\left[ \hat{p} - 2 \mu \partial_y \hat{v} \right]_{y=H} + i k \bar{U}'(H) \hat{\eta} - \rho g \hat{\eta}
- \gamma k^2 \hat{\eta} = 0
$$

$$
\left[ \partial_y \hat{u} + i k \hat{v}\right]_{y=H} +  \bar{U}''(H)
\hat{\eta}
= 0
$$


### Matricial formulation

At the end, we can write the equations in the form of a generalised
eigenvalue problem:
$$
- i \omega B \hat{q} = A \hat{q}
$$

The building of the matrices is done in function at the end of this program

%}

%{
## Parameters, definition of geometry and differential operators
%}

 
function [] = main()
close all; 
global y dy dyy Z I INT n H nu gamma g U Uy Uyy

% Physical parameters
H=1;        % the height of the domain, from -1 to 1
nu = 1e-3;  % viscosity
gamma = 0.; % surface tension
g = 1; % gravity
alpha = 1; % wavenumber (k)
Re = 1/nu;
Umax = 1;

% numerical parameters
n=100;      % the number of grid points

%{
## Definition of geometry, base flow and differential operators
%}

% differentiation and integration operators
[dy,dyy,wy,y] = dif1D('cheb',0,H,n); 
Z=zeros(n,n); I=eye(n); 
INT=([diff(y)',0]+[0,diff(y)'])/2; 

U = Umax*y/H.*(2-y/H); % velocity profile
Uy = Umax*2/H.*(1-y/H); % first derivative
Uyy = - Umax*2/H^2*ones(1,n); % second derivative


%{

## Resolution of the eigenvalue problem

%}

[s,UU] = EV_FreeSurf(alpha); 
omega = 1i*s;

%{
### Plotting the spectrum
%}

figure(1);
for ind=1:length(s)
  h=plot(real(omega(ind)),imag(omega(ind)),'*'); hold on
  set(h,'buttondownfcn',{@plotmode,UU(:,ind),omega(ind),alpha});
end
xlabel('\omega_r');
ylabel('\omega_i');
ylim([-2 ,max(imag(omega)) ] );
title({['Temporal spectrum for k = ',num2str(alpha), ', Re = ',num2str(Re)], 'Click on eigenvalues to see the eigenmodes'});
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','spectrum.png');

%{
![** Figure : Temporal spectrum of the free surface shear flow](FreeSurfaceShearFlow/spectrum.png)
%} 


plotmode([],[],UU(:,1),omega(1),alpha);

pause(0.1);
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','mode.png');
%{
![** Figure : Eigenmmode of the free surface shear flow for k=1](FreeSurfaceShearFlow/mode.png)
%} 

%{

## Loop over alpha

%}

figure(5);numev=5;%number of ev to plot
for alpha = .1:.1:2
    [s,UU] = EV_FreeSurf(alpha); 
    alpha
    s(1)
    omega = 1i*s;
    subplot(2,1,1);
    plot(alpha*ones(1,numev),real(omega(1:numev)),'*r');hold on;
    xlabel('k');ylabel('\omega_r');
    subplot(2,1,2);
    plot(alpha*ones(1,numev),imag(omega(1:numev)),'*r');hold on;
    xlabel('k');ylabel('\omega_i');
end
   

pause(0.1);
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','branches.png');
%{
![** Figure : Temporal branches of the free surface shear flow](FreeSurfaceShearFlow/branches.png)
%} 

end

%{
## Function performing the eigenvalue computation
%}

function [s,UU] = EV_FreeSurf(alpha)
global y dy dyy Z I INT n H nu gamma g U Uy Uyy

dx=1i*alpha*I; dxx=-alpha^2*I;
Delta=dxx+dyy;

% base flow
S=-diag(U)*dx+nu*Delta;

% the matrices
A=[S, -diag(Uy), -dx, Z(:,1); ...
   Z, S, -dy, Z(:,1); ...
   dx, dy, Z, Z(:,1); ...
   Z(1,:),I(n,:),Z(1,:),-U(n)*1i*alpha];

B = blkdiag(I,I,Z,1);

% Note that the last line (number 3n+1) is used to enforce the kinematic
% boundary condition


% boundary conditions
% lines 1 and n+1 (corresponding to point y=0 for u and v) are used to
% enforce no-slip boundary conditions at the bottom.
% lines n and 2n (corresponding to point y=H for u and v) are used to
% enforce the two dynamic conditions.


loc=[1 n n+1 2*n];  
B(loc,:)=0;  
A(1,:) = [I(1,:),Z(1,:),Z(1,:),0];
A(n+1,:) = [Z(1,:),I(1,:),Z(1,:),0];
A(n,:)   = [Z(n,:), 2*nu*dy(n,:), -I(n,:), g+gamma*alpha^2-1i*alpha*Uy(n)];
A(2*n,:) = [dy(n,:), 1i*alpha*I(n,:), Z(n,:),   Uyy(n)]; 

% computing eigenmodes 
[UU,S]=eig(A,B);

% sorting by ascending real part and removing large-eigenvalue spurious modes
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); UU=UU(:,o);
rem=abs(s)>100; s(rem)=[]; UU(:,rem)=[];

if(nargout==1) 
    s=s(1); 
    % this is a trick to allow to use the function as s=EV_FreeSurface(alpha). 
    % This allows to pass the function as a handle for fzero.
end

end % function EV_FreeSurf


function [] = plotmode(~,~,mode,omega,alpha)
    global y dy dyy
    figure(2);hold off;
    N = length(y);
    u = mode(1:N);
    v = mode(N+1:2*N);
    p = mode(2*N+1:3*N);
    eta = mode(3*N+1);
    vorticity = (dy*u)-1i*alpha*v;
    subplot(1,3,1);hold off;
    plot(real(u),y,'b-',imag(u),y,'b--');hold on;
    plot(real(v),y,'g-',imag(v),y,'g--');hold on;
    plot(real(p),y,'k-',imag(p),y,'k--');hold on;
    ylabel('y'); ylim([0,1.5]);
    legend({'Re(u)','Im(u)','Re(v)','Im(v)','Re(p)','Im(p)'},'Interpreter','latex')
    title('Structure of the eigenmode');
    % plot 2D reconstruction
    Lx=2*pi/alpha; Nx =30;  x=linspace(-Lx/2,Lx/2,Nx);
    yy = y;
    pp = 2*real(p*exp(1i*alpha*x));
    uu=2*real(u*exp(1i*alpha*x));
    vv=2*real(v*exp(1i*alpha*x));
    vorticityvorticity=2*real(vorticity*exp(1i*alpha*x));
    subplot(1,3,2:3); hold off;
    contourf(x,yy,vorticityvorticity,10); hold on; 
    quiver(x,yy,uu,vv,'k'); hold on;
    plot(x,y(N)+0.1*real(eta*exp(1i*alpha*x)),'r-','linewidth',3)
    ylim([0:1.5]);
    xlabel('x'); ylabel('y'); title({'Structure of the eigenmode (2D reconstruction)',['for k = ',num2str(alpha) , ' ; omega = ',num2str(omega)]});
end% function plotmode
