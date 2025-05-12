%{
# Stability of the plane Poiseuille flow in the COMPRESSIBLE CASE

This code is adapted from
[http://basilisk.fr/sandbox/easystab/Poiseuille_temporal_viscous.m](Poiseuille_temporal_viscous.m)
from the Easystab project. We consider the stability properties in the COMPRESSIBLE case,
 considering adiabatic perturbations.

The parameters describing the base flow are the Reynolds number $Re = U_m L
/ \nu$ and the Mach number $M = U_m / c_0$ where $U_m$ is the maximum
velocity in the chanel, $L$ its width and $c_0$ the speed of sound.

%}

function [] = Poiseuille_temporal_viscous_compressible_adiabatic % main program
global y dy dyy;
close all; 
set (0, 'defaultaxesfontsize', 14,'defaultLineLineWidth',2); 

k=1;    % the wave number
L=2;        % the height of the domain, from -1 to 1
Re=1e4;    % the Reynolds number
M = 0.1;
n=200;      % the number of grid points



%{
# Differentiation matrices
We need to compute the derivatives in $x$ and also in $y$. 
But in fact as is usually done for stability of parallel flows, 
we can do a Fourier transform in the direction where the system does not change, 
so here the numerical differentiation is done in $y$, and the differentiation in $x$ 
simply ammounts to multiplication with $ik$.
%}

[dy,dyy,wy,y] = dif1D('cheb',-1,2,n); 
Z=zeros(n,n); I=eye(n); 
dx=i*k*I; dxx=-k^2*I;
Delta=dxx+dyy;
INT=([diff(y)',0]+[0,diff(y)'])/2; 

% base flow
U=1-y.^2; Uy=-2*y;
S=-diag(U)*dx+Delta/Re;

%{
# Construction of the matrices

The problem is written as follows:

$$
\lambda  {\mathcal B} \, \hat{q} = {\mathcal A} \, \hat{q}
$$

with 
$$
{\mathcal B} = 
\left[
\begin{array}{ccc} 
1 & 0 & 0 \\ 
0 & 1 & 0 \\
0 & 0 & M^2  
\end{array} 
\right] 
$$

$$
{\mathcal A} = 
\left[
\begin{array}{ccc} 
-i k \bar{U} + Re^{-1} ( \partial_y^2 - k^2) & - \partial_y \bar{U} & - i k \\ 
0 & -i k \bar{U} + Re^{-1} ( \partial_y^2 - k^2) & - \partial_y \\
-i k  & -\partial_y  & -i k M^2 \bar{U} 
\end{array} 
\right] 
$$

Here we build the matrices $A$ and $B$. 
%}


A=[ ...
    S, -diag(Uy), -dx; ...
    Z, S, -dy; ...
    -dx, -dy, -M^2*diag(U)*dx ];
B=blkdiag(I,I,M^2*I);

%{
# Locations on the grid
For easily imposing the boundary condition, we first define the locations on the grid of the u velocity as well as for v and p. This step is not very important here, whereas it is important with more variables, and especially usefull in 2D.
%}



%{
# Boundary conditions
%}
% We have to modify the lines corresponding to first and last point for u and v.
% We first determine the corresponding indices:
loc=[1 n n+1 2*n ];  
% we then change the lines
III=eye(3*n);
C=III(loc,:);
A(loc,:)=C;
B(loc,:)=0;  

%{
# Computing the eigenmodes
We can use the *eig* function for DAE just as for dynamic systems. the difference is that we will get infinite eigenvalues. Thinking of the incompressible system as the limit of the compressible one, you realize that these infinite eigenvalues are the limit of the sound modes. We remove them, we sort the eigenvalues and eigenvectors.
%}

% computing eigenmodes 
[U,S]=eig(A,B);
% sort and remove spurious eigenvalues
lambda=diag(S);  [t,o]=sort(-real(lambda)); lambda=lambda(o); U=U(:,o);
rem=abs(lambda)>1000; lambda(rem)=[]; U(:,rem)=[];



figure(1);
for ind=1:length(lambda)
  h=plot(-imag(lambda(ind))/k,real(lambda(ind)),'ro'); hold on
  set(h,'buttondownfcn',{@plotmode,U(:,ind),lambda(ind),k});
end
xlabel('c_r');
ylabel('\lambda_r');
axis([-20,20,-.1,0.1]);
title({['Temporal spectrum for k = ',num2str(k), ', Re = ',num2str(Re)], 'Click on eigenvalues to see the eigenmodes'});
set(gcf,'paperpositionmode','auto');
print('-dpng','Spectrum.png');

% second figure with zoom in the range c_r = [0 - 1]
figure(4);
subplot(2,1,1); plot(-imag(lambda)/k,real(lambda),'ro'); hold on
plot([-20 20],[0 0],':k')
xlabel('c_r');ylabel('\lambda_r');axis([-20,20,-.5,0.02]);
subplot(2,1,2); plot(-imag(lambda)/k,real(lambda),'ro'); hold on
plot([-20 20],[0 0],':k')
xlabel('c_r');ylabel('\lambda_r');axis([0,1,-.5,0.1]);
print('-dpng','-r80','SpectrumPoisComp_withZoom.png');

%{
![Spectrum (two representations using different scales](Poiseuille_temporal_viscous_compressible_adiabatic/SpectrumPoisComp_withZoom.png)
%}



if ~isempty(getenv('OCTAVE_AUTORUN')) 
% to generate a figure for the most amplified mode when running on the server
  figure(2);
  plotmode(0,0,U(:,1),lambda(1),k)
  print('-dpng','-r80','Mode.png');
end  
%{
![Most amplified mode for this set of parameters 
(dashed lines indicate the location of the critical layer defined by $c_r=-\lambda_i/k=U(y_c)$ .](Poiseuille_temporal_viscous_compressible_adiabatic/Mode.png)
%}


end

function [] = plotmode(~,~,mode,lambda,k) % function to plot one mode
    global y dy dyy
    figure(2);hold off;
    N = length(y);
    u = mode(1:N);
    v = mode(N+1:2*N);
    p = mode(2*N+1:3*N);
    vorticity = (dy*u)-1i*k*v;
    subplot(1,3,1);hold off;
    plot(real(u),y,'b-',imag(u),y,'b--');hold on;
    plot(real(v),y,'g-',imag(v),y,'g--');hold on;
    plot(real(p),y,'k-',imag(p),y,'k--');hold on;
    ylabel('y'); 
    %legend({'Re(u)','Im(u)','Re(v)','Im(v)','Re(p)','Im(p)'})
    if (-imag(lambda)/k)<1&&(-imag(lambda)/k)>0
      ycrit = sqrt(1+imag(lambda)/k);
      plot([-1,1],[ycrit,ycrit],'k:',[-1,1],[-ycrit,-ycrit],'k:')
    end
    % plot 2D reconstruction
    Lx=2*pi/k; Nx =30;  x=linspace(-Lx/2,Lx/2,Nx);
    yy = y;
    pp = 2*real(p*exp(1i*k*x));
    uu=2*real(u*exp(1i*k*x));
    vv=2*real(v*exp(1i*k*x));
    vorticityvorticity=2*real(vorticity*exp(1i*k*x));
    subplot(1,3,2:3); hold off;
    contourf(x,yy,pp,10); hold on; yskip = round(N/Nx);
    quiver(x,yy(1:yskip:end),uu(1:yskip:end,:),vv(1:yskip:end,:),'k','linewidth',1); hold on;
    %axis equal;
    xlabel('x'); ylabel('y'); title({'Structure of the eigenmode (pressure)',['lambda = ',num2str(lambda)]});
end% function plotmode


