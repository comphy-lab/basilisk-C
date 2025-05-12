%{


*This document belongs to the [easystab](http://basilisk.fr/sandbox/easystab/README) project, please consult the main page of the project for explanations.*

*The present program was specifically designed as a support for the course ["Introductions to hydrodynamical instabilities"](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) for the "DET" Master's cursus in Toulouse. For the underlying theory please see [Lecture notes for chapter 7](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md#lecture-7-shear-flow-instabilities-ii-viscous-instabilities)*

*To use this program: click "raw page source" in the left column, copy-paste in the Matlab/Octave editor, and run.*

# Stability of the plane Poiseuille flow in primitive variables

In this code, we solve the linear eigenvalue problem for stability analysis of the Plane Poiseuille flow. 
This is the flow in a flat channel between two infinite plates, driven by a pressure gradient. Here is a sketch of the flow:

![Sketch of the Poiseuille flow.](poiseuill_uvp_sketch.png)

We perform the analysis using the primitive variables u,v and p the horizontal and vertical components
 of the velocity and the pressure. Thus we use the Navier-Stokes equations, linearized 
 about a parabolic velocity profile, plus as well the continuity equation. 
 Since there is no time derivative in the continuity equation for incompressible flow this 
 system becomes a *descriptor system*, also known as *differential algebraic equation*. 
 This is a nice formulation because there are no fourth order derivative as for the Orr-Sommerfeld equation.
 
The parameters are the wave number $k=2\pi/\lambda$ where $\lambda$ is the wavelength in 
the direction of the flow, the height of the domain $L$ which must be $2$, the Reynolds number 
and the number of grid points.
%}

function [] = Poiseuille_temporal_viscous % main program
global y dy dyy;
close all; 

k=.4;    % the wave number
L=2;        % the height of the domain, from -1 to 1
Re=100000;    % the Reynolds number
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
0 & 0 & 0  
\end{array} 
\right] 
$$

$$
{\mathcal A} = 
\left[
\begin{array}{ccc} 
-i k \bar{U} + Re^{-1} ( \partial_y^2 - k^2) & - \partial_y \bar{U} & - i k \\ 
0 & -i k \bar{U} + Re^{-1} ( \partial_y^2 - k^2) & - \partial_y \\
i k  & \partial_y  & 0 
\end{array} 
\right] 
$$

Here we build the dynamics matrices $A$ and $B$. 
%}


A=[ ...
    S, -diag(Uy), -dx; ...
    Z, S, -dy; ...
    dx, dy, Z];
B=blkdiag(I,I,Z);

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
ylabel('\omega_i');
axis([0,1,-.5,0.05]);grid on;
title({['Temporal spectrum for k = ',num2str(k), ', Re = ',num2str(Re)], 'Click on eigenvalues to see the eigenmodes'});
set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','Spectrum.png');
%{
![Marginal stability curve](Poiseuille_temporal_viscous/Spectrum.png)
%}



if (getenv('OCTAVE_AUTORUN')=='true') 
% to generate a figure for the most amplified mode when running on the server
  figure(2);
  plotmode(0,0,U(:,1),lambda(1),k)
  print('-dpng','-r100','Mode.png');
%{
![Most amplified mode for this set of parameters 
(dashed lines indicate the location of the critical layer defined by $c_r=-\lambda_i/k=U(y_c)$ .](Poiseuille_temporal_viscous/Mode.png)
%}
end  

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
    legend({'Re(u)','Im(u)','Re(v)','Im(v)','Re(p)','Im(p)'})
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
    contourf(x,yy,vorticityvorticity,10); hold on; 
    quiver(x,yy,uu,vv,'k'); hold on;
    %axis equal;
    xlabel('x'); ylabel('y'); title({'Structure of the eigenmode (vorticity)',['lambda = ',num2str(lambda)]});
end% function plotmode


