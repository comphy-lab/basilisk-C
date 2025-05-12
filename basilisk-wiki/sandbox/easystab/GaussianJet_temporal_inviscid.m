%{

*This document belongs to the [easystab](http://basilisk.fr/sandbox/easystab/README) project, check the main page of the project to understand the general philosophy of the project.*

*The present program was specifically designed as a support for the course ["Introductions to hydrodynamical instabilities"](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) for the "DET" Master's cursus in Toulouse. For the underlying theory please see [Lecture notes for chapter 6](http://basilisk.fr/sandbox/easystab/LectureNotes_Inviscid.md)*

*To use this program: click "raw page source" in the left column, copy-paste in the Matlab/Octave editor, and run.*

# Kelvin-Helmholtz instability of a shear layer 

(temporal analysis, viscous case, discretization based on Rayleigh Equation)

We study the stability of a parallel shear flow U(y) in the inviscid case.
We solve the Rayleigh equation, written as follows :

$$(\bar{U} - c) (\partial_y^2 - k^2) \hat \psi - \partial_y^2 \bar{U} \hat \psi = 0$$

Version modified on jan. 9 2020, D.F

%}

%function [] = main()
close all;
global y dy dyy w Z I U Uy Uyy
loopk = 1; % set to 0 to skipp the loop over k 
%{

### Derivation matrices

Here we wish to solve a problem in an infinite domain. 
We use Chebyshev discretization with stretching. See [differential_equation_infinitedomain.m]() to see how this works.

%}

% Numerical parameters
N=201;      % the number of grid points
discretization = 'chebInfAlg'; 
[dy,dyy,w,y] = dif1D(discretization,0,3,N,0.9999);
Z=zeros(N,N); I=eye(N); 

%{

### Base flow

The base flow is $U(y) = exp(-y^2)$.

We also need to compute its first and second derivatives:

%}

U=exp(-y.^2); 
Uy=-2*y.*exp(-y.^2);
Uyy = (-2+4*y.^2).*exp(-y.^2);

%{

### Eigenvalue computation

We compute the eigenvalues/eigenmodes with the function [KH_inviscid](#Function_KH_Inviscid), 
defined at the end of this program.
%}

%%
% Physical parameters
alpha=1.;    % the wave number

[c,UU] = KH_inviscid(alpha,N);
omega = c*alpha;

%{ 
### Plot the spectrum :
%}

figure(1);hold off;
for ind=1:length(omega)
  %%%% plot one eigenvalue
  h=plot(real(omega(ind)),imag(omega(ind)),'*'); hold on
  %%%%  assign the corresponding event on click 
  set(h,'buttondownfcn',{@plotmode,UU(:,ind),omega(ind),alpha});
end
xlabel('\omega_r');
ylabel('\omega_i');
title({['Temporal spectrum for k = ',num2str(alpha)], 'Click on eigenvalues to see the eigenmodes'});
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','KH_temporal_inviscid_spectrum.png');
%{
![** Figure : Temporal spectrum of a tanh shear layer](KH_temporal_inviscid/KH_temporal_inviscid_spectrum.png)
%} 

%{ 
### Plot the unstable eigenmode :
%}

plotmode([],[],UU(:,1),omega(1),alpha);

pause(0.1);
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','KH_temporal_inviscid_mode.png');
%{
![** Figure : Unstable eigenmmode of the tanh shear layer for k=0.5](KH_temporal_inviscid/KH_temporal_inviscid_mode.png)
%} 


%{

## Loop over k to draw the amplification rate as function of the wavenumber

%}
%%
if loopk
    alphatab = 0:.01:2.5;
    lambdatab = [];
    for alpha=alphatab
        [s,~] = KH_inviscid(alpha,N);
        lambdatab = [lambdatab alpha*s(1)];
        figure(4);
        plot(alphatab(1:length(lambdatab)),imag(lambdatab),'b');hold on;
        pause(0.01);
    end
    ylabel('\omega_i');
    xlabel('k');
    title('Growth rate as function of k');
    set(gcf,'paperpositionmode','auto');
    print('-dpng','-r80','KH_temporal_inviscid.png');

end%if loopk

%end%function main
%{

![** Figure : results for the Kelvin-Helmholtz instability (temporal) of a tanh shear layer](KH_temporal_inviscid/KH_temporal_inviscid.png)


## Function KH_inviscid
 
Here is the function to compute the eigenvalues/eigenmodes.

%}
%%
function [s,UU] = KH_inviscid(alpha,N)

global y dy dyy I U Uyy

%{ 

The Rayleigh equation is written under the form

$$
 c \left[ \partial_y^2 - k^2 \right] \hat \psi = 
\left[\bar{U} (\partial_y^2 - k^2) - \partial_y^2 \bar{U} \right] \hat \psi
$$

which is a generalized eigenvalue problem. After discretizing the operators, is can be written under the matricial form:

$$c B X  = A X$$ 

where $X$ is the discretized version of $\hat \psi$.


for more details on the theory please see [Lecture notes for chapter 6](http://basilisk.fr/sandbox/easystab/LectureNotes_Inviscid.md)
%} 

%differential operators
dx=1i*alpha*I; dxx=-alpha^2*I;
Delta=dxx+dyy;

% the matrices
B=Delta;
A=diag(U)*Delta-diag(Uyy);

%Boundary conditions : we use Dirichlet at the boundaries of the stretched domain
indBC=[1,N];
A(indBC,:)=I(indBC,:);
B(indBC,:)=0;  

% compute the eigenmodes 
[UU,S]=eig(A,B);

% sort the eigenvalues by increasing imaginary part
s=diag(S);  [~,o]=sort(-imag(s)); s=s(o); UU=UU(:,o);
rem=abs(s)>1000; s(rem)=[]; UU(:,rem)=[];

end%function KH_inviscid

%{

## Function plotmode
 
This function is used for plotting an eigemode when clicking on the corresponding eigenvalue in the "spectrum" plot.

%}


function [] = plotmode(~,~,psiM,omega,alpha)
    global y dy dyy
    Yrange = 5;
    figure(2);hold off;
    % plot psi(y)
    u=dy*psiM;
    v=-1i*alpha*psiM;
    vorticity = -alpha^2*psiM-dyy*psiM;
    subplot(1,3,1);hold off;
    plot(real(psiM),y,'r-',imag(psiM),y,'r--');hold on;
    plot(real(u),y,'b-',imag(u),y,'b--');hold on;
    plot(real(v),y,'g-',imag(v),y,'g--');hold on;
    ylabel('y'); ylim([-Yrange,Yrange]);
%legend({'$Re(\hat \psi)$','$Im(\hat \psi)$','$Re(\hat u)$','$Im(\hat u)$','$Re(\hat v)$','$Im(\hat v)$'},'Interpreter','latex');
    legend({'Re(\psi)','Im(\psi)','Re(u)','Im(u)','Re(v)','Im(v)'});
    title('Structure of the eigenmode');
    % plot 2D reconstruction
    Lx=2*pi/alpha; Nx =30;  x=linspace(-Lx/2,Lx/2,Nx);
    psiM(abs(y)>Yrange,:)=[];
    u(abs(y)>Yrange,:)=[];
    v(abs(y)>Yrange,:)=[];
    vorticity(abs(y)>Yrange,:)=[];
    yy = y;
    yy(abs(y)>Yrange)=[];
    psipsi = 2*real(psiM*exp(1i*alpha*x));
    uu=2*real(u*exp(1i*alpha*x));
    vv=2*real(v*exp(1i*alpha*x));
    vorticityvorticity=2*real(vorticity*exp(1i*alpha*x));
    N = length(psiM);
    sely=1:2:N;
    subplot(1,3,2:3); hold off;
    contourf(x,yy,uu,10); hold on; 
    quiver(x,yy(sely),uu(sely,:),vv(sely,:),'k'); hold on;
    %axis equal;
    xlabel('x'); ylabel('y'); title({'Structure of the eigenmode (2D reconstruction)',['for k = ',num2str(alpha) , ' ; omega = ',num2str(omega)]});
end% function plotmode

%{

# Exercices/Contributions

* Please check the convergence of the results by comparing with other discretization methods (e.g. Hermite)
* Please modify the program to treat the case of a 2D jet defined as $U(y) = 1/cosh(y)^N$ where $N$ is an integer 
  defining the 'steepness' of the profile (try for instance $N=1$ for a smooth profile and $N=10$ for a very steep profile). 
  Compare the latter case with theoretical results for a "top-hat" jet. 
  You may need to adapt the mesh parameters to have enough points in the region of the steep gradients ! 

%}