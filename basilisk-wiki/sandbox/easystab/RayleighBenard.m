%{
# Rayleigh-Benard instability

*This document belongs to the [easystab](http://basilisk.fr/sandbox/easystab/README) project, please consult the main page of the project for explanations.*

*The present program was specifically designed as a support for the course ["Introductions to hydrodynamical instabilities"](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) for the "DET" Master's cursus in Toulouse. For the underlying theory please see [Lecture notes for chapter 5](http://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md) *

*To use this program: click "raw page source" in the left column, copy-paste in the Matlab/Octave editor, and run.*


This program solves the linear stability problem for the Rayleigh-Bénard instability 
of a horizontal layer of fluid heated from below. This program is adapted from [/sandbox/easystab/rayleigh_benard.m](), which
solved the problem in the case of slip conditions ("Free-Free boundaries") and validated the
approach by comparing with analytical results which exist in this case. 
The present program considers the more physical case of no-slip conditions at the walls
("Rigid-Rigid boundaries").

The theory is presented in [lecture 5](http://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md) 
of the M2 DET course. The  problem is set in generalised eigenvalue form 
$\lambda B X = A X$. The construction of the matrices and resolution is done by [function RB(k,Ra,Pr)](#function-rb) defined at the end of this document.

%}


function [] = RayleighBenard()
global y dy dyy w Z I N

% parameters
N=61; % number of gridpoints

% differentiation matrices
[dy,dyy,w,y] = dif1D('cheb',0,1,N);
Z=zeros(N,N); I=eye(N); 


%{ 

## Validation of the code

We first validate the code by looking for the instability threshold.

According to the litterature, the neutral conditions are :

%}
Ra = 1708;
k = 3.117;

%{
These value correspond to a neutral mode whatever the Prandtl number. Here
we take : 
%}

Pr = 10.;

%{
We compute the leading eigenvalue using the function [FK](#function-fk)
defined at the bottom of this program:
%}

lambdamax = RB(k,Ra,Pr)

%{
This value is close to zero. If we want an even better evaluation of the
threshold we can use fzero to find the value of Ra where the eigenvalue
is exactly zero:
%}

Rac = fzero(@(Ra)(real(RB(k,Ra,Pr))),Ra)


%{    
## Spectrum, and structure of the leading eigenmode for Ra = 2000.

We now take a value of Ra above the threshold and solve the eigenvalue
problem for these parameters :

%}

Ra = 2000;
Pr = 10;
k = pi ;% this means that the wavelenghth is twice the plate spacing)

[s,U] = RB(k,Ra,Pr);

%{
We plot the spectrum in figure 1
%}

%set(0,'defaultAxesFontSize',16);
figure(1);
plot(imag(s),real(s),'go');
hold on; plot([-1, 1],[0,0],'k:');
ylabel('\lambda_r');ylim([-200, 50]);
xlabel('\lambda_i');xlim([-1,1]);
title(['Spectrum  for Ra =',num2str(Ra),' ; k = ',num2str(k),' ; Pr = ',num2str(Pr)]),
saveas(gcf,'RB_spectrum','svg');

%{ 
![Figure : Spectrum for Ra = 2000](RayleighBenard/RB_spectrum.svg)
%}


Mode = U(:,1);
figure(2);
ModeU = imag(Mode(1:N)); ModeV = real(Mode(N+1:2*N)); ModeP = real(Mode(2*N+1:3*N));ModeT = real(Mode(3*N+1:4*N));
plot(ModeU/max(abs(ModeU)),y,'r',ModeV/max(abs(ModeV)),y,'b',ModeP/max(abs(ModeP)),y,'g',ModeT/max(abs(ModeT)),y,'k');
ylabel('y');
legend({'Im( u(y) )','v(y)','p(y)','\theta(y)'});
title({['Leading Eigenmode for Ra =',num2str(Ra)], 'Plot of eigenfunctions' });
saveas(gcf,'RB_Mode','svg');

 
%{ 
![Figure : Most amplified eigenmode for Ra = 2000](RayleighBenard/RB_Mode.svg)
%}

%{
 
Reconstruct the 2D structure of the mode and plot

%}
 

Amp=0.2;step = 3;
Xarray = linspace(0,4,40);
Yarray = y(1:step:end);
[XX,YY] = meshgrid(Xarray,Yarray);
UU = -Amp*ModeU(1:step:end)*sin(Xarray*pi);
VV = Amp*ModeV(1:step:end)*cos(Xarray*pi);
TT = Amp*ModeT(1:step:end)*cos(Xarray*pi);

figure(20);
contourf(XX,YY,TT);hold on;
quiver(XX,YY,UU,VV,'k');hold off;
axis equal;
title({['Leading Eigenmode for Ra =',num2str(Ra)], '2D reconstruction' });
%saveas(gcf,'RB_mode2D_perturb','svg');
print('-dpng','-r120','RB_mode2D_perturb.png');
%{ 
![Figure : 2D reconstruction of leading eigenmode for Ra = 2000 (colors are for **temperature perturbations**)](RayleighBenard/RB_mode2D_perturb.png)
%}

TTtot = TT+(1-Yarray)*ones(1,40);
figure(21);
contourf(XX,YY,TTtot);hold on;
quiver(XX,YY,UU,VV,'k');hold off;
axis equal;
title({'2D reconstruction of eigenmode ',' superposed to temperature field at equilibrium' });
%saveas(gcf,'RB_mode2D_full','svg');
print('-dpng','-r120','RB_mode2D_full.png');
  
%{ 
![Figure : 2D reconstruction of leading eigenmode for Ra = 2000 (colors are for full temperature field)](RayleighBenard/RB_mode2D_full.png)
%}

pause(10);
disp('Next part of the program, parametric study, will start in 10 seconds...');

%{

## Parametric study
 
We have to characterize the variation of the leading $\lambda$ as function
of the two parameters $k$ and $Ra$ (we still consider $Pr = 10$).

### Computations
 
First we compute $\lambda(k)$ for several values of $Ra$ and plot the
results in figure 2.
  
%}

for Ra = [1500:250:2500];
    ktab = [1:.1:5];
    smaxtab = [];
    for k = ktab
        [s,U] = RB(k,Ra,Pr);
        smaxtab = [smaxtab real(s(1))];
    end
    figure(3);
    subplot(2,1,1);
    plot(ktab,smaxtab);hold on;
    pause(0.1);
end

plot(ktab,0*ktab,'k:');
xlabel('k');
ylabel('\lambda');
%title('Growth rate \lambda(k) for various values of Ra');
legend('Ra=1500','Ra=1750','Ra=2000','Ra=2250','Ra=2500');
legend('Location','SouthEast');
pause(0.1);
    
%{
We now build the marginal stability curve $Ra_c(k)$ corresponding to the location in the $[Ra,k]$-plane where the 
leading eigenvalue is exactly zero. 
    
For this we do a loop over k and look for the location where
$\lambda_{max}$ is exactly zero, using again fzero as done above.
Results are plotted in figure 2.    
%}
    
Ra = 3000;
Ratab=[];
ktab = 1.5:.1:5;
for k=ktab
    Ra = fzero(@(Ra)(real(RB(k,Ra,Pr))),Ra);
    Ratab= [Ratab Ra];
end
    
figure(3);
subplot(2,1,2);
plot(ktab,Ratab);
xlabel('k');ylabel('Ra');
ylim([1000,3000]);
%title('Neutral curve');
saveas(gcf,'RB_neutralcurve','svg');

### Results

%{ 
![Figure : amplification rate as function of $k$ for several values of $\lambda$ (top) ; boundary of the unstable domain in the (Ra/k) plane (bottom)](RayleighBenard/RB_neutralcurve.svg)
%}


end



%{

# FUNCTION RB

Here is the definition of the function RB which performs the eigenvalue
computation. Note that this function is designed so that it can be used
in two ways :
    
- [s,U] = RB(k,Ra,P) will return a vector s containing the 10 leading
eigenvalues and an array U containing (in column) the corresponding eigenvectors.
    
    
- lambdamax = RB(k,Ra,P) will return only one value corresponding to the leading eigenvalue. 
This is useful to use this function with fzero.
    
%}

function [s,U] = RB(k,Ra,Pr)
global y dy dyy w Z I N


% renaming the differentiation matrices
dx=1i*k*I; dxx=-k^2*I;
Delta=dxx+dyy;

%{
## System matrices


As explained in  the [lecture notes](http://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md#case-of-a-horizontal-cell-of-large-dimension), 
after nondimensionalization the system of equations can be written in a matrix form
$$
\lambda B q =A q
$$
with the matrices
$$
q=\left(\begin{array}{c}
u \\ v\\ p\\ \theta
\end{array}\right)
, \quad
A=\left(\begin{array}{cccc}
Pr \Delta&0&-\partial_x&0\\
0&Pr \Delta&-\partial_y& Pr \\
\partial_x&\partial_y&0&0\\
0& Ra &0&\Delta
\end{array}\right)
, \quad
B=\left(\begin{array}{cccc}1&0&0&0\\0&1&0&0\\0&0&0&0\\0&0&0&1\end{array}\right)
$$

NB in this program we use a slightly different nodimensionalization 
in which the velocities are resclaed by $\sqrt(Re)$; this is equivalent.


%}

% system matrices
A=[Pr*Delta, Z, -sqrt(Ra)*dx, Z; ...
   Z, Pr*Delta, -sqrt(Ra)*dy, Pr*sqrt(Ra)*I;  ...
   dx, dy, Z, Z;  ...
   Z, sqrt(Ra)*I, Z, Delta];
B=blkdiag(I,I,Z,I);

%{
## Boundary conditions

The natural conditions are that $u$ and $v$ are zero at the walls, 
and the temperature perturbation $\theta$ is also zero at the walls 
(the temperature is imposed at the wall, without perturbations). 

$$
\begin{array}{l}
u|_0=0\\
u|_L=0\\
v|_0=0\\
v|_L=0\\
\theta|_0=0\\
\theta|_L=0\\
\end{array}
$$
thus the boundary conditions are expressed $Cq=0$, with the constraint matrix

$$
C=\left(\begin{array}{cccc}
I|_0&0&0&0\\
I|_L&0&0&0\\
0&I|_L&0&0\\
0&I|_0&0&0\\
0&0&0&I|_L\\
0&0&0&I|_0\\
\end{array}\right)
$$

%}

% boundary conditions
II=eye(4*N); 
u0=1; uL=N; v0=N+1; vL=2*N; T0=3*N+1; TL=4*N;
loc=[u0,uL,v0,vL,T0,TL];
C(loc,:)=II(loc,:);
A(loc,:)=C(loc,:);
B(loc,:)=0; 

%{
## Computing eigenmodes

In the way the linear system $\lambda B q=Aq$ hzs been constucted, the matrix $B$ is not invertible, so we solve a generalized eigenvalue problem which will have infinite eigenvalues. We remove them form the results.  
%}

% computing eigenmodes
[U,S]=eig(A,B);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

if(nargout==1) 
    s=(s(1)); 
    % this is a trick to allow to use the function as s=RT(k,Ra,Pr). 
    % This allows to pass the function as a handle for fzero.
end

end

function [] = plotmode(~,~,mode,lambda,k) % function to plot one mode
    global y dy dyy
    figure(2);hold off;
    N = length(y);
    u = mode(1:N);
    v = mode(N+1:2*N);
    p = mode(2*N+1:3*N);
  %  p = p/max(p)*max(v); % renormalization
    T = mode(3*N+1:4*N);
    subplot(1,3,1);hold off;
    plot(real(u),y,'b-',imag(u),y,'b--');hold on;
    plot(real(v),y,'g-',imag(v),y,'g--');hold on;
    plot(real(p),y,'r-',imag(p),y,'r--');hold on;
    plot(real(T),y,'k-',imag(T),y,'k--');hold on;
    ylabel('y'); 
    legend({'Re(u)','Im(u)','Re(v)','Im(v)','Re(p)','Im(p)','Re(T)','Im(T)'})
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
    TT=2*real(T*exp(1i*k*x));
    subplot(1,3,2:3); hold off;
    contourf(x,yy,TT,10); hold on; 
    quiver(x,yy,uu,vv,'k'); hold on;
    %axis equal;
    xlabel('x'); ylabel('y'); title({'Structure of the eigenmode ( color : \theta ; vectors : (u,v)  )',['lambda = ',num2str(lambda)]});
end% function plotmode

%{


# Exercices/Contributions

* Please superpose the neutral curve onto the results of Drazin & Reid
(Fig. 2.2, page 53)
* Drazin & Reid (page 51) found that a second mode becomes unstable above a
critical Rayleigh number of 17610 for k=5.365. Please check these results
with the present code.
* Please show that the flow is always stable when the hot plate is at the top
* Please draw the velocity field and the temperature field corresponding to the five most unstable modes.
* Please write a code that does the advection of tracer particles by the velocity field of the first eigenmode in a stable case and in an unstable case
* Extension considering surface tension and a free surface, also called [Rayleigh Benard Marangoni convection](stab2014/rayleigh_benard_marangoni.m)
%}
