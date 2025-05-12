%{
# Rayleigh-Benard instability

This program solves the linear stability problem for the Rayleigh-Bénard instability 
of a horizontal layer of fluid heated from below. This program is adapted from [/sandbox/easystab/rayleigh_benard.m](), which
solved the problem in the case of slip conditions ("Free-Free boundaries") and validated the
approach by comparing with analytical results which exist in this case. 
The present program considers the more physical case of no-slip conditions at the walls
("Rigid-Rigid boundaries"). In this case no simple analytical solutions exist, but numerical solutions
are given in many books (Drazin & Reid, Chandrasekar, etc...)

%}


function [] = RayleighBenard()
set(0,'defaultAxesFontSize',18);
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

Pr = 10.

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

Rac = fzero(@(Ra)(RB(k,Ra,Pr)),Ra)


%{    
## Spectrum, and structure of the leading eigenmode for Ra = 2000.

We now take a value of Ra above the threshold and solve the eigenvalue
problem for these parameters :

%}

Ra = 2000;
k = pi ;% this means that the wavelenghth is twice the plate spacing)

[s,U] = RB(k,Ra,Pr);

%{
We plot the spectrum and the structure of the leading eigenmode in figure 1
%}

figure(1);
subplot(2,1,1);
plot(real(s),imag(s),'go');
xlabel('\lambda_r');
ylabel('\lambda_i');
title(['Spectrum for Ra =',num2str(Ra),' ; k = ',num2str(k),' ; Pr = ',num2str(Pr)]),

Mode = U(:,1);
subplot(2,1,2);
ModeU = imag(Mode(1:N)); ModeV = real(Mode(N+1:2*N)); ModeP = real(Mode(2*N+1:3*N));ModeT = real(Mode(3*N+1:4*N));
plot(y,ModeU/max(abs(ModeU)),'r',y,ModeV/max(abs(ModeV)),'b',y,ModeP/max(abs(ModeP)),'g',y,ModeT/max(abs(ModeT)),'k');
legend({'$Im(\hat{u})$','$\hat{v}$','$\hat{p}$','$\hat{\theta}$'},'Interpreter','Latex');
title(['Leading Eigenmode for Ra =',num2str(Ra)]);

 set(gcf,'paperpositionmode','auto');
 print('-dpng','-r100','RB_spectrumandmode.png');

%{
 
Reconstruct the 2D structure of the mode

%}
 
figure(20);
Amp=0.2;
Xarray = linspace(0,4,60);
Yarray = y;
[XX,YY] = meshgrid(Xarray,Yarray);
UU = Amp*ModeU*sin(Xarray*pi);
VV = Amp*ModeV*cos(Xarray*pi);
TT = Amp*ModeT*cos(Xarray*pi);
contourf(XX,YY,TT);hold on;
quiver(XX,YY,UU,VV,'k');hold off;
axis equal;
title('Rayleigh-Bénard mode : [u,v] (vectors) and T'' (colors) ')  

 
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
    figure(2);
    subplot(2,1,1);
    plot(ktab,smaxtab);hold on;
    pause(0.1);
end
    plot(ktab,0*ktab,'k:');
    xlabel('k');
    ylabel('\lambda');
    title('Growth rate \lambda(k) for various values of Ra');
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
    Ra = fzero(@(Ra)(RB(k,Ra,Pr)),Ra);
    Ratab= [Ratab Ra];
end
    
    figure(2);
    subplot(2,1,2);
    plot(ktab,Ratab);
    xlabel('k');ylabel('Ra_c(k)');
    ylim([1000,3000]);
    title('Neutral curve');
    
    set(gcf,'paperpositionmode','auto');
    print('-dpng','-r100','RB_neutralcurve.png');

### Results

%{ 
![The figure](david/RB_neutralcurve.png)
%}


end

%{ 
![The figure](david/RB_spectrumandmode.png)



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


As explained in  the [lecture notes](http://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md#case-of-a-horizontal-cell-of-large-dimension), the system of equations can be written in a matrix form
$$
Eq_t=Aq
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
%}

% system matrices
A=[Pr*Delta, Z, -dx, Z; ...
   Z, Pr*Delta, -dy, Pr*I;  ...
   dx, dy, Z, Z;  ...
   Z, Ra*I, Z, Delta];
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

For the linear system $Eq_t=Aq$ we have imposed the boundary conditions in a way that $E$ is not invertible, 
so we solve a generalized eigenvalue problem which will have infinite eigenvalues corresponding to the boundary 
conditions constraints. We remove them form the result and we plot only the five eigenvalues that have largest growth rate. 
These modes will be relevant for the stability analysis. here the wave speed of the mode is zeor, these are stationnary mode, 
they have no propagation in the $x$ direction, so we plot only the growth rate.
%}

% computing eigenmodes
[U,S]=eig(A,B);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

if(nargout==1) 
    s=real(s(1)); 
    % this is a trick to allow to use the function as s=RT(k,Ra,Pr). 
    % This allows to pass the function as a handle for fzero.
end

end

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
