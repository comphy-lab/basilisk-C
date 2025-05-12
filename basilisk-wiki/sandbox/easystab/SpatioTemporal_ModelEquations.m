function [] = main()

close all;
global model;
model = 'Crow';
R = 1;

%{

# Spatial stability analysis of a model equation

This program computes the spatial branches $k(\omega)$ for a few
model equations of the form $D(k,\omega,R)=0$ where $R$ is a control parameter modelling the effect of convection.

Typically, the flow is *convectivelty unstable* for $R<R_c$ where $R_c$ is some threshold depending upon the problem and becomes *absolutely unstable* for $R>R_c$. In the convective casd, spatial stability analysis leads to two well-separated kind of branches calledk^+$ and $k^-$, while in the absolute case this disctinction can no longer be made. This program illustrates this transition for some model equaions. 


Several cases are available :

* Kupfer equation (cf. Schmid & Henningson's book)
* Crow equation (cf. Fabre, Cossu & Jacquin 2000).
* Model equation for shear layer (cf. M2R-DET lecture 9 and exercice ; to be finished)
* Ginsburg-Landau equation ; to be finalized)

%}

%{
## Temporal analysis 
%}

ktab = [0.01:.002:2];
omegaguess = ktab(1)*1i;
for j = 1:length(ktab)
    k = ktab(j);
    omegatab(j) = myfzero(@(omega)(reldisp(k,omega,R)),omegaguess);
    omegaguess = omegatab(j);
end
omegaguess = ktab(1)*(-1i);
for j = 1:length(ktab)
    k = ktab(j);
    omegatab2(j) = myfzero(@(omega)(reldisp(k,omega,R)),omegaguess);
    omegaguess = omegatab2(j);
end


figure(1);
%subplot(2,1,1);
title(['Temporal stability analysis for the ',model,' Model']);hold on;
xlabel('k');ylabel('\omega_i');
plot(ktab,imag(omegatab),ktab,imag(omegatab2));
legend('unstable branch','stable branch');
%subplot(2,1,2);
%xlabel('k');ylabel('c_r');legend('unstable branch','stable branch');
%plot(ktab,real(omegatab)./ktab,ktab,real(omegatab2)./ktab);

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','Crow_temporal.png');
%{
![** Figure : Temporal stability results](SpatioTemporal_ModelEquations/Crow_temporal.png)
%} 


%{
# Spatial analysis 
%}

switch model
    
   
    case 'Crow'

figure(2);
R = 0.9;
omegatab = [.1:.01:5];
 kplus = spatialbranch(R,omegatab,(1-.1i)*omegatab(1)/(R*1.61));
 kmoins = flip(spatialbranch(R,flip(omegatab),-0.5-1.8i));
 subplot(2,2,1);
 title(['k(omega) for R = ' ,num2str(R)]);hold on;
 plot(omegatab,-imag(kplus),'r',omegatab,-imag(kmoins),'b');
 xlabel('\omega');ylabel('-k_i');legend('k+ branch','k- branch');
 subplot(2,2,3);
 %xlabel('k_r');ylabel('-k_i');hold on;
 %plot(real(kplus),-imag(kplus),'r',real(kmoins),-imag(kmoins),'b')%,real(kmoins2),-imag(kmoins2),'g');
 xlabel('\omega');ylabel('c_r');hold on;
 plot(omegatab,real(kplus)./omegatab,'r',omegatab,real(kmoins)./omegatab,'b');
 ylim([0 2]);
 legend('k+ branch','k- branch');
% 
 R = 1.1;
 kplus = spatialbranch(R,omegatab,(1-.1i)*omegatab(1)/(R*1.61)); 
 ks = 2.2600 + 0.0010i;
 kmoins = flip(spatialbranch(R,flip(omegatab),ks));
 subplot(2,2,2);
 title(['k(omega) for R = ' ,num2str(R)]);hold on;
 xlabel('\omega');ylabel('-k_i');
 plot(omegatab,-imag(kplus),'r',omegatab,-imag(kmoins),'b');
 legend('branch 1','branch 2');
 subplot(2,2,4);
 %xlabel('k_r');ylabel('-k_i');hold on;
 %plot(real(kplus),-imag(kplus),'r',real(kmoins),-imag(kmoins),'b');
 xlabel('\omega');ylabel('c_r');hold on;
 plot(omegatab,real(kplus)./omegatab,'r',omegatab,real(kmoins)./omegatab,'b');
 ylim([0 2]);
 legend('branch 1','branch 2');

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','Crow_spatial.png');
%{
![** Figure : Spatial stability results](SpatioTemporal_ModelEquations/Crow_spatial.png)
%} 


 case 'Kupfer'

figure(2);
R = 4;
kg = Starting(R);
omegatab = [0:.01:1]*50;
kplus = spatialbranch(R,omegatab,kg(1)); 
kmoins1 = spatialbranch(R,omegatab,kg(2));
kmoins2 = spatialbranch(R,omegatab,kg(3));
subplot(2,2,1);
title(['k(omega) for R = ' ,num2str(R)]);hold on;
xlabel('\omega');ylabel('-k_i');%legend('unstable branch','stable branch');
plot(omegatab,-imag(kplus),'r',omegatab,-imag(kmoins1),'b',omegatab,-imag(kmoins2),'g');
subplot(2,2,3);
xlabel('k_r');ylabel('-k_i');%legend('unstable branch','stable branch');hold on;
plot(real(kplus),-imag(kplus),'r',real(kmoins1),-imag(kmoins1),'b',real(kmoins2),-imag(kmoins2),'g');

R = .99;
kg = Starting(R);
omegatab = [0:.01:1];
kplus = spatialbranch(R,omegatab,kg(1)); 
kmoins1 = spatialbranch(R,omegatab,kg(2));
kmoins2 = spatialbranch(R,omegatab,kg(3));
subplot(2,2,2);
title(['k(omega) for R = ' ,num2str(R)]);hold on;
xlabel('\omega');ylabel('-k_i');%legend('unstable branch','stable branch');
plot(omegatab,-imag(kplus),'r',omegatab,-imag(kmoins1),'b',omegatab,-imag(kmoins2),'g');
subplot(2,2,4);
xlabel('k_r');ylabel('-k_i');legend('unstable branch','stable branch');hold on;
plot(real(kplus),-imag(kplus),'r',real(kmoins1),-imag(kmoins1),'b',real(kmoins2),-imag(kmoins2),'g');

end % switch


end %function main


function D = reldisp(k,omega,R)
global model

switch model
    
    case('Piecewise')
        
%{ 

### Piecewise shear layer (Godreche & Manneville)        
        
** Analysis

We write :
$$
\psi = \left\{ \begin{array}{ll} 
 B_1 e^{k y} , & y<-\delta\\
 A_0 e^{- k y} + B_0 e^{- k y}, \quad  & -\delta<y< \delta \\
 A_2 e^{-k y} , & y>\delta
\end{array}
\right.
$$

%} 

    D = 4*k^2*(omega/k-1)^2-R^2*((2*k-1)^2-exp(-4*k))+1e-3i;
    
     case('Kupfer')
    
    % model equation following Kupfer (1987); see Schmid & Henningson, sec. 7.2
         
    D = omega - ( (k-1i)^3/3 + 1i - k*R); 
    
     case('Crow')
    
    % model equation of the Crow instability ; cf. Fabre Cossu & Jacquin (2000)
    g = -log(1/3)-0.577+1/4;
    s = 1; %sign(real(k));
    nu = 1e-2;
    varpi = k^2/2*(log(2/k*s)+g);
    phi = s*k*besselk(1,s*k)+k^2*besselk(0,s*k);
    chi = s*k*besselk(1,s*k);
    Req = 1.61/R;
    D = (omega-k*Req)*(omega-k*Req+nu*1i*k^2) + (1-phi+varpi)*(1+chi-varpi);
   
    case('GL')
    sigma0 = 0.1;
    D = (omega-k)-R*1i*k^2;
    
    case('ShearLayerModel')
    
    D = (omega-k)-R*1i*k*(1-k);
    

end
end

function x = myfzero(F,guess)
% this function finds a zero (possibly complex) of the the function F.
x = guess;
A = 1;
j = 1;
epsilon = 1e-6;
while (A>1e-10)&&(j<100)
    Fx = feval(F,x);
    dFx = (feval(F,x+epsilon)-feval(F,x-epsilon))/(2*epsilon);
    x = x-Fx/dFx;
    A = abs(Fx);
    j = j+1;
end
if j==100
    warning('non convergence in Newton!')
end
end

function kbranch = spatialbranch(R,omegatab,kguess) 
for j = 1:length(omegatab)
    omega = omegatab(j);
    kbranch(j) = myfzero(@(k)(reldisp(k,omega,R)),kguess);
    kguess = kbranch(j);
end
end

function ks = Starting(R)
% starting points for the spatial branches for omega = 0 (Kupfer)
global model
switch model

    case 'Kupfer'
ks(1) = -(12i*R+(-12i)+4*sqrt(-4*R^3-9*R^2+18*R-9))^(1/3)/4 ...
      -R/(12i*R-12i+4*sqrt(-4*R^3-9*R^2+18*R-9))^(1/3)+1i+.5i*sqrt(3.D0)...
      *((12i*R-12i+4*sqrt(-4*R^3-9*R^2+18*R-9))^(1/3)/2-2*R/(12i*R+(-12i)...
      +4*sqrt(-4*R^3-9*R^2+18*R-9))^(1/3));

ks(2) =  -(12i*R+(-12i)+4*sqrt(-4*R^3-9*R^2+18*R-9))^(1/3)/4 ...
      -R/(12i*R-12i+4*sqrt(-4*R^3-9*R^2+18*R-9))^(1/3)+1i-.5i*sqrt(3.D0)...
      *((12i*R-12i+4*sqrt(-4*R^3-9*R^2+18*R-9))^(1/3)/2-2*R/(12i*R+(-12i)...
      +4*sqrt(-4*R^3-9*R^2+18*R-9))^(1/3));
  
ks(3) = (12i*R+(-12i)+4*sqrt(-4*R^3-9*R^2+18*R-9))^(1/3)/2 ...
     +2*R/(12i*R/(-12i)+4*sqrt(-4*R^3-9*R^2+18*R-9))^(1/3)+1i;
    
    case 'Crow'
     
       omega = .1
       ks(1) = omega*(1-.1i)/(R/1.61);       
       ks(2) = omega/(R/1.61)*(1-1i);
end

end


