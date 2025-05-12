%{
*This document belongs to the [easystab](http://basilisk.fr/sandbox/easystab/README) project, check the main page to understand the general philosophy of the project.*

*The present program was specifically designed as a support for the course ["Introductions to hydrodynamical instabilities"](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) for the "DET" Master's cursus in Toulouse. For the underlying theory please see [Lecture notes for chapter 7](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md#lecture-7-shear-flow-instabilities-ii-viscous-instabilities)*

*To use this program: click "raw page source" in the left column, copy-paste in the Matlab/Octave editor, and run.*



# Stability of the plane Poiseuille flow in primitive variables.

This code is adapted from [sandbox/easystab/poiseuille_uvp.m]().

The numerical treatment of the problem is the same as in this initial program,
but here the eigenvalue calculation is done inside a function
"EV_Poiseuille" to allow to perform loop over the parameters.

%}

%{
## Definition of geometry and differential operators
%}


clear all; clf;
global y dy dyy Z I INT n
L=2;        % the height of the domain, from -1 to 1
n=100;      % the number of grid points

% differentiation and integration operators
scale=-2/L;
[y,DM] = chebdif(n,2); 
dy=DM(:,:,1)*scale;    
dyy=DM(:,:,2)*scale^2;    
y=y/scale; 
Z=zeros(n,n); I=eye(n); 
INT=([diff(y)',0]+[0,diff(y)'])/2; 


%{

## Validation 

The valudation is done in the program [sandbox/easystab/poiseuille_uvp.m]()
by comparing with results from Drazin & Reid and is not reproduced here
%}


%{
## Curves lambda(k) and c(k) for various values of Re
%}

for Re = [5772 6000 6500 7000 8000 10000 20000 50000 100000] 
 alphatab = [0.5:0.05:1.5];
 stab = [];
  for alpha = alphatab
  s = EV_Poiseuille(Re,alpha); 
  stab = [stab s];
end
    
  figure(1);
   subplot(2,1,1);
   plot(alphatab,real(stab));hold on;
   grid on;
   subplot(2,1,2);
   plot(alphatab,-imag(stab)./alphatab); hold on;
   pause(0.1); 
    
end
    
   figure(1);subplot(2,1,1);
   title('Growth rate');
   xlabel('k');ylabel('\lambda_r');    figure(1);subplot(2,1,2);
   title('Phase velocity');
   xlabel('k');ylabel('c_r');
   legend('Re = 5772', 'Re=6000', 'Re=8000','Re=10 000','Re=20 000','Re= 50 000', 'Re = 100 000');
    
   set(gcf,'paperpositionmode','auto');
   print('-dpng','-r100','PlanePoiseuille_SigmaOmegaK.png');

%{
![Growth rate and phase velocity as function of k for various Re](/sandbox/easystab/david/PlanePoiseuille_SigmaOmegaK.png)
%}


%{

## Construction of the Marginal stability curve

%}
   
  Retab = [5772 5800 6000 6500 7000  8000 10000 15000 20000 35000 50000 75000 100000]; 
  % first point : from known result
  kinftab = 1.02;
  ksuptab = 1.02; 
  % second point : guess 
  Re = Retab(2)
  kinf = fzero(@(alpha)(real(EV_Poiseuille(Re,alpha))),0.97)
  ksup = fzero(@(alpha)(real(EV_Poiseuille(Re,alpha))),1.05)
  kinftab = [kinftab kinf];
  ksuptab = [ksuptab ksup];
  
  for Re = Retab(3:end)
     Re
     kinf = fzero(@(alpha)(real(EV_Poiseuille(Re,alpha))),kinftab(end))
     ksup = fzero(@(alpha)(real(EV_Poiseuille(Re,alpha))),ksuptab(end))
     kinftab = [kinftab kinf];
     ksuptab = [ksuptab ksup];
  end
  
  figure(3);
  plot(Retab,kinftab,'b',Retab,ksuptab,'b');
 xlabel('Re');ylabel('k');
 set(gcf,'paperpositionmode','auto');
 print('-dpng','-r100','PlanePoiseuille_NeutralCurve.png');
 
 %{
![Marginal stability curve](/sandbox/easystab/david/PlanePoiseuille_NeutralCurve.png)
%}
 
  
%{
## Function performing the eigenvalue computation
%}

function [s,U] = EV_Poiseuille(Re,alpha)

global y dy dyy Z I INT n

dx=i*alpha*I; dxx=-alpha^2*I;
Delta=dxx+dyy;

% base flow
U=1-y.^2; Uy=-2*y;
S=-diag(U)*dx+Delta/Re;

% the matrices
A=[ ...
    S, -diag(Uy), -dx; ...
    Z, S, -dy; ...
    dx, dy, Z];
E=blkdiag(I,I,Z);


% locations on the grid
u=1:n; v=n+1:2*n; p=2*n+1:2*n;

% boundary conditions
III=eye(3*n);
loc=[u(1) u(n) v(1) v(n) ];  
C=III(loc,:);
E(loc,:)=0;  A(loc,:)=C;



% computing eigenmodes 
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

if(nargout==1) 
    s=s(1); 
    % this is a trick to allow to use the function as s=EV_Poiseuille(Re,alpha). 
    % This allows to pass the function as a handle for fzero.
end

end % function EV_Poiseuille
