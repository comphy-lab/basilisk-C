%{
# Stability of the plane Poiseuille flow in primitive variables.

This code is adapted from [sandbox/easystab/poiseuille_uvp.m]() from the easystab project.

The numerical treatment of the problem is the same as in this initial program,
but here the eigenvalue calculation is done inside a function
"EV_Poiseuille" to allow to perform loop over the parameters.

%}

%{
## Definition of geometry and differential operators
%}

 
function [] = TS_PlanePoiseuille()
close all; 
global y dy dyy Z I INT n

% Physical parameters
L=2;        % the height of the domain, from -1 to 1
Re = 5e4;
alpha = .75;

% numerical parameters
n=100;      % the number of grid points
iloops = 1; % set to zero to skip the parametric study

% differentiation and integration operators
[dy,dyy,wy,y] = dif1D('cheb',-1,2,n); 

Z=zeros(n,n); I=eye(n); 
INT=([diff(y)',0]+[0,diff(y)'])/2; 


%{

## Resolution of the eigenvalue problem

The validation is done in the program [sandbox/easystab/poiseuille_uvp.m]()
by comparing with results from Drazin & Reid and is not reproduced here
%}

[s,UU] = EV_Poiseuille(Re,alpha); 
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
ylim([-.25 .05]);
xlim([0 1.1*alpha]);
title({['Temporal spectrum for k = ',num2str(alpha), ', Re = ',num2str(Re)], 'Click on eigenvalues to see the eigenmodes'});
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','spectrum.png');

%{
![** Figure : Temporal spectrum of the plane poiseuille flow](TS_PlanePoiseuille/spectrum.png)
%} 


plotmode([],[],UU(:,1),omega(1),alpha);

pause(0.1);
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','mode.png');
%{
![** Figure : Unstable eigenmmode of the plane poierull flow  for k=0.75](TS_PlanePoiseuille/mode.png)
%} 



%{
## Curves lambda(k) and c(k) for various values of Re
%}

if iloops
 for Re = [5772 6000 8000 10000 20000 50000 100000] 
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
   xlabel('k');ylabel('\lambda_r');    
   figure(1);subplot(2,1,2);
   title('Phase velocity');
   xlabel('k');ylabel('c_r');
   legend('Re = 5772', 'Re=6000', 'Re=8000','Re=10 000','Re=20 000','Re= 50 000', 'Re = 100 000');
    
   set(gcf,'paperpositionmode','auto');
   print('-dpng','-r100','PlanePoiseuille_SigmaOmegaK.png');

%{
![Growth rate and phase velocity as function of k for various Re](TS_PlanePoiseuille/PlanePoiseuille_SigmaOmegaK.png)
%}


%{

## Construction of the Marginal stability curve

%}
   
  Retab = [5772 5800 6000 8000 10000 15000 20000 35000 50000 75000 100000]; 
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
 title('Marginal stability curve of the plane Poiseuille flow')
 set(gcf,'paperpositionmode','auto');
 print('-dpng','-r100','PlanePoiseuille_NeutralCurve.png');
 
%{
![Marginal stability curve](TS_PlanePoiseuille/PlanePoiseuille_NeutralCurve.png)
%}
end
 
end

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


function [] = plotmode(~,~,mode,omega,alpha)
    global y dy dyy
    Yrange = 1;
    figure(2);hold off;
    N = length(y);
    u = mode(1:N);
    v = mode(N+1:2*N);
    p = mode(2*N+1:3*N);
    vorticity = (dy*u)-1i*alpha*v;
    subplot(1,3,1);hold off;
    plot(real(u),y,'b-',imag(u),y,'b--');hold on;
    plot(real(v),y,'g-',imag(v),y,'g--');hold on;
    plot(real(p),y,'k-',imag(p),y,'k--');hold on;
    ylabel('y'); ylim([-Yrange,Yrange]);
    legend({'$Re(\hat u)$','$Im(\hat u)$','$Re(\hat v)$','$Im(\hat v)$','$Re(\hat p)$','$Im(\hat p)$'},'Interpreter','latex')
    title('Structure of the eigenmode');
    % plot 2D reconstruction
    Lx=2*pi/alpha; Nx =30;  x=linspace(-Lx/2,Lx/2,Nx);
    p(abs(y)>Yrange,:)=[];
    u(abs(y)>Yrange,:)=[];
    v(abs(y)>Yrange,:)=[];
    vorticity(abs(y)>Yrange,:)=[];
    yy = y;
    yy(abs(y)>Yrange)=[];
    pp = 2*real(p*exp(1i*alpha*x));
    uu=2*real(u*exp(1i*alpha*x));
    vv=2*real(v*exp(1i*alpha*x));
    vorticityvorticity=2*real(vorticity*exp(1i*alpha*x));
    subplot(1,3,2:3); hold off;
    contourf(x,yy,vorticityvorticity,10); hold on; 
    quiver(x,yy,uu,vv,'k'); hold on;
    xlabel('x'); ylabel('y'); title({'Structure of the eigenmode (2D reconstruction)',['for k = ',num2str(alpha) , ' ; omega = ',num2str(omega)]});
end% function plotmode
