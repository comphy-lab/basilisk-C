%{

# Instabilities of a gaussian jet 
(temporal analysis, viscous case, uvp formulation)

(This code is adapted from [sandbox/easystab/KH_temporal_viscous.m](KH_temporal_viscous.m) from the easystab project).

The base flow is 
$$
U(y) = exp(-y^2)
$$


We look for solutions under eigenmode form :
$$
[u,v,p] = [\hat u(y), \hat v(y), \hat p(y)] e^{i k x} e^{-i \omega t}
$$ 

We will consider separately symmetrical and antisymmetrical modes and use a
semiinfinite domain.
WARNING : only antisymmetrical modes are validated

%}

%function [] = main()
close all;
global y D DD w Z I U Uy discretization



% numerical parameters
N=100;      % the number of grid points
loopk = 1; % set to 0 to skipp the loops over k and Re to build the curves 

%{

### Derivation matrices

Here we use Chebyshev discretization with stretching, considering a doubly
infinite domain

%}

discretization = 'chebInfAlg'; % 
[D,DD,w,y] = dif1D(discretization,0,3,N,0.99999);
 
Z=zeros(N,N); I=eye(N); 
Ndim=N;

% base flow
U=exp(-y.^2); 
Uy=-2*y.*exp(-y.^2); 

%%
%{
We compute the eigenvalues/eigenmodes with the function [KH](#Function KH), 
defined at the end of this program
%}
% physical parameters
alpha=.5;    % the wave number
Re=200;    % the Reynolds number
symmetry = -1; % 1-> symmetric modes ; -1 -> antisymmetric modes 
[s,UU] = KH(alpha,Re,N);
omega = 1i*s;

%{
### Plotting the spectrum
%}

figure(101);hold off;
for ind=1:length(s)
  h=plot(real(omega(ind)),imag(omega(ind)),'*'); hold on
  set(h,'buttondownfcn',{@plotmode,UU(:,ind),omega(ind),alpha});
end
xlabel('\omega_r');
ylabel('\omega_i');
ylim([-1 .25]);
xlim([-.5*alpha 1.5*alpha]);
title({['Temporal spectrum for k = ',num2str(alpha), ', Re = ',num2str(Re)], 'Click on eigenvalues to see the eigenmodes'});
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','KH_temporal_viscous_spectrum.png');

%{
![** Figure : Temporal spectrum of a tanh shear layer](KH_temporal_viscous/KH_temporal_viscous_spectrum.png)
%} 


plotmode([],[],UU(:,1),omega(1),alpha);

pause(0.1);
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','KH_temporal_viscous_mode.png');
%{
![** Figure : Unstable eigenmmode of the tanh shear layer for k=0.5](KH_temporal_viscous/KH_temporal_viscous_mode.png)
%} 



%%
%{

### Now we do loops over k and Re to plot the growth rate curves
$\omega_i(k)$

%}

if loopk
  for Re = 30  %[10 30 100 300 1000]
    alphatab = 0:.025:3;
    lambdatab = [];
    for alpha=alphatab
         smode = s(1);
        if(real(smode)<0)
            sfilter = s(-imag(s)/alpha>0.1);
            [~,ind]= max(real(sfilter));
            smode = sfilter(ind);
        end
        lambdatab = [lambdatab smode];
    end
    figure(3);hold on;
    subplot(2,1,1);
    plot(alphatab(1:length(lambdatab)),real(lambdatab));hold on;
    xlabel('k');ylabel('\omega_i');
    legend('Re=10','Re=30','Re=100','Re=300','Re=1000');
    legend('Location','South');
    subplot(2,1,2);
    plot(alphatab(1:length(lambdatab)),-imag(lambdatab)./alphatab(1:length(lambdatab)));hold on;
    xlabel('k');ylabel('\omega_r/k');
    hold on;
    pause(0.1);
  end

 figure(3);
    %xlabel('k');ylabel('\omega_i');
    %title('growth rate for several values of Re');

    set(gcf,'paperpositionmode','auto');
    print('-dpng','-r80','KH_temporal_viscous_curves.png');
    
end
%%  
%end
%{

![** Figure : results for the Kelvin-Helmholtz instability (temporal) of a tanh shear layer](KH_temporal_viscous/KH_temporal_viscous_curves.png)

%}

%{

#Function KH

%}

%%
function [s,UU] = KH(alpha,Re,N)
global y D DD w Z I U Uy discretization % to use these objects within the function

% renaming the differentiation matrices
dy=D; dyy=DD;
dx=1i*alpha*I; dxx=-alpha^2*I;
Delta=dxx+dyy;

%{
### System matrices
%}

% the matrices
S=-diag(U)*dx+Delta/Re;
A=[ ...
    S, -diag(Uy), -dx; ...
    Z, S, -dy; ...
    dx, dy, Z];
E=blkdiag(I,I,Z);

% Boundary conditions

    % Dirichlet conditions are used for fd, cheb, etc...
    III=eye(3*N);
    indBC=[1,N,N+1,2*N];
    C=III(indBC,:);
 %   if(symmetry==1)
 %       % change condition for ux to Neumann
 %       C(1,1:N) = dy(1,:);
 %   else    
 %       C(3,N+1:2*N) = dy(1,:); 
 %   end
    A(indBC,:)=C;
    E(indBC,:)=0;  

% computing eigenmodes 

[UU,S]=eig(A,E);

% sort the eigenvalues by decreasing real part and remove the spurious ones
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); UU=UU(:,o);
rem=abs(s)>1000; s(rem)=[]; UU(:,rem)=[];

end %function KH

function [] = plotmode(~,~,mode,omega,alpha)
    global y dy dyy
    Yrange = 5;
    figure(102);hold off;
    N = length(y);
    u = mode(1:N);
    v = mode(N+1:2*N);
    p = mode(2*N+1:3*N);
    %vorticity = (dy*u)-1i*k*v;
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
    %vorticity(abs(y)>Yrange,:)=[];
    yy = y;
    yy(abs(y)>Yrange)=[];
    pp = 2*real(p*exp(1i*alpha*x));
    uu=2*real(u*exp(1i*alpha*x));
    vv=2*real(v*exp(1i*alpha*x));
    %vorticityvorticity=2*real(vorticity*exp(1i*alpha*x));
    subplot(1,3,2:3); hold off;
    contourf(x,yy,pp,10); hold on; 
    quiver(x,yy,uu,vv,'k'); hold on;
    %axis equal;
    xlabel('x'); ylabel('y'); title({'Structure of the eigenmode (2D reconstruction)',['for k = ',num2str(alpha) , ' ; omega = ',num2str(omega)]});
end% function plotmode

%{

# Exercices/Contributions

* Please compare the results with the inviscid ones obtained using the program [KH_temporal_inviscid.m]()
* Please validate the results by comparing with the litterature (Drazin & Reid)
* Please try other discretisation methods, for instance Hermite ----> [sandbox/easystab/kelvin_helmholtz_hermite.m](kelvin_helmholtz_hermite.m)
* Please look at the structure of the adjoint eigenmode and compute the nonnormality factor.

%}