%{

# Kelvin-Helmholtz instability of a shear layer 
(Spatial analysis, viscous case, uvp formulation)

(This code is adapted from [sandbox/easystab/KH_temporal_viscous.m](KH_temporal_viscous.m) from the easystab project).

We start from the linearised Navier-Stokes for a parallel base flow defined
as $U(y)$.



$$
\begin{array}{l}
\rho \partial_t u= -U \partial_y u - U_y v -\partial_x p +\mu\Delta u\\
\rho \partial_t v= -U \partial_y v - \partial_y p +\mu \Delta v\\
\partial_x u+ \partial_y v=0.\\
\end{array}
$$

The base flow is 
$$
U(y) = (\Delta U) tanh(y) + U_m 
$$

Here $a = U_m / (\Delta U)$ is the convection parameter. $a > 1$
corresponds to a shear layer developing in the $+x$ direction while $|a|<1$ 
corresponds to a situation with a counter-flow in the $-x$ direction in the
region $y < 0$.



We look for solutions under eigenmode form :
$$
[u,v,p] = [\hat u(y), \hat v(y), \hat p(y)] e^{i k x} e^{-\omega t}
$$ 



We have seen that differentiation with respect to $x$ ammounts to multiplication by $i k$, and differentiation 
with respect to $t$ ammounts to multiplication by $-i\omega$, thus we have 
$$
\begin{array}{l}
-i \omega \hat{u}=  -i k U \hat{u} - U_y \hat{v} - \hat{p}+\mu (-k^2 \hat{u}+\hat{u}_{yy})\\
- i \omega  \hat{v}=-i k U \hat{u} -\hat{p}_y+\mu(-k^2 \hat{v}+\hat{v}_{yy}) \\
i k \hat{u}+\hat{v}_y=0\\
\end{array}
$$

When considering *spatial stability* formalism, this problem is not a standard linear eigenvalue problem because the eigenvalue $k$ appears
 quadratically (in the viscous terms).
To transform the problem into a regular eigenvalue one, we have to
introduce the auxiliary unknowns $\hat{u}_1 = k \hat{u}$ and $\hat{v}_1 = k \hat{v}$.

We can thus write the problem as follows:

$$ k B q =  A q$$
with 
$$
q = \left[ \begin{array}{c} \hat{u} \\ \hat{u}_1 \\  \hat{v} \\ \hat{v}_1
\\ \hat{p} \end{array} \right],
$$

$$
A = \left[ \begin{array}{ccccc} 
0 & 1 & 0 & 0 & 0 \\
i \omega + Re^{-1} \partial_y^2  & 0 & - \partial_y \overline{U} & 0 & 0 \\
0 & 0 & 0 & 1 & 0 \\
0 & 0 & i \omega + Re^{-1} \partial_y^2 & 0 & - \partial_y  \\
0 & i & \partial_y & 0 & 0 
 \end{array} \right]
$$

$$ 
B = \left[ \begin{array}{ccccc} 
1 & 0 & 0 & 0 & 0 \\
i \overline{U}  & Re^{-1}  & 0 &  0 & 0 \\
1 & 0 & 0 & 0 & 0 \\
0 & 0 & i \overline{U}  & Re^{-1} & 0 \\
0 & 0 & 0 & 0 & 0 
 \end{array} \right]
$$




%}

function [] = main()
clear all; close all;
global y D DD w Z I discretization
% 'global' allows to use these objects in the function as well
 
% Physical parameters
omega=1;  %  the frequency
Re=200;    %  the Reynolds number
a = 2;     %  the coflow parameter 

% numerical parameters
N=90;      % the number of grid points
discretization = 'chebInfAlg'; % you may try Hermite as well 
iloop = 1;

%{

### Derivation matrices

Here we use Chebyshev discretization with stretching. See
[differential_equation_infinitedomain.m]() to see how this works.

%}


[D,DD,w,y] = dif1D(discretization,0,3,N,0.9999);
% [D,DD,w,y] = dif1D('cheb',-10,20,N);
 
Z=zeros(N,N); I=eye(N); 
Ndim=N;



%{
We compute the eigenvalues/eigenmodes with the function [KH](#Function KH), 
defined at the end of this program
%}

[sphys,Uphys,s,UU] = KH(omega,Re,a);

%sphys = s;Uphys = UU;
%sort = criterion>1e-3|imag(s)>-0.05;
%sphys(sort) = []; Uphys(:,sort)=[];

%{
### Plotting the spectrum
%}



figure(1);
plot(real(sphys),-imag(sphys),'o'),hold on;
grid on
figure(1);
for ind=1:length(s)
  h=plot(real(s(ind)),-imag(s(ind)),'*'); hold on
  set(h,'buttondownfcn',{@plotmode,UU(:,ind),s(ind),omega,Re,a});
end
ylim([-.5,.2]);
if ~isempty(sphys); xlim([min(real(sphys))-.1,max(real(sphys))+.1]); end
title(['Spatial spectrum for \omega = ',num2str(omega),'; Re= ',num2str(Re),'; a= ',num2str(a)]);
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','KH_spatial_spectrum.png');

%{
![** Figure : Temporal spectrum of a tanh shear layer](KH_spatial_viscous/KH_spatial_spectrum.png)
%} 


plotmode([],[],Uphys(:,1),sphys(1),omega,Re,a);

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','KH_spatial_mode.png');

%{
![** Figure : Temporal spectrum of a tanh shear layer](KH_spatial_viscous/KH_spatial_mode.png)
%} 


%{

# Loop over omega

%}
if iloop
    omegatab = [0.05:0.05:3];
    for j=1:length(omegatab)
        omega = omegatab(j)
        [k,U] = KH(omega,Re,a);
        k
        if ~isempty(k)
            k(length(k)+1:10,1)=0;
        else
            k = zeros(10,1);
        end
        ktab(:,j) = k(1:10);
    end
    figure(3);
    for j=1:10
        plot(omegatab,-imag(ktab(j,:))); hold on;
    end
    title(['Spatial growth rate for  Re= ',num2str(Re),'; a= ',num2str(a)]);
    set(gcf,'paperpositionmode','auto');
    print('-dpng','-r80','KH_spatial_spatialbranch.png');
end%if iloop
%{
![** Figure : Temporal spectrum of a tanh shear layer](KH_spatial_viscous/KH_spatial_branch.png)
%} 



end

%{

#Function KH

%}

function [sphys,Uphys,s,UU] = KH(omega,Re,a)
global y D DD w Z I discretization 
N = length(y);

% base flow
%U = (1-1./cosh(y))+a; Uy = D*U; % wake 
U=tanh(y)+a; Uy=(1-tanh(y).^2); % shear layer

% renaming the differentiation matrices
dy=D; dyy=DD;
Re1 = Re^-1;

%{
### System matrices
%}

% the matrices
S=1i*omega*I+Re1*dyy;
A = [    S,    -diag(Uy),   Z,      Z,      Z;  ...
         Z,      S,        -dy,     Z,      Z;  ...
         Z,     dy,         Z,      Z,      Z;  ...
         Z,     Z,          Z,      I,      Z;  ...
         Z,     Z,          Z,      Z,      I;
    ];

iU = 1i*diag(U);
B=[ 
    iU,     Z,      1i*I,   Re1*I,      Z;      ...
    Z,      iU,     Z,      Z,          Re1*I;  ...
    -1i*I,  Z,      Z,      Z,          Z;       ...
    I,      Z,      Z,      Z,          Z;      ...
    Z,      I,      Z,      Z,          Z       ... 
    ];

% Boundary conditions
if(strcmp(discretization,'her')==1)
    indBC = [];
    % No boundary conditions are required with Hermite interpolation which
    % naturally assumes that the function tends to 0 far away.
else
    % Dirichlet conditions are used for fd, cheb, etc...
 %   III=eye(5*N);
 %   indBC=[1,N,N+1,2*N];
 %   C=III(indBC,:);
 %   A(indBC,:)=C;
 %   B(indBC,:)=0;  
end

% computing eigenmodes 

[UU,S]=eig(A,B);

% sort the eigenvalues by decreasing real part and remove the spurious ones
s=diag(S);  [~,o]=sort(imag(s)); s=s(o); UU=UU(:,o);
rem=(abs(s)>2)|(abs(s)<1e-4); s(rem)=[]; UU(:,rem)=[];

modeS = UU(1:N,:); % U-component of modes
a0=fft([modeS; flipud(modeS(2:N-1,:))]);        
a0=a0(1:N,:).*[0.5; ones(N-2,1); 0.5]/(N-1);   % a0 contains Chebyshev coefficients 
criterion = sum(abs(a0(N*9/10:N,:)).^2)/sum(abs(a0(1:N)).^2); 
criterion = criterion';
           % this criterion is the 'energy' contained in the 10% coefficients of highest order.
           % this criterion is usually largest for spurious modes than for physical ones
sphys = s;Uphys = UU;
rem = criterion>0.01|imag(s)>0.05|abs(real(s))<3e-2;
sphys(rem) = []; Uphys(:,rem)=[];           


%if nargin==1
%    sphys = sphys(1);
%end

end %function KH

function [] = plotmode(~,~,mode,k,omega,Re,a)
    global y dy dyy
    Yrange = 4;
    figure(2);hold off;
    N = length(y);
    u = mode(1:N);
    v = mode(N+1:2*N);
    p = mode(2*N+1:3*N);
%    vorticity = (dy*u)-1i*k*v;
    subplot(1,3,1);hold off;
    plot(real(u),y,'b-',imag(u),y,'b--');hold on;
    plot(real(v),y,'g-',imag(v),y,'g--');hold on;
    plot(real(p),y,'k-',imag(p),y,'k--');hold on;
    ylabel('y'); ylim([-Yrange,Yrange]);
    legend({'$Re(\hat u)$','$Im(\hat u)$','$Re(\hat v)$','$Im(\hat v)$','$Re(\hat p)$','$Im(\hat p)$'},'Interpreter','latex')
    title('Structure of the eigenmode');
    % plot 2D reconstruction
    Lx=min([abs(2*pi/real(k)),20]); 
    Nx =30;  
    x=linspace(-Lx/2,Lx/2,Nx);
    p(abs(y)>Yrange,:)=[];
    u(abs(y)>Yrange,:)=[];
    v(abs(y)>Yrange,:)=[];
 %   vorticity(abs(y)>Yrange,:)=[];
    yy = y;
    yy(abs(y)>Yrange)=[];
    pp = 2*real(p*exp(1i*k*x));
    uu=2*real(u*exp(1i*k*x));
    vv=2*real(v*exp(1i*k*x));
  %  vorticityvorticity=2*real(vorticity*exp(1i*alpha*x));
    subplot(1,3,2:3); hold off;
    contourf(x,yy,uu,10); hold on; 
    quiver(x,yy,uu,vv,'k'); hold on;
    xlabel('x'); ylabel('y'); title({'Structure of the eigenmode (2D reconstruction)',['for omega = ',num2str(omega) , ' ; k = ',num2str(k)]});
end% function plotmode


%{

# Exercices/Contributions

* Please ...
%}