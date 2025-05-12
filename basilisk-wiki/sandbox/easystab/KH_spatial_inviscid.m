%{

This program solves the spatial stability problem for a tanh shear layer
in the inviscid case.

References : Huerre & Monkewitz (1985) ; Godr�che & Manneville (1998)

NOTE : this program has been validated for $R = 0.5$ (strongly convective
regime). The references predict C/A transition at R = 1.315 but up to now
I could not obtain reliable results in this range. Must implement better
criteria to sort physical/numerical modes and/or better discretization
methods.

## Equations

The base flow is 
$$
U(y) = 1 + R  \, tanh(y)  
$$


We look for solutions under eigenmode form :
$$
[u',v',p'] = [\hat u(y), \hat v(y), \hat p(y)] e^{i k x} e^{-\omega t}
$$ 



We have seen that differentiation with respect to $x$ ammounts to multiplication by $i k$, and differentiation 
with respect to $t$ ammounts to multiplication by $-i\omega$, thus we have 
$$
\begin{array}{l}
-i \omega \hat{u}=  -i k \bar{U} \hat{u} - \bar{U}_y \hat{v} - i k \hat{p}\\
- i \omega  \hat{v}=-i k \bar{U}\hat{u} -\hat{p}_y \\
i k \hat{u}+\hat{v}_y=0\\
\end{array}
$$

We can thus write the problem as follows:

$$ k B q =  A q$$
with 
$$
q = \left[ \begin{array}{c} \hat{u} \\  \hat{v} 
\\ \hat{p} \end{array} \right],
$$

$$
A = \left[ \begin{array}{ccccc} 
i \omega   & \overline{U}_y & O  \\
0 & i \omega  & - \partial_y  \\
0 & \partial_y & 0 
 \end{array} \right]
$$

$$ 
B = \left[ \begin{array}{ccccc} 
i \bar{U} & 0 & i \\
0 & i \bar{U} & 0 \\
- i & 0 & 0 
 \end{array} \right]
$$




%}

function [] = main()
clear all; close all;
global y D DD w Z I discretization
% 'global' allows to use these objects in the function as well
 
% Physical parameters
omega=0.3;  %  the frequency
R = 0.5;     %  the coflow parameter 

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



[s,UU] = KH(omega,R);

%sphys = s;Uphys = UU;
%sort = criterion>1e-3|imag(s)>-0.05;
%sphys(sort) = []; Uphys(:,sort)=[];

%{
### Plotting the spectrum
%}

Re = Inf;

figure(1);
%plot(real(ss),-imag(sphys),'o'),hold on;
grid on
figure(1);
for ind=1:length(s)
  h=plot(real(s(ind)),-imag(s(ind)),'*'); hold on
  set(h,'buttondownfcn',{@plotmode,UU(:,ind),s(ind),omega,Re,R});
end
%ylim([-.5,.2]);
%if ~isempty(sphys); xlim([min(real(sphys))-.1,max(real(sphys))+.1]); end
title(['Spatial spectrum for \omega = ',num2str(omega),'; Re= ',num2str(Re),'; R= ',num2str(R)]);
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','spectrum.png');

%{
![** Figure : Temporal spectrum of a tanh shear layer](KH_spatial_inviscid/spectrum.png)
%} 


plotmode([],[],UU(:,1),s(1),omega,Re,R);

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','spatial_mode.png');

%{
![** Figure : Temporal spectrum of a tanh shear layer](KH_spatial_inviscid/spatial_mode.png)
%} 


%{

# Loop over omega

%}
if iloop
    pause(1);
    omegatab = [0.05:0.025:2];
    for j=1:length(omegatab)
        omega = omegatab(j)
        [k,U] = KH(omega,R);
        k
        if ~isempty(k)
            k(length(k)+1:10,1)=0;
        else
            k = zeros(10,1);
        end
        ktab(:,j) = k(1:10);
    end
    figure(3);
    subplot(2,1,1);
    for j=1:1
        plot(omegatab,-imag(ktab(j,:)),'*'); hold on;
    end
        title(['Spatial growth rate for  Re= ',num2str(Re),'; R= ',num2str(R)]);
    subplot(2,1,2);
    for j=1:1
        plot(omegatab,real(ktab(j,:)),'*'); hold on;
    end
    title(['k_r for  Re= ',num2str(Re),'; R= ',num2str(R)]);
    set(gcf,'paperpositionmode','auto');
    print('-dpng','-r80','spatialbranch.png');
end%if iloop
%{
![** Figure : Temporal spectrum of a tanh shear layer](KH_spatial_inviscid/spatialbranch.png)
%} 



end

%{

#Function KH

%}

function [s,UU] = KH(omega,R)
global y D DD w Z I discretization 
N = length(y);

% base flow
%U = (1-1./cosh(y))+a; Uy = D*U; % wake 
U=1+R*tanh(y); Uy=R*(1-tanh(y).^2); % shear layer

% renaming the differentiation matrices
dy=D; dyy=DD;
%Re1 = Re^-1;

%{
### System matrices
%}

% the matrices
S=1i*omega*I;
A = [   S,    -diag(Uy),   Z ;       ...
       Z,      S,         +dy ;      ...
       Z,     -dy,         Z       ...
    ];

iU = 1i*diag(U);
B=[ 
   [ iU,     Z,      -1i*I] ;      ...
   [ Z,      iU,     Z] ;        ...
   [ 1i*I,  Z,      Z ]        ...
    ];

% Boundary conditions
if(strcmp(discretization,'her')==1)
    indBC = [];
    % No boundary conditions are required with Hermite interpolation which
    % naturally assumes that the function tends to 0 far away.
else
    % Dirichlet conditions are used for fd, cheb, etc...
    III=eye(3*N);
    indBC=[1,N,N+1,2*N];
    C=III(indBC,:);
    A(indBC,:)=C;
    B(indBC,:)=0;  
end

% computing eigenmodes 

[UU,S]=eig(A,B);

% sort the eigenvalues by decreasing real part and remove the spurious ones
s=diag(S);  [~,o]=sort(imag(s)); s=s(o); UU=UU(:,o);
rem=(abs(s)>5)|(abs(s)<1e-4)|abs(real(s))<1e-1; s(rem)=[]; UU(:,rem)=[];

%modeS = UU(1:N,:); % U-component of modes
%a0=fft([modeS; flipud(modeS(2:N-1,:))]);        
%a0=a0(1:N,:).*[0.5; ones(N-2,1); 0.5]/(N-1);   % a0 contains Chebyshev coefficients 
%criterion = sum(abs(a0(N*9/10:N,:)).^2)/sum(abs(a0(1:N)).^2); 
%criterion = criterion';
           % this criterion is the 'energy' contained in the 10% coefficients of highest order.
           % this criterion is usually largest for spurious modes than for physical ones
%sphys = s;Uphys = UU;
%rem = criterion>0.01|imag(s)>0.05|abs(real(s))<3e-2;
%sphys(rem) = []; Uphys(:,rem)=[];           


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