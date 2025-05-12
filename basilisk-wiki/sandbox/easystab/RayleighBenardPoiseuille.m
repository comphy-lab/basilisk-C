%{
# Mixed-convection instabilities in Rayleigh-Benard-Poiseuille flow

*This document belongs to the [easystab](http://basilisk.fr/sandbox/easystab/README) project, please consult the main page of the project for explanations.*

*The present program was specifically designed as a support for the course ["Introductions to hydrodynamical instabilities"](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) for the "DET" Master's cursus in Toulouse. For the underlying theory please see [Lecture notes for chapter 7](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md#lecture-7-shear-flow-instabilities-ii-viscous-instabilities)*

*To use this program: click "raw page source" in the left column, copy-paste in the Matlab/Octave editor, and run.*


This program was written for the exam of February 2022.

%}


function [] = RayleighBenard()

set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesFontSize',13)

global y dy dyy w Z I N

% parameters
N=101; % number of gridpoints

% differentiation matrices
[dy,dyy,w,y] = dif1D('cheb',0,1,N);
Z=zeros(N,N); I=eye(N); 


%{    
## Spectrum, and structure of the leading eigenmode.

We now take a value of Ra above the threshold and solve the eigenvalue
problem for these parameters :

%}
%%
global Ra Re Pr
Ra = 2000;
k = 2 ;
Pr = 10.;
Re = 20000;

%Re = 0;
%Ra = 2000;
%k = 3;
%Pr = 10;

[s,U] = RBP(k,Ra,Pr,Re);

%{
We plot the spectrum in figure 1
%}

if (Re < 1000)
% assuming scaling for thermal
 figure(1);hold off;
 for ind=1:length(s)
   h=plot(-imag(s(ind)),real(s(ind)),'ro'); hold on
   set(h,'buttondownfcn',{@plotmodeR,U(:,ind),s(ind),k});
 end
 xlabel('\lambda_i'); 
 if (Re ==0&&Ra>0) 
     xlim([-1 1]);
 end
 ylabel('\lambda_r'); ylim([-400, 50+0*ceil(max(real(s)))]);
 grid on;
else
    
% scaling for hydro
figure(1);hold off;
for ind=1:length(s)
  h=plot(-imag(s(ind)/(Pr*Re*k)),real(s(ind)/(Pr*Re)),'ro'); hold on
  set(h,'buttondownfcn',{@plotmodeR,U(:,ind),s(ind),k});
end
xlabel('-\lambda_i/(Re Pr k)');
ylabel('\lambda_r/(Re Pr)'); ylim([-.5, .1]); grid on;
end

%title({['Temporal spectrum for k = ',num2str(k), ', Ra = ',num2str(Ra),' , Re = ',num2str(Re)], 'Click on eigenvalues to see the eigenmodes'});
title(['Temporal spectrum for k = ',num2str(k), ', Ra = ',num2str(Ra),' , Re = ',num2str(Re)]);
saveas(gcf,['RBP_spectrum_Ra',num2str(Ra),'_Re',num2str(Re),'_k',num2str(k)],'png');
saveas(gcf,'RB_spectrum','svg'); 
%{ 
![Figure : Spectrum for Ra = 2000](RayleighBenardPoiseuille/RB_spectrum.svg)
%}
%%
%{
 
Reconstruct the 2D structure of the mode and plot

%}


%if strcmp(getenv('OCTAVE_AUTORUN'),'true') 
% to generate a figure for the most amplified mode when running on the server
  figure(2);
  plotmodeR(0,0,U(:,1),s(1),k)
%  print('-dpng','-r100','Mode.png');
saveas(gcf,['Mode_Ra',num2str(Ra),'_Re',num2str(Re),'_k',num2str(k)],'png')
saveas(gcf,'Mode','svg')
  

%{ 
![Figure : 2D reconstruction of the most amplified leading eigenmode for Ra = 2000 (colors are for **temperature perturbations**)](RayleighBenardPoiseuille/Mode.svg)
%}  
  
%% 
%{
  
  ## Loop over Re for Ra = 2000
  
  %}
   
  
if 0
  %% loop over Re for Ra = 2000
  k = 3; RRa = 2000;
  Re_tab = [0:.1:5];
  s_tab = [];
  for RRe = Re_tab
      s = RBP(k,RRa,Pr,RRe);
      s_tab = [s_tab,s];
  end
  figure(20);
  subplot(2,1,1); plot(Re_tab,real(s_tab));
  title(['Eigenvalue as function of Re for Ra = ',num2str(RRa),' ; k = ',num2str(k)]);
  xlabel('Re');ylabel('\lambda_r');grid on;
  subplot(2,1,2); plot(Re_tab,-1/k*imag(s_tab)./(Re_tab*Pr));
  xlabel('Re');ylabel('-\lambda_i/( k Re Pr)');grid on;
  saveas(gcf,'EV_Re','png');
end 
  
if 0
  %% loop over Ra for Re = 20000
  k = 2; RRe = 20000;
  Ra_tab = [-3000:150:3000];
  s_tab = [];
  for RRa = Ra_tab
      s = RBP(k,RRa,Pr,RRe);
      s_tab = [s_tab,s];
  end
  figure(20);
  subplot(2,1,1); plot(Ra_tab,real(s_tab));
  title(['Eigenvalue as function of Ra for Re = ',num2str(RRe),' ; k = ',num2str(k)]);
  xlabel('Ra');ylabel('\lambda_r');grid on;
  subplot(2,1,2); plot(Ra_tab,-1/k*imag(s_tab)./(RRe*Pr));
  xlabel('Ra');ylabel('-\lambda_i/( k Re Pr)');grid on;
  saveas(gcf,'EV_Ra','png');
end 
  

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

function [s,U] = RBP(k,Ra,Pr,Re)
global y dy dyy w Z I N

U0 = 4*y.*(1-y);
dU0 = 4-8*y;

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

A=[Pr*Delta-Re*Pr*1i*k*diag(U0), -Re*Pr*diag(dU0), -dx, Z; ...
   Z, Pr*Delta-Re*Pr*1i*k*diag(U0), -dy, Pr*I;  ...
   dx, dy, Z, Z;  ...
   Z, Ra*I, Z, Delta-Re*Pr*1i*k*diag(U0)];
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
s=diag(S) ;
[t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1e20; s(rem)=[]; U(:,rem)=[];

if(nargout==1) 
    s=(s(1)); 
    % this is a trick to allow to use the function as s=RT(k,Ra,Pr). 
    % This allows to pass the function as a handle for fzero.
end

end

function [] = plotmodeR(~,~,mode,lambda,k) % function to plot one mode
    global y dy dyy
    global Ra Pr Re
    figure(2);hold off;
    N = length(y);
   
    %% normalisation basée sur v(0,H/2) doit etre réelle
    mode = mode/mode(N+(N+1)/2);
    u = mode(1:N)/max(abs(mode(1:2*N)));
    v = mode(N+1:2*N)/max(abs(mode(1:2*N)));
    p = mode(2*N+1:3*N)/max(abs(mode(2*N+1:3*N)));
    T = mode(3*N+1:4*N)/max(abs(mode(3*N+1:4*N)));
    subplot(1,3,1);hold off;
    plot(real(u),y,'b-',imag(u),y,'b--');hold on;
    plot(real(v),y,'g-',imag(v),y,'g--');hold on;
    plot(real(p),y,'r-',imag(p),y,'r--');hold on;
    plot(real(T),y,'k-',imag(T),y,'k--');hold on;
    ylabel('y'); 
    legend({'Re(u)','Im(u)','Re(v)','Im(v)','Re(p)','Im(p)','Re(T)','Im(T)'})
  %  if (-imag(lambda)/k)<1&&(-imag(lambda)/k)>0
  %    ycrit = sqrt(1+imag(lambda)/k);
  %    plot([-1,1],[ycrit,ycrit],'k:',[-1,1],[-ycrit,-ycrit],'k:')
  %  end
    % plot 2D reconstruction
    Lx=2*pi/k; Nx =30;  x=linspace(-Lx/2,Lx/2,Nx);
    yy = y;
    pp = 2*real(p*exp(1i*k*x));
    uu=2*real(u*exp(1i*k*x));
    vv=2*real(v*exp(1i*k*x));
    TT=2*real(T*exp(1i*k*x));
    if (Ra==0)
        TT = 0*TT;
    end
    subplot(1,3,2:3); hold off;
    contourf(x,yy,TT,10); hold on; 
    quiver(x,yy,uu,vv,'k'); hold on;
    %axis equal;
    xlabel('x'); ylabel('y'); 
%    title({'Structure of the eigenmode ',...
%          ['lambda = ',num2str(lambda),  ', Ra = ',num2str(Ra),' , Re = ',num2str(Re), ', k = ',num2str(k)],...
%          'color : \theta ; vectors : (u,v)'  });
    title({[ '\lambda = ',num2str(lambda)],...
          [' Ra = ',num2str(Ra),' , Re = ',num2str(Re), ', k = ',num2str(k)]});
    saveas(gcf,['Mode_Ra',num2str(Ra),'_Re',num2str(Re),'_k',num2str(k)],'png')

end% function plotmodeR

