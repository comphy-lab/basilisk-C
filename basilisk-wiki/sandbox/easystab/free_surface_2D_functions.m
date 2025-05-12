%{
# Free surface oscillations in 2D (with surface tension) 

Here is the version like [free_surface_2D.m](), but using [dif1D.m]() and [dif2D.m]().

%}


clear all; clf

% parameters
Nx=20; % gridpoints in x
Ny=20; % gridpoints in y 
Lx=1; % domain size in x
Ly=1; % domain size in y
pts=5; % number of points in finite difference stencils
sigma=1; % surface tension
rho=1; % fluid density

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,Lx,Nx);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,Ly,Ny);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% useful matrices
ZZ=zeros(NN+Nx,NN+Nx);    
zx=zeros(Nx,Nx);
II=eye(NN+Nx);            
ix=eye(Nx);
l.eta=((NN+1):(NN+Nx))';
l.phi=(1:NN)';

% system matrices
A=blkdiag(D.xx+D.yy,zx);
A(l.eta(2:end-1),l.phi)=D.y(l.top,:);
E=blkdiag(Z,ix);

% boundary conditions
loc=[l.top;l.left;l.bot;l.right;l.eta([1,2,Nx])];

c=[D.x([l.left;l.right],:)*II(l.phi,:); ...
   D.y(l.bot,:)*II(l.phi,:); ...
   d.x([1,Nx],:)*II(l.eta,:); ...
   d.wx*II(l.eta,:)];

Ca=[c; ...
    sigma*d.xx(2:end-1,:)*II(l.eta,:)];

Ce=[0*c; ...
    rho*II(l.top,:)]; 

A(loc,:)=Ca;
E(loc,:)=Ce;

% compute eigenmodes
disp('computing eigenmodes');
[U,S]=eig(full(A),full(E));
s=diag(S);  [t,o]=sort(abs(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
rem=imag(s)<0; s(rem)=[]; U(:,rem)=[];


% validation
subplot(1,4,1)
lambda=Lx./([1:15]/2); 
alpha=2*pi./lambda;

stheo=i*sqrt(sigma*alpha.^3.*tanh(alpha*Ly))';
plot(real(s),imag(s),'b.',real(stheo),imag(stheo),'ro');
axis([-1,1,-10,100]);
xlabel('real part of eigenvalue'); ylabel('imaginary part of eigenvalue');
title('spectra'); legend('numeric','theory','location','north')
grid on


% show velocity field and free surface
sel=[2,3,4];
for ind=1:length(sel);
    subplot(1,4,ind+1)
    
    % select eigenvector
    q=U(:,sel(ind)); 
    
    % extract u and v and reshape
    u=reshape(D.x*q(l.phi),Ny,Nx);
    v=reshape(D.y*q(l.phi),Ny,Nx);
    quiver(x,y,real(u),real(v),'b'); hold on 
    quiver(x,y,imag(u),imag(v),'r');
    
    % free surface
    e=q(l.eta);
    plot(x,real(e)+Ly,'b-',x,imag(e)+Ly,'r-')
    xlabel('x'); ylabel('y'); title(['mode ' num2str(sel(ind))]);
    axis equal; axis([0,Lx,0,2*Ly]);
    grid on


end
hold off

% print figure
set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','free_surface_2D.png');
