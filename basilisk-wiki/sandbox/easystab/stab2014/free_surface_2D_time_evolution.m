%{
# Time evolution of free surface oscillations in 2D (with surface tension) 

Here we use the code [free_surface_2D.m](/sandbox/easystab/free_surface_2D.m) and using the exponential in time, like in [free_surface_gravity_particles.m](/sandbox/easystab/free_surface_gravity_particles.m), to show the time evolution of the eigenmodes.

Dependency:

* [chebdif.m](/sandbox/easystab/free_surface_2D.m) for the Chebychev differentiation matrices

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

%1D differentiation matrices
scale=-2/Lx;
[x,DM] = chebdif(Nx,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=x/scale; x=x-x(1); 
intx=([diff(x)',0]+[0,diff(x)'])/2;

scale=2/Ly;
[y,DM] = chebdif(Ny,2); 
dy=DM(:,:,1)*scale;    
dyy=DM(:,:,2)*scale^2;    
y=y/scale; y=y-y(end); 
inty=([diff(y)',0]+[0,diff(y)'])/2;

% 2D differentiation matrices
Dx=kron(dx,eye(Ny));
Dxx=kron(dxx,eye(Ny));
Dy=kron(eye(Nx),dy);
Dyy=kron(eye(Nx),dyy);
[X,Y]=meshgrid(x,y);


% vectors for coordinate selections 
NN=Nx*Ny;
dom=reshape(1:NN,Ny,Nx);
top=dom(1,2:end-1); top=top(:);
bot=dom(end,2:end-1); bot=bot(:);
left=dom(:,1); left=left(:);
right=dom(:,end); right=right(:);
phi=(1:NN)';
eta=((NN+1):(NN+Nx))';

% useful matrices
Z=zeros(NN,NN); ZZ=zeros(NN+Nx,NN+Nx);    zx=zeros(Nx,Nx);
I=eye(NN);      II=eye(NN+Nx);            ix=eye(Nx);

% system matrices
A=blkdiag(Dxx+Dyy,zx);
A(eta(2:end-1),phi)=Dy(top,:);
E=blkdiag(Z,ix);

% boundary conditions
loc=[top;left;bot;right;eta([1,2,Nx])];

c=[Dx([left;right],:)*II(phi,:); ...
   Dy(bot,:)*II(phi,:); ...
   dx([1,Nx],:)*II(eta,:); ...
   %II(eta([1,Nx],:),:);...
   intx*II(eta,:)];

Ca=[c; ...
    sigma*dxx(2:end-1,:)*II(eta,:)];

Ce=[0*c; ...
    rho*II(top,:)]; 

A(loc,:)=Ca;
E(loc,:)=Ce;

% compute eigenmodes
disp('computing eigenmodes');
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(abs(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
rem=imag(s)<0; s(rem)=[]; U(:,rem)=[];

%{
# Validation

Just as for the example in [wave_like.m#validation](), we have chosen boundary conditions on $\eta$ so that by chance our numerical solution should be the same as a wavelike solution periodic in $x$. We do here a theory similar to that of [free_surface_gravity.m#validation]() but with surface tension (pressure jump through the interface) instead of gravity (pressure gradient in the fluid). And also as a difference we use velocity potential instead of stream-function (for a change...).

%}

% validation
figure(2)
subplot(1,4,1)
lambda=Lx./([1:15]/2); 
alpha=2*pi./lambda;

stheo=i*sqrt(sigma*alpha.^3.*tanh(alpha*Ly))';
plot(real(s),imag(s),'b.',real(stheo),imag(stheo),'ro');
axis([-1,1,-10,100]);
xlabel('real part of eigenvalue'); ylabel('imaginary part of eigenvalue');
title('spectra'); legend('numeric','theory','location','north')
grid on

%{
# Velocity fields

It is nice to see what the flow looks like for the different modes of oscillation of our fluid/interface system. We need to get the velocity field and also the connection with the interface. We plot both the real part (in blue) and the imaginary part (in red). 

You can recognize directly by eye that these configurations corresponds to stable oscillations because at the top of the domain, the velocity field is in the direction opposite to the deflection of the free surface: the flow pushes back the surface toward rest (but it will overshoot and go to the other side, so the fluid will push back again and so on). If this was unstable, the flow would instead push the interface away from its rest position.

%}

% show velocity field and free surface
sel=[2,3,4];
for ind=1:length(sel);
    subplot(1,4,ind+1)
    
    % select eigenvector
    q=U(:,sel(ind)); 
    
    % extract u and v and reshape
    u=reshape(Dx*q(phi),Ny,Nx);
    v=reshape(Dy*q(phi),Ny,Nx);
    quiver(x,y,real(u),real(v),'b'); hold on 
    quiver(x,y,imag(u),imag(v),'r');
    
    % free surface
    e=q(eta);
    plot(x,real(e)+Ly,'b-',x,imag(e)+Ly,'r-')
    xlabel('x'); ylabel('y'); title(['mode ' num2str(sel(ind))]);
    axis equal; axis([0,Lx,0,2*Ly]);
    grid on


end
hold off

%{
![The spectrum and the velocity field of the first three eigenmodes](/sandbox/easystab/free_surface_2D.png)

From this figure, you can see the validation of the code as well as the velocity field of the fluid. We can imagine the real motion of the fluid and the interface from the velocity vectors here, but it will be more interesting to see the real motion of the interface with the time evolution which is what we will do in the next part.

%}


nper=1;      % number of periods of oscillation
nt=100;       % number of time steps per period
alpha=2*pi/Lx; %wave number
modesel=3;  % which mode to animate

%select the eigenmode
top=dom(1,1:end); top=top(:);
q=U(:,modesel); 
lambda=s(modesel);

%scale mode amplitude 
q=0.05*q/abs(q(eta(1)));

tvec=linspace(0,nper*2*pi/abs(lambda),nper*nt); 
dt=tvec(2)-tvec(1); 
x=linspace(0,Lx,Nx);
qq=2*real(exp(lambda*tvec(1))*q);

% initialize tracer particles
[X,Y]=meshgrid(x,y);
[px,py]=meshgrid(linspace(0,Lx,Nx),linspace(0,Ly,Ny)); 
py=py.*(1+kron(ones(Ny,1),(2*real(exp(lambda*tvec(1))*q(eta))/Ly)'));
px=px(:);py=py(:);

figure(1)
filename = 'free_surface_time_evolution.gif';

for ind=1:nper*nt
    % expand mode to physical space 
    qq=2*real(exp(lambda*tvec(ind))*q);
    
    
    % plot free surface
    plot(x,Ly+qq(eta,:),'b-',x,0*x+Ly,'k--'); hold on
    
    plot(mod(px,Lx),py,'k.');hold on
    
    [px,py]=meshgrid(linspace(0,Lx,Nx),linspace(0,Ly,Ny)); 
    py=py.*(1+kron(ones(Ny,1),(2*real(exp(lambda*tvec(ind))*q(eta))/Ly)'));
    px=px(:);py=py(:);
    
    xlabel('x');    
    ylabel('y'); 
    axis equal; axis([0,Lx,0,1.3*Ly]); 
    grid off
    hold off
    drawnow

    frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if ind == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
%{
At this step, we can advect also particles using this time dependent velocity field to give a realistic impression of the fluid/interface motion. But we have a problem when using Taylor expansions to express the velocity and pressure at the position of the moving interface as a function of the values at the fixed position Ly.
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ongoing project to advect the particles%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     % compute velocity at position of particles
%     pu=interp2(X,Y,reshape(Dx*q(phi),[Ny,Nx]),px,py);
%     pv=interp2(X,Y,reshape(Dy*q(phi),[Ny,Nx]),px,py);
%     
%     %Taylor expansion
%     pu(py>Ly)=Dx*q(j)+q(eta,:)*Dy(Ny,:)*Dx*q(phi);
%     pv(py>Ly)=Dy*q(phi)+q(eta,:)*Dy(Ny,:)*Dy*q(phi);
%     
%     puu=2*real(exp(lambda*tvec(ind))*pu);
%     pvv=2*real(exp(lambda*tvec(ind))*pv);
%  
%     % advect particles
%     px=px+puu*dt;
%     py=py+pvv*dt;
   
    
end

%{

![Time evolution of free surface 2D (with surface tension)](/sandbox/easystab/stab2014/free_surface_time_evolution.gif)

From this animation, you can see the time evolution of the third mode of the interface of the free surface waves.
%}