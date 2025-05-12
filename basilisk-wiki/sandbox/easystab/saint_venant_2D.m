%{

%}
%%%% parameters 
Nx=40; % number of grid nodes in z
Ny=40; %number of grid nodes in r
Lx=15; % length in z of the domain [0,Lz]
Ly=15;
pts=5; % number of points in finite difference stencils
H=0.1; % depth of water
b=0.1; % viscosité
f=0.; % coriolis
g=1; % gravity
tmax=50;
dt=0.3;


% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('fp',-Lx/2,Lx,Nx,pts);
[d.y,d.yy,d.wy,y]=dif1D('fp',-Ly/2,Ly,Ny,pts);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);
II=blkdiag(I,I,I);
DDx=blkdiag(D.x,D.x,D.x);
DDy=blkdiag(D.y,D.y,D.y);

% system matrices
E=II;
A=[-b*I, f*I,-g*D.x; ...
    -f*I, -b*I, -g*D.y; ...
    -H*D.x, -H*D.y, Z];

% locations in the state
u=(1:NN)';
v=(NN+1:2*NN)';
eta=(2*NN+1:3*NN)';

% add corners in top and bot
l.top=[l.ctl; l.top; l.ctr];
l.bot=[l.cbl; l.bot; l.cbr];


% % boundary conditions
% dir=[u([l.left;l.right]);v([l.bot;l.top])];
% neux=eta([l.left;l.right]);
% neuy=eta([l.top;l.bot]);
% loc=[dir; neux; neuy ];
% 
% C=[II(dir,:); ...
%     DDx(neux,:); ...
%     DDy(neuy,:)];  
% 
% E(loc,:)=0; 
% A(loc,:)=C;
         
% march in time matrix 
M=(E-A*dt/2)\(E+A*dt/2);

% initial condition
eta0=exp(-(X-1).^2-(Y+0.5).^2);
q=[X(:)*0; X(:)*0; eta0(:)]; 

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

for ind=1:Nt    
    q=M*q; % one step forward
    
    % plotting
    subplot(1,2,1);
    quiver(X,Y,reshape(q(u),Ny,Nx),reshape(q(v),Ny,Nx)); 
    subplot(1,2,2);
    mesh(X,Y,reshape(q(eta),Ny,Nx)); %shading interp;
    axis([-Lx/2,Lx/2,-Ly/2,Ly/2,-1,1])
    caxis(0.2*[-1,1])
    xlabel('x'); ylabel('y');
    drawnow
end
%legend('position','velocity'); title('Vibrating string')