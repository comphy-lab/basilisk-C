%{

A code that looks very much like [jet_2D.m]() except for the fact that the top boundary is not a slip boundary. We build the inflow profile with a Blasius boundary layer profile to see how it behaves in 2D with the Navier-Stokes equations. The Blasius profile is computed in [blasiusf.m](), a function version of the code [blasius.m]().

See [jet_2D.m]() and then [venturi.m]() to learn how this works.

Dependency:

* [easypack.zip]()
* [blasiusf.m]()

%}

clear all; clf; format compact

Lyvec=25:25:150;
coco=jet(length(Lyvec));


for tre=1:length(Lyvec);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% parameters 
Ly=Lyvec(tre) % domain height
Re=50; % Inflow reynolds number
Nx=70; % number of grid nodes in x
Ny=3*Ly; %number of grid nodes in r
Lx=100; % domain length in x

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% preparing boundary conditions
NN=Nx*Ny;
l.u=(1:NN)'; l.v=l.u+NN; l.p=l.v+NN;
II=speye(3*NN);
D.lap=D.yy+D.xx;

neuploc=[l.ctl;l.ctr;l.ctr-Ny];  % where to impose the neumann condition on the pressure
p0loc=2*Ny; % where to impose zero pressure
dir=[l.left;l.top;l.bot]; % where to put Dirichley on u and w

loc=[l.u(dir); l.v(dir); l.p(p0loc); ...
    l.u(l.right); ...
    l.v(l.right); ...
    l.p(neuploc)];

C=[II([l.u([l.left;l.bot]);l.v(dir);l.p(p0loc)],:); ...     % Dirichlet on u,v,and p
   D.y(l.top,:)*II(l.u,:); ... % slip at top for u
   D.x(l.right,:)*II(l.u,:); ...   % Neuman on u at outflow
   D.x(l.right,:)*II(l.v,:); ...   % Neumann on v at outflow
   D.lap(neuploc,:)/Re*II(l.v,:)-D.x(neuploc,:)*II(l.p,:)]; % neuman constraint on pressure

% Compute Blasius velocity profile
[ybla,ubla]=blasiusf(Ly,100); %plot(ubla,ybla)

uprof=tanh(y/y(2));%plot(uprof,y,'.-');

% initial guess
V=zeros(NN,1);
U=uprof*ones(1,Nx);
P=-X/Re; P=P-P(p0loc); % pressure zero at p0loc
q0=[U(:);V(:);P(:)];

% Newton iterations
disp('Newton loop')
q=q0;
quit=0;count=0;
while ~quit     
 
    % the present solution and its derivatives
    U=q(l.u); V=q(l.v); P=q(l.p);
    Ux=D.x*U; Uy=D.y*U; 
    Vx=D.x*V; Vy=D.y*V;
    Px=D.x*P; Py=D.y*P;

    % nonlinear function
    f=[-U.*Ux-V.*Uy+D.lap*U/Re-Px; ...
       -U.*Vx-V.*Vy+D.lap*V/Re-Py; ...
       Ux+Vy];
    
    % Jacobian 
    A=[-(spd(U)*D.x+spd(Ux)+spd(V)*D.y)+D.lap/Re, -spd(Uy), -D.x; ...
       -spd(Vx),-(spd(U)*D.x+spd(V)*D.y+spd(Vy))+D.lap/Re, -D.y; ...
       D.x, D.y, Z];
     
    % Boundary conditions 
    f(loc)=C*(q-q0);
    A(loc,:)=C;

%     % plotting
%     subplot 311;
%     surf(X,Y,reshape(P-1,Ny,Nx)); view(2); shading interp; hold on
%     
%     sely=1:Ny; selx=1:6:Nx;
%     ww=reshape(U,Ny,Nx); vv=reshape(V,Ny,Nx); 
%     quiver(X(sely,selx),Y(sely,selx),ww(sely,selx),vv(sely,selx),'k');
%     axis([0,Lx,0,Ly]);
%     xlabel('z'); ylabel('r'); title('Pressure P'); grid off;hold off
%     
%     subplot 312;
%     surf(X,Y,reshape(U,Ny,Nx)); view(2); shading interp; 
%     xlabel('x'); ylabel('y'); title('horizontal velocity U'); grid off
%     
%     subplot 313;
%     surf(X,Y,reshape(V,Ny,Nx)); view(2); shading interp; 
%     xlabel('x'); ylabel('y'); title('vertical velocity V'); grid off
%     drawnow
    
    % convergence test  
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if count>50|res>1e5; disp('no convergence');break; end
    if res<1e-5; quit=1; disp('converged'); continue; end
    
    % Newton step
    q=q-A\f;   
    count=count+1;
end

% velocity profiles
%clf;
uu=reshape(U,Ny,Nx);
if tre==1; plot(ubla,ybla,'r'); hold on; end
for ind=10:10:Nx
plot(uu(:,ind),y/sqrt(x(ind)/Re),'color',coco(tre,:))
end
ylim([0,15])
grid on
drawnow

end

set(gcf,'paperpositionmode','auto')
%print('-dpng','-r80','boundary_layer_2D.png')
