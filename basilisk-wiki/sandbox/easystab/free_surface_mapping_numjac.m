%{

In this code, I test to have the domain mapping as one of the unknown of the nonlinear problem.
In this version, I use the numerical Jacobian. To see the use of the analytical Jacobian, please check out [free_surface_mapping.m]().

%}

clear all; clf;
disp('%%%%%%%%%')

% parameters
Lx=1;
Ly=1;
Nx=15;
Ny=15;
rho=1;
sig=1;
V=Lx*Ly*0.9; % the initial volume
U=2.01; % inflow velocity
method='fd';
delta=0.1; % continuation length
eps=1e-8; % for the numerical Jacobian

% differentiation
[d.x,d.xx,d.wx,x]=dif1D(method,0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D(method,0,Ly,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% useful
l.top=[l.ctl; l.top; l.ctr]; % include the corners to the top
l.bot=[l.cbl; l.bot; l.cbr]; % include the corners to the bot
l.phi=(1:NN)';
l.e=l.phi(end)+(1:Nx)'; % where we store the interface deformation
l.p0=l.e(end)+1; % where we store the reference pressure

n=NN+Nx+1; % numbers of degrees of freedom
ZZ=zeros(n,n);    zx=zeros(Nx,Nx);
II=eye(n);        ix=eye(Nx);

%initial guess
phi=U*X(:);
e=1+0*x;
p0=0;
q0=[phi; e; p0];
q=q0;


% Newton iterations
quit=0;count=0;
while ~quit
    
    
%{
# Numerical Jacobian

In most of the codes, I build the analytical Jacobian, building it using the differentiation matrices and the formula that I derived by hand. Here instead I build it numerically by finite differences.

%}
    % nonlinear function and numerical Jacobian
    A=zeros(n,n);
    for gre=0:n
        
        qq=q+((1:n)==(gre))'*eps;
        
        % variables and their derivatives
        e=qq(l.e);
        ex=d.x*e; exx=d.xx*e;
        a=1+ex.^2;
        E=ones(Ny,1)*e';E=E(:);
        Ex=ones(Ny,1)*ex'; Ex=Ex(:);
        Exx=ones(Ny,1)*exx'; Exx=Exx(:);
        
        Dp.x=D.x-diag(Y(:).*E.^-1.*Ex)*D.y;
        Dp.y=diag(E.^-1)*D.y;
        Dp.lap=D.xx ...
            +diag(-2*Y(:).*E.^-1.*Ex)*(D.x*D.y) ...
            +diag(Y(:).^2.*E.^-2.*Ex.^2+E.^-2)*D.yy ...
            +diag(Y(:).*(2*E.^-2.*Ex.^2-E.^-1.*Exx))*D.y;
        
        % the present solution and its derivatives
        Dp.xE=Dp.x(l.top,:);
        Dp.yE=Dp.y(l.top,:);
        Dp.lapE=Dp.lap(l.top,:);
        Dp.xyE=Dp.xE*D.y;
        
        phi=qq(l.phi); phixE=Dp.xE*phi; phiyE=Dp.yE*phi;
        p0=qq(l.p0);
%{
# Nonlinear function

%}
        % nonlinear function
        feps=[Dp.lap*phi; ... % potential flow
            -sig*(exx.*a.^-1.5)-(p0+rho/2*(U^2-(phixE.^2+phiyE.^2))); ... % pressure at interface
            d.wx*Ly*e-V];  % volume
        
        % boundary conditions
        loc=[l.top;l.left;l.bot;l.right;l.e([1,Nx])];
        
        C=[D.x([l.left;l.right],:)*II(l.phi,:); ...
            D.y(l.bot,:)*II(l.phi,:); ...
            ix([1,Nx],:)*II(l.e,:)];
        
        feps(loc)=[C*(qq-q0); ...
            phiyE-phixE.*ex];
        
        % store Jacobian
        if gre==0; f=feps; else
            A(:,gre)=(feps-f)/eps;
        end
    end
    
    % show solution
    Yp=Y.*reshape(E,Ny,Nx); % physical mesh
    plot(x,1+0*x,'k--');hold on
    u=Dp.x*phi; v=Dp.y*phi;
    surf(X,Yp,reshape(p0-rho/2*(u.^2+v.^2),Ny,Nx)-10,'facealpha',0.2); shading interp;
    quiver(X(:),Yp(:),u,v,'k');hold off
    axis([-0.1,Lx+0.1,-0.1,1.5*Ly]);
    title(V);
    drawnow
    
    % convergence test
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if count>50||res>1e5||isnan(res); disp('no convergence'); quitcon=1; break; end
    if res<1e-7; quit=1; disp('converged'); continue; end
    
    % Newton step
    q=q-A\f;
    count=count+1;
end


% print figure
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','free_surface_mapping.png')

%{

The figure showing the deformed domain with velocity field and pressure as color map. We are very close to the fold of the bifurcation branch, we see this by the curvature of the interface.

![The figure](free_surface_mapping.png)
%}