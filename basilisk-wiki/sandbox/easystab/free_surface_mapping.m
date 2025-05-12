%{

In this code, I test to have the domain mapping as one of the unknown of the nonlinear problem.
In this version, I use the analytical expression of the Jacobian (not so easy...). Please see [free_surface_mapping_numjac.m]() for a version with numerical Jacobian. 

If you would like to understand better the mapping, see [stretching_formula.m](), [jacobian_formula.m]() and domain_derivative_2D_mapping.m]().

Here I do in addition some Keller pseudoarclength continuation like in [meniscus_overturn_keller.m]() to follow the nonlinear branch. There is a nice fold of the branch, which corresponds to things we observe as well for a 1D model of the (axisymmetric) capillary Venturi in [capillary_venturi_continuation.m] (Another difference is that here I do the continuation on fluid volume instead of doing it on inflow velocity). 

%}

clear all; clf;
disp('%%%%%%%%%')

% parameters
Lx=1;
Ly=1;
Nx=40;
Ny=40;
rho=1;
sig=1;
V=Lx*Ly; % the desired volume
U=1.9; % inflow velocity
method='fd';
pts=5;
delta=0.1; % 

% computational-space differentiation
[d.x,d.xx,d.wx,x]=dif1D(method,0,Lx,Nx,pts);
[d.y,d.yy,d.wy,y]=dif1D(method,0,Ly,Ny,pts);
[Dc,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% useful
l.top=[l.ctl; l.top; l.ctr]; % include the corners to the top
l.bot=[l.bot;]; % include the corners to the bot
l.phi=(1:NN)'; % where we store the velocity potential
l.e=l.phi(end)+(1:Nx)'; % where we store the interface deformation
l.p0=l.e(end)+1; % where we store the reference pressure
l.V=l.p0+1;
l.phie=[l.phi; l.e];

n=NN+Nx+2; % numbers of degrees of freedom
ZZ=zeros(n,n);    zx=zeros(Nx,Nx);
II=eye(n);        ix=eye(Nx);
T=kron(speye(Nx,Nx),ones(Ny,1)); % transformation operator from 1D to 2D
yy=Y(:);

%initial guess
phi=U*X(:);
e=1+0*x;
p0=0;
q0=[phi; e; p0; V];
dir=[zeros(n-1,1);-1];   % initial direction
q=q0;

disp('Continuation loop')
ind=0;quitcon=0;
while ~quitcon&ind<200
    qprev=q;
    q=q+dir*delta; % new prediction of solution
    
    % Newton iterations
    quit=0;count=0;
    while ~quit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % the present solution and its computational_space derivatives
        phi=q(l.phi); phicy=Dc.y*phi; phicxy=Dc.x*phicy;  phicxx=Dc.xx*phi; phicyy=Dc.yy*phi;
        e=q(l.e); ex=d.x*e; exx=d.xx*e; E=T*e; Ex=T*ex; Exx=T*exx; a=(1+ex.^2);
        p0=q(l.p0);
        V=q(l.V);
        
%{

# Physical-space differentiation matrices

Just as in [stretching_formula.m]().
%}
        % physical-space differentation matrices
        D.lap=Dc.xx ...
            +spd(-2*yy.*E.^-1.*Ex)*(Dc.x*Dc.y) ...
            +spd(yy.^2.*E.^-2.*Ex.^2+E.^-2)*Dc.yy ...
            +spd(yy.*(2*E.^-2.*Ex.^2-E.^-1.*Exx))*Dc.y;
        D.x=Dc.x+spd(-yy.*E.^-1.*Ex)*Dc.y;
        D.y=spd(E.^-1)*Dc.y;
        
%{

# Jacobians of the differentiations

Differentiation is now a nonlinear operator! Because the stretching function is one of the unknowns...

%}
        % jacobians of the physical-space differentiation operators
        D.lapphie=-2*spd(phicxy.*yy)*T*(diag(e.^-1)*d.x-diag(e.^-2.*ex)) ...
            +spd(phicyy)*2*(spd(yy.^2)*T*(diag(e.^-2.*ex)*d.x-diag(e.^-3.*ex.^2))-T*diag(e.^-3)) ...
            +spd(phicy.*yy)*T*(diag(4*e.^-2.*ex)*d.x-diag(e.^-1)*d.xx+diag(-4*e.^-3.*ex.^2+e.^-2.*exx));
        D.xphie=spd(yy.*phicy)*T*(diag(-e.^-1)*d.x+diag(e.^-2.*ex));
        D.yphie=spd(phicy)*T*diag(-e.^-2);
        
        D.lapj=[D.lap,D.lapphie];
        D.xj=[D.x,D.xphie];
        D.yj=[D.y,D.yphie];
        
        % physical-space derivatives
        philap=D.lap*phi; phix=D.x*phi; phiy=D.y*phi;
        phixE=phix(l.top); phiyE=phiy(l.top);
        
        % nonlinear function
        f=[philap; ... % potential flow
            -sig*(exx.*a.^-1.5)-(p0+rho/2*(U^2-(phixE.^2+phiyE.^2))); ... % pressure at interface
            d.wx*Ly*e-V; ...
            dir'*(q-qprev)-delta];  % volume
        
%{
# The Jacobian of the system
%}
        % Jacobian
        Ae=-sig*(spd(a.^-1.5)*d.xx-spd(3*ex.*exx.*a.^-2.5)*d.x)*II(l.e,:) ...
            +rho*(spd(phixE)*D.xj(l.top,:)+spd(phiyE)*D.yj(l.top,:))*II(l.phie,:) ...
            -ones(Nx,1)*II(l.p0,:);
        
        A=[D.lapj*II(l.phie,:); ...
            Ae; ...
            d.wx*Ly*II(l.e,:)-II(l.V,:);
            dir'];
        
%{
# Boundary conditions

Almost all the boundary conditions are now nonlinear since most of them involve spatial derivatives on the 2D grid.

%}
        % boundary conditions
        neux=[l.left; l.right];
        neuy=l.bot; phi0loc=2+Ny;
        loc=[l.phi(phi0loc);neux;neuy;l.top;l.e([1,Nx])];
        
        C=[II(l.phi(phi0loc),:); ...
            D.xj(neux,:)*II(l.phie,:); ... % u at left and right
            D.yj(neuy,:)*II(l.phie,:); ... % v at bottom
            ix([1,Nx],:)*II(l.e,:); ... % free-surface pinned at start and end
            (D.yj(l.top,:)-spd(ex)*D.xj(l.top,:))*II(l.phie,:)-spd(phixE)*d.x*II(l.e,:)]; % interface is a streamline
        
        f(loc)=[phi(phi0loc); ...
            phix(neux)-U*(neux*0+1); ...
            phiy(neuy); ...
            e([1,Nx])-[1;1]; ...
            phiyE-phixE.*ex];
        A(loc,:)=C;
        
        % show solution
        subplot 121
        Yp=Y.*reshape(E,Ny,Nx); % physical mesh
        plot(x,1+0*x,'k--');hold on
        p=reshape(p0-rho/2*(phix.^2+phiy.^2),Ny,Nx);
        surf(X,Yp,p,'facealpha',0.4); shading interp;
        selx=(1:5:Nx); sely=1:1:Ny; u=reshape(phix,Ny,Nx);v=reshape(phiy,Ny,Nx);
        quiver(X(sely,selx),Yp(sely,selx),u(sely,selx),v(sely,selx),'k');hold off
        axis([-0.1,Lx+0.1,-0.1,1.5*Ly]);
        xlabel('x'); ylabel('y'); title('free surface and potential flow');
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % show nonlinear branch
    subplot 122
    plot(p0,V,'b.'); 
    xlabel('pressure'); ylabel('Volume'); title('bifurcation diagram');
    grid on; drawnow;hold on
    
    % New direction
    dir=A\[zeros(n-1,1);1]; % new direction
    dir=dir/norm(dir); % normalization
    
    % adjust the continuation jump
    if count<3; delta=1.2*delta; elseif count>3; delta=delta/1.2; end
    disp([num2str(ind) ' %%%%%%% Continuation V=' num2str(V) ' delta=' num2str(delta)]); 
    if ind==15; disp('computation finished'); break; end
    ind=ind+1;
end

% print figure
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','free_surface_mapping.png')

%{

The figure showing the deformed domain with velocity field and pressure as color map. We have allready passed the fold as you can see on the right subplot, so the volume after reaching a minimum is increasing again, and the curvature of the interface is starting to localize in the center.

![The figure](free_surface_mapping.png)
%}
