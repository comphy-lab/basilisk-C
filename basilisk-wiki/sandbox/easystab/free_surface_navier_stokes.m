%{

This is a viscous flow described with the Navier-Stokes equation. The top is a free surface that can deform accordingly to the pressure jump through the interface due to surface tension, and the normal and tangential stress in the flow. At the inflow we impose a profile that looks like a boundary layer.

%}

clear all; clf; format compact

%paramvec=linspace(5,1,7);
paramvec=linspace(1,0.85,10);
%paramvec=0.95

for tre=1:length(paramvec);
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    param=paramvec(tre)
    
    % parameters
    Re=500; % Reynolds number
    Nx=100;  % number of grid nodes in z
    Ny=50;  %number of grid nodes in r
    Lx=10;  % domain length
    Ly=1;   % domain height
    sig=5;  % surface tension
    V=Lx*Ly*param; % desired volume
    umax=1; % max inflow velocity
    mu=1/Re;
    
    % differentiation
    [d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
    [d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,5);
    [Dc,l,X,Y,Z,I,NN]=dif2D(d,x,y);
    
    % useful
    l.top=[l.ctl; l.top; l.ctr]; % include the corners to the top
    l.u=(1:NN)'; l.v=l.u+NN; l.p=l.v+NN; l.e=l.p(end)+(1:Nx)'; l.p0=l.e(end)+1; % location vectors
    l.ue=[l.u;l.e]; l.ve=[l.v;l.e]; l.pe=[l.p;l.e]; % compound location vectors
    
    p0loc=2*Ny; % where to impose pressure=p0
    t=l.top; % a shortcut because its often used
    
    n=3*NN+Nx+1; % numbers of degrees of freedom
    ZZ=spalloc(n,n,0);    zx=spalloc(Nx,Nx,0);
    II=speye(n);        ix=speye(Nx);
    Iu=II(l.u,:); Iv=II(l.v,:); Ip=II(l.p,:); Ie=II(l.e,:); Ip0=II(l.p0,:);
    T=kron(speye(Nx,Nx),ones(Ny,1)); % transformation operator from 1D x to 2D grid
    yy=Y(:);
    
    % initial guess
    v=zeros(NN,1);
    %u=umax*Y(:).*(2-Y(:));
    u=umax*tanh(Y/0.1);
    p=-mu*X*umax; p=p-p(p0loc);
    e=Ly+0*sin(pi*x);
    p0=0;
    q0=[u(:);v(:);p(:);e;p0];
    
    % Newton iterations
    disp('Newton loop'); tlu=0;
    if tre==1; q=q0; end
    quit=0;count=0;
    while ~quit
        tic
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % the present solution and its computational-space derivatives
        e=q(l.e); ex=d.x*e; exx=d.xx*e; E=T*e; Ex=T*ex; Exx=T*exx; a=1+ex.^2;
        u=q(l.u); v=q(l.v); p=q(l.p); p0=q(l.p0);
        ux=Dc.x*u; uy=Dc.y*u; uxy=Dc.y*ux; uyy=Dc.yy*u;
        vx=Dc.x*v; vy=Dc.y*v; vxy=Dc.y*vx; vyy=Dc.yy*v;
        px=Dc.x*p; py=Dc.y*p;
        
%{
# Physical space differentiation matrices

This is just like we did in [stretching_formula.m]().

%}
        % physical-space differentation matrices
        D.lap=Dc.xx ...
            +spd(-2*yy.*E.^-1.*Ex)*(Dc.x*Dc.y) ...
            +spd(yy.^2.*E.^-2.*Ex.^2+E.^-2)*Dc.yy ...
            +spd(yy.*(2*E.^-2.*Ex.^2-E.^-1.*Exx))*Dc.y;
        D.x=Dc.x+spd(-yy.*E.^-1.*Ex)*Dc.y;
        D.y=spd(E.^-1)*Dc.y;
        
%{

# Jacobians of the differentiation matrices

Since the domain deformation is one of the unknown of th system, the differentiation is no longer a linear operation (because it depends nonlinearly on the domain stretching $\eta$). We do just like we did in the test case [jacobian_formula.m](). 

%}
        % jacobians of the physical-space differentiation operators
        AA=-2*spd(yy)*T*(diag(e.^-1)*d.x-diag(e.^-2.*ex));
        BB=2*(spd(yy.^2)*T*(diag(e.^-2.*ex)*d.x-diag(e.^-3.*ex.^2))-T*diag(e.^-3));
        CC=spd(yy)*T*(diag(4*e.^-2.*ex)*d.x-diag(e.^-1)*d.xx+diag(-4*e.^-3.*ex.^2+e.^-2.*exx));
        DD=T*(diag(-e.^-1)*d.x+diag(e.^-2.*ex));
        EE=T*diag(-e.^-2);
        
        D.xue=spd(yy.*uy)*DD;  D.yue=spd(uy)*EE;  D.lapue=spd(uxy)*AA+spd(uyy)*BB+spd(uy)*CC;
        D.xve=spd(yy.*vy)*DD;  D.yve=spd(vy)*EE;  D.lapve=spd(vxy)*AA+spd(vyy)*BB+spd(vy)*CC;
        D.xpe=spd(yy.*py)*DD;  D.ype=spd(py)*EE;
        
        Iux=[D.x,D.xue]*II(l.ue,:);  Iuy=[D.y,D.yue]*II(l.ue,:);  Iulap=[D.lap,D.lapue]*II(l.ue,:);
        Ivx=[D.x,D.xve]*II(l.ve,:);  Ivy=[D.y,D.yve]*II(l.ve,:);  Ivlap=[D.lap,D.lapve]*II(l.ve,:);
        Ipx=[D.x,D.xpe]*II(l.pe,:);  Ipy=[D.y,D.ype]*II(l.pe,:);
        
        % physical-space derivatives
        ux=D.x*u; uy=D.y*u;
        vx=D.x*v; vy=D.y*v;
        px=D.x*p; py=D.y*p;
   
%{
# The nonlinear function

The first thing we do is to compute the expression of the normal and tangential stresses at the free surface. The stress tensor in a viscous Newtonian fluid is 
$$
\sigma=
\begin{pmatrix}
2\mu u_x-p & \mu(u_y+v_x) \\ 
\mu(u_y+v_x) & 2\mu v_y-p
\end{pmatrix}
$$
We are interested at the components of the stress along our free-surface. We first need to find the expresion of its normal and tangent vectors. We can say that the normal is along the gradient of a function for which the free-surface is an iso-line. For instance take
$$
F(x,y)=y-\eta(x)
$$
then its gradient, normalized to unit one is the normal to the interface
$$
\overrightarrow{n}=
\begin{pmatrix}
-\frac{\eta_x}{\sqrt{1+\eta_x^2}}\\
\frac{1}{\sqrt{1+\eta_x^2}}\\
\end{pmatrix}
$$
and the tangent vector is 
$$
\overrightarrow{t}=
\begin{pmatrix}
\frac{1}{\sqrt{1+\eta_x^2}}\\
\frac{\eta_x}{\sqrt{1+\eta_x^2}}\\
\end{pmatrix}
$$
Then, the normal stress is 
$$
\sigma_n=\overrightarrow{n}\sigma \overrightarrow{n}
=\frac{\eta_x^2(2\mu u_x-p)-2\mu\eta_x(u_y+v_x)+2\mu v_y-p}{1+\eta_x^2}
$$
and the tangential stress is
$$
\sigma_t=\overrightarrow{t}\sigma \overrightarrow{n}
=\mu\frac{(1-\eta_x^2)(u_y+v_x)+2\mu\eta_x(v_y-u_x)}{1+\eta_x^2}
$$
The conditions at the interface are that the normal stress has a jump equal to the pressure jump through a curved interface with surface tension, and the tangential stress is zero. 

We have then Navier-Stokes in the flow domain. And the last equation imposes the value of the fluid volume.
%}
        % nonlinear function
        stressn=-sig*(exx.*a.^-1.5)+(ex.^2.*(2*mu*ux(t)-p(t))-2*mu*ex.*(uy(t)+vx(t))+2*mu*vy(t)-p(t))./a; % normal stress
        stresst=mu*((1-ex.^2).*(uy(t)+vx(t))+2*ex.*(vy(t)-ux(t)))./a; % tangential stress
        
        f=[-(u.*ux+v.*uy)+D.lap*u/Re-px; ...
            -(u.*vx+v.*vy)+D.lap*v/Re-py; ...
            ux+vy; ...
            stressn; ... % jump in normal stress
            d.wx*Ly*e-V]; % volume
%{
# The Jacobian

We build here the expression of the Jacobian of $f$. 

%}
        % analytical Jacobian
        Astressn=-sig*(spd(a.^-1.5)*d.xx-spd(3*ex.*exx.*a.^-2.5)*d.x)*Ie ...
            +diag(2*mu*ex.^2./a)*Iux(t,:)-diag(2*mu*ex./a)*Iuy(t,:) ...
            -diag(2*mu*ex./a)*Ivx(t,:)+diag(2*mu./a)*Ivy(t,:) ...
            -Ip(t,:);
        
        Astresst=mu*(diag((1-ex.^2)./a)*Iuy(t,:)-diag(2*ex./a)*Iux(t,:) ...
            +diag((1-ex.^2)./a)*Ivx(t,:)+diag(2*ex./a)*Ivy(t,:));
        
        A=[-(spd(u)*Iux+spd(ux)*Iu+spd(v)*Iuy+spd(uy)*Iv)+Iulap/Re-Ipx; ...
            -(spd(u)*Ivx+spd(vx)*Iu+spd(v)*Ivy+spd(vy)*Iv)+Ivlap/Re-Ipy; ...
            Iux+Ivy; ...
            Astressn; ...
            d.wx*Ly*Ie];
%{
# Boundary conditions

They are: 

* Dirichlet for $u$ and $v$ at inflow (imposed inflow profile set by the initial guess) and bottom (no-slip)
* Neumann for $u$ and $v$ at outflow ($x$ derivative is zero) 
* The pressure must equal the reference pressure $p_0$ close to the bottom left corner
* The Neumann condition on pressure to avoid spurious pressure modes
* The height of the interface is imposed at $x=0$ and $x=L_x$.
* The free surface is a streamline, or said differently: no penetration through the interface:
$$
v-u\eta_x=0
$$
* The tangential stress is zero

%}
        % Boundary conditions
        dir=[l.u([l.left;l.bot]); l.v([l.left;l.bot])]; % where to put Dirichlet on u and v
        neuploc=[l.ctl;l.ctr;l.ctr-Ny];  % where to impose the neumann condition on the pressure
        
        loc=[dir; ...
            l.u(l.right); ...
            l.v(l.right); ...
            l.p(p0loc); ...
            l.p(neuploc); ...
            l.e([1,Nx]); ...
            l.u(t);
            l.v(t)];
        
        f(loc)=[q(dir)-q0(dir); ... % all Dirichlet
            ux(l.right); ... % Homogeneous Neuman on u at outflow
            vx(l.right); ...
            p(p0loc)-p0; % pressure equal to p0 at p0loc
            D.lap(neuploc,:)*v/Re-px(neuploc); ... % neumann constraint on pressure
            e(1)-Ly; ... % free-surface pinned at both ends
            e(Nx)-Ly; ...
            v(t)-u(t).*ex; ... % free-surface is a streamline
            stresst];
        
        C=[II(dir,:); ...     % all Dirichlet
            Iux(l.right,:); ...
            Ivx(l.right,:); ...
            Ip(p0loc,:)-Ip0;  % pressure equal to p0 at p0loc
            Ivlap(neuploc,:)/Re-Ipx(neuploc,:); ... % Neumann constraint on pressure
            ix(1,:)*Ie; ... % free-surface pinned at both ends
            ix(Nx,:)*Ie; ...
            Iv(t,:)-diag(ex)*Iu(t,:)-spd(u(t))*d.x*Ie; ... % free-surface is a streamline
            Astresst];
        
        A(loc,:)=C;
        tjac=toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % plotting
        Yp=Y.*reshape(E,Ny,Nx); % physical mesh
        ww=reshape(u,Ny,Nx); vv=reshape(v,Ny,Nx);
        vo=vx-uy; % vorticity
        surf(X,Yp,reshape(vo,Ny,Nx)-10,'facealpha',1); view(2); shading interp; hold on
        sely=1:2:Ny; selx=round(linspace(1,Nx,10));
        quiver(X(sely,selx),Yp(sely,selx),ww(sely,selx),vv(sely,selx),'k','MaxHeadSize',0.01);
        plot(x,Ly+0*x,'b--',Lx+y*0,y,'b--');
        axis([0,1.1*Lx,0,1.2*Ly]);
        xlabel('x'); ylabel('y'); title('Vorticity and velocity field'); grid off;hold off
        %caxis([-30,0]-10)
        drawnow
        
        % convergence test
        res=norm(f);
        disp([num2str(count) '  ' sprintf('%9.3g',res) '  tjac:' sprintf('%7.3g',tjac) '  tlu:' sprintf('%7.3g',tlu)]);
        if count>50||res>1e5; disp('no convergence');break; end
        if res<1e-5; quit=1; disp('converged'); continue; end
        
        % Newton step
        tic; q=q-A\f; tlu=toc;
        count=count+1;
    end
    
end
set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','free_surface_navier_stokes.png')

%{

Here is the final figure:

![The figure](free_surface_navier_stokes.png)

%}
