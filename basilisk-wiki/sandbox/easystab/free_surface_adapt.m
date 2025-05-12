%{

In [free_surface_2D.m]() we compute the linear oscillations of a free surface on top of a square of fluid, with a flat interface at rest. in [capillary_venturi_continuation.m] we compute the shape of the interface in a capillary bridge with throughflow with a 1D model for the fluid. Here we do something like that, but we want to get the nonlinear shape of the interface over a square, and vary the volume of fluid (so the interface will not be flat) and with a flow. The fluid dynamics is described using a potential flow.

The fact that the interface is not flat and in contact with a fluid domain makes that we need to adapt the shape of the fluid mesh to have as a boundary the free surface. This is done with [map2D.m]() and this is new here to have the mapping as a variable of the Newton iterations.
%}

clear all; clf;
disp('%%%%%%%%%')

% parameters
Lx=1;
Ly=1;
Nx=20;
Ny=20;
rho=1;
sig=1;
V=Lx*Ly*0.9; % the initial volume
U=0.1; % inflow velocity
method='fd';
delta=0.05; % continuation length

% differentiation
[d.x,d.xx,d.wx,x]=dif1D(method,0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D(method,0,Ly,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);
X0=X; Y0=Y; D0=D;

%{

# The unknowns

For the boundaries, we need to locate all the sides of the 2D fluid domain. We include the corners in the top and bottom. The variables are the velocity potential in the fluid domain
$$
\phi(x,y),
$$
the deformation of the interface
$$
\eta(x)
$$
and the reference pressure
$$
p_0,
$$
which is the pressure at a point where we nknow the velocity, for instance the bottom left corner (at this point, the horizontal velocity is imposed by the boundary conditions on $\phi$ and the vertical velocity is zero).

As unknowns, we have the velocity potential on the 2D grid, the values of the interface displacement at the top, and the reference pressure, thus

    Nx*Ny+Nx+2

unknown in our system.
%}
% useful matrices
l.top=[l.ctl; l.top; l.ctr]; % include the corners to the top
l.bot=[l.cbl; l.bot; l.cbr]; % include the corners to the bot
l.phi=(1:NN)';
l.e=l.phi(end)+(1:Nx)'; % where we store the interface deformation
l.p0=l.e(end)+1; % where we store the reference pressure

n=NN+Nx+1; % numbers of degrees of freedom
ZZ=zeros(n,n);    zx=zeros(Nx,Nx);
II=eye(n);        ix=eye(Nx);

%{
We build an initial gues that satisfies all the boundary conditions. Since the velocity components are
$$
u=\phi_x, v=\phi_y
$$
chosing the initial value of $\phi$ as below gives a uniform horizontal velocity $U$. The initial direction of continuation *dir* is just a reduction of the volume. if you want to increase the volution, put $-1$ instead for the last element of this vector.
%}
%initial guess
phi=U*X(:);
E=1+0*x;
e=0*x;
p0=0;
q0=[phi; e; p0];
q=q0; 

% Newton iterations
quit=0;count=0;
while ~quit
    
    %{
# Domain mapping

Here is one of the core element of this code: the fact that the Newton loop will adjust the shape $\eta(x)$ of the free surface and that the fluid domain will adjust accordingly such that the top grid cells fit with the free surface. For this we extract the actual guess for $\eta$ and we stretch vertically the domain using [map2D.m](). For this we have saved the differentiation matrices on the square domain *D0* and the original $y$ coordinate *Y0*.
    %}
    % the new domain
    e=q(l.e);
    E=E+e;
    %q(l.e)=0;
    Y=Y0.*repmat(E'/Ly,Ny,1);
    D=map2D(X,Y,D0);
    
    % the present solution and its derivatives
    D.xE=D.x(l.top,:);
    D.yE=D.y(l.top,:);
    D.xxE=D.xx(l.top,:);
    D.yyE=D.yy(l.top,:);
    D.xyE=D.xE*D.y;
    
    phi=q(l.phi); phixE=D.xE*phi; phiyE=D.yE*phi; phixxE=D.xxE*phi; phiyyE=D.yyE*phi; phixyE=D.xyE*phi;
    ex=d.x*e;  exx=d.xx*e;
    Ex=d.x*E; Exx=d.xx*E; a=1+Ex.^2;
    p0=q(l.p0);
    
    %{
# Nonlinear function and boundary conditions

Here is the nonlinear function that we would like to be zero. It starts with the kinematics of the velocity fluid in the mapped domain, it should satisfy the Poisson equation
$$
\phi_{xx}+\phi_{yy}=0.
$$
We then have the equation describing the pressure jump through the interface. The outer pressure is zero. The pressure jump through the curved interface is related to the surface tension $\sigma$ and the curvature of the interface, just like in [meniscus.m#the-shape-of-a-meniscus]()
$$
p_{in}=-\sigma \frac{\eta_{xx}}{(1+\eta_x^2)^{3/2}}
$$
This capillary pressure should be equal to the hydrodynamic pressure which we get from the velocity potential $\phi$ and the Bernoulli equation
$$
\rho \phi_t+\frac{\rho}{2}(\phi_x^2+\phi_y^2)+p=cste.
$$
In this equation, we have $\phi_t=0$ since we search a steady solution, and for the constant we consider the bottom left corner where the velocity is known because of the boundary condition
$$
cste=P_0+\frac{\rho}{2}U^2
$$
which gives
$$
p=P_0+\frac{\rho}{2}U^2-\frac{\rho}{2}(\phi_x^2+\phi_y^2)
$$
this thus gives the equation for $\eta$
$$
-\sigma \frac{\eta_{xx}}{(1+\eta_x^2)^{3/2}}-\left[P_0+\frac{\rho}{2}U^2-\frac{\rho}{2}(\phi_x^2|_\eta+\phi_y^2|_\eta)\right]=0.
$$
This should be imposed at every gridpoint of the free surface $\eta$, thus *Nx* equations. Please note that it is the pressure at the position of the interface that comes into play, this is why we have specified that the velocity field should be considered at the position $\eta$
$$
\phi_x^2|_\eta+\phi_y^2|_\eta
$$

We then have a scalar equation to impose the total volume of fluid, which constrains the values of $\eta$
$$
\int_0^{L_x} \eta dx-V=0
$$
where $V$ is the desired volume of fluid. This scalar equation is placed in the line of the system corresponding to the reference pressure $P_0$.

For the boundary conditions we impose no penetration at the left, right and bottom boundaries
$$
\begin{pmatrix}{l}
\phi_x|_{x=0}\\
\phi_x|_{x=L_x}\\
\phi_y|_{y=0}\\
\end{pmatrix}=0
$$
we then impose the the free surface is pinned at the top of the inflow and the top of the outflow
$$
\begin{pmatrix}{l}
\eta|_{x=0}\\
\eta|_{x=L_x}\\
\end{pmatrix}=0
$$
then we have a nonlinear boundary condition telling that there is no penetration of the fluid through the free-surface. The advection equation for a free surface by a velocity field $u=\phi_x, v=\phi_y$ is
$$
\eta_t+\phi_x|_\eta \eta_x=\phi_y|_\eta.
$$
Recall here that the free-surface is advected by the velocity at the psition of the free-surface, thus the position $|_\eta$. We look for a steady state so the time derivative is zero thus
$$
\phi_y|_\eta-\phi_x|_\eta \eta_x=0
$$
which should be satisfied at every point of the top boundary, thus *Nx* equations.

Our problem is thus fully determined by the two functions $f(q)$ and $c(q)$ as follows
$$
f(q)
=
f
\begin{pmatrix}
\phi\\ \eta \\ P_0
\end{pmatrix}
=
\begin{pmatrix}{l}
\phi_{xx}+\phi_{yy}\\
-\sigma \frac{\eta_{xx}}{(1+\eta_x^2)^{3/2}}-\left[P_0+\frac{\rho}{2}(U^2-(\phi_x^2|_\eta+\phi_y^2|_\eta))\right]\\
L_xL_y+\int_0^{L_x} \eta dx-V\\
\end{pmatrix}=0,
$$
$$
c(q)=
\begin{pmatrix}{l}
\phi_x|_{x=0}\\
\phi_x|_{x=L_x}\\
\phi_y|_{y=0}\\
\eta|_{x=0}\\
\eta|_{x=L_x}\\
\phi_y|_\eta-\phi_x|_\eta \eta_x\\
\end{pmatrix}=0
$$

# Unknown free-surface

We have a nonlinear equation on a domain that is unknown. We will first do a simplification of the unknown domain shape by saying that the true value of $\eta$ can be written $\eta=E+e$ and assuming that $e$ is small. This is the *flattening* of the boundary conditions at the free surface (see [pedagogy#systems-with-free-surfaces]() and [domain_derivative_1D_adapt.m]()). Here $E$ is the actual known guess for the free-surface position, and $e$ is a small deviation.
This gives us a nonlinear function on a known domain defined by the free-surface position $y=E(x)$. We solve this nonlinear function using Newto iterations. When a solution is found, with velocity potential $\phi$, volume $V$, pressure $P_0$ and interface displacement $e$, we adapt the domain by updating its shape $E_{n+1}=E_n+e$ and solve again the nonlinear function for a new velocity field, pressure, volume and interface displacement. We eventually converge to the full solution by iterating this way.

This is a *splitting* of the nonlinearity in two steps: the geometrical nonlinearity on one hand and the flow nonlinearity on the other hand.

Introducing the decomposition of $\eta\rightarrow E+e$ and linearizing with respect to $e$ we get

$$
f\begin{pmatrix}
\phi\\ e \\ P_0
\end{pmatrix}
=
\begin{pmatrix}{l}
\phi_{xx}+\phi_{yy}\\
-\sigma \left[\frac{E_{xx}}{(1+E_x^2)^{3/2}}+\frac{e_{xx}}{(1+E_x^2)^{3/2}}-\frac{3E_{xx}E_x e_x}{(1+E_x^2)^{5/2}}\right]
-\left[P_0+\frac{\rho}{2}(U^2-
(\phi_x^2|_E+\phi_y^2|_E+2(\phi_x|_E\phi_{xy}|_E+\phi_y|_E\phi_{yy}|_E)e))\right]\\
\int_0^{L_x} E dx +\int_0^{L_x} e dx-V\\
\end{pmatrix}=0,
$$

$$
c(q)=
\begin{pmatrix}{l}
\phi_x|_{x=0}\\
\phi_x|_{x=L_x}\\
\phi_y|_{y=0}\\
e|_{x=0}\\
e|_{x=L_x}\\
\phi_y|_E-\phi_x|_E E_x+(\phi_{yy}|_E-\phi_{xy}|_E E_x)e-\phi_x|_Ee_x \\
\end{pmatrix}=0
$$
    %}
    % nonlinear function
    fe=-sig*(Exx.*a.^-1.5+exx.*a.^-1.5-3*Exx.*Ex.*ex.*a.^-2.5) ...
        -(p0+rho/2*(U^2-( ...
        phixE.^2+phiyE.^2 ...
        +2*(phixE.*phixyE+phiyE.*phiyyE).*e))); ...
        
    f=[D.xx*phi+D.yy*phi; ... % potential flow
        fe; ...                    % pressure at interface
        d.wx*E+d.wx*e-V];  % volume
    
    %{
# Jacobian

We build now the jacobian of these functions, we have
$$
f(q+\tilde{q})\approx f(q)+A\tilde{q}
$$
with
$$
A\tilde{q}=\begin{pmatrix}
\partial_{xx}+\partial_{yy} & 0 & 0 \\
A_{\eta\phi} & A_{\eta\eta} & -1 \\
0 & \int_0^{L_x} dx & 0 \\
\end{pmatrix}
\begin{pmatrix}
\tilde{\phi}\\ \tilde{\eta}\\ \tilde{P_0}
\end{pmatrix}
$$
with
$$
A_{\eta\phi}=\rho\left[
(\phi_x|_E+e\phi_{xy}|_E)\partial_x|_E
+(\phi_y|_E+e\phi_{yy}|_E)\partial_y|_E
+e\phi_x\partial_{xy}|_E+e\phi_y|_E\partial_{yy}|_E
\right]
$$
and
$$
A_{\eta\eta}=-\sigma \left[ \frac{1}{(1+E_x^2)^{3/2}}\partial_{xx} - \frac{3E_{xx}E_x}{(1+E_x^2)^{5/2}}\partial_x\right]
+\rho(\phi_x|_E\phi_{xy}|_E+\phi_y|_E\phi_{yy}|_E)
$$

We do the same thing for the boundary conditions
$$
c(q+\tilde{q})\approx c(q)+C\tilde{q}
$$
with
$$
C\tilde{q}=\begin{pmatrix}{ccc}
\partial_x|_{x=0} & 0 & 0 \\
\partial_x|_{x=L_x} & 0 & 0 \\
\partial_y|_{y=0} & 0 & 0 \\
0 & I|_{x=0} & 0 \\
0 & I|_{x=L_x} & 0 \\
C_\phi & C_\eta & 0
\end{pmatrix}
\begin{pmatrix}
\tilde{\phi}\\ \tilde{\eta}\\ \tilde{P_0}
\end{pmatrix}
$$
with
$$
C_\phi=\partial_y|_E-(E_x+e_x)\partial_x|_E+e\partial_{yy}|_E-E_xe\partial_{xy}|_E
$$
$$
C_\eta=\phi_{yy}-\phi_{xy}E_x-\phi_x\partial_x
$$
        %}
        % analytical jacobian
        Aephi=rho*(diag(phixE+e.*phixyE)*D.xE ...
            +diag(phiyE+e.*phiyyE)*D.yE ...
            +diag(e.*phixE)*D.xyE ...
            +diag(e.*phiyE)*D.yyE);
        
        Aee=-sig*(diag(a.^-1.5)*d.xx-diag(3*Ex.*Exx.*a.^-2.5)*d.x) ...
            +rho*diag(phixE.*phixyE+phiyE.*phiyyE);
        
        A=[D.xx+D.yy, ZZ(l.phi,l.e), ZZ(l.phi,l.p0); ...
            Aephi, Aee, -ones(Nx,1); ...
            ZZ(l.p0,l.phi), d.wx, 0];
        
        % boundary conditions
        loc=[l.top;l.left;l.bot;l.right;l.e([1,Nx])];
        
        C=[D.x([l.left;l.right],:)*II(l.phi,:); ...
            D.y(l.bot,:)*II(l.phi,:); ...
            ix([1,Nx],:)*II(l.e,:)];
        
        cnl=phiyE-phixE.*Ex ...
            +(phiyyE-phixyE.*Ex).*e ...
            -phixE.*ex;
        
        Cnl=[D.yE-diag(Ex+ex)*D.xE+diag(e)*D.yyE-diag(Ex.*e)*D.xyE, ...
            diag(phiyyE-phixyE.*Ex)-diag(phixE)*d.x,zeros(Nx,1)];
        
        f(loc)=[C*(q-q0);cnl];
        A(loc,:)=[C; Cnl];
        
        % show solution
        plot(x,E,'b.-',x,E+e,'r.-',x,1+0*x,'k--');hold on
        u=D.x*phi; v=D.y*phi;
        surf(X,Y,reshape(p0-rho/2*(u.^2+v.^2),Ny,Nx)-10,'facealpha',0.2); shading interp;
        quiver(X(:),Y(:),u,v,'k');hold off
        axis([-0.1,Lx+0.1,-0.1,1.5*Ly]);
        title(V); drawnow
        
        % convergence test
        res=norm(f,inf);
        disp([num2str(count) '  ' num2str(res)]);
        if count>50||res>1e5||isnan(res); disp('no convergence'); quitcon=1; break; end
        if res<1e-7; quit=1; disp('converged'); continue; end
        
        % Newton step
        q=q-A\f;
        count=count+1;
end

% print figure
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','free_surface_nonlinear.png')

%{

The color map represents the pressure field:

![The figure](free_surface_nonlinear.png)

# Exercices/Contributions

* Please
* Please



%}