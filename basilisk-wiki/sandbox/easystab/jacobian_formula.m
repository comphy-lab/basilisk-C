%{

Validation: Comparing the expression of the analytical Jacobian and the numerical Jacobian for the functions on a mapped domain. We are only interested in the part of the Jacobian where the expression of the stretching plays a role.
This is a follow-up of the code [stretching_formula.m]() and to be used in [domain_derivative_mapping.m]().

%}

disp('%%%%%%%')
clear all; clf;

% parameters
Nx=10;           % number of grid points in x
Ny=10;            % number of grid points in y
Lx=1;      % the length in x of the computational domain
Ly=1;       % the length in x of the computational domain
eps=1e-8;   % Small perturbation for computing the numerical Jacobian
method='cheb'; % discretization method
pts=5; % number of points in the stencil in case of finite differences

% differentiation
[d.x,d.xx,d.wx,x]=dif1D(method,0,Lx,Nx,pts);
[d.y,d.yy,d.wy,y]=dif1D(method,0,Ly,Ny,pts);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% An operator to transform from the 1D x grid to the 2D grid
T=kron(speye(Nx,Nx),ones(Ny,1)); 
yy=Y(:);

l.phi=(1:NN)';
l.e=NN+[1:Nx]';
n=NN+Nx;

% base state
e=1+0.2*sin(pi*x/Lx);
phi=cos(pi*X).*cos(pi*Y);
q=[phi(:); e(:)];

% Computing Jacobian with respect to e
A=zeros(3*NN,Nx);
for gre=0:Nx
    
    % Perturbation of one degree of freedom
    if gre>0;
        qq=q+((1:n)==(l.e(gre)))'*eps;
    else
        qq=q;
    end
    
    % computational domain derivatives
    phi=qq(l.phi); phiy=D.y*phi;  phixx=D.xx*phi; phiyy=D.yy*phi; phixy=D.x*phiy;
    e=qq(l.e); ex=d.x*e; exx=d.xx*e; E=T*e; Ex=T*ex; Exx=T*exx;
    
%{
# The physical space differentiation matrices

Here instead of using the general function [map2D.m]() to get the differentiation matrices on the deformed domain, we build them from the analytical mapping as derived and tested in [stretching_formula.m#the-mapping](). 

%}
    % Physical space diferentiation matrices
    Dp.lap=D.xx ...
        +spd(-2*yy.*E.^-1.*Ex)*(D.x*D.y) ...
        +spd(yy.^2.*E.^-2.*Ex.^2+E.^-2)*D.yy ...
        +spd(yy.*(2*E.^-2.*Ex.^2-E.^-1.*Exx))*D.y;
    Dp.x=D.x+spd(-yy.*E.^-1.*Ex)*D.y;
    Dp.y=spd(E.^-1)*D.y;
    
    % Compute the physical space derivatives
    phiplap=Dp.lap*phi;
    phipx=Dp.x*phi;
    phipy=Dp.y*phi;
    
    % nonlinear function
    feps=[phiplap; phipx; phipy];
    
    % store Jacobian
    if gre==0;

        % The analytical expressions for the Jacobian
        f=feps;
        
%{

# The analytical expression of the Jacobian

When we just had to get the Jacobian of derivatives it was easy because these are linear operators, thus the Jacobian is the same as the operator itself. Now it is a bit more complicated because the expression of the derivative depends on the domain stretching $\eta$, and thus the derivatives are nonlinear functions of $\eta$^[cf jerome's notebook c107p130-131]. In the following we denote
$$
\partial_\bar{x}^{\eta}
$$
the physical-space differentiation operator based on a domain stretched by $\eta$. And we are interested in how this thing depends on $\eta$ for a small perturbation $\tilde{\eta}$.

# First derivative in $y$

We have
$$
\partial^\eta_\bar{y}\phi=\eta^{-1}\phi_y
$$
thus considering the perturbations of $\phi$ and $\eta$ we have
$$
\partial^{\eta+\tilde{\eta}}_\bar{y}(\phi+\tilde{\phi})
=(\eta+\tilde{\eta})^{-1}(\phi_y+\tilde{\phi}_y)\approx
\eta^{-1}\phi_y+\eta^{-1}\partial_y\tilde{\phi}-\eta^{-2}\phi_y\tilde{\eta}
$$
which we can write in matrix form as
$$
\partial^{\eta+\tilde{\eta}}_\bar{y}(\phi+\tilde{\phi})
\approx
\partial^\eta_\bar{y}\phi+
\begin{pmatrix} 
\partial^\eta_\bar{y} &-\eta^{-2}\phi_y\end{pmatrix}
\begin{pmatrix}
\tilde{\phi}\\
\tilde{\eta}
\end{pmatrix}
$$

# First derivative in $x$

We have
$$
\partial^\eta_\bar{x}\phi=\phi_x-y\eta^{-1}\eta_x\phi_y
$$
then considering the perturbations of $\phi$ and $\eta$ we have
$$
\partial^{\eta+\tilde{\eta}}_\bar{x}(\phi+\tilde{\phi})
\approx
\partial^{\eta}_\bar{x}\phi+
\begin{pmatrix} 
\partial^\eta_\bar{x} &-y\phi_y(\eta^{-1}\partial_x-\eta^{-2}\eta_x)
\end{pmatrix}
\begin{pmatrix}
\tilde{\phi}\\
\tilde{\eta}
\end{pmatrix}
$$


# For the laplacian

Now for the physical-space Laplacian 
$$
\bar{\Delta}^{\eta}=\partial^\eta_{\bar{xx}}\partial^\eta_{\bar{yy}}
$$
we have 
$$
\bar{\Delta}^{\eta+\tilde{\eta}}(\phi+\tilde{\phi})
=\bar{\Delta}_\eta \phi+ 
\begin{pmatrix}
\bar{\Delta}^{\eta} & \phi_{xx}A+\phi_{yy}B+\phi_{xx}C
\end{pmatrix}
\begin{pmatrix}
\tilde{\phi}\\
\tilde{\eta}
\end{pmatrix}
$$
with $A, B$ and $C$ the Jacobians of the functions of $\eta$
$$
\begin{array}{l}
a(\eta)=-2y\eta^{-1}\eta_x\\
b(\eta)=y^2\eta^{-2}\eta_x^2+\eta^{-2}\\
c(\eta)=y(2\eta^{-2}\eta_x^2-\eta^{-1}\eta_{xx}
\end{array}
$$
After some algebra, we find
$$\begin{array}{l}
A=-2y(\eta^{-1}\partial_x-\eta^{-2}\eta_x)\\
B=2y^2\eta^{-2}(\eta_x\partial_x-\eta^{-1}\eta_x^2)-2\eta^{-3}\\
C=y[ 4\eta^{-2} (\eta_x\partial_x-\eta^{-1}\eta^2_{x})-\eta^{-1}(\partial_{xx}-\eta^{-1}\eta_{xx})]
\end{array}
$$

Now the difficulty with these expressions is that $y$ is defined for everypoint of the grid, but $\eta$ is only define on a 1D grid, so we need to use an operator that multiplies the size *Nx* vector of $\eta$ to give the grid equivalent of $\eta$. This is the operator *T* built above.

%}
        AA=-2*spd(yy)*T*(diag(e.^-1)*d.x-diag(e.^-2.*ex));
        BB=2*(spd(yy.^2)*T*(diag(e.^-2.*ex)*d.x-diag(e.^-3.*ex.^2))-T*diag(e.^-3));
        CC=spd(yy)*T*(diag(4*e.^-2.*ex)*d.x-diag(e.^-1)*d.xx+diag(-4*e.^-3.*ex.^2+e.^-2.*exx));
        
        Dp.lapphie=spd(phixy)*AA+spd(phiyy)*BB+spd(phiy)*CC;
        Dp.xphie=spd(yy.*phiy)*T*(diag(-e.^-1)*d.x+diag(e.^-2.*ex));
        Dp.yphie=spd(phiy)*T*diag(-e.^-2);
        
        
    else
        % One column of the numerical Jacobian
        A(:,gre)=(feps-f)/eps;
    end
end

% Validation: compare numerical and analytical Jacobians
err_lap=norm(A(1:NN,:)-Dp.lapphie)
err_x=norm(A(NN+(1:NN),:)-Dp.xphie)
err_y=norm(A(2*NN+(1:NN),:)-Dp.yphie)

% Show the function phi and e
Yp=reshape(E.*Y(:),Ny,Nx);
mesh(X,Yp,reshape(phi,Ny,Nx));
hold on
plot(x,e,'k-',[0 Lx],[Ly Ly],'k--');
view(103,66);
xlabel('x');ylabel('y');zlabel('h');




set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','jacobian_formula.png')


%{

This gives the reassuring screen output

    err_lap =
       9.5728e-04
    err_x =
       9.9049e-06
    err_y =
       1.0995e-06

# The figure

The figure just shows the reference state $\phi$ and $\eta$ used to test the analytical expressions of the Jacobian. 

![](jacobian_formula.png)

# Exercices/Contributions

* Please
* Please

%}
