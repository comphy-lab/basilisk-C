%{

This code is a validation of the stretching formula for [free_surface_mapping.m]().
In this code, we do not use the numerical mapping [map2D.m](): we want the analytical formula because the stretching is itself one of the unknown of the system.
%}




clear all; clf; format compact
disp('%%%%%%%%%')

% parameters
Lx=1;
Ly=1;
Nx=20;
Ny=20;
method='cheb';

% differentiation
[d.x,d.xx,d.wx,x]=dif1D(method,0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D(method,0,Ly,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

%{
# Stretching function

We introduce a stretching of the mesh
$$
\bar{x}=x, \bar{y}=\eta(x)y
$$
where $x$ and $y$ are the rectangular computational coordinates, and $\bar{x}$ and $\bar{y}$ are the physical coordinates.
In this code, everything that is related to the physical domain has subscript *p*.

Also, in the doce, we have te variable *e* which is a vector with *Nx* elements, and *E* which is an array with $N_y\timesN_x$ elements, where we copy *e* for every grid point of the mesh.
%}
% stretching function eta
e=1-0.3*sin(x*pi);
ex=d.x*e; exx=d.xx*e;
E=ones(Ny,1)*e';E=E(:);
Ex=ones(Ny,1)*ex'; Ex=Ex(:);
Exx=ones(Ny,1)*exx'; Exx=Exx(:);

% physical mesh
Xp=X;
Yp=Y.*reshape(E,Ny,Nx);

%{
# The mapping

The first useful thing is to get the expression of the derivative of the computational coordinates with respect to the physical coordinates
$$
\begin{array}{l}
x_\bar{x}=1\\
y_\bar{x}=-y\eta^{-1}\eta_x\\
x_\bar{y}=0\\
y_\bar{y}=\eta^{-1}\\
\end{array}
$$

Then we use the composition of derivatives to get the expression of the differentiation in physical space in terms of the computational differentiation matrices
$$
\begin{array}{l}
\phi(x,y)_{\bar{x}}=\phi_x x_\bar{x}+\phi_y y_\bar{x}=\eta^{-1}\phi_y\\
\phi(x,y)_{\bar{y}}=\phi_x x_\bar{y}+\phi_y y_\bar{y}=\phi_x-y\eta^{-1}\eta_x\phi_y\\
\end{array}
$$
We do the same operations for the expression of the Laplacian 
$$
\phi_{\bar{x}\bar{x}}+\phi_{\bar{y}\bar{y}}=(\phi_\bar{x})_\bar{x}+(\phi_\bar{y})_\bar{y}
$$
and we get
$$
\phi_{\bar{x}\bar{x}}+\phi_{\bar{y}\bar{y}}=
\begin{array}{l}
\phi_{xx}(1)\\
+\phi_{xy}(-2y\eta^{-1}\eta_x)\\
+\phi_{yy}(y^2\eta^{-2}\eta_x^2+\eta^{-2})\\
+\phi_y(y[2\eta^{-2}\eta_x^2-\eta^{-1}\eta_{xx}]))
\end{array}
$$
%}
% Build physical space differentiation matrices
Dp.x=D.x ...
    -diag(Y(:).*E.^-1.*Ex)*D.y;
Dp.y=diag(E.^-1)*D.y;

Dp.lap=D.xx ...
    +diag(-2*Y(:).*E.^-1.*Ex)*(D.x*D.y) ...
    +diag(Y(:).^2.*E.^-2.*Ex.^2+E.^-2)*D.yy ...
    +diag(Y(:).*(2*E.^-2.*Ex.^2-E.^-1.*Exx))*D.y;

% Validation
phi=cos(pi*Xp).*sin(pi*Yp);
phix=-pi*sin(pi*Xp).*sin(pi*Yp);
phixx=-pi^2*cos(pi*Xp).*sin(pi*Yp);
phiy=pi*cos(pi*Xp).*cos(pi*Yp);
phiyy=-pi^2*cos(pi*Xp).*sin(pi*Yp);

mesh(Xp,Yp,phi); hold on; 
plot3(x,Ly*e,0*x,'k-');
view(112,38); hold off
xlabel('x'); ylabel('y'); zlabel('phi'); title('phi');
legend('phi','e')

err_x=norm(phix(:)-Dp.x*phi(:))
err_y=norm(phiy(:)-Dp.y*phi(:))
err_lap=norm(phixx(:)+phiyy(:)-Dp.lap*phi(:))

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','stretching_formula.png');

%{
Which gives the reassuring screen output

    err_x =
       5.1937e-08
    err_y =
       3.8849e-13
    err_lap =
       1.8580e-05

and the figure 
![Figure](stretching_formula.png)
%}










