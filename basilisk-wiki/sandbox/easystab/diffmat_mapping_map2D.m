%{
This code is just like [diffmat_mapping.m]() but here we use the function [map2D.m]() instead of explicitely computing the differentiation matrices for the mapped domain.

%}
clear all; clf

% parameters and flags
Nx=21; % gridpoints in x 
Ny=22; % gridpoints in x  
Lx=pi; % domain size in x
Ly=pi; % domain size in y

% 1D and 2D differentiation matrices
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,Lx,Nx,3);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,Ly,Ny,3);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

%{
# The mapping

Here we use the function [map2D.m]() to compute the differentiation matrices that account for the mapping. The inputs are the coordinate arrays $X$ and $Y$ once mapped, and the original differentiation matrices in the structure *D*. The output is the new *D* that works for the mapped domain. 

%}

etay=1+0.1*x; 
etax=1-0.1*cos(y);  
X=X.*repmat(etax,1,Nx); 
Y=Y.*repmat(etay',Ny,1);   
D=map2D(X,Y,D);

% test the mapping
f=cos(X).*sin(Y);
fx=-sin(X).*sin(Y); fxx=-cos(X).*sin(Y);
fy=cos(X).*cos(Y);  fyy=-cos(X).*sin(Y);

% quantifying the error
ex=norm(fx(:)-D.x*f(:))
ey=norm(fy(:)-D.y*f(:))
exx=norm(fxx(:)-D.xx*f(:))
eyy=norm(fyy(:)-D.yy*f(:))


%{
# Size of the approximation error

The last commands give this screen output, approximation error close to machine accuracy.

* ex = 4.2286e-11
* ey = 2.6175e-10
* exx = 8.6893e-10
* eyy = 1.8709e-08

# Links

Please see [diffmat_mapping.m#links]() for links to other codes.

# Exercices/Contributions

* Please 
* Please
* ...

%}
