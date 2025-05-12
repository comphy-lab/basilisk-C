%{
# 3D differentiation matrices

This can be usefull, and this is very similar to what is done in [diffmat_2D_dif2D.m](), but with a third direction *z*. 

%}

clear all; format compact

% parameters and flags
n=20;   % to change easily the size of the problem
Nx=n+1; % gridpoints in x 
Ny=n+2; % gridpoints in y  
Nz=n+3; % gridpoints in z

Lx=2*pi % domain size in x
Ly=2*pi % domain size in y
Lz=2*pi; % % domain size in z

N=Nx*Ny*Nz % number of degrees of freedom

% 1D and 2D differentiation matrices
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,Lx,Nx,3);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,Ly,Ny,3);
[d.z,d.zz,d.wz,z]=dif1D('cheb',0,Lz,Nz,3);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

%{

# 3D differentiation

Now that the 1D and 2D differentiation matrices are built using [dif1D.m]() and [dif2D.m](), we can procedd to build the 3D differentiation matrices and the associated grid. We follow the same ideas of using the *kron* operator. For a general discussion on dimensionality, please see [pedagogy#1D,-2D-and-3D]().

%}
% 3D differentiation matrices
DD.x=kron(speye(Nz),D.x);
DD.xx=kron(speye(Nz),D.xx);
DD.y=kron(speye(Nz),D.y);
DD.yy=kron(speye(Nz),D.yy);
DD.z=kron(d.z,speye(Nx*Ny));
DD.zz=kron(d.zz,speye(Nx*Ny));

% Build the grid itself with meshgrid
[X,Y,Z]=meshgrid(x,y,z);

%{
# Validation
As we usually do, we compare exact derivatives and numerical derivatives.
%}

% Analytical derivatives
f=cos(X).*sin(Y).*cos(Z);
fx=-sin(X).*sin(Y).*cos(Z);
fxx=-cos(X).*sin(Y).*cos(Z);
fy=cos(X).*cos(Y).*cos(Z);
fyy=-cos(X).*sin(Y).*cos(Z);
fz=-cos(X).*sin(Y).*sin(Z);
fzz=-cos(X).*sin(Y).*cos(Z);

% Numerical derivatives
fX=reshape(DD.x*f(:),Ny,Nx,Nz);
fXX=reshape(DD.xx*f(:),Ny,Nx,Nz);
fY=reshape(DD.y*f(:),Ny,Nx,Nz);
fYY=reshape(DD.yy*f(:),Ny,Nx,Nz);
fZ=reshape(DD.z*f(:),Ny,Nx,Nz);
fZZ=reshape(DD.zz*f(:),Ny,Nx,Nz);

% Validation
err_x=max(max(max(abs(fx-fX))))
err_xx=max(max(max(abs(fxx-fXX))))
err_y=max(max(max(abs(fy-fY))))
err_yy=max(max(max(abs(fyy-fYY))))
err_z=max(max(max(abs(fz-fZ))))
err_zz=max(max(max(abs(fzz-fZZ))))

% Loop to show the function f
for ind=1:Nz
    mesh(X(:,:,ind),Y(:,:,ind),f(:,:,ind));
    axis([0 Lx 0 Ly -1 1]);
    xlabel('x');ylabel('y');zlabel('f');
    title(['z=' num2str(z(ind))]);
    drawnow;pause(0.1)
end

%{
And here is the (comforting) screen output telling the difference between the analytical derivative and the numerical derivative:

     err_x =
        2.2843e-14
     err_xx =
        4.2755e-13
     err_y =
        1.6542e-14
     err_yy =
        9.5714e-13
     err_z =
        2.0789e-14
     err_zz =
        6.7557e-13
# Links

In [poisson3D.m]() we use the 3D differentiation matrices to solve a Poisson problem in 3D. 

Since I did not take care of computing the integration weights in 3D and since I do not do very much things in 3D for now, I did not write a code dif3D.m which would be the equivalent in 3D of [dif1D.m]() and [dif2D.m]().


# Exercices/Suggestions

* Please
* Please
* Please


%}
