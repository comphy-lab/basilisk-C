%{
# 3D differentiation matrices

We have [diffmat_2D.m](),so we can add the direction z to change it into 3d.And i test it with an exponential function.And i will compare the error with the trigonometric function.

%}
clear all; format compact

%%%% parameters and flags
n=20;
Nx=n+1; % gridpoints in x 
Ny=n+2; % gridpoints in y  
Nz=n+3; % gridpoints in z

Lx=2*pi % domain size in x
Ly=2*pi % domain size in y
Lz=2*pi; % % domain size in z

N=Nx*Ny*Nz

%1D differentiation matrices
scale=-2/Lx;
[x,DM] = chebdif(Nx,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 

scale=-2/Ly;
[y,DM] = chebdif(Ny,2); 
dy=DM(:,:,1)*scale;    
dyy=DM(:,:,2)*scale^2;    
y=(y-1)/scale; 

scale=-2/Lz;
[z,DM] = chebdif(Nz,2); 
dz=DM(:,:,1)*scale;    
dzz=DM(:,:,2)*scale^2;    
z=(z-1)/scale; 

% 2D differentiation matrices
Dx=kron(dx,speye(Ny));
Dxx=kron(dxx,speye(Ny));
Dy=kron(speye(Nx),dy);
Dyy=kron(speye(Nx),dyy);
%{
Here we just add a direction z,and we are now in 3d
%}
% 3D differentiation matrices,here we just add a direction z,and we are now in 3d
DDx=kron(speye(Nz),Dx);
DDxx=kron(speye(Nz),Dxx);
DDy=kron(speye(Nz),Dy);
DDyy=kron(speye(Nz),Dyy);
DDz=kron(dz,speye(Nx*Ny));
DDzz=kron(dzz,speye(Nx*Ny));

[X,Y,Z]=meshgrid(x,y,z);
%{
This is the derivation maths,and i have changed a trigonometric funcion to exponential function
%}
% Analytical derivatives
f=exp(X).*sin(Y).*sin(Z);
fx=exp(X).*sin(Y).*sin(Z)
fxx=exp(X).*sin(Y).*sin(Z);
fy=exp(X).*cos(Y).*sin(Z);
fyy=-exp(X).*sin(Y).*sin(Z);
fz=exp(X).*sin(Y).*cos(Z);
fzz=-exp(X).*sin(Y).*sin(Z);

% Numerical derivatives
fX=reshape(DDx*f(:),Ny,Nx,Nz);
fXX=reshape(DDxx*f(:),Ny,Nx,Nz);
fY=reshape(DDy*f(:),Ny,Nx,Nz);
fYY=reshape(DDyy*f(:),Ny,Nx,Nz);
fZ=reshape(DDz*f(:),Ny,Nx,Nz);
fZZ=reshape(DDzz*f(:),Ny,Nx,Nz);
%{
Here is the code for calculating the error
%}
% Validation
err_x=max(max(max(abs(fx-fX))))
err_xx=max(max(max(abs(fxx-fXX))))
err_y=max(max(max(abs(fy-fY))))
err_yy=max(max(max(abs(fyy-fYY))))
err_z=max(max(max(abs(fz-fZ))))
err_zz=max(max(max(abs(fzz-fZZ))))
%{
Here is the code for showing the graph
%}
% Loop to show the function f
for ind=1:Nz
    mesh(X(:,:,ind),Y(:,:,ind),f(:,:,ind));
    axis([0 Lx 0 Ly -1 1]);
    xlabel('x');ylabel('y');zlabel('f');
    title(['z=' num2str(z(ind))]);
    drawnow;pause(0.1)
end
%This code is to print a picture in right size
print('-djpeg','-r80','3d.jpg');
%{
# 3D the graph
This is graph of derivation 3d of the test funtion
![alt text](/3d.jpg)


%}
%{
# The new error
I find that the error of exponential function is a little bigger than trigonometric funtion

err_x =
  8.2991e-012
  
err_xx =
  3.9483e-010
  
err_y =
  8.9813e-012
  
err_yy =
  4.9109e-010
  
err_z =
  1.1937e-011
  
err_zz =
  5.0538e-010

