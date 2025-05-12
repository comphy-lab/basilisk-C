%{

This code is a combination of [poisson2D.m]() and [diffmat_mapping_map2D.m](): we solve the Poisson problem on a stretched mesh (non-rectangular). Please see [README#differentiation-with-a-non-rectangular-mesh]() for more examples on using stretched meshes. For instance [diffmat_mapping.m]() and [diffmat_mapping_map2D.m]().

%}

clear all; clf; clc

% parameters and flags
Nx=20; % gridpoints in x 
Ny=15; % gridpoints in x  
Lx=1; % domain size in x
Ly=1; % domain size in y

% 1D and 2D differentiation matrices
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,Lx,Nx,3);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,Ly,Ny,3);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% mesh mapping
etay=1+0.3*x; 
etax=1-0.3*cos(y);  
X=X.*repmat(etax,1,Nx); 
Y=Y.*repmat(etay',Ny,1);   
D=map2D(X,Y,D);

% System matrix
A=D.xx+D.yy;

% Forcing
k=2; ell=1;
b=-pi^2*(k^2+ell^2)*sin(pi*k*X).*sin(pi*ell*Y);
b=b(:);
solexact=sin(pi*k*X).*sin(pi*ell*Y);

%{
# Boundary conditions
Here it is a little bit different from [poisson2D.m](), since the exact solution that we use is zeros on the boundaries of a square domain of size 1. Here the domain is not square, so instead of imposing a homogeneous boundary condition, we impose that on the boundaries of the domain, the computed solution has to be equal to the value of the exact solution.
%}

% boundary conditions
loc=[l.top; l.bot; l.left; l.right];
A(loc,:)=I(loc,:);
b(loc)=solexact(loc);

% solving the linear system
f=A\b;

% plotting the result
subplot(1,2,1);
mesh(X,Y,reshape(f,Ny,Nx)); view(2)
xlabel('x'); ylabel('y'); zlabel('f')
title('Poisson problem');

subplot(1,2,2);
mesh(X,Y,reshape(f,Ny,Nx)-solexact);
xlabel('x'); ylabel('y'); zlabel('f')
title('The error');

% printing the figure
set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','poisson2D_mapping.png')

% display of the error
err=norm(f-solexact(:))

%{
![Poisson solution and error](poisson2D_mapping.png)

# Exercices/contributions

* Please 
* Please
* ...

%}
