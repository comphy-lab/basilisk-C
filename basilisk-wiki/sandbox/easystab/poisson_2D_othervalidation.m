%{
# Poisson problem in 2D
Here we solve a Poisson problem in 2D with another forcing solution to validate the [first case](http://basilisk.fr/sandbox/easystab/poisson2D.m). The boundary conditions are the same. 
%}

clear all; clf

%%%% parameters and flags
Nx=40; % gridpoints in x 
Ny=40; % gridpoints in x  
Lx=1; % domain size in x
Ly=1; % domain size in y

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

% 2D differentiation matrices
Dx=kron(dx,eye(Ny));
u=Dx; v=eye(Ny);
Dxx=kron(dxx,eye(Ny));
Dy=kron(eye(Nx),dy);
Dyy=kron(eye(Nx),dyy);
[X,Y]=meshgrid(x,y);

%{
# Extracting the boundary cell locations
The most important aspect of this example is to show how to treat the boundary conditions in 2D when we do the transformation from array representatio to the vector representatio of the unknown. The first thing to do is to build the variable `dom`, which tells the number of every cells in the vector representation, but transformed to the array shape. This way we have an array with the same shape as the physical variables, but which gives us the location of every point in the final vector representation. Then we can extract the `top`, `bottom`, `left`and `right`location vector from the row and lines of `dom`. Note that here we chose that the corners belong to the `top` and `bottom` and not to the `left` and `right`. Sometimes it may be better to do the opposite, but this will depend on the kind of boundary conditions you want to implement. 
%}
% locations at the boundary
dom=reshape(1:Nx*Ny,Ny,Nx);
top=dom(1,1:end); top=top(:); 
bot=dom(end,1:end); bot=bot(:); 
left=dom(2:end-1,1); left=left(:); 
right=dom(2:end-1,end); right=right(:); 

% System matrix
A=Dxx+Dyy;

%{
# The exact solution
We do a Poisson problem with homogeneous boundary conditions at all boundaries (the value of `f` is zero), but with a forcing term `b`
$$
\Delta f=b
$$
and the forcing is
$$
b = -2\pi^2 cos(\pi k X) sin(\pi l Y)^2) sin(\pi l y)^2  -  2 \pi^2 cos(\pi k Y) sin(\pi l X)^2
$$
and the exact solution is 
$$
f = sin(\pi X)^2 sin(\pi Y)^2
$$
This solution is easy to check by doing a Fourier transform. We suppose we have an harmonic solution
$$
f(x,y)=\hat{f}\exp(i\pi kx+i\pi ly)+\text{Complex conjugate}
$$
which is fine as well as the boundary conditions that we impose fit with this (which is the case here). We have 
$$f_x=-i\pi k \hat{f} \exp(...)+C.C.$$ 
and 
$$
f_x=-i\pi l \hat{f} \exp(...)+C.C,$$ 
thus
$$
(-\pi^2k^2-\pi^2l^2)\hat{f}\exp(...)=-\pi^2(k^2+l^2)\exp(...)
$$
thus $\hat{f}=1$.
%}

% New forcing case
k=2; l=1;
b=(-2*pi^2*cos(pi*k*X).*sin(pi*l*Y).^2)-2*pi^2*cos(pi*k*Y).*sin(pi*l*X).^2; 
b=b(:);

%{
# Imposing the boundary conditions
To impose the boundary conditions, we replace the lines of $A$ with the lines of the identity matrix corresponding to the elements on the boundaries, and we replace the boundary elements of the forcing term with zeros. Doing this we effectively impose that the values of $f$ at the boundary cells must be zero.
%}

% boundary conditions
II=eye(Nx*Ny); ZZ=zeros(Nx*Ny,Nx*Ny);
loc=[top; bot; left; right];
A(loc,:)=II(loc,:);
b(loc)=0;
test=kron(b,ones(1,(Nx*Ny)));


% solving the linear system
f=A\b;

% plotting the result
subplot(1,2,1);
mesh(X,Y,reshape(f,Ny,Nx));
xlabel('x'); ylabel('y'); zlabel('f')
title('Poisson problem');

subplot(1,2,2);
solexact=(sin(pi*X).^2).*(sin(pi*Y).^2);
mesh(X,Y,reshape(f,Ny,Nx)+solexact);
xlabel('x'); ylabel('y'); zlabel('f')
title('The error');
