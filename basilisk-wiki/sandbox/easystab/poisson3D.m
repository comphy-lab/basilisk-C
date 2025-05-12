%{
# Poisson problem in 3D

Just as for [poisson2D.m]() but here in 3D. To see how we build the differentiation matrices for a 3D grid, please see [diffmat_3D.m](). 
%}

clear all; format compact; clf

% parameters
n=10;
Nx=n+1; % gridpoints in x 
Ny=n+2; % gridpoints in y  
Nz=n+3; % gridpoints in z

Lx=1; % domain size in x
Ly=1; % domain size in y
Lz=1; % % domain size in z

N=Nx*Ny*Nz % number of degrees of freedom

% 1D and 2D differentiation matrices
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,Lx,Nx,3);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,Ly,Ny,3);
[d.z,d.zz,d.wz,z]=dif1D('cheb',0,Lz,Nz,3);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% 3D differentiation matrices
DD.x=kron(speye(Nz),D.x);
DD.xx=kron(speye(Nz),D.xx);
DD.y=kron(speye(Nz),D.y);
DD.yy=kron(speye(Nz),D.yy);
DD.z=kron(d.z,speye(Nx*Ny));
DD.zz=kron(d.zz,speye(Nx*Ny));

[X,Y,Z]=meshgrid(x,y,z);

%{
# Extracting the boundary cell locations

%}

% locations at the boundaries
dom=reshape(1:N,Ny,Nx,Nz);
top=dom(end,:,:); top=top(:); 
bot=dom(1,:,:); bot=bot(:); 
left=dom(2:end-1,1,:); left=left(:); 
right=dom(2:end-1,end,:); right=right(:); 
front=dom(2:end-1,2:end-1,1); front=front(:); 
back=dom(2:end-1,2:end-1,end); back=back(:); 

% Show locations
subplot(1,2,1);
plot3(X(front),Y(front),Z(front),'b*'); hold on
plot3(X(left),Y(left),Z(left),'m*'); 
plot3(X(bot),Y(bot),Z(bot),'k*'); 
axis([0 Lx 0 Ly 0 Lz]);box on; grid on
xlabel('x');ylabel('y');zlabel('z');
legend('front','left','bot')

% System matrix
A=DD.xx+DD.yy+DD.zz;

% forcing
k=1; l=1; m=1;
b=-pi^2*(k^2+l^2+m^2)*sin(pi*k*X).*sin(pi*l*Y).*sin(pi*m*Z);
b=b(:);

%{
# Imposing the boundary conditions
%}

% boundary conditions
II=speye(N); ZZ=spalloc(N,N,0);
loc=[top; bot; left; right; front; back];
A(loc,:)=II(loc,:);
b(loc)=0;

% solving the linear system
f=A\b;
f=reshape(f,Ny,Nx,Nz);

% Validation
solexact=sin(pi*k*X).*sin(pi*l*Y).*sin(pi*m*Z);
err=max(max(max(abs(f-solexact))))


% Loop to show the function f
subplot(1,2,2);
for ind=1:Nz
    mesh(X(:,:,ind),Y(:,:,ind),f(:,:,ind));
    axis([0 Lx 0 Ly -1 1]);
    xlabel('x');ylabel('y');zlabel('f');
    title(['z=' num2str(z(ind))]);
    drawnow;pause(0.1)
end

%{

Here is the (comforting) screen output telling the difference between the analytical solution and the numerical solution:

    err =
       2.3601e-10

And the figure:

![The figure](poisson3D.png)

# Exercices/contributions

* Please 
* Please
* Please

%}
