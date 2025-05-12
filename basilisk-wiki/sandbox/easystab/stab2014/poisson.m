%{
# Theory for solving the equation of poisson
After solving the problem of derivation 3d,we know that the summation of second derivation 3D in the three direction is exactly a laplacian.So we can use the code of 3d derivation to solve the equation of poisson
I will solve a problem poisson with the result exact:f=cos(k*X*pi).*cos(l*Y*pi).*cos(m*Z*pi)
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

% 3D differentiation matrices
DDx=kron(speye(Nz),Dx);
DDxx=kron(speye(Nz),Dxx);
DDy=kron(speye(Nz),Dy);
DDyy=kron(speye(Nz),Dyy);
DDz=kron(dz,speye(Nx*Ny));
DDzz=kron(dzz,speye(Nx*Ny));

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
%{
Here A=DDXX+DDyy+DDzz is the derivation 3D,and it is exactly a Laplacian,so we can use the code de 3D to solve the problem de poisson

 %}
% System matrix,this is the theory maths from derivation of 3d to a Laplacian
A=DDxx+DDyy+DDzz;

%{
# I change the test function and the boundary condition
I have changed the test function from sine to cosine,the original function sin(pi*k*X).*sin(pi*l*Y).*sin(pi*m*Z) is from [poisson.m](http://basilisk.fr/sandbox/easystab/poisson3D.m).And
i change the the boundary conditon by using identity matrix instead of sparse matrix,and i also change the value of 'b' in boudary condition to move the position of graph,and this method come from the boundary condition of [poisson.m](http://basilisk.fr/sandbox/easystab/poisson3D.m),i just add a direction z.
 %}
% forcing
k=1; l=1; m=1;
b=-pi^2*(k^2+l^2+m^2)*cos(pi*k*X).*cos(pi*l*Y).*cos(pi*m*Z);
b=b(:);
% boundary conditions
II=eye(Nx*Ny*Nz);
loc=[top; bot; left; right; front; back];
A(loc,:)=II(loc,:);
b(loc)=0.5;

% solving the linear system
f=A\b;
f=reshape(f,Ny,Nx,Nz);

% Validation
solexact=cos(pi*k*X).*cos(pi*l*Y).*cos(pi*m*Z);
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
print('-djpeg','-r80','poisson.jpg');
%{ 
The new error is 1.5, it's much bigger than the error of sine function( err=2.3601e-10)

N =
        1716
err =
     1.5
   ![poission](/poisson.jpg)
 