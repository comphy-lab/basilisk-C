%{
#Introduction
If you need any explanations on the basics of the code please check this page : [../diffmat_2D.m]()
Here we are going to test to differents points on this code. First, we need to verify the second and crossed derivatives and realize a convergence study.
The loop on the indice 'i' is used for the convergence study. We stopped at 40 because we found out that if we were to go further the error would start to increase due to high sensibility of the method used for the derivatives.
%}

clear all; clf

for i=3:40
%%%% parameters and flags
Nx=i; % gridpoints in x 
Ny=i; % gridpoints in x  
Lx=2*pi; % domain size in x
Ly=pi; % domain size in y

%{
#1D Matrices
%}

%1D differentiation matrices
scale=-2/Lx;
[x,DM] = chebdif(Nx,1); 
dx=DM(:,:,1)*scale;    
x=(x-1)/scale; 

scale=-2/Ly;
[y,DM] = chebdif(Ny,1); 
dy=DM(:,:,1)*scale;    
y=(y-1)/scale; 

%{
# 2D matrices
%}

% 2D differentiation matrices
Dx=kron(dx,eye(Ny));
Dy=kron(eye(Nx),dy);
[X,Y]=meshgrid(x,y);

%{
# Test
%}

% Analytical derivatives
f=cos(X).*sin(Y);
fx=-sin(X).*sin(Y);
fy=cos(X).*cos(Y);
fxx=-cos(X).*sin(Y);
fyy=-cos(X).*sin(Y);
fxy=-sin(X).*cos(Y);


%{
# Numerical derivative
%}

% Nuerical derivatives
fX=reshape(Dx*f(:),Ny,Nx);
fY=reshape(Dy*f(:),Ny,Nx);
fXX=reshape(Dx*Dx*f(:),Ny,Nx);
fYY=reshape(Dy*Dy*f(:),Ny,Nx);
fYX=reshape(Dy*Dx*f(:),Ny,Nx);
fXY=reshape(Dx*Dy*f(:),Ny,Nx);

%{
Here we save the infinite norm on the whole domain and we plot it to see the evolution of the error with the size of the meshgrid.
%}

er=fXY-fxy;
maxval=max(er);
maxxval=max(maxval);

figure(3)
hold on
plot(i, maxxval,'r*');set(gca,'Yscale','log');
set(gca,'Yscale','log');
legend('error with the size of the meshgrid')
xlabel('mesgrid size'); ylabel('error'); title('Visualisation of the error for the y-second derivatives');
%{
![Evolution of the error with the size of the meshgrid](/sandbox/easystab/stab2014/diffmat_2D_graphanderror.png)

The fact that the error re-increase after 20 points is due to the method used to calculate the derivatives which is very sensitive to the truncation error. We can see that the error dimish at approximately the same speed as the exponential because we are studying functions that are very easy to derive.
%}
%{
# Structure of the matrices
 %}
end

% showing the structure of Dx and Dy
figure(1)
subplot(1,2,1); spy(Dx); title('Struture of the Dx matrix');
subplot(1,2,2); spy(Dy); title('Structure of the Dy matrix');

%{
![The structure of the matrices](/sandbox/easystab/stab2014/diffmat_2D_spy.png)

# Validation
 %}

% results
figure(2)
subplot(3,3,1);mesh(X,Y,fX); xlabel('x'); ylabel('y'); title('Dx*f');
subplot(3,3,2);mesh(X,Y,fY); xlabel('x'); ylabel('y'); title('Dy*f');
subplot(3,3,3);mesh(X,Y,fx-fX); xlabel('x'); ylabel('y'); title('fx-Dx*f'); 
subplot(3,3,4);mesh(X,Y,fy-fY); xlabel('x'); ylabel('y'); title('fy-Dy*f');

subplot(3,3,5);mesh(X,Y,fxx-fXX); xlabel('x'); ylabel('y'); title('fxx-Dx*Dx*f');
subplot(3,3,6);mesh(X,Y,fyy-fYY); xlabel('x'); ylabel('y'); title('fyy-Dy*Dy*f');
subplot(3,3,7);mesh(X,Y,fxy-fXY); xlabel('x'); ylabel('y'); title('fxy-Dx*Dy*f');
subplot(3,3,8);mesh(X,Y,fxy-fYX); xlabel('x'); ylabel('y'); title('fyx-Dy*Dx*f');

%{
![Visualisation of the graph's functions and the error](diffmat_2Dgraphanderror.png)

%}
%{
The difference in the order of the error between the second derivatives in x and y is due to the difference of the wavelength if we consider the x or the y axis.

return to [../diffmat_2D.m]()

You can find the rest of our contributions here : [./alois&antoine]()
%}
