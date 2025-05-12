%{
# Convergence of 2D differentiation matrices

We want to do a convergence study with grid resolution to validate [diffmat_2D.m](/sandbox/easystab/diffmat_2D.m).

So we need to try a range of *Ny* rows and *Nx* columns to show the
convergence of 2D differentiation matrices.
%}
clear all; clf
Nmin=2;Nmax=60; % gridpoints minimum and maximum

%{
We put a double loop with the range of *Ny* and *Nx* between of each
minimum and maximum.
%}
for Ny=Nmin:Nmax
    for Nx=Nmin:Nmax

%%%% parameters and flags  
Lx=2*pi; % domain size in x
if Nx==Ny
    hx(1,Nx)=(Lx/Nx)^1;
    hx(2,Nx)=Nx*Ny;
end
Ly=pi; % domain size in y

%{
We don't change this part of code from original code.
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

% 2D differentiation matrices
Dx=kron(dx,eye(Ny));
Dy=kron(eye(Nx),dy);
[X,Y]=meshgrid(x,y);

% Analytical derivatives
f=cos(X).*sin(Y);
fx=-sin(X).*sin(Y);
fy=cos(X).*cos(Y);

% Numerical derivatives
fX=reshape(Dx*f(:),Ny,Nx);
fY=reshape(Dy*f(:),Ny,Nx);

%{
# Spectral error

To do a convergence study, we calculate the spectral error. It's the norm
two of difference between the analytical and numerical derivatives.
We calculate, in first, the spectral error for x and for y. Next, the
spectral error for grid resolution.
%}

%calcul of spectral error
a=norm(fx-fX,2);
b=norm(fy-fY,2);

%place the spectral error in a matrices
Cx(Nx,Ny)=a;    %spectral error in x for different Nx
Cy(Nx,Ny)=b;    %spectral error in y for different Ny
Ct(Nx,Ny)=sqrt(a^2+b^2);    %spectral error for couple x,y for different Nx/Ny
    end
end

%{
# Study convergence

To validate the convergence of solution for 2D differentiation matrices, we
trace the spectral error for grid resolution.
To compare different grids, we put different factors ahead Ny.
%}

%factor to multipliate Ny
fact1=1;
fact2=2/3;
fact3=1/3;

%loop to calculate the error in function of Nx and Ny
for c=Nmin:Nmax
    Cx1D(c)=Ct(c,round(c*fact1));    
    xerr(c)=c*round(c*fact1);
    Cx1Dbis(c)=Ct(c,round(c*fact2));
    xerrbis(c)=c*round(c*fact2);
    Cx1Dbis2(c)=Ct(c,round(c*fact3));
    xerrbis2(c)=c*round(c*fact3);
end

%calculate 
X_min=min([min(xerr) min(xerrbis) min(xerrbis2)]);
X_max=min([max(xerr) max(xerrbis) max(xerrbis2)]);

%{
# Drawing the convergence

We can now plot the convergence study for different couple of Nx/Ny. 
%}

%draw the errors for differents configuration of Nx/Ny
figure(1)
loglog(xerr,Cx1D,'*',xerrbis,Cx1Dbis,'*',xerrbis2,Cx1Dbis2,'*',hx(2,:),hx(1,:),'-')
xlim([X_min X_max]);
xlabel('Number of points of grid');ylabel('Spectral error')
title('Comparison of convergence study with grid resolution ');legend('Ny=Nx','Ny=0.7*Nx','Ny=0.3*Nx','h1');
set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','diffmat_2D_convergence_grid_resolution')

%{
![The study convergence of spectral error for different couple of Nx/Ny](/sandbox/easystab/stab2014/diffmat_2D_convergence_grid_resolution.png)

We can see that solutions converge. For small Ny against Nx (factor between 0.25 and 0.4), the spectral
error is decrease more slowly. For bigger Ny (factor between 0.4 and 1), the spectral
error is decrease more quickly. We need to put a factor enough to create a good grid
resolution.
At final, the error aims to 10e-14, this error is numerical error.

We can see that the error conduct like order 1 until around 200 points
for Ny near Nx. When Ny become small against Nx, the number of points
increase, around 400 points for Ny=0.3Nx for example.
So more the grid is of square shape, more the error will increase by order
quickly.
%}

%{
We can now validate the convergence for 2D differentiation matrices!
%}

%{
# Contributor's page
Link to page of contributor [Fabien](/sandbox/easystab/stab2014/fabien.m)
%}
