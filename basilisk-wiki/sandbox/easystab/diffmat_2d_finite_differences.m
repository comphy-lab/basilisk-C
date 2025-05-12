%{
# 2D differentiation matrices

In this code, we show and test the differentiation matrices for a 2D domain, based of the code written in this article [diffmat_2D.m](). This time, we will use the finite differences method in comparaison to the method of chebychev
%}
clear all; 
clf;
close ALL;

%%%% parameters and flags
Nx=50; % gridpoints in x 
Ny=50; % gridpoints in x  
%{
About the domain siez on both axes, it is necessary to have a periodical domain. Because we are actually using periodicals boundaries.
%}
Lx=2*pi % domain size in x
Ly=2*pi % domain size in y

%{
# 1D matrices
Here we just build the 1D differentiation matrice as we use to do. We build the dx and dy with different sizes and number of points.


##1D differentiation matrices using the method of chebyshev
%}
scale=-2/Lx;
[x,DM] = chebdif(Nx,1); 
dx=DM(:,:,1)*scale;    
x=(x-1)/scale; 

scale=-2/Ly;
[y,DM] = chebdif(Ny,1); 
dy=DM(:,:,1)*scale;    
y=(y-1)/scale; 

%{
##1D differentiation matrices using the method of finite differences and
periodicals boundaries
%}

% the grid
xx=linspace(0,Lx,Nx+1)';xx=xx(1:end-1); 
hx=xx(2)-xx(1); % the grid size for the x axe
yy=linspace(0,Ly,Ny+1)';yy=yy(1:end-1)
hy=yy(2)-yy(1); % the grid size for the y axe

% 1d diffmat for the x axe

dfdx=zeros(Nx,Nx);
dfdx(1,[2,Nx])=[1/2,-1/2]/hx; 
dfdx(Nx,[1,Nx-1])=[1/2,-1/2]/hx;
for ind=2:Nx-1
    dfdx(ind,ind-1:ind+1)=[-1/2, 0, 1/2]/hx;
end

% 1d diffmat for the y axe
dfdy=zeros(Ny,Ny);
dfdy(1,[2,Ny])=[1/2,-1/2]/hy; 
dfdy(Ny,[1,Ny-1])=[1/2,-1/2]/hy;
for ind=2:Ny-1
    dfdy(ind,ind-1:ind+1)=[-1/2, 0, 1/2]/hy;
end


%{
# 2D matrices
We build *Dx* and *Dy* out of *dx* and *dy* and with identity matrices of the good size. Remember that Octave is case-sensitive, that is `D`can be a different variable from `d`. We had already the 1D meshes `x` and `y`, and we now as well build mesh arrays using the function `meshgrid`. This functino does something close to `kron`and is very convenient for building arrays from functions of `x`and `y`, in a much easier way than one may think at first. Here too, `X`and `Y` are the 2D equivalents of `x`and `y`. 
And we also build *DFDx* and *DFDy* out of *dfdx* and *dfdy* and with identity matrices of the good size.
%}

% 2D differentiation matrices
Dx=kron(dx,eye(Ny));
Dy=kron(eye(Nx),dy);
[X,Y]=meshgrid(x,y);

DFDx=kron(dfdx,eye(Ny));
DFDy=kron(eye(Nx),dfdy);
[XX,YY]=meshgrid(xx,yy);


%{
# Test
To test ouuntitled.figr differentiation matrices, we will compare the analytical and numerical derivatives for trigonometric functions. You see here how convenient is `meshgrid`. Be aware here in uilding `f` `fx` and `fy`that we use the `.*`multiplication which is the element by element multiplication, as opposed to the `*`multiplication which by default is the matrix multiplication.
%}

% Analytical derivatives

f=cos(X).*sin(Y);
fx=-sin(X).*sin(Y);
fy=cos(X).*cos(Y);

%now the analytical derivatives for the finite differences grid
fdf=cos(XX).*sin(YY);
fdfx=-sin(XX).*sin(YY);
fdfy=cos(XX).*cos(YY);
%{
# Numerical derivative
The derivative are simply obtained by multiplication of the vector representation of $f$. This representation is obtained by saying that we want all the elements of the array *f*. This is a shortcut notation to transform an aray into a vector. Then we do a *reshape* in order to come back to the array representation for the plotting.
%}

% Numerical derivatives
fX=reshape(Dx*f(:),Ny,Nx)
fY=reshape(Dy*f(:),Ny,Nx)
fdfX=reshape(DFDx*fdf(:),Ny,Nx)
fdfY=reshape(DFDy*fdf(:),Ny,Nx)

%{
# Structure of the matrices
Using the *spy* command we display in a figure the sparsity structure of the two differentiation matrices. We see that they are very different as sketched in the figure above. Since Chebychev differentation matrices in 1D are full, the *Dy* is block diagonal with *dy* stacked on the diagonal, whereas *Dx* is banded.
 %}

% showing the structure of Dx, Dy, DFDx and DFDy
set(gcf,'paperpositionmode','auto');
figure(1)
subplot(2,2,1); spy(Dx); title('Dx');
subplot(2,2,2); spy(Dy); title('Dy');
subplot(2,2,3); spy(DFDx); title('DFDx');
subplot(2,2,4); spy(DFDy); title('DFDy');
print('-dpng','-r100','diffmat_cheb_dfd.png'); % save the figure
%{
![The structure of the matrices](/sandbox/easystab/diffmat_cheb_dfd.png)

# Comparison
We show the shape of the numerical derivatives and the erreor between the exact derivative and the numerical derivative. We have chosen to take very few gridpoints  because otherwise the error would be close to machine accuracy.
 %}

% results
figure(2)
subplot(2,4,1);mesh(X,Y,fX); xlabel('x'); ylabel('y'); title('Dx*f');
subplot(2,4,2);mesh(XX,YY,fdfX); xlabel('x'); ylabel('y'); title('DFDx*fdf');
subplot(2,4,3);mesh(X,Y,fY); xlabel('x'); ylabel('y'); title('Dy*f');
subplot(2,4,4);mesh(XX,YY,fdfY); xlabel('x'); ylabel('y'); title('DFDy*fdf');
subplot(2,4,5);mesh(X,Y,fx-fX); xlabel('x'); ylabel('y'); title('fx-Dx*f'); 
subplot(2,4,6);mesh(XX,YY,fdfx-fdfX); xlabel('x'); ylabel('y'); title('fdfx-DFDx*f');
subplot(2,4,7);mesh(X,Y,fy-fY); xlabel('x'); ylabel('y'); title('fy-Dy*f');
subplot(2,4,8);mesh(XX,YY,fdfy-fdfY); xlabel('x'); ylabel('y'); title('fdfy-DFDy*f');

print('-dpng','-r120','diffmat_cheb_dfd_err.png'); % save the figure

%{
The derivate is similar on the x axe with the method of chebichev and the method of
finite differences. And the first one is more accuracy. 
On the opposite, there is a difference between the two method on the y axe
%}
%{
![The comparison of exact and numerical derivatives](/sandbox/easystab/diffmat_cheb_dfd_err.png)

![The comparison of exact and numerical errors](/sandbox/easystab/diffmat_cheb_dfd_err_2.png)

return to [diffmat_2D.m]()

%}