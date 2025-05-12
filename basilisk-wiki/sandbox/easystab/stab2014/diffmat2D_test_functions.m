%{
# TEST FUNCTIONS FOR [diffmat_2D.m](/sandbox/easystab/diffmat_2D.m)

This code is based upon [diffmat_2D.m](diffmat_2D.m). 

We simply use other test functions (in diffmat_2D.m, the function used is f = cos(X).*sin(Y)) for the validation of the derivative. We calculate de difference between the numerical solution and theory. You should bear in mind to do the calculations with the highest number of grid points possible in order to have the lowest error and hence an almost exact soution ! Here we choose Nx = Ny = 100. You will see that if you choose a lower grid points number, the code will not work for complicated functions.

When you use functions containing for example logarithm or 1/X, in other words functions that are not defined in a certain point like 0, you need to change X and Y which you use to create you grid (see below).
%}

clear all; clf

%%%% parameters and flags
Nx=100; % gridpoints in x 
Ny=100; % gridpoints in x  
Lx=2*pi; % domain size in x
Ly=pi; % domain size in y

%{
# 1D matrices

Here we just build the 1D differentiation matrice. We build "dx" (differentiation matrix for the variable x) and "dy" (differentiation matrix for the variable y) with different sizes and number of points.
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

We buid *Dx* and *Dy* out of *dx* and *dy* and with identity matrices of the good size. Remember that Octave is case-sensitive, that is `D`can be a different variable from `d`. We had already the 1D meshes `x` and `y`, and we now as well build mesh arrays using the function `meshgrid`. This functino does something close to `kron`and is very convenient for building arrays from functions of `x`and `y`, in a much easier way than one may think at first. Here too, `X`and `Y` are the 2D equivalents of `x`and `y`. 
%}

% 2D differentiation matrices
Dx=kron(dx,eye(Ny));
Dy=kron(eye(Nx),dy);
[X,Y]=meshgrid(x,y);

%{
# Test

To test our differentiation matrices, we will compare the analytical and numerical derivatives for sevral functions : polynomial functions, trigonometric functions, etc. You see here how convenient is `meshgrid`. We are building here `f` `fx` and `fy`that we use the `.*`multiplication which is the element by element multiplication, as opposed to the `*`multiplication which by default is the matrix multiplication.
%}

% Analytical derivatives

% If you'd like to test another of these functions, uncomment them and run
% the code. 

% polynomial functions

% f=(X.^5).*(Y.^6) ;
% fx=5*(X.^4).*(Y.^6);
% fy=6*(X.^5).*(Y.^5);

% trigonometric functions

%f=(sin(X).^2).*(sin(Y).^3);
%fx=2*cos(X).*sin(X).*(sin(Y).^3);
%fy=(sin(X).^2).*(3*cos(Y).*(sin(Y).^2));

% exponential and log-log functions

% X and Y must not contain any 0s because logarithm for instance is not
% defined in 0. Moreover, you can't divide by 0. So, we just add 1 to both X and Y in order no to have 0s in the grid. 
X=X+1;  
Y=Y+1;

f=exp(X).*log(Y);
fx=exp(X).*log(Y);
fy=exp(X)./Y;

% Numerical derivatives
fX=reshape(Dx*f(:),Ny,Nx);
fY=reshape(Dy*f(:),Ny,Nx);

% showing the structure of Dx and Dy
figure(1)
subplot(1,2,1); spy(Dx); title('Dx');
subplot(1,2,2); spy(Dy); title('Dy');

% results
figure(2)
subplot(2,2,1);mesh(X,Y,fX); xlabel('x'); ylabel('y'); title('Dx*f');
subplot(2,2,2);mesh(X,Y,fY); xlabel('x'); ylabel('y'); title('Dy*f');
subplot(2,2,3);mesh(X,Y,fx-fX); xlabel('x'); ylabel('y'); title('fx-Dx*f'); 
subplot(2,2,4);mesh(X,Y,fy-fY); xlabel('x'); ylabel('y'); title('fy-Dy*f');

%{
# Figures 

What you see below are the first second x and y-derivatives of the chosen function followed by  the difference between it's the analytical x and y-derivative and the numeerical solution. The defference is minimal, close to zero, proving the efficieny and the precision of the numerical method.

![A polynmial function](diffmat_2D_test_functions_poly.png)

![A trigonometric function](diffmat_2D_test_functions_trig.png)

![A function with logarithm](diffmat_2D_test_functions_log.png)

%}
