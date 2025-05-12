%{

This code is a demonstration of the spatial integration in 2D and with a mapped mesh. To learn more about mapped meshes, please see [README#differentiation-with-a-non-rectangular-mesh](). We test the integration in 1D in [integration.m]().

**Dependency**

* [easypack.zip]() 


%}


clear all; clf
format compact 
format long

% parameters
Lx=2*pi;
Ly=2*pi;
Nx=20;
Ny=20;

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,Lx,Nx);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,5);
[D,l]=dif2D(d);

%{
# Integration in 1D

The variable *d* is a structure. The integration weights in $x$ are in *d.wx*, and the integration weights in *y* are in *d.wy*. The grid vectors are *x* and *y*, and are column vectors. The integration weights are row-vectors such that the matrix multiplication with column vectors is the scalar product: the product of each elements and the sum of them
$$
\left(\begin{array}{cccc}
w_1 & w_2 & \dots & w_N
\end{array}\right)
\left(\begin{array}{c}
f_1 \\ f_2 \\ \vdots \\ f_N
\end{array}\right)
=
w_1f_1+w_2f_2+\dots+w_Nf_N
$$

# Stretching in *y*
Here we test the 2D integration when there is a stretching of the mesh. We do two cases of simple stetching, the first one we stretch in y and the second we stretch in x. For the stretching, we simply scale down the grid by building the vector *eta* as long as the *x* grid.
%}

% stretching in y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');disp('stretching in y');
[X,Y]=meshgrid(x,y);
eta=0*x+1+0.5*x.*(x-Lx)*4/Lx^2; 
Y=Y.*repmat(eta',Ny,1);  
Dm=map2D(X,Y,D);

%{
We chose to integrate something very simple: a constant function *f* of value 1, such that the surface integral of *f* is equal to the total area of the mesh.
%}
% the function to integrate
f=0*X+1;
subplot(2,2,1);
mesh(X,Y,f); view(2); axis tight
xlabel('x');ylabel('y');title('mapped domain')

%{
Then to perform the integration is simple, we multiply element by element the integration weights in *D.w* by the elements of *f* and we sum. The theoretical surface is he surface of the square $L_xL_y$ minus the surface of he parabola defined by *eta*.
%}
% surface integral
itheo=2*Lx^2/3
inum=sum(sum(Dm.w.*(0*X+1)))
err=abs(itheo-inum)

%{
For the integration in *y*, we multiply element by element the integration weights *D.wx* with the function *f*, then we sum on the direction on which we integrate. Since here we integrate along x, the sum is done along the dimension 2 of the arrays. This value should be equal of Ly*eta.
%}
% integration in y
disp('Integration in y');
iy=sum(Dm.wy.*f,1);
subplot(2,2,2);
plot(x,iy,'b',x,Ly*eta,'r.')
xlabel('x');ylabel('integral in y');
axis([0,Lx,0,Ly]); grid on



%{
# Stretching in x
We chose something simpler for the stertching in $x$, we just displace the mesh to the right with the same function eta. The total surface is thus unchanged from $L_xL_y$.
%}
% stretching in x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');disp('stretching in x');
[X,Y]=meshgrid(x,y);
eta=1+y.*(y-Ly)*4/Ly^2; 
X=X-repmat(eta,1,Nx);  
Dm=map2D(X,Y,D);

% the function to integrate
f=0*X+1;
subplot(2,2,3);
mesh(X,Y,f); view(2); axis tight
xlabel('x');ylabel('y');title('mapped domain')

% surface integral
itheo=d.wy*(eta*Lx)
inum=sum(sum(Dm.w.*(0*X+1)))
err=itheo-inum

% integration in x
ix=sum(Dm.wx.*f,2);
subplot(2,2,4);
plot(y,ix,'b',y,Lx,'r.')
xlabel('y');ylabel('integral in x');
xlim([0,Ly]); grid on

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','integration_2D.png')

%{
# The figure

![The shape of the domains and the integrals](integration_2D.png)

# Exercices/Contributions

* Please test other types of stretchings
* Please integrate without stretching
* Please integrate a non-constant function

%}
