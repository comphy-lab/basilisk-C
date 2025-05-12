%{
Just like for [multi_domain_1D.m]() and [multi_domain_2D.m]() but with an interesting domain (that we could not do otherwise than using multi-domains). Please check these other codes for the detailed explanations.

%}

clear all; clf; 

% parameters 
nx=20; % number of grid nodes in x
ny=20; % number of grid nodes in y
method='cheb'; % which discretization to use
pts=9;  % how mnay points in the finite difference stencil

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D(method,0,1,nx,pts);
[d.y,d.yy,d.wy,y]=dif1D(method,0,1,ny,pts);
[D1,l1,X1,Y1,Z1,I1,NN1]=dif2D(d,x,y);

[d.x,d.xx,d.wx,x]=dif1D(method,1,2,nx,pts);
[d.y,d.yy,d.wy,y]=dif1D(method,0,1,ny,pts);
[D2,l2,X2,Y2,Z2,I2,NN2]=dif2D(d,x,y);

[d.x,d.xx,d.wx,x]=dif1D(method,0,1,nx,pts);
[d.y,d.yy,d.wy,y]=dif1D(method,1,2,ny,pts);
[D3,l3,X3,Y3,Z3,I3,NN3]=dif2D(d,x,y);

D.x=blkdiag(D1.x,D2.x,D3.x);
D.y=blkdiag(D1.y,D2.y,D3.y);
D.xx=blkdiag(D1.xx,D2.xx,D3.xx);
D.yy=blkdiag(D1.yy,D2.yy,D3.yy);
I=blkdiag(I1,I2,I3);
Z=blkdiag(Z1,Z2,Z3);

NN=NN1+NN2+NN3;
f1=(1:NN1)'; 
f2=(f1(end)+(1:NN2))';
f3=(f2(end)+(1:NN3))';

X=[X1(:);X2(:);X3(:)];
Y=[Y1(:);Y2(:);Y3(:)];


%{
We solve simply
$$
f_{xx}+f_{yy}=0
$$
%}
% System matrices
A=D.xx+D.yy;
b=0*X;

%{
For validation, we use an analytical solution of the Poisson problem with Dirichlet boundary conditions.
$$
f(x,y)=3x^2y-y^3
$$
%}
% the analytical solution
ftheo=3*X.^2.*Y-Y.^3;

%{
To get the connections, it is good to draw a sketch of how the domains are located with respect to each others.
%}
% Boundary conditions
dir=[f1([l1.left; l1.bot]); ...
    f2([l2.right; l2.top; l2.bot]); ...
    f3([l3.left; l3.top; l3.right])];
loc=[dir; f1(l1.right);f2(l2.left);f1(l1.top);f3(l3.bot)];
C=[I(dir,:); ...
    I(f1(l1.right),:)-I(f2(l2.left),:); ...
    D.x(f1(l1.right),:)-D.x(f2(l2.left),:); ...
    I(f1(l1.top),:)-I(f3(l3.bot),:); ...
    D.y(f1(l1.top),:)-D.y(f3(l3.bot),:)];

A(loc,:)=C;
b(loc)=C*ftheo;

% solve system
f=A\b;

% validation
err=norm(f-ftheo,2)

mesh(X1,Y1,reshape(f(f1),ny,nx),'edgecolor','b');
hold on
mesh(X2,Y2,reshape(f(f2),ny,nx),'edgecolor','r');
mesh(X3,Y3,reshape(f(f3),ny,nx),'edgecolor','m');
view(72,50)

xlabel('x');ylabel('y');
legend('domain 1','domain 2','domain 3');

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','multi_domain_2D_tri2.png');

%{

Here is the figure produced by the code as it is:

![](multi_domain_2D_tri2.png)

And here is an other test with 
$$
f_{xx}+f{yy}=-1
$$
and homogeneous boundary conditions.

![](multi_domain_2D_tri.png)


%}