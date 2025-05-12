%{
This is like in [multi_domain_1D.m]() but in 2D. And not yet with an interesting domain to allow comparison with an analytical solution. Please see [multi_domain_2D_tri.m]() for a complex geometry.
%}

clear all; clf; 

% parameters 
nx=20; % number of grid nodes in x
ny=20; % number of grid nodes in y

% en 2D
% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,0.5,nx);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,1,ny);
[D1,l1,X1,Y1,Z1,I1,NN1]=dif2D(d,x,y);

[d.x,d.xx,d.wx,x]=dif1D('cheb',0.5,0.5,nx);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,1,ny);
[D2,l2,X2,Y2,Z2,I2,NN2]=dif2D(d,x,y);

D.x=blkdiag(D1.x,D2.x);
D.y=blkdiag(D1.y,D2.y);
D.xx=blkdiag(D1.xx,D2.xx);
D.yy=blkdiag(D1.yy,D2.yy);
I=blkdiag(I1,I2);
Z=blkdiag(Z1,Z2);

NN=NN1+NN2;
f1=(1:NN1)'; 
f2=(NN1+(1:NN2))';
X=[X1(:);X2(:)];
Y=[Y1(:);Y2(:)];


% System matrices
A=D.xx+D.yy;
k=1; l=1;
b=-pi^2*(k^2+l^2)*sin(pi*k*X).*sin(pi*l*Y);

%{
Like in 1D, we impose the usual boundary conditions on the borders of the domain, and we also need to impose the conditions of connection between the two domains. We impose that the value of the solution on common edge are the same and we impose that the normal derivative is the same on common edges.
%}
% Boundary conditions
dir=[f1([l1.left; l1.top; l1.bot]); ...
    f2([l2.right; l2.top; l2.bot])];
loc=[dir; f1(l1.right);f2(l2.left)];
C=[I(dir,:); ...
    I(f1(l1.right),:)-I(f2(l2.left),:); ...
    D.x(f1(l1.right),:)-D.x(f2(l2.left),:)];

A(loc,:)=C;
b(loc)=0;

% solve system
f=A\b;

%{
The analytical solution is the same as for [poisson2D.m]().
%}
% validation
ftheo=sin(pi*k*X).*sin(pi*l*Y);
err=norm(f-ftheo,2)

mesh(X1,Y1,reshape(f(f1),ny,nx),'edgecolor','b');
hold on
mesh(X2,Y2,reshape(f(f2),ny,nx),'edgecolor','r');
xlabel('x');ylabel('y');
legend('domain 1','domain 2');

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','multi_domain_2D.png');

%{

![](multi_domain_2D.png)

%}