%{
This is a code to test the idea of patching several domains to build more complex geometries. Here is just a test to check how it works in 1D, the more interesting case will be in 2D. Please check out [multi_domain_2D.m]() and [multi_domain_2D_tri.m]().

We solve a Poisson problem
$$
f_{xx}=-1
$$
with boundary conditions
$$
f(0)=0, f(2)=0
$$
%}

clear all; clf; 

% parameters 
nx=10; % number of grid nodes in x

% differentiation 
[d1.x,d1.xx,d.wx,x1]=dif1D('fd',0,1,nx,5);
[d2.x,d2.xx,d2.wx,x2]=dif1D('fd',1,1,nx,5);

%{
From this point, we build the differentiation matrices for the total system.
%}
% combining
d.x=blkdiag(d1.x,d2.x);
d.xx=blkdiag(d1.xx,d2.xx);

% useful
N=2*nx;
f1=(1:nx)'; 
f2=(nx+(1:nx))';
I=eye(N);
Z=zeros(N,N);

% System matrices
A=d.xx;
b=-ones(N,1);

%{
As boundary conditions, we simply enforce the usual ones at the boundary, but we also need to impose the connection between the two domains. We enforce that at the common point, the value of $f$ should be the same on the two domains, but also the first derivative should be the same.
%}
% Boundary conditions
dir=[f1(1); f2(nx)];
loc=[dir; f1(nx); f2(1)];
C=[I(dir,:); ...
    I(f1(nx),:)-I(f2(1),:); ...
    d.x(f1(nx),:)-d.x(f2(1),:)];

A(loc,:)=C;
b(loc)=0;

% solve system
f=A\b;

%{
For the validation, we compare with the analytical solution
$$
f_{theo}=x(2-x)/2.
$$
%}
% validation
ftheo=x.*(2-x)/2;

plot(x1,f(f1),'b',x2,f(f2),'r',x,ftheo,'mo')
legend('domain 1','domain 2','theory');
xlabel('x');ylabel('f');

err=norm(f-ftheo,2)

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','multi_domain_1D.png');


%{

![](multi_domain_1D.png)

%}