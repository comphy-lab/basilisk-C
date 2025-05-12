%{

The poisson equation, once more, is used to test a multi-domain with a rectangular mesh mapped into a hollow disc. 
This is a preparation to do Navier-Stokes around the cylinder and get the Von Karman street of vortices.


Dependency;

* [easypack.zip]()
* [Lspacing.m]() that makes the stretching in the radial direction

%}

clear all; clf; format compact

% parameters 
Rmax = 1;
Rmin = 0.1;
Lref = 2*Rmin; % Diameter of inlet
Ntheta=50; % number of grid nodes in x
Nrho=30; %number of grid nodes in y
pts=5; % number of points in finite difference stencils
method='fd';

% DOMAIN 1
% differentiation 
[d.x,d.xx,d.wx,x]=dif1D(method,0,1,Ntheta,pts);
[d.y,d.yy,d.wy,y]=dif1D(method,0,1,Nrho,pts);
[D1,l1,~,~,Z1,I1,NN1]=dif2D(d,x,y);
% Meshing
[Theta, Rho] = meshgrid(linspace(-pi/2,pi/2,Ntheta),Rmin + Lspacing(Nrho,2,Rmax-Rmin));
[X2, Y2] = pol2cart(Theta,Rho);
X1 = -X2; Y1=Y2;
D1=map2D(X1,Y1,D1);

% DOMAIN 2
% differentiation
[d.x,d.xx,d.wx,x]=dif1D('fd',0,1,Ntheta,pts);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,1,Nrho,pts);
[D2,l2,~,~,Z2,I2,NN2]=dif2D(d,x,y);
D2=map2D(X2,Y2,D2);

D.x=blkdiag(D1.x,D2.x);
D.y=blkdiag(D1.y,D2.y);
D.xx=blkdiag(D1.xx,D2.xx);
D.yy=blkdiag(D1.yy,D2.yy);
Z=blkdiag(Z1,Z2);
II = blkdiag(I1,I2);
X=[X1(:);X2(:)];
Y=[Y1(:);Y2(:)];

NN=NN1+NN2;

D.lap=D.yy+D.xx; % if cylindrical +ym1*D.y;
A = D.lap;

K = 1; L = 1;
b=-pi^2*(K^2+L^2)*sin(pi*K*X).*sin(pi*L*Y);
ftheo=sin(pi*K*X).*sin(pi*L*Y);

%{
# Boundary conditions
%}

%%%% preparing boundary conditions
f1=(1:NN1)';  
f2=f1+NN2;

% Boundary Conditions
dir1=[l1.top;l1.bot;l1.cbl;l1.cbr;l1.ctl;l1.ctr]; % where to put Dirichlet on u1 and v1
lnk1 = [l1.left;l1.right]; % Link
loc1=[f1(dir1); % Dirichlet
      f1(lnk1); % Link
    ];
dir2=[l2.cbl;l2.bot;l2.cbr;l2.top;l2.ctl;l2.ctr]; % where to put Dirichlet on u2 and v2
lnk2 = [l2.left;l2.right]; % Link
loc2=[f2(dir2); % Dirichlet
      f2(lnk2)]; % Link

lnk = [lnk1;lnk2];
loc = [loc1;loc2];
dir = [dir1;dir2];

    
C=[II([f1(dir1);f2(dir2)],:);% Dirichlet on f
    II(f1(lnk1),:)-II(f2(lnk2),:); % Link Continuity
    D.x(f1(lnk1),:)-D.x(f2(lnk2),:)]; % Link Tangency on f

A(loc,:)=C;
b(loc)=C*ftheo; % same boundary conditions as the theoretical solution

% solve system
f=A\b;

%plot
mesh(X1,Y1,reshape(f(f1),Nrho,Ntheta),'edgecolor','b');
hold on
mesh(X2,Y2,reshape(f(f2),Nrho,Ntheta),'edgecolor','r');
xlabel('x');ylabel('y'); view(-10,60); 
legend('domain 1','domain 2');

% validation
err=norm(f-ftheo)

set(gcf,'paperposition','auto')
print('-dnpg','-r80','poisson_2D_hollow_disc.png');

%{
The screen output for the validation

    err =
       2.6102e-04

And the figure

![The figure](poisson_2D_hollow_disc.png)

%}
