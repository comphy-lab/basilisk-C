%{
# The thermal entrance problem

%}

% thermal entrance

clear all; clf; format compact

%%%% parameters 
Pe=10000; % Peclet number
Nx=101; % number of grid nodes in z
Ny=40; %number of grid nodes in y
Lx=30; % length in x of the domain [0,Lz]
Ly=1; % Radius of the domain
pts=5; % number of points in finite difference stencils
xch=Lx/4; % location of x where there is a temperature change
lch=Lx/100; % length in xfor the change of the temperature

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,Ly,Ny);
[D,l,X,Y]=dif2D(d,x,y);
NN=Nx*Ny; Z=spalloc(NN,NN,0); I=speye(NN); 

% System matrices 
A=(D.yy+diag(1./Y(:))*D.y+D.xx)/Pe-diag(1-Y(:).^2)*D.x;
b=zeros(NN,1);

% initial guess
T0=0*X; T0=(1-tanh((X-xch)/lch))/2; 
%mesh(X,Y,reshape(T0,Ny,Nx)); break

% boundary conditions
dir=[l.left;l.top;l.ctl;l.cbl;l.ctr];
neux=[l.right;l.cbr];
neuy=l.bot;
loc=[dir;neux;neuy];
C=[I(dir,:); D.x(neux,:); D.y(neuy,:)];
A(loc,:)=C;
b(loc,:)=C*T0(:);

% solve system
T=A\b;

% show solution
subplot(2,1,1);
T=reshape(T,Ny,Nx);
surf(X,Y,T); shading interp; view(2)
xlabel('x');ylabel('r');title('temperature distribution')

% validation
subplot(2,1,2);
indvec=1:10:Nx;
co=jet(length(indvec));
% rescale the profiles in y
for ind=1:length(indvec);
    lo=indvec(ind);
plot(T(:,lo),Pe^(1/3)*(1-y)/x(lo)^(1/3),'color',co(ind,:));
hold on
end

% the analytical solution for the profile (Leveque)
yy=linspace(0,4,100);
f=gammainc(2*yy.^3/9,1/3);
plot(f,yy,'k','linewidth',2);
xlabel('temperature');ylabel('scaled y')
title('self-similar solution and numerical');
hold off
axis([0,1,0,3])

