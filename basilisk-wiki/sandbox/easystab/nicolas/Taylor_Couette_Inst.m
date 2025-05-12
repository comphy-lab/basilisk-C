function [] = main()

%{
# Stability of Couette-Taylor Flow in primitive variables
This code was inspired by the following codes: [Stability of Couette flow]() and [stability of poiseuille flow](). 
%}


clear all; clf;
if exist("~/sandbox/easystab",'dir')
   addpath("~/sandbox/easystab");
end

% phyzical parameters
nu=1    % viscosity
m=1;    % the wave number in theta
k=1;    % the wave number in z
L=2;        % the height of the domain, from -1 to 1
n=100;      % the number of grid points
R1=1; % radius of the first cylinder
R2=4.035/3.55; %radius of the second cylinder
omega1=100/(R2-R1); % speed rotation of the first cylinder (int)
omega2=0/(R2-R1); % speed rotation of the second cylinder (ext)
Re=(omega1*R1*(R2-R1))/nu % the Reynolds Number of the inner cylinder 
Re2=(omega2*R2*(R2-R1))/nu % Reynolds number of the outer cylinder

% differentiation and integration
[D,DD,w,r]=dif1D('cheb',R1,R2-R1,n,3);
Z=zeros(n,n) ; I=eye(n) ;


% renaming the differentiation matrices
I=eye(n); Z=zeros(n,n);
dr=D; drr=DD;
dtheta=(i*m)*I; dtheta2=-m^2*I;
dz=i*k*I; dzz=-k^2*I;
L= (drr - dtheta2./r.^2 - dzz + dr./r); 
Lprime = L - diag(1/r.^2);
Lter = diag(2i*m/r.^2);



% base flow

gamma = ((omega2*R2^2)-(omega1*R1^2))/(R2^2-R1^2);
beta=(omega1-omega2)/(R1^(-2)-R2^(-2));
Vbase=gamma.*r + beta./r ;
DVbase = dr*Vbase;


% the matrices
A=[ ...
    diag(Vbase./r.*dtheta) + Lprime/Re, diag(-2*Vbase./r) - Lter/Re,  Z,  -dr; ...
    diag(Vbase./r + DVbase) + Lter/Re,  diag(Vbase./r.*dtheta)  + Lprime/Re,  Z,  -dtheta./r; ...
    Z,  Z,  diag(Vbase./r.*dtheta) + L/Re, -dz; ...
    1+r.*dr , dtheta , dz , Z ];

    
    
E=blkdiag(I,I,I,Z);

% locations on the grid
u=1:n; v=u+n; w=v+n; p=w+n;

% boundary conditions
III=eye(4*n);
DDD=blkdiag(dr,dr,dr,dr);

loc=[u(1) u(n) v(1) v(n) w(1) w(n) ];  
C=III(loc,:);
E(loc,:)=0;  A(loc,:)=C;

%{
# Validation of the results
 (here we have to include a comparison with a relevant reference)
%}

% computing eigenmodes and plotting
figure(1);
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
plot(imag(s)/k,real(s),'ro')
grid on
xlabel('wave speed'); ylabel('exponential growth rate')

set(gcf,'paperpositionmode','auto');
print('-dsvg','-r80','Spectrum2.svg');

%{
 ![**Figure :** Representation of the Spctrum for m=1,k=1,R2/R1 = 0.8.  ](Taylor_Couette_Inst/Spectrum2.svg)
%}


figure(2);
eigenmode = U(:,1); % first eigenmode
plot(real(eigenmode(u)),r,'r-',imag(eigenmode(u)),r,'r--',real(eigenmode(v)),r,'b-',imag(eigenmode(v)),r,'b--');
xlabel('u,v');ylabel('r');
legend('Re(u)','Im(u)','Re(v)','Im(v)');
legend('Location','East');
grid on;

set(gcf,'paperpositionmode','auto');
print('-dsvg','PlanePoiseuille_Mode_Re10000_k1.svg');




%%
ur = real(U(1:n,1));
ui = imag(U(1:n,1));
ntheta = 50;
theta = linspace(0,2*pi,ntheta);
cosmtheta=cos(m*theta);
sinmtheta=sin(m*theta);
u2D = ur.*cosmtheta;
[R,THETA] = meshgrid(r,theta);
figure;contourf(THETA,R,u2D'); colorbar;

end