%{
# Wave-like pertubations

This code is adapted from wave_like.m from the easystab project.

Modified by  D. Fabre to experiment the following things :

- Psi formulation compared to u-v-p formulation,

- adherence conditions versus no-slip conditions.

%}




%{
# Main program :
%}


clear all; clf;


% parameters
N=50;      % number of gridpoints
L=1;        % Fluid height in y
rho=1;      % fluid density
mu=1;    % fuid viscosity
alpha=4 ;    % wavenumber in x
g=1;        % gravity


[D,DD,wy,y]=dif1D('cheb',0,L,N,3);
u=1:N; v=u+N; p=v+N;



% validation
figure(1);
alphavec=linspace(0,10,100);
betavec=2*pi./[2*L,L,2/3*L,1/2*L];
for ind=1:length(betavec)
    stheo=-mu*(alphavec.^2+betavec(ind).^2);
    plot(alphavec,stheo,'r-','DisplayName', 'theory'); hold on
end

for alpha = 1:10
[s,U] = Wavelike_UVP(mu,alpha,L,rho,N);
[spsi,Upsi] = Wavelike_psi(mu,alpha,L,rho,N);

plot(alpha,s(1:4),'b*','DisplayName', 'UVP');
plot(alpha,spsi(1:4),'go','DisplayName', 'Psi');
xlabel('alpha');ylabel('exponential growth rate'); title('validation')
end


 hr = findobj('Color','r');
 hb = findobj('Color','b');
 hg = findobj('Color','g');
 hv = [hr(1) hb(1) hg(1)];
legend(hv);
%legend('analytical','UVP','psi');
grid on

set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','wave_like_UVP_Psi.png');


% show the velocity field of the first three eigenvectors
figure(2);
Lx=2*pi/alpha;
x=linspace(0,Lx,20);
for ind=1:4
    subplot(2,2,ind);
    q=U(:,ind);
    qphys=2*real(q*exp(i*alpha*x));
    quiver(x,y(1:2:N),qphys(u(1:2:N),:),qphys(v(1:2:N),:));
    axis equal; axis([0,Lx,0,L]);xlabel('x');ylabel('y'); 
    title(['mode' num2str(ind)]);
end

set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','wave_like_Modes.png');




%{


![The figure](wave_like_UVP_PSI.png)

Here are the eigenvalues, compared to the theory.

![The figure](wave_like_UVP_Modes.png)

On the figure, we show the velocity field for the four least stable eigenmodes. 
All off course have the same wavelength $2\pi/\alpha$ in $x$, but we see that as we go to more decaying modes, they have more and more oscillations in the wall-normal direction. The first mode corresponds to two counter-rotating vortices occupying the whole height between the walls. The second mode is four counter-rotating vortices, since now from $y=0$ to $L$, space is made for two vortices, and so on for the following modes.



%}






%{

# Function for U-V-P formulation :

%}



function [s,U] = Wavelike_UVP(mu,alpha,L,rho,N);

% 1D differentiation matrices
[D,DD,wy,y]=dif1D('cheb',0,L,N,3);
I=eye(N); Z=zeros(N,N);

% renaming the matrices
dy=D; dyy=DD;
dx=i*alpha*I; dxx=-alpha^2*I;
Delta=dxx+dyy;

% location vectors
u=1:N; v=u+N; p=v+N;

% System matrices
A=[mu*Delta, Z, -dx; ...
   Z, mu*Delta, -dy; ...
   dx, dy, Z];

E=[rho*I, Z, Z; ...
   Z, rho*I, Z; ...
   Z, Z, Z];

% boundary conditions
loc=[u(1),u(N),v(1),v(N)];
II=eye(3*N); DD=blkdiag(dy,dy,dy);
C=[DD([u(1),u(N)],:); II([v(1),v(N)],:)];

E(loc,:)=0;  
A(loc,:)=C;


% compute eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
end


%{

# Function for Psi formulation (Orr-Sommerfeld) :

%}

function [s,U] = Wavelike_psi(mu,alpha,L,rho,N);

% differentiation matrices 
[y,DM] = chebdif(N,4); 
y = L*y/2;
d.y=2/L*DM(:,:,1);    
d.yy=(2/L)^2*DM(:,:,2);    
d.yyyy=(2/L)^4*DM(:,:,4);    

I=eye(N);Z=zeros(N,N);
d.x=i*alpha*I;

% base flow and derivatives
u0=0;
up=0;
upp=0;

% laplacian
k2=alpha^2;
lap=(d.yy-k2*I); 
lap2=(d.yyyy-2*k2*d.yy+k2^2*I); 

% system matrices
LOS=mu*lap2;     

A=[LOS]; 
E=[lap];

% boundary conditions
psi=1:N;% location vectors

loc=[1,2,N-1,N];
C=[I([1,N],:) ; d.yy([1,N],:) ];

A(loc,:)=C;
E(loc,:)=0;
    
% computing eigenmodes 
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
end
