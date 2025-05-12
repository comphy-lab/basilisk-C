%{
# The Orr-Sommerfeld and Squire equation

Here we do just as in [poiseuille_uvp.m](), except that instead of using the primitive variables formulation $(u,v,p)$, we use the Orr-Sommerfeld and Squire equations, with the variables $v$ and $w$. We test it on the Poiseuille flow.
%}

clear all; clf

N=200; % number of grid points
alpha=1; % wavenumber in x
beta=0; % wavenumber in z
Re=10000; % Reynolds number

% differentiation matrices 
[y,DM] = chebdif(N,4); 
d.y=DM(:,:,1);    
d.yy=DM(:,:,2);    
d.yyyy=DM(:,:,4);    

I=eye(N);Z=zeros(N,N);II=blkdiag(I,I);
d.x=i*alpha*I;
d.z=i*beta*I;

% base flow and derivatives
u=1-y.^2;
up=-2*y;
upp=-2;

% laplacian
k2=alpha^2+beta^2;
lap=(d.yy-k2*I); 
lap2=(d.yyyy-2*k2*d.yy+k2^2*I); 

% system matrices
LOS=-diag(u)*lap*d.x+diag(upp)*d.x+lap2/Re;     
LSQ=-diag(u)*d.x+lap/Re;
LC= -diag(up)*d.z;

A=[LOS,Z;LC,LSQ]; 
E=[lap,I];

%{
# Boundary conditions

They are: zero at the walls for $v$ and $w$ (no-slip) and zero derivative at the walls for $v$ (that comes from the no-slip condition on $u$). The contraint matrix is thus
$$
C=\left(\begin{array}{c}
I|_{-1}&0\\
I|_{1}&0\\
0&I|_{-1}\\
0&I|_{1}\\
\partial_y|_{-1}&0\\
\partial_y|_{1}&0\\
\end{array}\right)
$$
such that the boundary conditions are implemented
$$
Cq=0
$$
where $q$ is the state vector
$$
q=\left(\begin{array}{c}
v\\
w
\end{array}\right)
$$
%}
% boundary conditions
v=1:N; w=v+N; % location vectors
loc=[v([1,2,N-1,N]),w([1,N])];
C=[II([v([1,N]),w([1,N])],:); ... % Dirichlet on v and w
    d.y([1,N],:), Z([1,N],:)];    % Neumann on v

A(loc,:)=C;
E(loc,:)=0;
    
% computing eigenmodes 
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

% loading the scanned figure for comparison and plotting
a=imread('poiseuille_spectra.png'); ss=size(a);
xx=linspace(0,1,ss(2));
yy=linspace(0,-1,ss(1));
image(xx,yy,a); axis xy
hold on

plot(-imag(s)/alpha,real(s),'ro')
grid on
axis([0,1,-1,0.1]);
xlabel('wave speed'); ylabel('exponential growth rate')

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','orr_sommerfeld_squire.png');

%{
![The spectrum](orr_sommerfeld_squire.png)

# Exercices/Contributions

* Please write down the equations and how to obtain them in the comments
* Please test it for a case with non-zero $\beta$ (an oblique wave)
* Please test it for the Blasius boundary layer
* Please do a quantitative comparison of the accuracy with the primitive variables formulation [poiseuille_uvp.m]().


%}

