%{
# Eigenmodes of the diffusion equation 2D
Here we modified the [../diffmat_2D.m]() code to visualise at the Eigenmodes of the diffusion problem in 2D. You can find the code for the 1D problem here : [../diffusion_eigenmodes.m]()
%}

clear all; clf

%%%% parameters and flags
Nx=30; % gridpoints in x 
Ny=30; % gridpoints in x  
Lx=2*pi % domain size in x
Ly=pi % domain size in y

%{
# 1D matrices
%}

%1D differentiation matrices
scale=-2/Lx;
[x,DM] = chebdif(Nx,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 

scale=-2/Ly;
[y,DM] = chebdif(Ny,2); 
dy=DM(:,:,1)*scale;  
dyy=DM(:,:,2)*scale^2;    
y=(y-1)/scale; 

%{
# 2D matrices
%}

% 2D differentiation matrices
Dx=kron(dx,eye(Ny));
Dy=kron(eye(Nx),dy);
Dxx=kron(dxx,eye(Ny));
Dyy=kron(eye(Nx),dyy);

[X,Y]=meshgrid(x,y);

%{

%}


%{
# Test

Now that we have built the differenciation matrices we are going to bring to light the first nine eigen modes without the march in time. Indeed, given that we are not trying to solve the problem here the march in time isn't necessary.

%}


I=eye(Nx*Ny);
mu=1;

A=mu*(Dxx+Dyy);
E=eye(Nx*Ny);

% locations at the boundary
dom=reshape(1:Nx*Ny,Ny,Nx);
top=dom(1,1:end); top=top(:); 
bot=dom(end,1:end); bot=bot(:); 
left=dom(2:end-1,1); left=left(:); 
right=dom(2:end-1,end); right=right(:); 

%Boundary conditions
loc=[top; bot; right; left];
C=I(loc,:);
A(loc,:)=C;
E(loc,:)=0;

% computing eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); 
s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

%Visualisation of the natural modes
for ind=1:9
    subplot(3,3,ind);
surf(X,Y,reshape(U(:,ind),Ny,Nx)); shading interp; view(2); axis tight
title(s(ind))
end

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','diffusion_eigenmodes_2D.png');

%{
![Visualisation of the first nine eigenmodes for a vibrating membrane](../diffusion_eigenmodes_2D.png)
 
%}