%{
A code that looks very much like [venturi.m]() but in 2D and without the domain mapping. At the left boundary we impose a horizontal  jet flow, which difuses due to viscosity all the way to a Poiseuille flow at the outflow boundary. 

%}

clear all; clf; format compact

% parameters 
Re=100; % reynolds number
Nx=100; % number of grid nodes in z
Ny=50; %number of grid nodes in r
Lx=2; % domain length
Ly=1; % domain height

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,1,Ny);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% preparing boundary conditions
NN=Nx*Ny;
l.u=(1:NN)'; l.v=l.u+NN; l.p=l.v+NN;
II=speye(3*NN);
D.lap=D.yy+D.xx;

neuploc=[l.ctl;l.ctr;l.ctr-Ny];  % where to impose the neumann condition on the pressure
p0loc=2*Ny; % where to impose zero pressure
dir=[l.left;l.top;l.bot]; % where to put Dirichley on u and w

loc=[l.u(dir); l.v(dir); l.p(p0loc); ...
    l.u(l.right); ...
    l.v(l.right); ...
    l.p(neuploc)];

C=[II([l.u(dir);l.v(dir);l.p(p0loc)],:); ...     % Dirichlet on u,v,and p
   D.x(l.right,:)*II(l.u,:); ...   % Neuman on v at outflow
   D.x(l.right,:)*II(l.v,:); ...   % Neumann on u at outflow
   D.lap(neuploc,:)/Re*II(l.v,:)-D.x(neuploc,:)*II(l.p,:)]; % neuman constraint on pressure

% initial guess
V=zeros(NN,1);
U=exp(-((Y-Ly/2)/(Ly/8)).^2);
P=-X/Re; P=P-P(p0loc); % pressure zero at p0loc
q0=[U(:);V(:);P(:)];

% Newton iterations
disp('Newton loop')
q=q0;
quit=0;count=0;
while ~quit     
 
    % the present solution and its derivatives
    U=q(l.u); V=q(l.v); P=q(l.p);
    Ux=D.x*U; Uy=D.y*U; 
    Vx=D.x*V; Vy=D.y*V;
    Px=D.x*P; Py=D.y*P;

%{
# The nonlinear function and the Jacobian

This is the core of the code. The function for which we look for a root is
$$
f(q)=f\begin{pmatrix}
u\\v\\p
\end{pmatrix}
=
\begin{pmatrix}
-uu_x-vu_y-p_x+\Delta u/Re \\
-uv_x-vv_y-p_y+\Delta v/Re \\
u_x+v_y
\end{pmatrix}=0
$$
with the Laplacian $\Delta=\partial_{xx}+\partial_{yy}$. For the Newton method we need the Jacobian, we do a small perturbation
$$
f(q+\tilde{q})\approx f(q)+A\tilde{q}
$$
with the Jacobian
$$
A=\begin{pmatrix}
-u\partial_x-u_x-v\partial_y+\Delta/Re& -u_y& -\partial_x\\
       -v_x&-u\partial_x-v\partial_y-v_y+\Delta/Re& -\partial_y\\
       \partial_x& \partial_y&0
\end{pmatrix}
$$
%}
    % nonlinear function
    f=[-U.*Ux-V.*Uy+D.lap*U/Re-Px; ...
       -U.*Vx-V.*Vy+D.lap*V/Re-Py; ...
       Ux+Vy];
    
    % Jacobian 
    A=[-(spd(U)*D.x+spd(Ux)+spd(V)*D.y)+D.lap/Re, -spd(Uy), -D.x; ...
       -spd(Vx),-(spd(U)*D.x+spd(V)*D.y+spd(Vy))+D.lap/Re, -D.y; ...
       D.x, D.y, Z];
     
    % Boundary conditions 
    f(loc)=C*(q-q0);
    A(loc,:)=C;

    % plotting
    subplot 311;
    surf(X,Y,reshape(P-1,Ny,Nx)); view(2); shading interp; hold on
    
    sely=1:Ny; selx=1:6:Nx;
    ww=reshape(U,Ny,Nx); vv=reshape(V,Ny,Nx); 
    quiver(X(sely,selx),Y(sely,selx),ww(sely,selx),vv(sely,selx),'k');
    axis([0,Lx,0,1]);
    xlabel('z'); ylabel('r'); title('Pressure P'); grid off;hold off
    
    subplot 312;
    surf(X,Y,reshape(U,Ny,Nx)); view(2); shading interp; 
    xlabel('x'); ylabel('y'); title('horizontal velocity U'); grid off
    
    subplot 313;
    surf(X,Y,reshape(V,Ny,Nx)); view(2); shading interp; 
    xlabel('x'); ylabel('y'); title('vertical velocity V'); grid off
    drawnow
    
    % convergence test  
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if count>50|res>1e5; disp('no convergence');break; end
    if res<1e-5; quit=1; disp('converged'); continue; end
    
    % Newton step
    q=q-A\f;   
    count=count+1;
end


set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','jet_2D.png')

%{

The screen output for the convergence:

    Newton loop
    0  24.5262
    1  3.0485
    2  0.24814
    3  0.00075455
    4  7.4988e-10
    converged

And the figure, where we see that there are two recirculating bubbles above and below the jet:

![The figure](jet_2D.png)

%}