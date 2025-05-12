%{

This troubleshooting example is the follow-up of [troubleshooting_jacobian_formula.m](). Here we want to show you how to chose and add constraints on the ressure to ensure that the Jacobian is invertible for Navier-Stokes. This will be useful for you if you have encountered this message when doing your Newton iterations for Navier-Stokes:

    Warning: Matrix is close to singular or badly scaled.
             Results may be inaccurate. RCOND = 2.585247e-17. 
    > In troubleshooting_jacobian_invertibility at 75

This is based on the example of a jet flow [jet_2D.m]().

%}

clear all; clf; format compact

% parameters
Re=10; % reynolds number
Nx=16; % number of grid nodes in z
Ny=15; %number of grid nodes in r
Lx=2; % domain length
Ly=1; % domain height
method='fd';

% differentiation
[d.x,d.xx,d.wx,x]=dif1D(method,0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D(method,0,1,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);
D.lap=D.yy+D.xx;

% useful
l.u=(1:NN)'; l.v=l.u+NN; l.p=l.v+NN; % locations
n=3*NN; % The total number of degrees of freedom
II=speye(n);

% initial guess
U=exp(-((Y-Ly/2)/(Ly/8)).^2);
V=0*X;
P=0*X; 
q0=[U(:);V(:);P(:)];
q=q0;

% the present state and its derivatives
U=q(l.u); V=q(l.v); P=q(l.p);
Ux=D.x*U; Uy=D.y*U;
Vx=D.x*V; Vy=D.y*V;
Px=D.x*P; Py=D.y*P;

% nonlinear function
f=[-U.*Ux-V.*Uy+D.lap*U/Re-Px; ...
    -U.*Vx-V.*Vy+D.lap*V/Re-Py; ...
    Ux+Vy];

% Jacobian
A=[-(spd(U)*D.x+spd(Ux)+spd(V)*D.y)+D.lap/Re, -spd(Uy), -D.x; ...
    -spd(Vx),-(spd(U)*D.x+spd(V)*D.y+spd(Vy))+D.lap/Re, -D.y; ...
    D.x, D.y, Z];

%{
# Constraints

A matrix is not invertible if some of the lines can be built as a sum of other lines. It means that you have less equations than unkowns, and this means that you did not build your Jacobian properly. In other words when you do your Nwton step

    q=q-A\f

you are solving the linear system of equations

    Ax=f

for the unknown $x$, and that there are several possible choices for $x$.

When you look at a matrix, you should remember that the way it looks depends on the basis you have chosen. Here the default basis is a set of functions with zero everywhere except at one of the grid points. To see better what the matrix is, you can chose another basis: the basis of its eigenvectors. In this basis the matrix is diagonal with its eigenvalues on the diagonal. If some of the eigenvalues are zero, then your matrix is not invertible. 

If an eigenvalue $s_i$ is zero, it means that $A$ times its associated eigenvector $q_i$ is zero. One example of such possible eigenvector is a field with zero $u$, zero $v$ but constant and nonzero $p$. Indeed the pressure comes in the Navier-Stokes equations only through its derivative. So here the $x$ and the $y$ derivatives will be zero. To avoid this mode, we impose that the pressure somewhere on the grid be 0. This is a point Dirichlet bboundary condition on the pressure. You can chose where to impose this constraint on the grid. Here I do it in the bottom-left corner.

We can impose yet some other constraints, by writing this equation
$$
\Delta v/Re-p_x=0
$$
We impose this at gridpoints *neuploc*.

%}
% boundary conditions
neuploc=[l.ctl+Ny];  % where to impose the neumann condition on the pressure
p0loc=l.cbl+Ny; % where to impose p=0
dir=[l.left;l.top;l.bot;l.right]; % where to put Dirichlet on u and w

%{
For the present test, we remove all the contraints on the pressure 
%}
% remove the constraints to test their effect
neuploc=[];
p0loc=[];

loc=[l.u(dir); 
    l.v(dir); 
    l.p(p0loc); ...
    l.p(neuploc) ...
    ];

C=[II(l.u(dir),:); ...
    II(l.v(dir),:); ... % Dirichlet on u,v
    II(l.p(p0loc),:); ...     % Dirichlet on p
    D.lap(neuploc,:)/Re*II(l.v,:)-D.x(neuploc,:)*II(l.p,:); ... % neuman constraint on pressure
    ]; 

f(loc)=C*(q-q0);
A(loc,:)=C;

%{
Now we have built the Jacobian and imposed the boundary conditions. We can test wether it is invertible. If not, you will get the warning message

    Warning: Matrix is close to singular or badly scaled.
             Results may be inaccurate. RCOND = 2.585247e-17. 
    > In troubleshooting_jacobian_invertibility at 75

%}
% test Newton step
q=q-A\f;

%{
Here we compute the eigenvalues and the eigenvectors. You can do this with a coarse grid because this is costly. And then we only keep the small eigenvalues (smaller in absolute value than $1e^{-3}$.
%}
% compute eigenvalues (S) and eigenvectors (EV) of A
[EV,S]=eig(full(A));
s=diag(S); sel=abs(s)<1e-3; s=s(sel); EV=EV(:,sel);

%{
Here we display on screen the small eigenvalues and for each of them, we display the norm of $u, v$ and $p$, then show the associated pressure and the derivatives of the pressure on the figure.
%}
s
for num=1:length(s)
    disp(' ');
    
    U=EV(l.u,num); norm_u=norm(U)
    V=EV(l.v,num); norm_v=norm(V)
    P=EV(l.p,num); norm_p=norm(P)

    subplot(3,2,num)
    surf(X,Y,reshape(real(P),Ny,Nx));
    xlabel('x'); ylabel('y'); zlabel('pressure');
    title(['mode ' num2str(num)]);

    subplot(3,2,num+2)
    surf(X,Y,reshape(real(D.x*P),Ny,Nx));
    xlabel('x'); ylabel('y'); zlabel('pressure');
    title(['Px of mode ' num2str(num)]);

    subplot(3,2,num+4)
    surf(X,Y,reshape(real(D.y*P),Ny,Nx));
    xlabel('x'); ylabel('y'); zlabel('pressure');
    title(['Py of mode ' num2str(num)]);
end


set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','troubleshooting_jacobian_invertibility.png')

%{
The screen output is:

    s =
       1.0e-12 *
       -0.0646
       -0.2388
 
    norm_u =
       1.1758e-14
    norm_v =
       5.5416e-15
    norm_p =
        1.0000
 
    norm_u =
       4.1408e-14
    norm_v =
       2.1974e-14
    norm_p =
         1
# The figure

The pressure comes into the equations only through its $x$ and $y$ derivatives and only in the gridpoints where Navier-Stokes is enforced. Here at the boundary points we impose the boundary conditions, so the pressure don't come on this points. The pressure fields that are displayed here have their derivative zero everywhere except at the boundary points.

![The figure](troubleshooting_jacobian_invertibility.png)

%}
