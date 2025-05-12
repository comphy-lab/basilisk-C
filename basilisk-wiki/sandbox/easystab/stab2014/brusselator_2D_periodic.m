%{
# Reaction-diffusion equations

We adapt a C code, [brusselator.c](http://basilisk.fr/src/examples/brusselator.c), solving coupled reaction-diffusion equations, for 2 conpounds $U$ and $V$ :
$$
\partial_t U = \nabla^2 U + k(ka - (kb + 1)U + U^2 V)
$$
$$
\partial_t V = DD \nabla^2 V  + k(kb U - U^2 V)
$$

%}

clear all; clf

%{
# Parameters

As the calculations are long, and the convergence depending on parameters,
we recommend to keep a spatial step of 0.5 and time step around 0.1

We keep the same fluid parameters as in [brusselator.c](http://basilisk.fr/src/examples/brusselator.c)
We only play on the value of µ.
%}
% parameters
Nx=40; % number of gridpoints in x
Ny=40; % number of gridpoints in y 
Lx=20; % domain length x
Ly=20; % domain length y
dt=0.1; % time step
tmax=270; % final time

mu=[7.8,4,3.3]; % values of main parameter µ
mut=length(mu);

kt=100./((mu-0.3).^4); % time factor correction to adapt to µ (low µ requires more calculation time)


% secondary parameters, not modified in this study
k=1;ka=4.5;DD=8;        
nu=sqrt(1/DD);
kbcrit=sqrt(1+ka*nu);

%{
# System matrices

We use periodical differentiation matrices, available in this [package](../easypack.zip), so we won't have to
set boundary conditions.

$$
\left(
\begin{array}{cc}
I & 0 \\
0 & I \\
\end{array}\right)
\left(\begin{array}{c}
U_t \\
V_t \\
\end{array}\right)=
\left(
\begin{array}{cc}
-k (kb+1) I+(Dxx+Dyy) & 0 \\
k.kb.I & DD (Dxx+Dyy)\\
\end{array}\right)
\left(\begin{array}{c}
U \\
V \\
\end{array}\right)+
\left(\begin{array}{c}
k.ka\\
0 \\
\end{array}\right)+
\left(\begin{array}{c}
k U^2 V\\
-k U^2 V \\
\end{array}\right)
$$

If we note 
$$
q=\left(
\begin{array}{c}
U \\
V \\
\end{array}\right)
$$

we have a system $$ E q_t=A q + b + c(q)$$
%}

% differentiation matrices
[d.x,d.xx,d.wx,x]=dif1D('fp',0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D('fp',0,Ly,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

II=blkdiag(I,I);
E=II;
u=1:NN;u=u(:);
v=NN+1:2*NN; v=v(:);



for mui=1:mut % Loop on the different values of µ

    kb=kbcrit*(1+mu(mui)); % parameter depending on the value of µ
    
% system matrices

A=[-k*(kb+1)*I+(D.xx+D.yy) , Z  ;...
    k*kb*I , DD*(D.xx+D.yy)];
b=[k*ka*ones(Nx*Ny,1);zeros(Nx*Ny,1)];

%{
# Initial conditions
The balance state is for U=ka and V=kb/ka.
We set the initial conditions to this state, with a random noise on V
concentration.

%}

% initial condition

uinit=ka*ones(Nx*Ny,1);
vinit=kb/ka+0.001*(rand(Nx*Ny,1)-0.5); %initial noise
q=[uinit(:);vinit(:)];

%{
# March in time

We have to integrate
$$
\int_t^{t+dt} E q_t dt=\int_t^{t+dt} (Aq+b+c(q)) dt
$$ 
We use an implicit method, which tends to be more stable than
Crank-Nicolson method in this case. Thus the integral becomes :
$$
E(q(t+dt)-q(t))\approx [A q(t+dt)+b+c(q(t))] dt
$$ 
$$
q(t+dt)\approx (E-Adt)^{-1}[Eq(t)+b dt+c(q(t)) dt]
$$

Please note that the computation of the 3 cases is quite long, especially for the last one, and may take about 3 min with an i7 on a 40x40 grid.
%}

% march in time matrix 
Mm=E-A*dt; % implicit method shows better convergence than crank nicolson

BB=Mm\(b*dt);

% marching loop
tvec=dt:dt:tmax;
Nt=length(tvec);


for ind=1:(Nt*kt(mui))
    
    kuuv=k*q(u).*q(u).*q(v); 
    C=[kuuv;-kuuv]; % non linear term
    
    q=BB+gmres(Mm,q+(C*dt),1000,1e-7); % one step forward, gmres is used to fasten computing.

    %plotting
    subplot(1,3,mui)
    surf(X,Y,reshape(q(u),Ny,Nx));
    view(2);
    axis([0,Lx-Lx/Nx,0,Ly-Ly/Ny]);
    shading interp;
    drawnow
    
end
title(['U concentration with mu=' num2str(mu(mui))]);
xlabel('x'); ylabel('y');

end

%{
# Results and validation

As in [brusselator.c](http://basilisk.fr/src/examples/brusselator.c), we
get some critical values of µ, for which the final concentration in U
shows patterns :

![](brusselator2d.png)

However, these critical values of µ are very different from the ones obtained in the
C program. Moreover, there seems to be a dependance on the spatial and time
step. 
The C program uses more complex and robust calculation method, as well as
fastest. This could explain these differences, and the difficulties i
encountered to make the calculation stable.
Any contribution on this point is welcomed.

%}