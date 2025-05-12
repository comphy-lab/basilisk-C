%{
# Vibrating string

We use the code of [the vibrating string](../vibrating_string.m) and adapt the matrix to add one string and one mass.

To run the code you may need Chebychev differentiation matrices, available [here](../chebdif.m)

%}

clear all; clf

% parameters
N=100; % number of gridpoints per string
L=10; % strings length
c1=0.3; % wave velocity string 1
c2=1; % wave velocity string 2
dt=0.1; % time step
tmax=60; % final time
mass=4; % mass
t1=0.2; %  string 1 tension
t2=1; %  string 2 tension

% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 
Z=zeros(N,N); ZZ=zeros(2*N,2*N);
I=eye(N); II=eye(2*N);I4=eye(4*N+2);

%{
# System matrices

each strings follows the same equation, $v_t=c^2 f_{xx}$ and the mass follows its own equation, $m v_t=-t_1 (\partial_x f_1)|_L + t_2 (\partial_x f_2)|_0$. With, for each string and the mass, $f_t=v$
.Thus the new array representation for the system is 
$$
Eq_t=Aq
$$
with the state matrices $E$ and $A$ and the state $q$ defined as
$$
\left(
\begin{array}{cccccc}
I &  &  &  &  & \\
 & I &  &  &  & \\
 &  & I &  &  & \\
 &  &  & I &  & \\
 &  &  &  & m & \\
 &  &  &  &  & 1\\
\end{array}\right)
\left(\begin{array}{c}
v1_t \\
f1_t \\
v2_t \\
f2_t \\
vm_t \\
fm_t \\
\end{array}\right)=
\left(
\begin{array}{cc}
0 & c_{1}^2\partial_{xx} & 0 & 0 & 0 & 0 \\
I & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & c_{2}^2\partial_{xx} & 0 & 0\\
0 & 0 & I & 0 & 0 & 0\\
0 & -t_1 \partial_x|_L & 0 & t_2 \partial_x|_0 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0\\
\end{array}\right)
\left(\begin{array}{c}
v1 \\
f1 \\
v2 \\
f2 \\
vm \\
fm \\
\end{array}\right)
$$
%}

% locations in the state
v1=1:N;
f1=N+1:2*N;
v2=2*N+1:3*N;
f2=3*N+1:4*N;
vm=4*N+1;
fm=4*N+2;

for i=1:3 % A loop is used to operate the validation on iterations 2 and 3.
          
    if i==2
       c1=1;c2=1;t1=1;t2=1;mass=0; %specific parameters for null mass validation
    elseif i==3
       c1=1;c2=1;t1=1;t2=1;mass=1000000; %specific parameters for infinite mass validation
    end

% system matrices
A1=[Z,c1^2*dxx; I, Z]; A2=[Z,c2^2*dxx; I, Z];
E=eye(fm);
E(vm,vm)=mass;
A=[A1,ZZ,ZZ(:,[1,2]) ; ...
    ZZ,A2,ZZ(:,[1,2]) ; ...
    Z(1,:),-t1*dx(N,:),Z(1,:),t2*dx(1,:),0,0 ; ...
    ZZ(1,:),ZZ(1,:),1,0];

%{
# Boundary conditions
We consider that the two strings are fixed at their extremity, $f_1|_0=0=f_2|_L$, and they are connected to the mass, $f_1|_L=f_2|_0=f_m$.

The boundary condition can thus be expressed $Cq=0$
with the state $q$ and the constraint matrix $C$ 
$$
Cq=0=\begin{pmatrix}
0&I|_0&0&0&0&0\\
0&I|_L&0&0&0&-1\\
0&0&0&I|_0&0&-1\\
0&0&0&I|_L&0&0\\
\end{pmatrix}
\begin{pmatrix}
v_1\\f_1\\v_2\\f_2\\vm\\fm
\end{pmatrix}, 
$$

%}

% boundary conditions
loc=[f1(1),f1(N),f2(1),f2(N)];
E(loc,:)=0; 
A(loc,:)=I4(loc,:);
A(f1(N),fm)=-1;
A(f2(1),fm)=-1;

%{
# Initial condition
At initial time we say here that the string is initially deformed as a sinus (this satisfies the boundary conditions), and that the velocity is zero.
%}

% initial condition
if i==1
    q=[zeros(N,1) ; sin(pi*x/L) ; zeros(N,1) ; sin(pi*x/L) ; 0 ; 0]; 
elseif i ==2
    q=[zeros(N,1) ; sin(pi*x/L/2) ; zeros(N,1) ; sin(pi*(x+L)/L/2) ; 0 ; 0];% initial conditions for null mass
elseif i==3
    q=[zeros(N,1) ; sin(pi*x/L) ; zeros(N,1) ; sin(pi*x/L) ; 0 ; 0];% initial conditions for infinite mass
end

%{
# March in time
%}

% march in time matrix 
M=(E-A*dt/2)\(E+A*dt/2);

% unified parameters
f=[f1,f2];
v=[v1,v2];
x2=[x;L+x];

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
for ind=1:Nt    
    q=M*q; % one step forward
    
    % plotting
    
    if i==1 % plotting general case
    plot(x2,q(f),'b',x2,q(v),'r--',L,q(fm),'blacko');    
    axis([0,2*L,-2,2]);xlabel('x');ylabel('f');
    title('strings and mass movements');
    
    elseif i==2 % plotting validation for null mass
    figure(2);subplot(211);
    plot(x2,q(f),'.b',x,cos(c1*pi/L/2*ind/Nt*tmax)*sin(x*pi/L/2),'r',x+L,cos(c1*pi/L/2*ind/Nt*tmax)*sin((x+L)*pi/L/2),'r');
    title('2 identical strings and a null mass');
    axis([0,2*L,-1.2,1.2]);
    
    else % ploting validaiton for infinite mass
    figure(2);subplot(212);
    plot(x2,q(f),'.b',x,cos(c1*pi/L*ind/Nt*tmax)*sin(x*pi/L),'r',x+L,cos(c1*pi/L*ind/Nt*tmax)*sin(x*pi/L),'r');
    title('2 identical strings and an infinite mass');
    axis([0,2*L,-1.2,1.2]);
    end
    drawnow
    
end

% Legends are displayed after calculation as they strongly slow down dynamical graphics
    if i==1
        figure(1);legend('position','speed');
    elseif i==2
        figure(2);subplot(211);legend('numerical solution','theorical solution','location','ne');
    else
        legend('numerical solution','theorical solution','location','s');
    end
end



%{
Here is the final result of the strings position (blue) and speed (red) with 2 different strings

![](two_strings.gif)


# Validation

Null mass
:   We first test the code with a null mass and 2 identical strings ( velocity and
    tension ). We expect to find the same behaviour than a single 2L length
    string. 
    We know the theorical solution for mode 1 with homogenous Dirichlet conditions and adapted initial conditions :
    $$
    f(x,t)=sin(\pi x/L)*cos(ct \pi /L)
    $$
    We plot the solution at tmax below
    
Infinite mass
:   Now we test with an infinite mass.
    It is supposed to behave as a fixed point, thus we should see 2 independant L length strings vibrating.
    In mode 1 with adapted boundary and initial conditions, each one should follow this equation
    $$
    f(x,t)=sin(\pi x/2L)*cos(ct \pi /2L)
    $$
    Plotting at tmax :


![](validation.png)



 

%}