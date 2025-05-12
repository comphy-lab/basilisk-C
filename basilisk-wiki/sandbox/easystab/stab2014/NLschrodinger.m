%{
# NONLINEAR SCHRÖDINGER'S EQUATION FOR WATER WAVES

We are going to solve the non-linear form of Schrödinger's equation for water waves using the march in time that we have learnt from  [vibrating_string.m](/sandbox/easystab/vibrating_string.m). Here is the equation :
$$iU_t + U_{xx} = |U|^{alpha}U 
$$

We use the differentiation matrices to turn this equation into a matrix system as seen in [diffmat.m](/sandbox/easystab/diffmat.m).

%}

clear all; clf

% parameters
N=100; % number of gridpoints
L=10; % domain length
U=0.1; % wave velocity
dt=0.001; % time step
tmax=5; % final time
x0=L/2; %x-coordinate at initial time
l0=0.5; length width

% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;
dxx=DM(:,:,1)*scale^2;        
x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 

%{
# System matrices

Here we turn the Schrödinger's equation into a system matrix. The equations is :
$$
iU_t + U_{xx} = |U|^{alpha}U
$$ 
This is an equation with a single derivative, we can write this equation in the following way :
$$
Eq_{t} = Aq + F
$$ 
F, being the nonlinear term stored in a vector.

We can put an identity matrix on the left of the equation, in this way, we can impose the
boundary condition more easliy.

So the equation becomes :
$$
\begin{pmatrix}
I & Z \\
Z & I  \\
\end{pmatrix}
\left(\begin{array}{l}
U_{t} \\
U_{xt} \\
\end{array}\right)
=
\begin{pmatrix}
Z & -D_{xx} \\
I & Z  \\
\end{pmatrix}
\left(\begin{array}{l}
U \\
U_{x} \\
\end{array}\right)
+
\left(\begin{array}{l}
Nonlinear Term \\
\end{array}\right)
$$

The nonlinear term is : 
$$ |U|^2*U $$

%}

% system matrices
E=i*I;
A=-0.5*dxx; 

%{
# Boundary conditions
We apply homogeneous Dirichlet conditions to the boundaries.

We believe we are having issues with our boundary conditions bacause the solution diverges.

%}

% boundary conditions
loc=[1;N];
E(loc,:)=0;
A(loc,:)=I(loc,:);

%{
# March in time

Here we use the same methode with [vibrating_string.m](/sandbox/easystab/vibrating_string.m) to solve the equation.
%}

% march in time matrix 
Mm=(E-A*dt/2);
M=Mm\(E+A*dt/2);

%{
# Initial condition

We choose an initial condition here but you can generate it randomly. 
%}

% initial condition
q=1*[x*0+exp(-((x-x0)/l0).^2)]; 

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

for ind=1:Nt   
    nl=0.1*abs(q).^2.*q;
    nl(loc)=0;
    q=M*q+Mm\nl*dt; % one step forward
    
    % plotting
    plot(x,real(q),'b',x,imag(q),'r');    
    axis([0,L,-2,2])    
    drawnow
end
legend('position'); title('NLS')
xlabel('x'); ylabel('q');

%{
Eventhough our solution diverges, we present it to you (right before it diverges) : 

![NLS_solution](divergence_nl_schrodinger.png)

%}

%}
