%{
# March in time of the vibrating string: $f_t+Uf_x=0$
In this code, we do the marching in time of the advection equation $f_{t}+Uf_{x}=0$ with just one boundary condition at the left end by changing the code given [vibrating_string.m](). 

%}
clear all; clf

% parameters
N=100; % number of gridpoints
L=10; % domain length
U=2; % wave velocity
dt=0.01; % time step
tmax=10; % final time
x0=L/8; %x-coordinate at initial time
l0=0.5; %length width

% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;        
x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 

%{
# System matrices

Here we build the matrices that give the discretization of the equations. The equations is
$$
f_{t}+Uf_{x}=0
$$ 
This is an equation with a single derivative, we can write this equation in
the following way
$$
f_{t}=-UD_xf
$$ 

We can put an identity matrix on the left of the equation, in this way, we can impose the
boundary condition more easliy.

So the equation becomes :
$$
If_{t}=-UD_xf
$$ 

Now we have a linear system with the form
$$
Ef_{t}=Af
$$ 
%}


% system matrices
E=I;
A=-U*dx; 

% boundary conditions
E(1,:)=0;
A(1,:)=I(1,:);
         
%{
# March in time

Here we use the same method as [vibrating_string.m]().
%}

% march in time matrix 
M=(E-A*dt/2)\(E+A*dt/2);

%{
# Initial condition
We need to describe the state of the system at initial time. We say here that the string is initially deformed as a bell curve, and that the velocity is zero. From that state we then perform a loop to repeatedly advance the state of a short forward step in time. We use drawnow to show the evolution of the simulation as a movie when running the code. We store for validation the string position at the midle of the domain, to do this without worrying about the the grid points are, we interpolate 
f
 with the function interp1.
%}

% initial condition
q=[x*0+exp(-((x-x0)/l0).^2)]; 


% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

for ind=1:Nt    
    q=M*q; % one step forward
    e(ind)=interp1(x,q,L/2); % store center point
    
    % plotting
    subplot(1,2,1);
    plot(x,q,'b');    
    axis([0,L,-2,2])    
    drawnow
end
legend('position'); title('Vibrating string')
xlabel('x'); ylabel('f');

%{
# Validation
In order to verify the solution, we draw the position of the point in the middle of L at different times.
%}

% time evolution of central point
subplot(1,2,2);
tt=linspace(0,tmax,500);
etheo=exp(-((L/2-U*tt-x0)/l0).^2);


plot(tvec,e,'b.-',tt,etheo,'r-');
legend('numerical','theory');
xlabel('time'); ylabel('f(L/2)');
%{

<center>
![validation](/advection_1d.gif)
![validation](/sandbox/easystab/stab2014/vibrating2.png)
</center>
%}
%{
# ZHAO's contribution
\ If you are interseted in advection 2D, I have also a programme [advection_2D.m](/sandbox/easystab/advection_2D.m)

\ You can also link to my contribution page [zhao.m](/sandbox/easystab/stab2014/zhao.m)
%}