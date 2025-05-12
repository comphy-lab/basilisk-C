
clear all; clf

% parameters
N=200; % number of gridpoints
L=70; % domain length
c=1; % wave velocity
dt=0.2; % time step
tmax=50; % final time
U=0.5;   % speed value
nu=0.2   %value of the kinematic viscosity
orig=0.5; % time in the past when the dirac started

% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 
II=eye(2*N);


% system matrices
E=I;
A=nu*dxx-dx*U;

% locations in the state
v=1:N;
f=N+1:2*N;

%{
# Boundary conditions
The original system has one variable and two time derivatives, so we should apply two boundary conditions. We will simply say that the position of the string at its both ends should be zero: two attachment points. This is a homogeneous Dirichlet condition applied on the first and last point of the state vector for the state position.
%}

% boundary conditions
E(1,:)=0;
A(1,:)=I(1,:);
 E(N,:)=0;
 A(N,:)=I(N,:);
         


% march in time matrix 
M=(E-A*dt/2)\(E+A*dt/2);

%{
# Initial condition
We need to describe the state of the system at initial time. We say here that the string is initially deformed as a dirac diffused at a time t0+dt (this satisfies the boundary conditions), and that the velocity is zero. From that state we then perform a loop to repeatedly advance the state of a short forward step in time. We use *drawnow* to show the evolution of the simulation as a movie when running the code. We store for validation the string position at the midle of the domain, to do this without worrying about the the grid points are, we interpolate $f$ with the function *interp1*. 
%}

% initial condition
t0=orig;
q=(exp(-(((x-L/4).^2)/(4*nu*t0))))/(2*sqrt(pi*nu*t0));
qO=(exp(-(((x-L/4).^2)/(4*nu*t0))))/(2*sqrt(pi*nu*t0));

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

for ind=1:Nt    
    t=ind*dt;
    q=M*q; % one step forward, numerical solution
    e(ind)=interp1(x,q,L/2); % store center point
    
    % plotting and validation
    qff=(exp(-(((x-U*t-L/4).^2)/(4*nu*(t+orig)))))/(2*sqrt(pi*nu*(t+orig))); %analytical solution
    qf=exp(-((x-t*U-L/4)/0.5).^2); % only the advection of the initial condition
    subplot(1,2,1);  plot(x,q,'b',x,qff,'r--',x,qO,'g'); 
  
    axis([0,L,-2,2])    
    
    subplot(1,2,2);
    plot(t,norm((q-qff)/max(q),2),'b*'); grid on;hold on
    
    drawnow
end
legend('position'); title('Vibrating string')
xlabel('x'); ylabel('f');





