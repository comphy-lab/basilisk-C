%{
# Vibrating string : One wavelength in one domain length

We adapted the program so the user can choose the mode 'n'. Here is the mode 3
   %}
% Modified by Phan & Boodoo


clear all; clf

% parameters
N=50; % number of gridpoints
L=2*pi; % domain length
c=1; % wave velocity
dt=0.02; % time step
tmax=20; % final time
n=3; % mode number

% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 
II=eye(2*N);

% system matrices
E=[I,Z; Z,I];
A=[Z,I; c^2*dxx,Z]; 

% locations in the state
v=1:N;
f=N+1:2*N;

% boundary conditions
loc=[f(1),f(N)];
E(loc,:)=0; 
A(loc,:)=II(loc,:);

% march in time matrix 
M=(E-A*dt/2)\(E+A*dt/2);

% initial condition
q=[zeros(N,1); sin(2*n*pi*x/L)]; 

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

for ind=1:Nt    
    q=M*q; % one step forward
    e(ind)=interp1(x,q(f),L/(4*n)); % store center point
    
    % plotting
    subplot(1,2,1);
    plot(x,q(f),'b',x,q(v),'r--');    
    axis([0,L,-2,2])    
    drawnow
end
legend('position','velocity'); title('Vibrating string')
xlabel('x'); ylabel('f');

% time evolution of central point
subplot(1,2,2);
Ttheo=L/n*c;
tt=linspace(0,tmax,500);
etheo=cos(2*pi*tt/Ttheo);
plot(tvec,e,'b.-',tt,etheo,'r-');
legend('numerical','theory');
xlabel('time'); ylabel(['f(L/',num2str(4*n),')']);
