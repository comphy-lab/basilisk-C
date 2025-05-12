
clear all; clf

% parameters
N=100; % number of gridpoints
L=10; % domain length
mu=0.01; % 
dt=0.05; % time step
x0=L/2; %x-coordinate at initial time
l0=0.5; length width

% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;   
dxx=DM(:,:,2)*scale^2;   

x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 

% system matrices
E=I;
A=mu*dxx; 

% boundary conditions
loc=[1;N];
E(loc,:)=0;
A(loc,:)=[I(1,:)-I(N,:); dx(1,:)-dx(N,:)]; 

% march in time matrix 
Mm=(E-A*dt/2);
M=Mm\(E+A*dt/2);

% initial condition
%q=1*0+0.1*[x*0+exp(-((x-x0)/l0).^2)]; 
q=0.1*sin(2*pi*x/L);

% marching loop 
quit=0;
while ~quit
    nl=0.1*q.*(1-q);
    nl(loc)=0;
    
    qnext=M*q+Mm\nl*dt; % one step forward
    
    % stop when nothing moves
    norm(qnext-q,2)<0.01
    q=qnext;
    
    
    % plotting
    %subplot(1,2,1);
    plot(x,q,'b');    
    axis([0,L,-2,2])    
    drawnow
    
end
legend('position'); title('Vibrating string')
xlabel('x'); ylabel('f');
