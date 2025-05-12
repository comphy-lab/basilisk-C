%{
# vibrating string with dissipation
%}
clear all; clf

% parameters
N=100; % number of gridpoints
L=2*pi; % domain length
c=1; % wave velocity
dt=0.5; % time step
tmax=70; % final time
coef=0.1;
rho=1.3
% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 
II=eye(2*N);

alpha=-I*coef;

% system matrices
E=[I,Z; Z,I];
A=[alpha,c^2*dxx; I, Z]; 

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
q=[zeros(N,1); sin(pi*x/L)]; 

% marching loop
tvec=1:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

s=(-coef+sqrt(coef^2-4*c^2*(pi/L)^2))/4;
filename = 'dissipation.gif'
for ind=1:Nt 
    q=M*q; % one step forward
    e(ind)=interp1(x,q(f),L/2); % store center point
    Ua=sin(pi*x/L)*exp(s*ind); %analytical solution
    e1(ind)=interp1(x,Ua,L/2); % store center point
    
 %     uncomment this part to see the evolution of the code
     % plotting
     figure(1);
     plot(x,q(f),'b',x,Ua,'g');
     legend('numerical','theory');
     title('evolution of the vibrating string with dissipation');
     axis([0,L,-2,2])
     drawnow
           
     frame = getframe(1);
     im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if ind == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05);
      end 
      
    
     
end
legend('position','analytic'); title('Vibrating string')
xlabel('x'); ylabel('f');



% time evolution of central point
figure(2)
plot(tvec,e,'b.-',tvec,e1,'r-');
legend('numerical','theory');
title('evolution of the center point')
xlabel('time'); ylabel('f(L/2)');


%{
![evolution of the vibrating string with dissipation](dissipation.gif)


![evolution of the center point](centr-point.jpg)

%}