# Convergence study of the time step with a backward Euler scheme
 %{
Done by Alexandre Bilczewski and Fadil Boodoo. 
In this work, we do a convergence study of the time step of the code vibrating_string.m
The scheme used for the march in time is Backward Euler, using the code of sandbox/easystab/stab2014/vibrating_string_forward_backward_euler.m
%}

clear all; clf

% parameters
N=101; % number of gridpoints
L=2*pi; % domain length
c=1; % wave velocity
% dt=0.5; % time step
tmax=20; % final time

%Loop in time step
for i=1:30
   dt=i/10 ;

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
A=[Z,c^2*dxx; I, Z]; 

% locations in the state
v=1:N;
f=N+1:2*N;

% boundary conditions
loc=[f(1),f(N)];
A(loc,:)=0; 
E(loc,:)=II(loc,:);

% march in time matrix 
d=(E-A*dt);
M=inv(d)*E;


% initial condition
q=[zeros(N,1); sin(pi*x/L)]; 


% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

for ind=1:Nt    
    q=M*q; % one step forward
    e(ind)=q(f((N+1)/2)); % store center point
    
end

% time evolution of central point
subplot(1,2,1);
Ttheo=2*L/c;
tt=linspace(0,tmax,500);
etheo=cos(2*pi*tt/Ttheo);



% calculation of the error
etheoNt=cos(2*pi*tvec/Ttheo);
eetheo = etheoNt-e.' ;
dtstock(i)=dt;
error(i)=max(eetheo);

end

% plot the error
xx2=dtstock;
plot(dtstock,error,'b.-');
title('error in function of the time step');
subplot(1,2,2);

% comparison of the numerical error and the theorical error
xx=1./dtstock;
loglog(xx,error,'b.-',xx,xx.^-1,'r'); 
legend('numerical error','theorical error order 1'); title('Convergence of the theorical error and the numerical error in loglog')
grid on
xlabel('dt'); ylabel('error');

%{
![convergence time with a backward euler scheme](/vibrating_string_convergence_time_euler_backward.png)
