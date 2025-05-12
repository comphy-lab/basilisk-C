%{
Comparison between the second derivative computed by multiplication twice by the first derivative matrix or once by the second derivative matrix
%}
clear all; clf

% parameters
L=2*pi; % domain length
N=30; % number of points
dt=0.02; % time step
tmax=10; % final time

% the grid
x=linspace(0,L,N)';
h=x(2)-x(1); % the grid size

% first derivative
D=zeros(N,N);
D(1,1:3)=[-3/2, 2, -1/2]/h;
for ind=2:N-1
    D(ind,ind-1:ind+1)=[-1/2, 0, 1/2]/h;
end
D(end,end-2:end)=[1/2, -2, 3/2]/h;

% second derivative
DD=zeros(N,N);
DD(1,1:3)=[1, -2, 1]/h^2;
for ind=2:N-1
    DD(ind,ind-1:ind+1)=[1, -2, 1]/h^2;
end
DD(end,end-2:end)=[1, -2, 1]/h^2;

tvec=dt:dt:tmax;
Nt=length(tvec);

for ind=1:Nt;
    t=ind*dt;
    ftheo=-sin(x);
    f1=D*D*sin(x);
    f2=DD*sin(x);
    err1=abs(f1-ftheo);
    err2=abs(f2-ftheo);
    
    % plotting
    subplot(1,2,1);
    plot(x,ftheo,'g',x,f1,'b',x,f2,'r'); 
    axis([0,L,-1,1]);
    xlabel('x');ylabel('f');
    legend('-sin(x)','D*D*sin(x)','DD*sin(x)');
    
    subplot(1,2,2);
    plot(x,err1,'b',x,err2,'r'); 
    axis([0,L,0,0.25]);
    xlabel('x');ylabel('error');
    legend('f1-ftheo','f2-ftheo');
    drawnow;
end
%{
![](/diffmat_test.jpg)
%}