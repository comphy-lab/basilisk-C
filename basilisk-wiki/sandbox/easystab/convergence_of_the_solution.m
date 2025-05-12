%{
# Convergence of the solution
%}

clear all; clf

% parameters
L=2*pi; % domain length
N=5; % number of points

% the grid
for z=1:100
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

I=eye(N);
DD(1,:)=I(1,:);
DD(N,:)=I(N,:);


b=sin (x); b(1)=0 ;b(N)=0 % boundary condition


f=DD\b; % calcul of the solution

error(z)=sum(abs(f+sin(x)))%error


N=N+1

end
figure(1)   % compare the theorical solution and the numirical solution
plot (f,'b')
hold on
plot (-sin(x),'r')

figure(2)
Nvec=5:104; 
errortheo=1./(Nvec.^1);
errortheo=errortheo/errortheo(end)*error(end);
plot(Nvec,error,'b.-',Nvec,errortheo,'r--') %compare the theorical error and the numerical error
legend('error','errortheo')
title('convergence of the theorical error and the numerical error')
xlabel('N');
ylabel('error')

figure(3)
loglog(Nvec,error,'b.-',Nvec,errortheo,'r--') %compare the theorical error and the numirical error in loglog
legend('error','errortheo')
title('convergence of the theorical error and the numerical error in loglog')
xlabel('N');
ylabel('error')

%{
![The results](/sandbox/easystab/convergence)
%}