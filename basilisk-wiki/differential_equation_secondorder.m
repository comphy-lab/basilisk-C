
%{
# Let's solve the differential equation af'+ bf" = c

Let's set g = f' ==> bg' + ag = c
Then, the general solution is g(x) = Aexp{\frac{-ax}{b}}
With the particular solution g = constant <=> g = c/a, we obtain the exact
solution :
$$
f = {\frac{-bA}{a}}exp{\frac{-ax}{b}} + c{\frac{x}{a}} + B, A and B being constants
$$

For the solving, we will choose a = b = c = 1

This code is based upon the differentiation matrix code diffmat.m
%}

for N=5:1:20 % Number of points
    
% Parameters
L=2*pi; % Domain length

% The grid
x=linspace(0,L,N)';
h=x(2)-x(1); % The grid size

% First derivative
D=zeros(N,N);
D(1,1:3)=[-3/2, 2, -1/2]/h;
for ind=2:N-1
    D(ind,ind-1:ind+1)=[-1/2, 0, 1/2]/h;
end
D(end,end-2:end)=[1/2, -2, 3/2]/h;

% Second derivative
DD=zeros(N,N);
DD(1,1:3)=[1, -2, 1]/h^2;
for ind=2:N-1
    DD(ind,ind-1:ind+1)=[1, -2, 1]/h^2;
end
DD(end,end-2:end)=[1, -2, 1]/h^2;

I=eye(N);

A = D+DD;
A([1,N],:)= I([1,N],:);
c=ones(N,1);

% Boundary Conditions
c(1)=-1; c(N)=-exp(-2*pi)+2*pi;

% Solve the system
f2=A\c;

% Convergence
e=sum(abs(f2-(-exp(-x)+x)));

figure(1)
loglog(N,e,'r*')
hold on
 
end
grid on
hold on 
title('Error vs Number of points'); xlabel('N'); ylabel('error')
print('-djpeg','fig1.jpg'); % save the figure

% Comparision between the numerical solution and theory
%{
![Here is the figure](/sandbox/easystab/differential_equation_secondorder)
%}
figure(2)
plot(x,f2,'b-',x,-exp(-x)+x,'r.-')
title('Comparision between the numerical solution and theory'); xlabel('x'); xlabel('f')
legend('numerical','theory')
grid on
print('-djpeg','fig2.jpg'); % save the figure





%{
Code by : Larry, Nathasha, Charline & Antoine
%}