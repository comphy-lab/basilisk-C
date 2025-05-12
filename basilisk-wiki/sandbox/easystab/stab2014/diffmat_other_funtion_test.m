%{
# Presentation of the differentiation matrix
If you want to know more about  differentiation matrix, please go to [./sandbox/easystab/diffmat.m](./sandbox/easystab/diffmat.m)
#
In this program we test other function to see the performence of differentiation matrix
%}

clear all; clf

% parameters
L=5; % domain length
N=15; % number of points
% the grid
x=linspace(0.1,L,N)';
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



%{
![Comparison of the numerical and exact derivative for the first and second derivative of a polynomial function](diffmat_otherfunction.jpg)

Now that the differentiation matrices are built and stored in D and DD, we
can test them for a polynomial function.
%}

% test of the derivatives
f=3*x.^(3)+4*x.^(2)+x;
fp=9*x.^(2)+8*x+1;
fpp=18*x+8;
plot(x,f,'b.-',x,fp,'b.-',x,D*f,'r*',x,fpp,'g.-',x,DD*f,'m*');
legend('f=3x^3+4x^2+x','9x²+8x+1','D*f','fpp=18x+8','DD*f')
print('-djpeg','diffmat_otherfunction.jpg'); % save the figure

% second test of the derivatives
f1=exp(4*x+1)+log(2*x.^(2))+1;
f1p=4*exp(4*x+1)+2./(x);
f1pp=16*exp(4*x+1)-2./(x.^(2));

figure(2)
semilogy(x,f1,'b.-',x,f1p,'y.-',x,D*f1,'r*',x,f1pp,'g.-',x,DD*f1,'m*');
grid on;
print('-djpeg','diffmat_otherfunction2.jpg'); % save the figure
legend('f1=exp(4*x+1)+log(2*x.^(2))+1','f1p=4*exp(4*x+1)+2./(x)','D*f1','fpp=16*exp(4*x+1)-2./(x.^(2))','DD*f1')

%{
![Comparison of the numerical and exact derivative for the first and second derivative of a exponential function](diffmat_otherfunction2.jpg)


%}