%{
Like the [diffmat](../diffmat.m) code, the principle here is to determine the matrix differentiation to have the third derivate of a function. The discretization of the function look like that : 

$$u'''=\frac{-u_{i-2}+2u_{i-1}+0u_{i}-2u_{i+1}+u_{i+2}}{2h^3}$$

We see that is a centered stencil and to have the first and last part of the matrix we have to use whether a upwind stencil whether a backwind stencil. Thus we have :

$$u'''=\frac{-u_{i}+3u_{i+1}-3u_{i+2}+u_{i+3}}{h^3}$$

The second and before last line of the matrix are here to refine the accuracy of the differentiation matrix.
%}

clear all; clf

% parameters
L=2*pi; % domain length
N=30; % number of points

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

% third derivative
DDD=zeros(N,N);
DDD(1,1:4)=[-1, 3, -3, 1]/h^3;
DDD(2,1:4)=[1, -1, 1, -1]/(40*h^3);
for ind=3:N-2
    DDD(ind,ind-2:ind+2)=[-1/2, 1, 0, -1, 1/2]/h^3;
end
DDD(end-1,end-3:end)=[1, -1, 1, -1]/(40*h^3);
DDD(end,end-3:end)=[-1, 3, -3, 1]/h^3;

plot(x,sin(x),'k.-',x,DDD*cos(x),'g.--',x,-cos(x),'b.-',x,DDD*sin(x),'r.--');
legend('-cos','DDD*cos','sin','DDD*sin')

%{
![](/Third_order.jpg)
%}