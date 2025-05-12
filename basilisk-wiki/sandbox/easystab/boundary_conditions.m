%{
# How to impose boundary conditions

This code is built upon [diffmat.m]() where the differentiation matrices
are built. Hee we solve a linear non-homogeneous differential 
equation in 1D with non-homogeneous boundary conditions
%}

clear all; clf

% parameters
L=2*pi; % domain length
N=25; % number of points

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


%{
# The differential equation

We solve this
$$
f_{xx}=1
$$
with the boundary conditions 
$$
f(1)=0, f(L)=1
$$

%}

% I build the linear differential equation
A=D;

% second order derivative
A=DD;

% boundary conditions
I=eye(N);
A([1,N],:)=I([1,N],:);
b=1+zeros(N,1); b([1,N])=[0,0];

% solve the system
f=A\b;

% plotting
plot(x,f,'b.-',x,x.^2/2-pi*x,'r-');
xlabel('x');ylabel('f');
legend('numerical','theory')
xlim([0,L]); grid on

set(gcf,'paperpositionmode','auto')
print('-dpng','-r75','boundary_conditions.png')
%{


# Exercices/Contributions

* Please put the two boundary conditions on the same side of the domain
(for instance, impose the value of $f$ and its derivative at $x=0$)
* Please solve $f_{xx}=\sin(x)$
* Please check the convergence with the number of gridpoints
* Please use Chebychev differentiation matrix [chebdif.m]() and compare with the finite
differences
* Please use a non-constant coefficient i your differential equation
$a(x)f_{xx}=0$ (first find the way to code that using differentiation
matrices)
%}