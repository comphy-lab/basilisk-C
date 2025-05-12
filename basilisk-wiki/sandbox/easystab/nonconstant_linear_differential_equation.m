clear all; 
clf

%{ 
# Resolution of non constant linear differential equations

This program go along with

1 [Classical differential equation](../classical_differential_equation.m) which is the base

3 [Burgers](../burgers_global.m) which is a very good code for all technics in non-linear functions

The resolution is exactly the same as the constant linear equation, we just
replace the coefficient $\tau$ by a matrix in wich $\tau (t)$ is on the
eye.

Here we consider the equation

$$ \frac{\partial y}{\partial t}+t y=0 $$
With $y(0)=1$, $\frac{\partial y}{\partial t}(0)=0$,
The solution is $y=e^{-x^2/2}$
%}

%%%%%%%%%%%%%%%%%%%%%%PARAMETERS OF THE MODEL%%%%%%%%%%%%%%%%%%%
% numerical parameters
Tmax=3; % domain length
res=zeros(10,1);
resvect=linspace(100,1000,10)';
for a=1:10
    N=100*a
 % number of points
t=linspace(0,Tmax,N)'; %DON'T CHANGE THIS

%vector of the excitation, choose an expression for b
b=zeros(N,1); 

%Nonconstant vectors
matt=diag(t); %We put the coefficient t into a matrix

%Boundary condition
type1= 1 ; %1 dirichlet, 2 Neumann, 3 f''
place1= 0 ; % in the physical dimension (we calculate the indice after)
val1=1;
type2= 2 ; % 1 dirichlet, 2 Neumann, 3 f''
place2= 0 ; % in the physical dimension (we calculate the indice after)
val2=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=t(2)-t(1); % the grid size
indice1=floor(place1/h)+1;
indice2=floor(place2/h)+1;

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

%%%%%%%%%%%%%%%%%%%%%%%Creation of the matrix%%%%%%%%%%%%%%%%%%%%%%%%%%
I=eye(N);
A=matt+D; %The Identity is now replaced by matt (for matricetime)


%%%%%%%%%%%%%%%%%%%%%%%%%%initial conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if type1==1, A(1,:)=I(indice1,:) ; end
 if type1==2, A(1,:)=D(indice1,:) ; end
 if type1==3, A(1,:)=DD(indice1,:); end
 
 if type2==1, A(N,:)=I(indice2,:) ; end
 if type2==2, A(N,:)=D(indice2,:) ; end
 if type2==3, A(N,:)=DD(indice2,:); end
 
 %{Attribution of the values}%
 b([1,N])=[val1,val2];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%Resolution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U=A\b;

%%%%%%%%%%%%%%%%%%%Theorical solution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=exp(-0.5*t.^2);


res(a)=log(norm(U-u)/N)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TRACE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(resvect,res,'B')
xlabel('number of points');ylabel('power of ten of the error');title('maximum of the error with the theory');
        set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','Nocnstanterror.png')


figure(2)
plot(t,U,'R*')
hold on
plot(t,u,'B.')
xlabel('t');ylabel('y');title('Theory and numerical results for dy/dt+ty=0');
hold off
        set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','Noconstantresult.png')


%{
# Figures

Here is the result for 1000 points (although 20 are enough) :
![alt text](/Noconstantresult.png)

Here is the behaviour of the exponential of the error :

![alt text](/Nocnstanterror.png)

# Improve !
You may improve it by :

1. Doing some classic case : legendre polynoms, Gamma function...

2. Using this technic on bigger code, like Brusselator with different coefficient. You may observe the different pattern if the coefficients change in the space.
%}

