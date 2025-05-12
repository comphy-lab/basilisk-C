clear all; clf

%{ 
Edition du 24 Janvier 2015
Le programme semble a peu pret tourner pour la solution analytique, sauf
pour le régime critique $\delta=0$. L'erreur provient sans doute des
conditions aux limites qui sont mal calculées, mais j'ai (PV) revérifié les
formules et cela semble cohérent. La partie analytique est décalée de 0.1,
c'est étrange.


L'equation differentielle est mal résolue numériquement. 

Un systeme général de conditions aux limites numériques a été proposé mais je ne sais
pas si il est vraiment exact dans la théorie. Il faudrait le verifier sur
la théorie au moins.

Il y a encore pas mal de choses qui sont codés de manières grossière,
nottament :
_Le nom des variable, qui est parfois long pour pas grand chose
_La gestion des places des conditions aux limites, il y a des divisions par
h et un +1 qui n'est pas tres "propre" en terme de code.
_L'ecriture de la gestion des conditions initiales, ainsi que du delta est
extremement lourde. Je connais mal Matlab et donc je ne sais pas comment le
coder de manière plus élégante.
_Faire un beau module de tracé serai chouette aussi
_Une fois que la solution numérique sera correcte, il faudra faire une
étude de convergence ainsi qu'une étude de la précision en fonction de l'amortissement (pour le fun). 
%}


%{ 
# Resolution of classical linear differential equations

In mechanics or electrokinetic physics, we often face equations like 
$$
\frac{1}{\omega_0^2} \frac{\delta^2 U}{\delta t^2}+ \tau \frac{\delta U}{\delta t}+ U =
f(t)
$$

(For example, a RLC dipole or an association of a mass, a spring and a
valve)

Theory on theses equation is well-known for f(t)=U_0, takin a solution with
the form $e^{rt} wich imply to resolve a second order system

$$ \Delta= \tau^2 - \frac{4}{ \omega_0^2}$$

if $$\Delta > 0 $$ then $$U(t)=A e^{\lambda_1 t}+B e^{\lambda_2 t}$$
if $$\Delta = 0 $$ then $$U(t)=(Ax+B) e^{\lambda t}$$
if $$\Delta < 0 $$ then $$U(t)=(A cos(\Omega t)+B sin(\Omega
t))e^{\Lambda t} $$

with $\lambda_1$ and $\lambda_2$ solution of the polynom, $\lambda$ unique
solution, and \Lambda, \Omega real and imaginary part of the solution.

A and B are determined by the boundary conditions. (We solve the system by
calculating coefficients and inverse the matrix 

Here we use differentiation matrix to solve this equation with different valuers of $\omega_0$ and $\tau$, for a
constant sollicitation and any kind of sollicitation.

We create ours parameters 
%}

%%%%%%%%%%%%%%%%%%%%%%PARAMETERS OF THE MODEL%%%%%%%%%%%%%%%%%%%
% numerical parameters
Tmax=10; % domain length
N=100; % number of points

% physical parameters
omega0=pi;
tau=0.1;

%vector of the excitation
b=1+zeros(N,1); 

%Boundary condition
type1= 1 ; %1 dirichlet, 2 Neumann, 3 f''
place1= 0 ; % in time dimension
val1=1;
type2= 2 ; % 1 dirichlet, 2 Neumann, 3 f''
place2= 0 ; % in dimension
val2=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the grid
t=linspace(0,Tmax,N)';
h=t(2)-t(1); % the grid size
place1=floor(place1/h)+1;
place2=floor(place2/h)+1;

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
We have 
$$ \frac{1}{\omega_0^2} \frac{\delta^2 U}{\delta t^2}+ \frac{1}{\tau \omega_0^2}  \frac{\delta U}{\delta t}+ U =
f(t)$$

So it's equivalent to
 $$\frac{1}{\omega_0^2} DD*U+ \frac{1}{\tau \omega_0^2} D*U+U=f(t)$$

Or 
 $$U(t)=\frac{f(t)}{1+\frac{1}{\omega_0^2} DD+ \frac{1}{\tau \omega_0^2}
 D}$$
%}

%%%%%%%%%%%%%%%%%%%%%%%Creation of the matrix%%%%%%%%%%%%%%%%%%%%%%%%%%
A=1+tau*D+(omega0^(-2))*DD;
I=eye(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%initial conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
We need two conditions to solve the problem }
Condition on the value of one point : A(a,:)=I(a,:)
Condition on the first derivate of one point : D(a,:)=I(a,:)
Condition on the second derivate of one point : DD(a,:)=I(a,:)

%}
 if type1==1, A(place1,:)=I(place1,:) ; end
 if type1==2, A(place1,:)=D(place1,:) ; end
 if type1==3, A(place2,:)=DD(place1,:); end
 
 if type2==1, A(place2,:)=I(place2,:) ; end
 if type2==2, A(place2,:)=D(place2,:) ; end
 if type2==3, A(place2,:)=DD(place2,:); end
 
%{
 2 conditions on the same point :
You do the first condition as before and the second on the last line of
your matrix A.
%}
 if place1==place2 
  if type1==1, A(place1,:)=I(place1,:) ; end
  if type1==2, A(place1,:)=D(place1,:) ; end
  if type1==3, A(place2,:)=DD(place1,:); end
  
  if type2==1, A(N,:)=I(place2,:) ; end
  if type2==2, A(N,:)=D(place2,:) ; end
  if type2==3, A(N,:)=DD(place2,:); end
 end
 
 %{Attribution of the values}%
 b([place1,place2])=[val1,val2];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%Resolution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U=A\b;

%%%%%%%%%%%%%%%%%%%%%%%%%%Analytical solution%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
You can find a good summary of the theory here : http://www.physagreg.fr/electrocinetique-3-rlc.php
First we calculate the different physical pertinent values. 
%}

alpha=tau*omega0/2
Q=1/(tau*omega0)
lambda=alpha*omega0
Delta=lambda^2-omega0^2

%{
You got a different equation for each domain of $\Delta$. For each part, we
estimate the root of the equation, and the values of fonction, first and
second derivate at the boundary conditions. then we solve the system
aA+bB=value1 cA+dB=value2, a b c d depends of the boundary condition chosen.
%}
if Delta>0
    lambda1=-lambda+sqrt(lambda^2-omega0^2);
    lambda2=-lambda-sqrt(lambda^2-omega0^2);
    
    if type1==1
    a=exp(lambda1*place1*h);
    b=exp(lambda2*place1*h);
    elseif type1==2
    a=lambda1*exp(lambda1*place1*h);
    b=lambda2*exp(lambda2*place1*h);     
    elseif type1==3
    a=lambda1^2*exp(lambda1*place1*h);
    b=lambda2^2*exp(lambda2*place1*h);    
    end    
    
    if type2==1
    c=exp(lambda1*place2*h);
    d=exp(lambda2*place2*h);        
    elseif type2==2
    c=lambda1*exp(lambda1*place2*h);
    d=lambda2*exp(lambda2*place2*h);        
    elseif type3==3
    c=lambda1^2*exp(lambda1*place2*h);
    d=lambda2^2*exp(lambda2*place2*h);           
    end
    
    Sol=[a,b;c,d]\[val1,val2]';
    A=Sol(1);
    B=Sol(2);
    u=A*exp(lambda1*t)+B*exp(lambda2*t);
end

if Delta==0
    lambda=-lambda;
    
    if type1==1
    a=place1*exp(lambda*place1*h);  
    b=exp(lambda*place1*h);   
    elseif type1==2
    a=(1+lambda*place1*h)*exp(lambda*place1*h);    
    b=lambda*exp(lambda*place1*h); 
    elseif type1==3
    a=(2*lambda+lambda^2*place1*h)*exp(lambda*place1*h);    
    b=lambda^2*exp(lambda*place1*h);
    end    
    
    if type2==1
    c=place2*exp(lambda*place2*h);      
    d=exp(lambda*place2*h);        
    elseif type2==2
    c=(1+lambda*place2*h)*exp(lambda*place2*h);     
    d=lambda*exp(lambda*place2*h);       
    elseif type3==2
    c=(2*lambda+lambda^2*place2*h)*exp(lambda*place2*h);      
    d=lambda^2*exp(lambda*place2*h);       
    end
    
    Sol=[a,b;c,d]\[val1,val2]';
    A=Sol(1);
    B=Sol(2);
    u=(A*t+B).*exp(lambda*t);
end

if Delta<0
    lambda=-lambda;
    omega=omega0*sqrt(1-alpha^2); 
    
    if type1==1
    a=cos(omega*place1*h)*exp(lambda*place1*h);
    b=sin(omega*place1*h)*exp(lambda*place1*h);
    elseif type1==2
    a=(-omega*sin(omega*place1*h)+lambda*cos(omega*place1*h))*exp(lambda*place1*h);
    b=(omega*cos(omega*place1*h)+lambda*sin(omega*place1*h))*exp(lambda*place1*h);
    elseif type1==3
    a=((lambda^2-omega^2)*cos(omega*place1*h)-2*lambda*omega*sin(omega*place1*h))*exp(lambda*place1*h);
    b=((lambda^2-omega^2)*sin(omega*place1*h)+2*lambda*omega*cos(omega*place1*h))*exp(lambda*place1*h);
    end    
    
    if type2==1
    c=cos(omega*place2*h)*exp(lambda*place2*h);
    d=sin(omega*place2*h)*exp(lambda*place2*h);        
    elseif type2==2
    c=(-omega*sin(omega*place2*h)+lambda*cos(omega*place2*h))*exp(lambda*place2*h);
    d=(omega*cos(omega*place2*h)+lambda*sin(omega*place2*h))*exp(lambda*place2*h);  
    elseif type3==3
    c=((lambda^2-omega^2)*cos(omega*place2*h)-2*lambda*omega*sin(omega*place2*h))*exp(lambda*place2*h);       
    d=((lambda^2-omega^2)*sin(omega*place2*h)+2*lambda*omega*cos(omega*place2*h))*exp(lambda*place2*h);
    end
    
    Sol=[a,b;c,d]\[val1,val2]';
    A=Sol(1);
    B=Sol(2);
    u=(A*cos(omega*t)+B*sin(omega*t)).*exp(lambda*t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TRACE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(t,U,'.');
hold on
plot(t,u);
hold off
title(sprintf('tau= %0.01f omega0= %0.01f, limite1 (%0.0f, %0.01f, %0.01f) limite 2(%0.0f ,%0.01f ,%0.01f)', tau, omega0,type1,place1*h,val1,type2,place2*h,val2))
xlabel('temps')
ylabel('amplitude')


