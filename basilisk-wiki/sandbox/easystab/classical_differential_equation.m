clear all; 
clf

%{ 
This program go along with : 

1 [Non-constant differential equation](../nonconstant_linear_differential_equation.m) which is the same with non constant coefficients

2 [Differential equation](../differential_equation.m) which is a less general version

3 [Burgers](../burgers_global.m) which is a very good code for all technics in non-linear functions


# Resolution of classical linear differential equations

In mechanics or electrokinetic physics, we often face equations like :
$$
\frac{1}{\omega_0^2} \frac{\delta^2 U}{\delta t^2}+ \tau \frac{\delta U}{\delta t}+ U =
f(t)
$$

(For example, an RLC dipole or an association of a mass, a spring and a
valve)

Theory on theses equation is well-known. for Differential equations, the
solution is with the form $$ \sum_0^n A_i e^{r_i t}$$ n is the order of the
differential equation (here , 2),$ r_i$ a root of the equation
$$\frac{1}{\omega^2} r^2 + \lambda r + 1 = 0$$

So we solve :

$$ \Delta= \tau^2 - \frac{4}{ \omega_0^2}$$

if $\Delta > 0 $ then $U(t)=A e^{\lambda_1 t}+B e^{\lambda_2 t}$

if $\Delta = 0 $ then $U(t)=(At+B) e^{\lambda t}$

if $\Delta < 0 $ then $U(t)=(A cos(\Omega t)+B sin(\Omega
t))e^{\Lambda t} $

with $\lambda_1$ and $\lambda_2$ solution of the polynom, $\lambda$ unique
solution, and \Lambda, \Omega real and imaginary part of the solution.

A and B are determined by the boundary conditions. (We solve the system by
calculating coefficients and inverse the matrix)

Here we use differentiation matrix to solve this equation with different valuers of $\omega_0$ and $\tau$, for a
constant sollicitation.
%}

%%%%%%%%%%%%%%%%%%%%%%PARAMETERS OF THE MODEL%%%%%%%%%%%%%%%%%%%
% numerical parameters
Tmax=50; % domain length
N=1000; % number of points
t=linspace(0,Tmax,N)'; %DON'T CHANGE THIS
omega0=1; 
tau=1; %the two parameters of the equation

%TYPE OF SOLLICITATION
b=zeros(N,1); solli=0; %No sollicitation, solli allow you to chose what kind of plot you will have
%b=zeros(N,1); solli=0; %No sollicitation
%b=42*ones(N,1), solli=1; %step
%b=sin(t), solli=1; %sinusoidal sollicitation
%b=exp(-t)+1, solli=1;

%BOUNDARY CONDITIONS
type1= 1 ; %1 dirichlet, 2 Neumann, 3 f''
place1= 1 ; % in the physical dimension (we calculate the indice after)
val1=1;
type2= 2 ; % 1 dirichlet, 2 Neumann, 3 f''
place2= 3 ; % in the physical dimension (we calculate the indice after)
val2=0;

for z=1:20 %This loop is made for plotting three different cases, if you use the program remove it
    
if z==1
omega0=1;
tau=0.6;
b=zeros(N,1); solli=0; %No sollicitation
type1= 1 ; %1 dirichlet, 2 Neumann, 3 f''
place1= 1 ; % in the physical dimension (we calculate the indice after)
val1=1;
type2= 2 ; % 1 dirichlet, 2 Neumann, 3 f''
place2= 3 ; % in the physical dimension (we calculate the indice after)
val2=0;
end

if z==2
omega0=1;
tau=1;
b=sin(t); solli=1; %sinusoidal sollicitation
type1= 1 ; %1 dirichlet, 2 Neumann, 3 f''
place1= 1 ; % in the physical dimension (we calculate the indice after)
val1=1;
type2= 1 ; % 1 dirichlet, 2 Neumann, 3 f''
place2= 3 ; % in the physical dimension (we calculate the indice after)
val2=0;
end

if z>2
omega0=1;
tau=1/(0.2*z+0.5);
b=zeros(N,1); solli=0; %No sollicitation
type1= 1 ; %1 dirichlet, 2 Neumann, 3 f''
place1= 1 ; % in the physical dimension (we calculate the indice after)
val1=1;
type2= 1 ; % 1 dirichlet, 2 Neumann, 3 f''
place2= 3 ; % in the physical dimension (we calculate the indice after)
val2=0;
end

%{
We create the numerical tool used for the solution
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the grid

h=t(2)-t(1); % the grid size

indice1=floor(place1/h)+1;
indice2=floor(place2/h)+1;
loc=[indice1;indice2];

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
 $$U(t)=\frac{f(t)}{1+\frac{1}{\omega_0^2} DD+\tau 
 D}$$
%}

%%%%%%%%%%%%%%%%%%%%%%%Creation of the matrix%%%%%%%%%%%%%%%%%%%%%%%%%%
I=eye(N);
A=I+tau*D+(omega0^(-2))*DD;

%%%%%%%%%%%%%%%%%%%%%%%%%%initial conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
We need two conditions to solve the problem.
We can use condition on the value (Dirichlet), on the derivate (Neumann) or
on the second derivate. We always need N equations for N unknown
parameters. So we need to replace 2 equations from the system A. We replace
the first one and the last one to assure the continuity of the solution. 

Condition on the value of one point : A(1,:)=I(a,:)
Condition on the first derivate of one point : D(1,:)=I(a,:)
Condition on the second derivate of one point : DD(1,:)=I(a,:)
Replace 1 by N for the other solution, a is the index of the condition.
Replace in b the value of the condition at the first and last place.

%}
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

%%%%%%%%%%%%%%%%%%%%%%%%%%Analytical solution%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
You can find a good summary of the theory here : http://www.physagreg.fr/electrocinetique-3-rlc.php
First we calculate the different physical pertinent values. 
%}

alpha=tau*omega0/2;
Q=1/(tau*omega0);
lambda=alpha*omega0;
Delta=lambda^2-omega0^2;

%{
You got a different equation for each domain of $\Delta$. For each part, we
estimate the root of the equation, and the values of fonction, first and
second derivate at the boundary conditions. then we solve the system
eE+bB=value1 cE+dB=value2, e b c d depends of the boundary condition chosen.
We took e and E to be sure there was no confusion between the different
notations used.
%}
if Delta>0
    lambda1=-lambda+sqrt(lambda^2-omega0^2);
    lambda2=-lambda-sqrt(lambda^2-omega0^2);
    
    if type1==1
    a=exp(lambda1*place1);
    e=exp(lambda2*place1);
    elseif type1==2
    a=lambda1*exp(lambda1*place1);
    e=lambda2*exp(lambda2*place1);     
    elseif type1==3
    a=lambda1^2*exp(lambda1*place1);
    e=lambda2^2*exp(lambda2*place1);    
    end    
    
    if type2==1
    c=exp(lambda1*place2);
    d=exp(lambda2*place2);        
    elseif type2==2
    c=lambda1*exp(lambda1*place2);
    d=lambda2*exp(lambda2*place2);        
    elseif type3==3
    c=lambda1^2*exp(lambda1*place2);
    d=lambda2^2*exp(lambda2*place2);           
    end
    
    Sol=[a,e;c,d]\[val1,val2]';
    E=Sol(1);
    B=Sol(2);
    u=E*exp(lambda1*t)+B*exp(lambda2*t);
end

if Delta==0
    lambda=-lambda;
    
    if type1==1
    a=place1*exp(lambda*place1*h);  
    e=exp(lambda*place1);   
    elseif type1==2
    a=(1+lambda*place1)*exp(lambda*place1);    
    e=lambda*exp(lambda*place1); 
    elseif type1==3
    a=(2*lambda+lambda^2*place1)*exp(lambda*place1);    
    e=lambda^2*exp(lambda*place1);
    end    
    
    if type2==1
    c=place2*exp(lambda*place2);      
    d=exp(lambda*place2);        
    elseif type2==2
    c=(1+lambda*place2)*exp(lambda*place2);     
    d=lambda*exp(lambda*place2);       
    elseif type3==2
    c=(2*lambda+lambda^2*place2)*exp(lambda*place2);      
    d=lambda^2*exp(lambda*place2);       
    end
    
    Sol=[a,e;c,d]\[val1,val2]';
    E=Sol(1);
    B=Sol(2);
    u=(E*t+B).*exp(lambda*t);
end

if Delta<0
    lambda=-lambda;
    omega=omega0*sqrt(1-alpha^2); 
    
    if type1==1
    a=cos(omega*place1)*exp(lambda*place1);
    e=sin(omega*place1)*exp(lambda*place1);
    elseif type1==2
    a=(-omega*sin(omega*place1)+lambda*cos(omega*place1))*exp(lambda*place1);
    e=(omega*cos(omega*place1)+lambda*sin(omega*place1))*exp(lambda*place1);
    elseif type1==3
    a=((lambda^2-omega^2)*cos(omega*place1)-2*lambda*omega*sin(omega*place1))*exp(lambda*place1);
    e=((lambda^2-omega^2)*sin(omega*place1)+2*lambda*omega*cos(omega*place1))*exp(lambda*place1);
    end    
    
    if type2==1
    c=cos(omega*place2)*exp(lambda*place2);
    d=sin(omega*place2)*exp(lambda*place2);        
    elseif type2==2
    c=(-omega*sin(omega*place2)+lambda*cos(omega*place2))*exp(lambda*place2);
    d=(omega*cos(omega*place2)+lambda*sin(omega*place2))*exp(lambda*place2);  
    elseif type3==3
    c=((lambda^2-omega^2)*cos(omega*place2)-2*lambda*omega*sin(omega*place2))*exp(lambda*place2);       
    d=((lambda^2-omega^2)*sin(omega*place2)+2*lambda*omega*cos(omega*place2))*exp(lambda*place2);
    end
   
    Sol=[a,e;c,d]\[val1,val2]';
    E=Sol(1);
    B=Sol(2);
    u=(E*cos(omega*t)+B*sin(omega*t)).*exp(lambda*t);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TRACE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if z==1
 figure(1)    
plot(t,U,'.');
hold on
plot(t,u,'K');
plot(t(indice1),u(indice1),'R.','markersize',20)
plot(t(indice2),u(indice2),'R.','markersize',20)
legend('numerical','theory','condition place')
title(sprintf('N= %0.f ,tau= %0.01f, omega0= %0.01f, limite1 (%0.0f, %0.01f, %0.01f) limite 2(%0.0f ,%0.01f ,%0.01f)',N, tau, omega0,type1,place1,val1,type2,place2,val2))
hold off
  set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','nosollicitation.png')

end

if z==2
figure(2)    
plot(t,U,'.');
hold on
plot(t,b,'K');
plot(t(indice1),u(indice1),'R.','markersize',20)
plot(t(indice2),u(indice2),'R.','markersize',20)
legend('answer','sollicitation','condition place')
title(sprintf('N= %0.f ,tau= %0.01f, omega0= %0.01f, limite1 (%0.0f, %0.01f, %0.01f) limite 2(%0.0f ,%0.01f ,%0.01f)',N, tau, omega0,type1,place1,val1,type2,place2,val2))
hold off
  set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','sollicitation.png')
end

if z>2
figure(3)
plot(t,u);
hold on
end

end


title('shapes of solution for the same boundarie conditions, different Q')
hold off
  set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','theoreticalsolution.png')

%{
# Figures
Here are the form for no sollicitations and different Damp coefficient

![alt text](/theoreticalsolution.png)


Here is the comparison between theory and our program 

![alt text](/nosollicitation.png)

There is a non-homogenous dirichlet condition on the first red point, a homogenous Neumann condition on the second point


Here is an example of sollicitation

![alt text](/sollicitation.png)

This is a very usefull information to combine with Bode diagramm, into having easily every information you need about electrocinetic ou mechanical system studied in preparatory classes

# Improve !

You may improve this program with :

1. Generalizing it to higher order

2. Trace Bode diagrams

3. Do a differential equation of the sollicitation, wich allow you to study a lot of mechanical system or filters

4. As this was my first program, some more powerful tools were developped afterward (like easypack). It could be useful to rewrite the program with

