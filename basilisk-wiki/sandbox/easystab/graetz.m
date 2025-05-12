%{
# Graetz solution

This is the thermal entrance problem. A classical problem in thermal flow: internal foreced convection. We have a Poiseuille flow of a given temperature in a pipe which wall is at that temperature, and at $x=0$ there is a sudden change of the pipe wall to another temperature. We are interested in the description of the way the temperature of the fluid will adjust
 to that of the pipe. We assume the solution in the form of a series of functions of separated variables, and applying the equation to each term of this series gives a set of eigenvalue/eigenvector solutions. This is what we solve here.


%}

format compact
clear all; clf
% if octave
%setenv("GNUTERM","X11")
clf

%%%% parameters
Ny=100; %number of grid nodes in y
Ly=1;   % domain size in y

%%%% chebychev in y
scale=-Ly/2;
[y,DM] = chebdif(Ny,4);
y=(y-1)*scale;
dy=DM(:,:,1)/scale; 
dyy=DM(:,:,2)/scale^2;
I=eye(Ny);

%{
$$ u_{Pois}\frac{\partial T}{\partial x}   =   ( 
\frac{  1}{  r}\frac{\partial }{\partial r}r\frac{\partial T}{\partial r})
$$
with $T(x,1)=0$ and $\partial_r T(x,0)=0$.
in $x=0$, $T(0,r)=1$


We look at the solution in term of a series of a function of $x$ times a function of $r$; it is straightforward to see that 
by the $x$ function is an exponential, which must decrease
$$T =\Sigma e^{-\alpha_i^2 x } c_i\theta_i(r)$$ 
and $u_{Pois}=(1-r^2)$, 
so coefficients $c_i$ will be obtained latter by $1= \Sigma_i c_i e^{-\alpha_i^2 x }  \theta_i(r)$

The problem to solve is 
$$  
s (1-r^2)\theta_i(r)=(\frac{\partial^2 }{\partial r^2} + \frac{  1}{  r}\frac{\partial }{\partial r})\theta_i(r) 
$$
were $s=-\alpha_i^2$. 
%}

% System matrices
A=dyy+diag(1./y)*dy;
E=-diag(1-y.^2);

%{
We can rewrite this system of equations in a matrix form
$$
s_i E\theta_i =A \theta_i
$$


$(\frac{\partial^2 }{\partial r^2} + \frac{  1}{  r}\frac{\partial }{\partial r})$
is A, and $-(1-r^2) I $ is E 

   

 
attention 
$y(0)=0$ and $y(Ny)=1$

$\theta_i(r=1)=0$ and $\theta_i'(  r=0)=0$

r=0 is y=0 so position 1 and r=1 is y=1 position Ny

and with the boundary conditions $\theta_i(1)=0$ and $\partial_r \theta_i(0)=0$.
%}

% boundary conditions 
loc=[1,Ny];
C=[dy(1,:); I(Ny,:)];
E(loc,:)=0;
A(loc,:)=C;   

%{
# compute eigenvalues
%}

[U,S]=eig(A,E);
s=diag(S);
rem=abs(s)>1e5; s(rem)=[];U(:,rem)=[];
[t,o]=sort(real(s)); s=s(o);U=U(:,o);

%{
Note that the eigen value have been sorted.

# show Eigen functions
%}

plot(y,U(:,1:5),'-');legend('1','2','3','4','5');

%figure
%plot(y,dy*U(:,1:5),'-');legend('1','2','3','4','5');
% figure
% plot(y,dyy*U(:,1:5),'-');legend('1','2','3','4','5');

%{
#Results
 
Here the five first values of $\alpha_i$   

    2.7044 
    6.6790
    10.6734
    14.6711
    18.6699
 
%}
disp('the first eigen values ');
sqrt(s(1:5))


%{
# Links



# Exercises
* Do the same in plane 2D
* Do the same for the unsteady heat equation 
$$ \frac{\partial T}{\partial t}   =   ( 
\frac{  1}{  r}\frac{\partial }{\partial r}r\frac{\partial T}{\partial r})
$$
with $T(t,1)=0$ and $\partial_r T(t,0)=0$.
in $t=0$, $T(0,r)=1$

 
#Bibliography

* A. Leontiev (1985) ”Théorie des échanges de chaleur et de masse”, ed. MIR.

%}
