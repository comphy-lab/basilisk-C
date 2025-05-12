%{

In this code we demonstrate the integration in 1D of a function. To see this in 2D please see [integration_2D.m](). To see how the integration weights are coded, please see [dif1D.m]().

We integrate $f(x)=x^2$ from 0 to 1, the exact value is,1/3.
%}

clear all; clf

% parameters
L=1; % domain length
N=15; % number of points

%{
# Trapezoidal integration
We start with a simple equispaced grid $x$ with grid spacing $h$. A basic way to integrate is to approximate the function between each grid nodes with a straight line. The integral is thus 
$$
\int_{x_1}^{x_N} f dx =h\left(\frac{f_1}{2}+f_2+f_3+ \dots +f_{N-1}+\frac{f_N}{2}\right).
$$

![trapezoidal integration](trapz.png)

We continue with the idea of building linear operators to do the operations in easystab, instead of manipulating the individual elements of the variables. Thus we build the line vector of the integration weights
$$
w=h\begin{pmatrix}\frac{1}{2}, & 1, & 1, & \dots & 1, & \frac{1}{2} \end{pmatrix}
$$
and the integral is thus computed by the vector product w*f.
%}

% equispaced grid
x=linspace(0,L,N)'; % the grid
f=x.^2; % the function
h=x(2)-x(1); % grid spacing

wx=h*[0.5,ones(1,N-2),0.5];
intf=wx*f; % integration

errtrap=abs(intf-1/3) % the error


%{
# Integration for the Chebychev grid
When using Chebychev polynomial we have a grid which is not equispaced, and since the function is apprximated with polynomilas we may reach a higher accuracy than using the naive trapezoidal integration. The work was done by Clenshaw and Curtis and the integration weights are coded in [clencurt.m](). 
%}
% Chebychev mesh and Clenshaw-Curtis integration weights
[dx,dxx,wx,x]=dif1D('cheb',0,L,N,3);

f=x.^2; % the function
intf=wx*f; % integration

errcheb=abs(intf-1/3) % the error


%{
# Validation
And here is the screen output for the approximation error. We see that with just 15 grid cells we get quite a nice accuracy with the Chebychev grid.

    errtrap =
       8.5034e-04
    errcheb =
       5.5511e-17
   
# Exercices/contributions

* Please
* Please
* ...

%}
