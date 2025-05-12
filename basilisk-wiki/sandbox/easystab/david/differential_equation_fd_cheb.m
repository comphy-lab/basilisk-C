%{

# Solving a linear differential equation : comparison between finite-differences and chebychev discretizations.

The code solves the differential equation with boundary conditions:
$$
\frac{d^2 F}{d x^2} = 1 ; \quad F(0) = 0 ; \quad
\left.\frac{dF}{dx}\right|_L = 1.
$$
Whose theoretical solution is : $F(x) = x^2/2 + (1-L) x$.



This code is an improvement of the previous code
[differential_equation_chebychev.m](http://basilisk.fr/sandbox/easystab/differential_equation_chebychev.m), but
uses dif1D.m to build the matrices.
In addition it also demonstrate the way to compute an integral of a function using the "weight" vector.


%}

    clear all; clf

    % parameters
    L=2*pi; % domain length
    N=25; % number of points

%{ 
## Solution with finite elements
%}
    %% building the differentiation matrices with dif1D:
    [dxfd,dxxfd,wxfd,xfd]=dif1D('fds',0,L,N);
    Z=zeros(N,N); I=eye(N);
  
     %% I build the linear differential equation
    Afd=dxxfd;
    
    %% boundary conditions
    Afd([1,N],:)=[I(1,:); dxfd(N,:)]; 
    b=1*ones(N,1); b([1,N])=[0,1];
   
 %% solve the system
    f=Afd\b;
   
%% plotting the results and comparing with theory
    subplot(1,2,1)
    plot(xfd,f,'g-+',xfd,xfd.^2/2+(1-L)*xfd,'r-');
    xlabel('x');ylabel('f');
    legend('numerical','theory')
    xlim([0,L]); grid on
    title('Finite Differences');
    
%{ 
## Solution with Chebychev
%}
    %% building the differentiation matrices with dif1D:
   [dxcheb,dxxcheb,wxcheb,xcheb]=dif1D('cheb',0,L,N);
    Acheb=dxxcheb;
    
    %% boundary conditions are treated exactly in the same way 
    Acheb([1,N],:)=[I(1,:); dxcheb(N,:)]; 
    b=1*ones(N,1); b([1,N])=[0,1];
    
    %% solving the problem
    g=Acheb\b;

    %% plotting the Chebychev differentiation way
    subplot(1,2,2)
    plot(xcheb,g,'b-+',xcheb,xcheb.^2/2+(1-L)*xcheb,'r-');    
    xlabel('x');ylabel('f');
    legend('numerical','theory')
    xlim([0,L]); grid on
    title('Chebychev');
    
   

    set(gcf,'paperpositionmode','auto')
    print('-dpng','-r75','differential_equation_fd_cheb.png')
    
 %{ 
## Computing an integral
    
    We want to compute the integral $I = \int_0^L f(x) dx$ where $f(x)$ is 
    the solution of the previous problem.
    
    Theoretical solution is $I = L^3/6+(1-L)*L^2/2$
 %}

    disp(' Computing the integral :') 
    Itheo = L^3/6+(1-L)*L^2/2
    Ifd = wxfd*f
    Icheb = wxcheb*g
    %}
    
    
    
    %{
    ![**Figure :** The results. Note that the Chebyshev disctetization results in a clustering of the point along the boundaries of the interval.  ](differential_equation_fd_cheb.png)
    %}