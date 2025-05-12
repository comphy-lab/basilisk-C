%{ 
# Differential equation : $f_{xx}=\sin(x)$

Solve the system : DD*f=sin(x) with the two boundary conditions f(0)=0,fx(L)=1

For solve it, we use the code build upon [differential_equation.m](./differential_equation.m).
%}
 clear all; clf

 
    % parameters
    L=2*pi; % domain length
    N=50; % number of points

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

    % second order derivative
    A=DD;
    
%{ 
We need to change boundary conditions. We put sin(x) ofr b and we change the place of the first derivative, we place D(1,:) instead D(N,:).
%}
    % boundary conditions
    I=eye(N);
    A([1,N],:)=[I(1,:); D(1,:)];
    b=sin(x); 

    % solve the system
    f=A\b;

%{ 
We also need to change the analytical solution to compare numerical and analytical solutions.
%}

    % plotting
    plot(x,f,'b.-',x,-sin(x)+x,'r.');
    xlabel('x');ylabel('f');
    legend('numerical','analytical')
    xlim([0,L]); grid on

    set(gcf,'paperpositionmode','auto')
    print('-dpng','-r80','differential_equation_sin.png')
    
    



%{
![Comparison between numerical and analytical solution for $f_{xx}=\sin(x)$](differential_equation_sin.png)

We can validate the differential equation for fxx=sin(x)
%}

%{
# Contributor's page
Link to page of contributor [Fabien](/sandbox/easystab/stab2014/fabien.m)
%}