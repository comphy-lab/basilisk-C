   %{

# Solving a linear second order differential equation
   %}

    clear all; clf

    % We are going to resolve a linear second order differential equation : 
    % Loop on number of points to highlight the evolution of the error's norm.
    
    L=10;
    for N=3:100
    % We Create the grid
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


    % We define the coefficients of the equation :

    a=1;
    b=1;
    c=0*ones(N,1);

    % We factorize f in the equation to create the matrix A, which is necessary to resolve the system

    A=a*D+b*DD;
    I=eye(N);
    A([1,N],:)=I([1,N],:);
    c([1,N])=[1,0];
    f=A\c;
    
    e(N)=norm(f-exp(-x),2);
    
    end
    
    % Plotting the results : 
    
    figure(1)
    subplot(1,3,1)
    plot(x,f,'b.-',x,exp(-x),'r.')
    title('Numerical and Analytical solutions')
    legend('numerical','theory')
    subplot(1,3,2)
    plot(x,f-exp(-x))
    title('Difference between analytical and numerical')
    subplot(1,3,3)
    plot(1:N,e)
    title('Norm error evolution')


    set(gcf,'paperpositionmode','auto')
    print('-dpng','-r80','differential_equation_secondorder.png')
    
    
%{
![Here is the figure](differential_equation_secondorder.png)
%}

