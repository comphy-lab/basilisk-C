%{ 
# Differential equation - Two boundaries conditions on the same side

We solve this differential equation
$$
f_{xx}= 1
$$
with the boundary conditions 
$$
f(0)=0, ~f_x(0)=1
$$

    %}
% Solve the system : DD*f=b with the two boundaries conditions on the same
% side. f(0)=0 (Dirichlet) & f'(0)=0 (Neumann)

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

    % boundary conditions
    I=eye(N);
    A([1,N],:)=[I(1,:); D(1,:)];
    b=1+zeros(N,1); b([1,N])=[0,0];

    % solve the system
    f=A\b;

    % plotting
    plot(x,f,'b.-',x,x.^2/2,'r.');
    xlabel('x');ylabel('f');
    legend('numerical','analytical')
    xlim([0,L]); grid on

    set(gcf,'paperpositionmode','auto')
    print('-dpng','-r80','differential_equation_sameside.png')
    
    



%{
![Here is the figure](/sandbox/easystab/differential_equation_sameside)
%}

