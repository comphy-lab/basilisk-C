%{
# Testing the domain geometry derivative

This is basically the same as [domain_derivative_1D.m]() except that here we adapt the domain accordingly to the solution at each Newton iteration. This means that it converges to the true solution of the problem instead of approximating the function by its Taylor expansion outside the computational domain.

%}

clear all; clf;

% parameters
N=50;           % number of grid points
p=-1.21;         % desired slope at L
L=1;            % initial guess of the domain size
method='fd';    % the differentiation

% useful
[d.x,d.xx,d.wx,x]=dif1D(method,0,L,N,5);
Z=zeros(N,N);
I=eye(N);

%initial guess
sol0=[x.*(1-x);0];
sol=sol0;

% Newton iterations
quit=0;count=0;
while ~quit

    % the present solution and its derivatives
    l=sol(N+1); h=sol(1:N); hx=d.x*h;  hxx=d.xx*h;  
    
    % Show present solution
    xx=linspace(1-2*abs(l),1+2*abs(l),20);
    plot(x,h,'b-',L+l,0,'ro',xx,h(N)+(xx-1)*hx(N),'r',x(N),h(N),'b.', ...
        [L,L],[-0.5,0.5],'k--',[0,1.5*L],[0,0],'k--');
    xlim([0,1.5*L]); ylim([-0.5,0.5]);grid on; drawnow;
    pause(1)  ;%  if count==1; break; end

%{
# Domain  adaptation

Here is the core of this code and the difference with [domain_derivative_1D.m](). The previous Newton iteration has provided us with a guess for the domain size by computing the value of $\ell$, so the domain size should now be $L+\ell$. Here we trust this value, and reconstruct the problem for this domain size and do another Newton step. For this we need to rebuild the grid $x$ and the differentiation matrices. 
%}
    % addapt the domain
    L=L+l; 
%{
We need to interpolate the previous solution $h$ on the new grid. The value of $h$ outside the grid of the previous iteration is approximated by a linear Taylor expansion, because this is what we have assumed to get the solution.
%}
    xx=linspace(1.0001*L,2*L,100)';
    hh=h(N)+(xx-L)*hx(N); % linear extrapolation outside of the domain

    [d.x,d.xx,d.wx,x]=dif1D(method,0,L,N,5); % New differentiation
    h=interp1([x; xx],[h; hh],x); % interpolation to the new grid
%{
Now we have transformed the domain so the best approximation for the domain size is now $L$ thus we reset $\ell=0$ for the next Newton iteration
%}
    sol(1:N)=h; % update h
    sol(N+1)=0; % set l to zero

%{
And then we resume as in [domain_derivative_1D.m]()
%}
    % the present solution and its derivatives (on the new domain)
    l=sol(N+1); h=sol(1:N); hx=d.x*h;  hxx=d.xx*h; 
    
    % nonlinear function
    f=[hxx+2; hx(N)+l*hxx(N)-p]; 

    % analytical jacobian
    A=[d.xx, Z(:,1); ...
       d.x(N,:)+l*d.xx(N,:), hxx(N)];

    % Boundary conditions
    loc=[1 N];
    f(loc)=[h(1); h(N)+l*hx(N)];
    A(loc,:)=[I(1,:), 0; ...
              I(N,:)+l*d.x(N,:),hx(N)];
    
    % convergence test
    res=norm(f,inf);
    disp([num2str(count) '  ' num2str(res)]);
    if count>50|res>1e5; disp('no convergence'); break; end
    if res<1e-9; quit=1; disp('converged'); continue; end

    % Newton step
    sol=sol-A\f;
    count=count+1;
end

set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','domain_derivative_1D_addapt.png')

%{
# The figure
Here we see the final result of the computation

![](domain_derivative_1D_addapt.png)

# Exercices/Contributions

* Please 

%}