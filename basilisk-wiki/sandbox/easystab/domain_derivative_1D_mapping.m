%{
# Testing the domain geometry derivative

In this code, this is the same as in [domain_derivative_1D_adapt.m]() but here we use a domain mapping to define the shape of the domain inside the equations instead of rebuilding the domain at each iteration.

%}

clear all; clf;

% parameters
N=50;           % number of grid points
p=-1.1;         % desired slope at L
method='cheb';    % the differentiation

% useful
[d.x,d.xx,d.wx,x]=dif1D(method,0,1,N,5);
Z=zeros(N,N);
I=eye(N);

%initial guess
sol0=[x.*(1-x);1];
sol=sol0;

% Newton iterations
quit=0;count=0;
while ~quit

    % the present solution and its derivatives
    l=sol(N+1); h=sol(1:N); hx=d.x*h;  hxx=d.xx*h;  
    
%{
# The nonlinear function

The computational domain is $x\in[0,1]$ and the physical domain is $\bar{x}\in[0,L]$, so we define a mapping $\bar{x}=Lx$. We get the derivatives of the physical variables from the computational derivatives like this
$$
\frac{\partial h}{\partial \bar{x}}=
\frac{\partial h}{\partial x}\frac{\partial x}{\partial \bar{x}}
$$
and we have from the mapping
$$
\frac{\partial x}{\partial \bar{x}}=\frac{1}{L}.
$$
we thus have
$$
h_\bar{x}=h{x}/L, h_{\bar{x}\bar{x}}=h_{xx}/L^2,  
$$
which we replace in our equations. We thus have the nonlinear function
$$
f(q)=
f\begin{pmatrix}
h\\L
\end{pmatrix}=
\begin{pmatrix}
h_{xx}/L^2+2\\
h_x|_1/L-p
\end{pmatrix}=0
$$
%}
    % nonlinear function
    f=[hxx/l^2+2; hx(N)/l-p]; 

%{
# The Jacobian

We have the actual guess $q$ for the solution, and we say that $q+\tilde{q}$ is the true solution, doing a linear approximation of $f$ arround $q$ (assuming $\tilde{q}$ is small)
$$
f(q+\tilde{q})=
f\begin{pmatrix}
h+\tilde{h}\\L+\tilde{L}
\end{pmatrix}\approx
\begin{pmatrix}
h_{xx}/L^2+2-2h_{xx}\tilde{L}/L^3+\tilde{h}_{xx}/L^2\\
h_x|_1/L-p-h_x/L^2\tilde{L}+\tilde{h}_x/L
\end{pmatrix}=
f(q)+A\tilde{q}
$$
with the expression for the Jacobian
$$
A=
\begin{pmatrix}
\partial_{xx}/L^2 & -2h_{xx}/L^3 \\
\partial_x|_1/L & -h_x|_1/L^2
\end{pmatrix}
$$
%}
% analytical jacobian
    A=[d.xx/l^2, -2*hxx/l^3; ...
       d.x(N,:)/l, -hx(N)/l^2];

%{
# Boundary conditions

With this formulation, the boundary conditions are easy; simply homogeneous Dirichlet at the last grid point:
%}
    % Boundary conditions
    loc=[1 N];
    f(loc)=[h(1); h(N)];
    A(loc,:)=[I(1,:), 0; ...
              I(N,:), 0];
    
    % Show present solution
    plot(l*x,h,'b-',l*x(N),h(N),'b.',[l,l],[-0.5,0.5],'k--',[0,1.5*l],[0,0],'k--');
    xlim([0,1.5*l]); ylim([-0.5,0.5]);grid on; drawnow;
    pause(1);
    
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
print('-dpng','-r80','domain_derivative_1D_mapping.png')

%{

![](domain_derivative_1D_mapping.png)


# Exercices/Contributions

* Please 
* Please 

%}