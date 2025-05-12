%{
# Testing the domain geometry derivative

Here we use a method close to the "flattening" used typically for the stability of flows with a free-surface, like for instance [free_surface_gravity.m](). 

The problem we awant to solve is
$$
h_{xx}+2=0
$$
on the domain, with boundary conditions
$$
h(0)=0, h(L)=0, h_x(L)=p
$$
so I have a second derivative problem for which I should impose two boundary conditions. But I want to impose one more condition, the value of the slope at the right boundary. So I need to let something free. Thus I let free the position of the right boundary $L$. Thus the right boundary is at $x=L+\ell$.

%}

clear all; clf;

% parameters
N=50;           % number of grid points
p=-1.2;         % desired slope at L

% useful
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,1,N);
Z=zeros(N,N);
I=eye(N);

%initial guess
sol0=[x.*(1-x);0];
sol=sol0;

% Newton iterations
quit=0;count=0;
while ~quit

    % the present solution and its derivatives
    h=sol(1:N); hx=d.x*h;  hxx=d.xx*h; 
    l=sol(N+1);
    
%{
# The nonlinear function

The computational domain is fixed on $[0,1]$, and the function is modeled outside of the domain by using a linear extrapolation. Thus the function that should be zero is
$$
f(q)=f\begin{pmatrix}h\\\ell \end{pmatrix}=
\begin{pmatrix}
h_{xx}+2\\
h|_L+\ell h_x|_L-p
\end{pmatrix}=0.
$$
%}
% nonlinear function
    f=[hxx+2; hx(N)+l*hxx(N)-p]; 

%{
# The Jacobian
The actual guess for the solution of the problem is $q$, and we search the perturbation $\tilde{q}$ such that $q+\tilde{q}$ is closer to the actual solution. We approximate this function at first order
$$
f(q+\tilde{q})=f\begin{pmatrix}h+\tilde{h}\\\ell+\tilde{\ell} \end{pmatrix}=
\begin{pmatrix}
h_{xx}+\tilde{h}_{xx}+2\\
h|_L+\tilde{h}|_L+(\ell+\tilde{\ell}) (h_x|_L+\tilde{h}_x|_L)-p
\end{pmatrix}\approx
f(q)+A\tilde{q}=0
$$
with
$$
A=\begin{pmatrix}{cc}
\partial_{xx} & 0\\
\partial_x|_L+\ell\partial_{xx}|_L & h_xx|_L
\end{pmatrix}
$$
%}
    % analytical jacobian
    A=[d.xx, Z(:,1); ...
       d.x(N,:)+l*d.xx(N,:), hxx(N)];

%{
# Boundary conditions

We have nonlinear boundary conditions because of the extrapolation
$$
\begin{array}
h(0)=0\\
h(L+\ell)\approx h|_L+\ell h_x|_L=0\\
\end{array}
$$
and the associated Jacobian lines after linearization are
$$
\begin{pmatrix}{cc}
I|_0 & 0\\
I|_L+\ell\partial_x|_L & h_x|_L
\end{pmatrix}
$$
%}
% Boundary conditions
    loc=[1 N];
    f(loc)=[h(1); h(N)+l*hx(N)];
    A(loc,:)=[I(1,:), 0; ...
              I(N,:)+l*d.x(N,:),hx(N)];
    
    % Show present solution
    x2=linspace(1-2*abs(l),1+2*abs(l),20);
    plot(x,h,'b-',1+l,0,'ro',x2,h(N)+(x2-1)*hx(N),'r',x(N),h(N),'b.',[1,1],[-0.5,0.5],'k--',[0,1.5],[0,0],'k--');
    xlim([0,1.5]);grid on; drawnow;
    
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
print('-dpng','-r80','domain_derivative_1D.png')

%{
# The figure

Here the solution $h$ is shown in blue and the linear extrapolation is shown in red.

![](domain_derivative_1D.png)

# Exercices/Contributions

* Please
* Please

%}