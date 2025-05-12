%{

In Easystab, we adopt a monolythic description the systems we study. For a nonlinear system for instance, everything is described in a single vector $f$ that depends nonlinearly upon all the variables (the unknowns). This is very powerful since we can use directly the methods of resolutions like for instance the Newton iterations and nonlinear continuation. On the other hand we have to deal with the Jacobians: build them and manipulate them. This page is a pedagogical approach to checking the expression of your Jacobian. 

The correct expression of the Jacobian is needed for the Newton iterations to find roots of nonlinear systems. If your Newton does not converge, it is often because you Jacobian is wrong. Getting the correct expression for the jacobians is not always straightforward, especially for a system with several equations and several interacting unknowns. This gets especially difficult when some of the variables are on a 2D grid (fluid velocity and pressure) and some other are one a 1D grid (free-surface deformation) and some others are scalars (reference pressure) like in [free_surface_navier_stokes.m]().

This means that you need to have a way to progressivelly test your code for the Jacobian. In this example, I show how to verify that the expression of the Jacobian $A$ is coherent with the expression of the nonlinear function $f$. For this I build both the numerical Jacobian by numerically perturbing $f$ one degree of freedom at a time, then comparing to my expression.

The numerical Jacobian is nice because it is obtained directly from the expression of the nonlinear function $f$: you don't need to derive formulas for the analytical Jacobian where ou can make mistakes that are hard to find. You are sure that the Jacobian fits with $f$ so you have a good chance that your Newton iterations will converge. On the other hand the computation is very long, whereas it is quick to build the matrix expression of the analytical Jacobian you have derived. For testing your formulas, you can take a coarse grid once for all and check that everything is OK, then remove the computation of the numerical Jacobian.

Another important aspect of Jacobians is that they should be invertible for the Newton method, this is discussed in [troubleshooting_jacobian_invertibility.m]().

Here we use as example the Navier-Stokes equation in a 2D domain from [jet_2D.m](). We take a coarse grid because computing the numerical Jacobian takes some time and also so that the sparsity structure of the Jacobian is more clear.

%}

clear all; clf; format compact

% parameters
Re=10; % reynolds number
Nx=20; % number of grid nodes in z
Ny=20; %number of grid nodes in r
Lx=2; % domain length
Ly=1; % domain height

% differentiation
[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,1,Ny);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);
D.lap=D.yy+D.xx;

% useful
l.u=(1:NN)'; l.v=l.u+NN; l.p=l.v+NN;
II=speye(3*NN);
n=3*NN; % The total number of degrees of freedom

% preparing boundary conditions
neuploc=[l.ctl;l.ctr;l.ctr-Ny];  % where to impose the neumann condition on the pressure
p0loc=2*Ny; % where to impose p=0
dir=[l.left;l.top;l.bot]; % where to put Dirichley on u and w

% initial guess
V=zeros(NN,1);
U=exp(-((Y-Ly/2)/(Ly/8)).^2);
P=-X/Re; P=P-P(p0loc); % pressure zero at p0loc
q0=[U(:);V(:);P(:)];
q=q0;

%{
# Perturbing the nonlinear function

The jacobian comes from a linear approximation of how $f$ depends on its variables close to a reference state $q$. The peturbation is called $\tilde{q}$ and should be small:
$$
f(q+\tilde{q})\approx f(q)+A\tide{q}
$$
We can build one by one the columns of the Jacobian by chosing simply $\tilde{q}$ full of zeros except at the row $i$ where we give it the small value $\varepsilon$. This small value should be small so that $f$ will effectively vary linearly but not too small to avoid round-off errors. We take Typically $\varepsilon=1e^{-8}$, half way (logarithmically) between 1 and machine accuracy. 

We then have 
$$
\tilde{q}_i=\varepsilon\begin{pmatrix}
0 \\ 0 \\ \vdots \\ 1 \\ \vdots \\ 0 \\ 0
\end{pmatrix}
$$
so that $A\tilde{q}_i=\varepsilon A_i
$$
where $A_i$ is the column number $i$ of $A$. Thus the column is
$$
A_i=\frac{f(q+\tilde{q}_i)-f(q)}{\varepsilon}.
$$
This means that we can get the Jacobian by $n+1$ evaluations of the nonlinear function $f$: one to get $f(q)$ and then $n$ times to get each column of $A$. 
%}

% loop to perturb independently each degree of fredom
eps=1e-8;
Anum=zeros(n,n);
for gre=0:n; 
%{
Below is just a compact and tricky formulation to add $\varepsilon$ to the line number *gre* of *qq*. When *gre* is zero at the first iteration of the loop, $q$ is not perturbed, this is to get the reference $f(q)$ which only need to be computed once. 
%}
    qq=q+((1:n)==(gre))'*eps;
    
    % the present solution and its derivatives
    U=qq(l.u); V=qq(l.v); P=qq(l.p);
    Ux=D.x*U; Uy=D.y*U;
    Vx=D.x*V; Vy=D.y*V;
    Px=D.x*P; Py=D.y*P;
    
    % nonlinear function
    feps=[-U.*Ux-V.*Uy+D.lap*U/Re-Px; ...
        -U.*Vx-V.*Vy+D.lap*V/Re-Py; ...
        Ux+Vy];
    
    % boundary conditions
    loc=[l.u(dir); l.v(dir); l.p(p0loc); ...
        l.u(l.right); ...
        l.v(l.right); ...
        l.p(neuploc)];
    
    C=[II([l.u(dir);l.v(dir);l.p(p0loc)],:); ...     % Dirichlet on u,v,and p
        D.x(l.right,:)*II(l.u,:); ...   % Neuman on v at outflow
        D.x(l.right,:)*II(l.v,:); ...   % Neumann on u at outflow
        D.lap(neuploc,:)/Re*II(l.v,:)-D.x(neuploc,:)*II(l.p,:)]; % neuman constraint on pressure
    
    feps(loc)=C*(qq-q0);
    
    
    if gre==0;
%{
# The analytical jacobian

At the first iteration of the loop, $q$ is not perturbed, this is the good time to store the reference $f(q)$ and build the analytical jacobian using all the derivatives of the parts of $q$. This expression should be the same as for [jet_2D.m]().
%}
        % This is the case without perturbations
        f=feps;
        Aana=[-(spd(U)*D.x+spd(Ux)+spd(V)*D.y)+D.lap/Re, spd(Uy), -D.x; ...
            -spd(Vx),-(spd(U)*D.x+spd(V)*D.y+spd(Vy))+D.lap/Re, -D.y; ...
            D.x, D.y, Z];
        Aana(loc,:)=C;
        
    else
%{
# The numerical Jacobian

We have computed $f_\varepsilon$ so we can store the corresponding column of the Jacobian using simple *finite differences*.
%}
        % We store one column of the Jacobian
        Anum(:,gre)=(feps-f)/eps;
    end
end

%{
# Validation

We first look at the sparsity structure of the numerical Jacobian. The Navier-Stokes equations is like that
$$
f\begin{pmatrix}
u\\v\\p
\end{pmatrix}
=
f\begin{pmatrix}
f_u(u,v,p)\\f_v(u,v,p)\\f_p(u,v,p)
\end{pmatrix}
$$
so a small perturbation gives
$$
f(q+\tilde{q})=
f\begin{pmatrix}
f_u(u+\tilde{u},v+\tilde{v},p+\tilde{p})\\f_v(u+\tilde{u},v+\tilde{v},p+\tilde{p})\\f_p(u+\tilde{u},v+\tilde{v},p+\tilde{p})
\end{pmatrix}
\approx
f\begin{pmatrix}
u\\v\\p
\end{pmatrix}
+
\begin{pmatrix}
A_{uu}&A_{uv}&A_{up}\\
A_{vu}&A_{vv}&A_{vp}\\
A_{pu}&A_{pv}&A_{pp}\\
\end{pmatrix}
\begin{pmatrix}
\tilde{u}\\\tilde{v}\\\tilde{p}
\end{pmatrix}
$$
%}

% validation
subplot 121 
spy(Anum,'b.'); hold on
v=1:n;
plot(v,v*0+l.u(end),'k--',v,v*0+l.v(end),'k--')
plot(v*0+l.u(end),v,'k--',v*0+l.v(end),v,'k--')
xlabel('columns'); ylabel('rows');
title('Nonzero elements of the numerical Jacobian')

subplot 122
spy(abs(Aana-Anum)>1e-4,'r.'); hold on
plot(v,v*0+l.u(end),'k--',v,v*0+l.v(end),'k--')
plot(v*0+l.u(end),v,'k--',v*0+l.v(end),v,'k--')
xlabel('columns'); ylabel('rows');
title('Numerical-Analytical mismatch')



set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','troubleshooting_jacobian_formula.png')

%{


# The figure

On the left figure below, we use the function *spy* to show the sparsity structure of the Jacobian, it draws a dot at every nonzero element. We have drawn dashed lines to emphasize the $3\times 3$ block structure of the Jacobian.

The right figure shows the nonzero elements of the difference between the numerical and analytical Jacobian. We chose to say that they are different when they differ from more than $1e^{-4}$. You see that there is a problem with the $A_{uv}$ block. 

You can also see that this mismatch does not have the shape of a differentiation matrix: only the diagonal terms of the block are in mismatch. Also, with some experience you will be able to distinguish the sparsity structure of derivatives in $x$ or $y$. See for instance [diffmat_2D.m]() for the sparsity structure of 2D Chebychev differentiation matrices. Here we use finite differences in $x$ and Chebychev in $y$. This will help you to quickly find by comparison which terms of your analytical Jacobian have troubles. You can already learn a lot by comparing the expression of the Jacobian above and its sparsity structure below.

![The figure](troubleshooting_jacobian_formula.png)


# Exercices/Contributions

* Please find the bug in the $A_{uv}$ block of the expression of the Jacobian.
* Please do this for another system.
%}
