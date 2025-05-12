%{
# The Blasius boundary layer

In this code, we solve the Blasius equations to get the boundary layer velocity profile. We use the Newton iterations as we used for the shape of the meniscus [meniscus.m](), by using the analytical Jacobian of the nonlinear function.
%}

clear all; clf;

% parameters
L=10; % box height
N=100; % number of gridpoints

%{
For the differentiation matrices, since we as well need the third derivative, we build directly the matrices from [chebdif.m]() as in [diffmat.m]() instead of using the function [dif1D.m]().
%}

% differentiation matrices
scale=-2/L;
[y,DM] = chebdif(N,3); 
dy=DM(:,:,1)*scale;    
dyy=DM(:,:,2)*scale^2;    
dyyy=DM(:,:,3)*scale^3;
y=(y-1)/scale;  	
I=eye(N); Z=zeros(N,N);

%{

# Initial guess

To start the computations and for validation, we use the order 4 Polhausen velocity profile for the boundary layer. The idea of Pohlhausen is to use an order 4 polynomial that satisfies no-slip at the wall and reaches the free-stream velocity $1$ with a zero slope at height $y=\delta_{99}$. The extra boundary condition is obtained from $\partial_y^2u=0$ at the wall.
This $\delta_{99}$ is the "total" height of the velocity defect of the profile. 
The rescaled coordinate $\eta$ is
$\eta=y/\delta_{99}$
The polynomial is
$$
\begin{array}{l}
u0=1-(1-\eta)^3(1+\eta), \eta<=1\\
u0=1, \eta>1\\
\end{array}
$$  
To find the evolution of $\delta_{99}$ we write the 
  Von K\'arm\'an equation.
  $$\frac{\partial \delta_2 u_e^2}{\partial x} + \delta_1 u_e \frac{d u_e}{d x} = \frac{\partial u}{\partial y}|_0$$
  where the displacement thickness is 
$$\delta_1=\int_0^\infty (1-u)dy$$
  and momentum thickness  
$$\delta_2=\int_0^\infty u(1-u)dy,$$
they are linked and linked to $\delta_{99}$ thanks to the integration of the profile, this gives 
$\delta_2= 37 \delta_{99}/75$ and $\delta_2= \delta_1/H$ with $H=189/74$. The shear at the wall is computed from the profile and depends of the thickness: 
 $\frac{\partial u}{\partial y}|_0= f_2 H /\delta_1$ with $f_2 =74/315$.
the VK equation is then (here as function of $\delta_1$)
 $$\frac{\partial }{\partial x} (\frac{\delta_1}{H})  = \frac{f_2 H}{\delta_1}$$
the solution (with $\sqrt{2 f_2}H= 9\sqrt{7/185} \simeq 1.7507$) is then 
$$\delta_1 = \sqrt{2 f_2}H x^{1/2}.$$
The shear is $\tau=\frac{f_2 H}{\delta_1} = (1/3)\sqrt{37/35} x^{-1/2}\simeq 0.34 x^{-1/2}$,
and finaly the result of the calculation of Polhausen estimation of the boundary layer total thickness gives
$$
\delta_{99}=(6\sqrt{35/37})x^{1/2}
$$
  
%}

% initial guess from order 4 Polhausen profile
scalepol=6*sqrt(35/37);
eta=y/scalepol; deta=dy*scalepol; % rescaled vertical coordinate
u0=1-(1-eta).^3.*(1+eta);
u0(eta>1)=1;

%{
but now, the variable of the Blasius equation $g$ is the integral of the velocity profile $u$ 
$$
u=g_y
$$
so we need to integrate the Polhausen profile to get our initial guess. We do this numerically because this is a nice example of the power of differentiation matrices. The equation above is a differential equation with one boundary condition at the wall $g(y=0)=0$, see [differential_equation.m]() for details
%}

A=deta; A(1,:)=I(1,:); u0(1)=0; % set up integration 
g0=A\u0; % compute integral

sol=g0;

% Newton iterations
quit=0;count=0;
while ~quit
    
    % the present solution and its derivatives
    g=sol; gy=dy*g; gyy=dyy*g; gyyy=dyyy*g;
    
%{
# The Blasius equation

Here is the nonlinear equation that we wish to solve
$$
0=g_{yyy}+gg_{yy}
$$
and the analytical Jacobian is
$$
A=\partial_{yyy}+g\partial_{yy}+g_{yy}I
$$
%}

    % nonlinear function
    f=2*gyyy+g.*gyy;
    
    % analytical jacobian
    A=2*dyyy+diag(g)*dyy+diag(gyy)*I;

%{
# Boundary conditions
The boundary conditions are homogeneous Dirichlet and Neuman at the wall and Neuman 1 at the top. The constraint matrix for the boundary conditions is thus
$$
C=\left(\begin{array}{c}
I|_0\\
\partial_y|_0\\
\partial_y|_L
\end{array}\right)
$$
so that the homogeneous and nonhomogeneous boundary conditions are expressed
$$
Cg=\left(\begin{array}{c}
0\\
0\\
1
\end{array}\right)
$$
%}
    % Boundary conditions
    loc=[1,2,N];
    C=[I(1,:); dy(1,:); dy(N,:)];
    f(loc)=C*g-[0;0;1];
    A([1,2,N],:)=C;
    
    % convergence test
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if count>50|res>1e5; disp('no convergence'); break; end
    if res<1e-5; quit=1; disp('converged'); continue; end
    
    % Newton step
    sol=sol-A\f;
    count=count+1;
end

%{
To recover the velocity profile from $g$ is
$$
u=g_y
$$
%}

% the velocity profile
u=dy*g;

%{
 
# Result and Plot

Here we plot the velocity profile, and we compare with the Polhausen profile

%}
% Showing the result
plot(u,y,'b-',u0,y,'r--'); xlabel('u'); ylabel('y'); ylim([0,y(end)]);
title('Blasius boundary layer');
legend('numerical','Polhausen');
grid on

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','blasius.png');

%{
We would like to compare some typical quantities of the profile to compare to the Polhausen one, first the displacement thickness and then the shear at the wall:
%}

% Compute the displacement thickness 
w=([diff(y)',0]+[0,diff(y)'])/2; % integration weights

d1num=w*(1-u)
d1pol=36/120*scalepol

shearnum=gyy(1)
shearpol=1/3*sqrt(37/35)

%{
![The velocity profile](blasius.png)

# Exercices/Contributions

* Please find a way to do a validation of this result.
* The Falkner-scan boundary-layer (with a positive or negative pressure gradient along the longitudinal direction: the boundary layer on a flat plate with an angle to the free-stream)-------------> [falkner_skan.m]()
* The stagnation point solution -------------> [hiemenz.m]()
* And the same thing using the Keller arclength continuation to get the fold of the branch----------> [falkner_skan_continuation.m]().

# Links
* the Blasius solution [blasius.m]() $\beta=0$
* the Falkner Skan solution [falkner_skan.m]() for one value of $\beta<0$ (with separation)
* the Hiemenz solution [hiemenz.m]() for $\beta=1$
* the Falkner Skan solution [falkner_skan_continuation.m]() for all the $\beta$
 
# bibliography
* H. Schlichting "boundary layer theory" 6th edition Mc Graw Hill (Chapter X approximate methods for steady solutions, page 206) 
* PY Lagrée "Boundary Layers" Master 2 UPMC  [http://www.lmm.jussieu.fr/~lagree/COURS/CISM/blasius_CISM.pdf]()
* PY Lagrée "Transferts thermiques et massiques dans les fluides" ENSTA   [http://www.lmm.jussieu.fr/~lagree/COURS/ENSTA/C6int.ENSTA.pdf]()



%}