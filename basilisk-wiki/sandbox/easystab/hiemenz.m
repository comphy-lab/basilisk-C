%{
# The Hiemenz boundary layer

This is a particular case of a Falkner-Skan boundary layer corresponding to $\beta=1$, that is the boundary layer arround a stagnation point. In this case, there is a easy solution ofthe Pohlhausen profile that we use as an initial guess for Newton and as a comparison for the final result.

See [blasius.m]() and [falkner_skan.m]() for other details.
We use the polynomial closure 
at order 4,  so that
$$\frac{u(\eta)}{\bar u_e}=(2 \eta - 2 \eta^3 +\eta^4) +  \frac{1}{6}\Lambda (\eta- 3 \eta^2 + 3 \eta^3  -\eta^4)$$
with 
$\Lambda= \delta^2 d\bar u_e/d\bar x$
in the Hiemenz case $\bar u_e=\bar  x$ then $\Lambda= \delta^2$, 
The Von K\'arm\'an equation
$$
 \frac{d}{d \bar x} (\bar u_e^2\tilde  \delta_{2})
 +    {\tilde \delta_{1}} {\bar u_e}  \frac{d \bar u_{e}}{d \bar x}   
  =
 \frac{u'(0)}{\tilde \delta} \bar u_{e},  \;\;\mbox{ reads 
  }\;\;
2 \tilde  \delta_{2}
 +    {\tilde \delta_{1}}  
  =
 \frac{u'(0)}{\tilde \delta },  $$
which is an equation 
were 
$\delta_1/\delta=(36 -\Lambda)/120$ and $\delta_2/\delta= 37/315 - \Lambda/945 - (\Lambda^2)/9072,$
and $u'(0) = (2 + \Lambda/6)$
and remember that $\Lambda= \delta^2$, we then substitute in VK:
$$
\delta \left(\frac{3}{10}-\frac{\delta^2}{120}\right)-\frac{\frac{\delta^2}{6}+2}{\delta}+2 \delta
   \left(-\frac{\delta^4}{9072}-\frac{\delta^2}{945}+\frac{37}{315}\right)=0$$
we
solve and find numericaly 
$\delta=2.65562$
this gives 
$\Lambda=7.05$ and $\delta_1=0.640617$, and $H=2.30809$ and $\tau= \frac{u'(0)}{\tilde \delta }=1.1957$


  

**Dependency**

* [chebdif.m]()

%}

clear all; clf;
%setenv("GNUTERM","X11")

% parameters
L=10; % box height
N=100; % number of gridpoints
beta=1; % pressure gradient parameter

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

% initial guess from order 4 Pohlhausen profile
d99=2.65562;
lambda=d99^2;
eta=y/d99; deta=dy*d99; % rescaled vertical coordinate
u0=1-(1-eta).^3.*(1+(1-lambda/6)*eta);
u0(eta>1)=1;
A=deta; A(1,:)=I(1,:); u0(1)=0; % set up integration 
g0=A\u0; % compute integral

sol=g0;


% Newton iterations
quit=0;count=0;
while ~quit
    
    % the present solution
    g=sol; gy=dy*g; gyy=dyy*g; gyyy=dyyy*g;
    
%{
# The FS equation

Here is the nonlinear equation that we wish to solve
$$
0=g_{yyy}+gg_{yy} + \beta (1 - g_y^2)
$$
and the analytical Jacobian is
$$
A=\partial_{yyy}+g\partial_{yy}+g_{yy}I -2 \beta g_y \partial_{y}
$$
%}

    % nonlinear function
    f=gyyy+g.*gyy+beta*(1-gy.^2);
    
    % analytical jacobian
    A=dyyy+diag(g)*dyy+diag(gyy)-2*beta*diag(gy)*dy;

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
    A(loc,:)=C;
    
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
% the solution
u=dy*g;

% Compute the displacement thickness and normalize it to unity
INT=([diff(y)',0]+[0,diff(y)'])/2; % integration weights
mt=INT*(1-u);

% Showing the result
plot(u,y,'b-',u0,y,'r--'); xlabel('u'); ylabel('y'); ylim([0,y(end)]);
title('Falkner-skan boundary layer');
legend('Falkner-Skan','Polhausen');
grid on

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','hiemenz.png');

% Compute the displacement thickness 
w=([diff(y)',0]+[0,diff(y)'])/2; % integration weights

d1num=w*(1-u)
d1pol=(36-lambda)/120*d99

shearnum=gyy(1)
shearpol=(2+lambda/6)/d99



%{
# Result and Plot

The   Hiemenz flow $f''' +  ff'' +  (1 - f'^2)=0$ as solution $f''(0)=1.2325$ (compare to 1.1957 for Pohlhausen4), the displacement thickness
is $\int (1-f')d \eta=0.6479$ (compare to 0.640617 for Pohlhausen4).


![The velocity profile](/sandbox/easystab/hiemenz.png)
 

# Links
* the Blasius solution [blasius.m]() $\beta=0$
* the Falkner Skan solution [falkner_skan.m]() for one value of $\beta<0$ (with separation)
* the Hiemenz solution [hiemenz.m]() for $\beta=1$
* the Falkner Skan solution [falkner_skan_continuation.m]() for all the $\beta$
 
# Exercices/Contributions

* Please check that the axi Hiemenz flow $f''' + 2  ff'' +  (1 - f'^2)=0$ as solution $f''(0)=1.31194$ (in 2D 1.2325), the displacement thickness
is $\int (1-f')d \eta=0.568902$ (compare to 0.6479 in 2D).
* Please check the convergent case $\beta \rightarrow \infty$ so   $f'''+1-f'^2=0.$ $f''_0=1.15,$ and $\int_0^\infty (1-f')d\eta=0.779 $




# Bibliography
see [blasius.m]()


%}