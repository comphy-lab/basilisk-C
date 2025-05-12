%{
# The Falkner Skan boundary layer

In this code, we solve the Falkner Skan equation to get the boundary layer velocity profile. 
We use the Newton iterations as we used for the shape of the meniscus, by using the analytical Jacobian of the nonlinear function.

This is very similar to [blasius.m]() but with an additional term for the pressure gradient.
%}

clear all; clf;
%setenv("GNUTERM","X11")

% parameters
L=12; % box height
N=100; % number of gridpoints
beta=-0.1; % pressure gradient parameter

% differentiation matrices
scale=-2/L;
[y,DM] = chebdif(N,3); 
dy=DM(:,:,1)*scale;    
dyy=DM(:,:,2)*scale^2;    
dyyy=DM(:,:,3)*scale^3;
y=(y-1)/scale;  	
I=eye(N); Z=zeros(N,N);

% initial guess from order 4 Polhausen profile
d99=6*sqrt(35/37);
eta=y/d99; deta=dy*d99; % rescaled vertical coordinate
u0=1-(1-eta).^3.*(1+eta); u0(eta>1)=1;
A=deta; A(1,:)=I(1,:); u0(1)=0; % set up integration 
g0=A\u0; % compute integral
sol=g0

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
%%y=y/mt;

% Showing the result
plot(u,y,'b-',u0,y,'r--'); xlabel('u'); ylabel('y'); ylim([0,y(end)]);
title('Falkner-skan boundary layer');
legend('numerical','Pohlhausen');
grid on

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','falkner_skan.png');


%{
# Result and Plot

![The velocity profile](/sandbox/easystab/falkner_skan.png)

# Exercices/Contributions

* Please validate this result with Pohlhausen profile. => [http://basilisk.fr/sandbox/easystab/stab2014/falkner_skan_validation]()
* Please add the Cooke transverse velocity profile for the boundary layer on a swept-wing
* Please validate the behavior for large $\beta$ (which is the self similar flow in a convergent) => [http://basilisk.fr/sandbox/easystab/stab2014/falkner_skan_validation]()


# Links
* the Blasius solution [blasius.m]() $\beta=0$
* the Falkner Skan solution [falkner_skan.m]() for one value of $\beta<0$ (with separation)
* the Hiemenz solution [hiemenz.m]() for $\beta=1$
* the Falkner Skan solution [falkner_skan_continuation.m]() for all the $\beta$ 
  
# Bibliography
see [blasius.m]()
 
 
%}