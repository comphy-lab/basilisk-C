%{

# Model for the boundary layer in curved geometry

We seek here to build a model based on the viscous/acceleration growth/decay of a boundary layer to compare to [venturi.m]().

The model used is the Interactive Boundary Layer Theory with integral method for the boundary layer.  It consists to solve the  conservation of flux in the ideal fluid layer:
$$
u_e(R-y_w-\delta_1)^2=1
$$
 It is as if the radius of the pipe, the initial one minus the shape of the stenosis: ($R-y_w$) was reduced because of the boundary layer thickness $\delta_1$.
And the second equation is the growth of the boundary layer thickness because of the effect of the pipe geometry and the variation of the ideal fluid velocity $u_e$, this is the Von K'arm'an boundary integral layer equation: 
$$
\frac{\partial }{\partial x} (\frac{\delta_1}{H})+\frac{\delta_1}{u_e}(1+\frac{2}{H})\frac{\partial u_e }{\partial x} =\frac{f_2H}{\delta_1 u_e}
$$

Note that in classical Boundary layer theory $u_e$ is given by $u_e(R-y_w)^2=1$, here we have a strong coupling.


 
$H$ and $f_2$ are two fixed parameters (for now).


%}

% if octave
%setenv("GNUTERM","X11")

clear all; clf;

% parameters
L=.5;             % domain length
N=200;            % number of grid points
H=2.59;           % H shape parameter
f2=0.22;          % f2 friction parameter
d0=0.01;          % thickness at inflow
R=1;              % rayon du tuyau
xb=.2;            %position of the bump
wb=0.05 ;         % width of the bump
alpha=0.3;        % degree of closure 

% Differentiation and integration
[D,DD,wx,x]=dif1D('cheb',0,L,N,3);
Z=zeros(N,N); I=eye(N); 
    
% differentiation and integration
scale=-2/L;
[x,DM] = chebdif(N,2); 
D=DM(:,:,1)*scale;    
DD=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 

% shape of the pipe is R-yw
yw=alpha*exp(-((x-xb)/wb).^2);

%initial guess
sol=[1+0*x;d0+0*1.7*x.^(.5)];
    
% Newton iterations
disp('Newton loop')
quit=0;count=0;
while ~quit     

    % the present solution and its derivatives
    u=sol(1:N); d=sol(N+1:2*N);
    ux=D*u; dx=D*d;

%{

# The nonlinear function

The first equation is the conservation of flux in the ideal fluid layer $u_e$:
$$
u_e(R-y_w-\delta_1)^2-1=0
$$
This is not a differential equation, but we keep it explicitely so that the expression of the Jacobian will be easier (we could have removed $u_e$ from the variable by saying 
$$
u_e=1/(R-y_w-\delta_1)^2
$$
and replaced that in the second equation)

And the second equation is the growth of the boundary layer thickness because of the effect of the pipe geometry and the variation of the centerline velocity $u_e$ (and as well viscosity):
$$
\frac{\partial }{\partial x} (\frac{\delta_1}{H})+\frac{\delta_1}{u_e}(1+\frac{2}{H})\frac{\partial u_e }{\partial x} -\frac{f_2H}{\delta_1 u_e}=0
$$
Note that the shear stress is $\tau = u_e\frac{f_2H}{\delta_1 }$ it is an important quantity 

%}
    % shear stress
    tau = (f2*H)*u./d ;
    
    % flux, nonlinear VK function
    f=[u.*(R-yw-d).^2-1; ...
       dx/H+(1+2/H)*d.*ux./u-tau./u./u];     

%{

# The analytical Jacobian

We get it by replacing in the equations $u_e$ by $u_e+\tilde{u_e}$ and the same for the thickness $\delta_1+\tilde{\delta_1}$. The Jacobian is the linear operator acting on the order one terms of $\tilde{u_e}$ and $\tilde{\delta_1}$. 

%}
    % analytical jacobian
    A=[diag((R-d-yw).^2),   diag(-2*u.*(R-d-yw)); ...
       (1+2/H)*(diag(1./u)*D-diag(ux./u.^2))+f2*H*diag(1./(d.*u.^2)), ...
       D/H+(1+2/H)*diag(ux./u)+f2*H*diag(1./(d.^2.*u))];
       
%{

# Boundary conditions

We have two unknowns but just one derivative, so there is one boundary condition to impose. We set at inflow ($x=0$) the value of the boundary layer thickness to $\delta_1(x=0)=\delta_{10}$.

%}

    % Boundary conditions
    f(1)=d(1)-d0;
    A(1,:)=[Z(1,:),I(1,:)];

    % convergence test  
    plot(x,u,'b-',x,d,'r-',x,R-yw,'b-',x,tau/10,'k--'); axis([0,L,0,4*R]);grid on; drawnow
    
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if count>666|res>1e5; disp('no convergence');conv=0;break; end
    if res<15; quit=1; disp('converged'); conv=1;continue; end
    
    % Newton step
    sol=sol-A\f;   
    count=count+1;
end

%{

# Plot of the results

we plot the ideal fluid velocity, the boundary layer thickness and the wall shear stress

%}

xlabel('x'); ylabel('y');title('Boundary layer');
legend('Ideal Fluid velocity','BL thickness','shear stress/10','geometry');
set(gcf,'paperpositionmode','auto');
print('-dpng','-r75','boundary_layer_pyl.png')

%{

![Evolution of the ideal fluid velocity and the boundary layer thickness](boundary_layer_pyl.png)

# Exercices/Contributions

* Please use a closure of $H(\delta_1^2 du_e/dx)$ and $f_2(H)$ and modify the code accordingly
* Please compare to [venturi.m]().


# Links
* the Blasius solution [blasius.m]() $\beta=0$
* the Falkner Skan solution [falkner_skan.m]() for one value of $\beta<0$ (with separation)
* the Hiemenz solution [hiemenz.m]() for $\beta=1$
* the Falkner Skan solution [falkner_skan_continuation.m]() for all the $\beta$
* the Navier Stokes and RNSP solution [venturi.m]().

# Bibliography

* [http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/lagreelorthois05.pdf]() 

* [http://www.lmm.jussieu.fr/~lagree/TEXTES/GLOTTE/1190.pdf]() 

%}
