%{
# The shape of a capillary bridge

Here we use the Keller pseudo arc-length continuation metho, just as we did in [meniscus_overturn_keller.m]() but for the shape of an axisymetric capillary bridge. The capillary bridge is the little bridge of liquid that you have between your fingers when they are wet, you touch them then you move them away. This is also the little bridge of water that sticks together the grains of sand of a sand castle on the beach. 

Here the bridge is held between two cylinders of radius $R$; with a liquid of surace tension $\sigma$, with a length $L$. We start the continuation with a straight bridge: a straight cylinder of radius $R$. We then progressivelly reduce the volume of the bridge by something that would be in an experiment, sucking the fluid from one of the ends of the bridge, for instance if the bridge is held between a cylinder and a small pipe connected to a syringe.

%}

clear all; clf; format compact

% parameters
L=2;             % domain length
n=100;           % number of grid points
rho=1;           % density
sig=1;           % surface tension
g=0;             % gravity
R=1;             % radius 
delta=0.2;      % continuation length
v0=pi*R^2*L;     % initial volume

% Differentiation matrices and integration
[D,DD,wx,x]=dif1D('cheb',0,L,n,3);
Z=zeros(n,n); I=eye(n); 

%{
# Initial guess

Initially, the bridge is a cylinder of constant radius $h=R$, and the volume is equal to that of the cylinder. For the continuation, the initial direction is just to decrease the volmue (thus the -1 at the place of the irection vector corresponding the where the volume is stored in the state of the system. if you want to increase the volume, replace the -1 by a +1).
%}

%initial guess
h0=x*0+R;
sol=[h0; v0];
dir=[zeros(n,1);-1];   % initial direction

% Continuation loop
ind=1;quitcon=0;
while ~quitcon 
    
   solprev=sol;
   sol=sol+dir*delta; % extrapolation of solution
    
   % Newton iterations
   quit=0;count=0;
   while ~quit     
    
      h=sol(1:n); v=sol(n+1); 
      hp=D*h;  hpp=DD*h; a=(1+hp.^2);
        
       %{
       # The nonlinear function and Jacobian
      
      The first thing to be computed is the value of the inner pressure due to the curvature of the interface. The interface is axisymmetric, so there are two components of the curvature $\kappa$  
      $$ 
      \kappa=\frac{1}{R_1}+\frac{1}{R_2}=\frac{1}{h(1+h_x^2)^{1/2}}-\frac{h_{xx}}{(1+h_x^2)^{3/2}}
      $$
      And because of the Lalace pressure jump due to the surface tension through a curved interface, the inner pressure should be $p=\sigma \kappa$. 
      
      Then the nonlinear function *f* is expressed simply as follows, that the pressure gradient along $x$ should be equal to the weight. This way we can as well add the effect of the gravity on the shape of the capillary bridge, provided that the gravity is aligned with the axis of the cylindrical bridge.
      
      To get the analytical Jacobian we need to differentiate the expression of the nonlinear function. In *f* the pressure comes through its $x$ derivative, and the expression of the Jacobian is quite lengthy, so we calculate analytically the Jacobian for the expresion of the pressure, and then we diferentiate this Jacobian numerically using a differentiation matrix. In this sense, the jacobian is not cmopletely analytical, let's say that it is semi-analytical.
      
      Here we should show the steps of the derivation of the Jacobian of the pressure. 
       %}
      
      % nonlinear function
      p=h.^-1.*a.^-0.5-hpp.*a.^-1.5;
      f=[-sig/rho*D*p-g; ...
         dir'*(sol-solprev)-delta];
    
      % Jacobian
      P=diag(-h.^-2.*a.^-0.5) ...
       +diag(3*hp.*hpp.*a.^-2.5-h.^-1.*hp.*a.^-1.5)*D ...
       +diag(-a.^-1.5)*DD; 
   
      A=[-sig/rho*D*P, Z(:,1); ...
         dir'];
%{
# Boundary conditions

We have two obvious boundary conditions, to impose that the radius of the bridge is $R$ at $x=0$ and at $x=L$. Then we see that we have a single variable, *h*, with two derivatives in its equation, but there is in fact a third derivative because the pressure comes through its derivative, so we need to impose in fact a third boundary condition. The most obvious thing is to impose the volume of the bridge which would be free otherwise. The volume of the bridge the sum of thin slices of thickness *dx* and radius *h*:
$$
v-\int_{x=0}^L \pi h^2 dx=0
$$
Thus the Jacobian of this "boundary" condition which is rather better understood as a "constraint" is obtained by saying that $v$ is perturbed as $v+\tilde{v}$ and $h$ is perturbed as $h+\tilde{h}$ where $\tilde{v}$ and $\tilde{h}$ are very small
$$
\begin{array}{l}
v+\tilde{v}-\int_{x=0}^L \pi (h+\tilde{h})^2 dx=0 \\
v+\tilde{v}-\int_{x=0}^L \pi (h^2+2h\tilde{h}+\tilde{h}^2) dx=0 \\
v+\int_{x=0}^L \pi h^2dx +\tilde{v}-\int_{x=0}^L \pi 2h\tilde{h} dx=0 
\end{array}
$$ 
where in the last expression we have neglected the nonlinear terms in $h$ (and the terms in $v$ are all linear).

Then you may ask yourself: "how do we do the integral in our matrices. This is simply done using the integration vector *wx* that we have built above, at the time of building the differentiation matrices, so
$$
v-\int_{x=0}^L \pi h^2 dx
$$
is coded as *v-pi*wx*h.^2*, and
$$
\tilde{v}-\int_{x=0}^L \pi 2h\tilde{h} dx
$$
is coded as *(-2*pi*wx.*h')*h*.
%}

      % Boundary conditions
      cloc=[1,n,n-1];
      f(cloc)=[h(1)-R; h(n)-R; v-pi*wx*h.^2];
      A(cloc,:)=[ ...
         I(1,:), 0; ...
         I(n,:), 0; ...
         -2*pi*wx.*h', 1];

      % Test of convergence
      res=norm(f);
      disp([num2str(count) ' ' num2str(res)]);
      if count>20||res>1e5; disp('-> no convergence'); quitcon=1;return; end
      if res<5e-6; disp('----> converged'); quit=1; continue; end
    
      % New solution
      sol=sol-A\f;
    
      count=count+1;
   end
   
%{
# Validation

Now we have a solution showing the shape of the liquid bridge for a given volume. There is no analytical solution that I know for this, but we can compare our results to results of the litterature. I scanned the figure from the paper "Capillary surfaces: stability from families of equilibria with application to the liquid bridge", by BJ Lowry and P Steen, Proceedings: mathematical and physical sciences, vol 449 N 1937, pages 411-439. This figures show the branch of the bifurcation diagram using the bridge inner pressure and the bridge volume. 

This branch is quite interesting because it has a fold bifurcation: if you reduce progressively the volume of the bridge, you reach a point where it is no longer possible to continue, and if you want to follow the nonlinear solution, you will have to start increasing again the volume. This is why the bridge just breaks if you progressively pump out its fluid.
%}
   
   % Showing the shape of the bridge
   subplot(2,1,1);
   plot(x,h);axis equal; axis([0,L,0,1.5*R]); 
   xlabel('x'); ylabel('r'); title('shape of the bridge');
   
   % drawing image from litterature for comparison
   subplot(2,1,2);
   if ind==1;
      aa=imread('bifurl2.png');s=size(aa);
      xx=linspace(-1,3,s(2)); yy=linspace(4,0,s(1));
      image(xx,yy,aa); axis xy; hold on
   end
   plot(p(1),v/v0,'b.'); 
   xlabel('pressure'); ylabel('Volume/V0'); title('Bifurcation branch');
   axis([0.5,2,0,1.1]);
   drawnow
    
   
   % new direction
   dir=A\[Z(:,1);1]; % the new direction
   dir=dir/norm(dir); % normalization
   ind=ind+1;
end

% print the figure
print('-dpng','-r100','capillary_bridge_continuation.png')

%{
![alt text](/capillary_bridge_continuation.png)

# Exercices/Contributions

* Please show what happens when the volume is larger than the volume of the cylinder
* Please draw a state space of the existence of the capillary bridge when varrying the length of the bridge (its aspect ratio) (below a given volume the bridge cannot exist as showed by the fold of the bifurcation branch)
* Please do the derivation of the Jacobian of the expression for the pressure (this is not coding, this is documenting and calculating with paper and pen)
* Please show that the continuation loops stops when the radius at the center of the bridge becomes zero. And show that in this limit the solution corresponds to two spheres: one on the left and one on the right, that touch eachother just by one point. This is another analytical solution that we could put on the bifurcation branch (since these are spheres you know their volume and you know the inner pressure)
* Please show that when the bridge becomes long, you may have more branches with more than two spheres at the end of the branch
* Please investigate what happens when the length of the bridge becomes longer than its end perimeter (related to the Rayleigh-Plateau instability)
* Please test and validate how the bridge behaves when there is gravity (and find some results in the litterature)
* Please investigate how the domain of existence of the bridge depends on the gravity: gravity increases or decreases the zone of existence?

%}
