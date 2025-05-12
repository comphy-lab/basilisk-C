%{
# The capillary Venturi

The capillary Venturi is a capillary bridge with a troughflow. See [capillary_bridge_continuation.m]() for the code for the nonlinear bifurcation branch of the simple capillary bridge, the code you see here is mostly built upon that one.

here the difference is that there is a flow through it. Experimentally you can realise this system by having a capillary bridge between the end of two pipes, and you blow through the left pipe exactly the same flow as you blow through the right pipe. If the two fluxes are not equal, then the volume of the bridge will increase or decrease. A simple way is to use a syringe-pump and use two syringes, one which is being emptied, connected to the left pipe, and one which is filled, connected to the right pipe.

This flow is interesting because it is a simpler representation of what happen through necks in retracting liquid ligaments, like for instance in "Recoil of a liquid filament: escape from pinch-off through creation of a vortex ring" by Hoepffner & Paré, Journal of Fluid Mechanics, vol 734 (2013) pages 183-197.

The only parameter in addition to [capillary_bridge_continuation.m]() is $U$ the velocity through the bridge, and here the Keller continuation is performed on $U$ instead on the volume of the bridge. This volume is kept constant.
%}

clear all; clf; format compact

% parameters
L=2;             % domain length
n=100;           % number of grid points
rho=1;           % density
sig=1;           % surface tension
g=0;             % gravity
R=1;             % radius 
delta=0.5;       % continuation length
v=0.99*pi*R^2*L; % volume
U=0;             % initial velocity
nu=0.001;        % viscosity

% Differentiation matrices and integration
[D,DD,wx,x]=dif1D('cheb',0,L,n,3);
Z=zeros(n,n); I=eye(n); 

%{
# Initial guess

We start with a straight cylinder with constant velocity. Now in the system we have three unknowns, the bridge radius as a function of $x$, the axial velocity $u$ as a function of $x$, and a scalar $U$ for the continuation.
%}

%initial guess
u0=x*0+U;
h0=x*0+R;
sol=[u0; h0; U];
dir=[zeros(2*n,1);1];   % initial direction

% Continuation loop
ind=1;quitcon=0;
while ~quitcon 
    
   solprev=sol;
   sol=sol+dir*delta; % extrapolation of solution
    
   % Newton iterations
   quit=0;count=0;
   while ~quit     
    
      u=sol(1:n); h=sol(n+1:2*n); U=sol(2*n+1); 
      hp=D*h;  hpp=DD*h; a=(1+hp.^2);
      up=D*u;  upp=DD*u;
        
%{
# The system

This is an axisymetric flow and the velocity profile should have some radial velocity and as the axial velocity, should be a functino of $x$ and $r$. Instead of this, we use a long wave approximation moe or less similar to the one leading to the Saint-Venant equation, so that we get a 1D equation for $h$ and $u$. This is done in the paper "Eggers, Jens & Dupont, Todd F. 1994 Drop formation in a one-dimensional approximation of the Navier-Stokes equation. J. Fluid Mech. 262, 205–221."

The equations are
$$
\begin{array}{l}
 -uu_x-\frac{p_x}{\rho}+\frac{3\nu(h^2u_x)}{h^2}=0 \\
\sigma \left[ \frac{1}{h(h+h_{x}^2)^{1/2}}-\frac{h_{xx}}{(1+h_x^2)^{3/2}} \right]=p \\
-uh_x-\frac{1}{2}u_xh=0
\end{array}
$$
with $\nu$ the kinematic viscosity.

%}

      % the nonlinear function
      p=h.^-1.*a.^-0.5-hpp.*a.^-1.5;
      
      f=[-u.*up-sig/rho*D*p+3*nu*(D*(h.^2.*up)).*h.^-2-g; ...
         -u.*hp-0.5*up.*h; ...
         dir'*(sol-solprev)-delta];
    
      % analytical Jacobian
      P=diag(-h.^-2.*a.^-0.5) ...
       +diag(3*hp.*hpp.*a.^-2.5-h.^-1.*hp.*a.^-1.5)*D ...
       +diag(-a.^-1.5)*DD; 
 
      A=[-diag(up)-diag(u)*D+diag(6*nu./h.*hp)*D+3*nu*DD, -sig/rho*D*P-diag(6*nu.*up.*hp./h.^2)+diag(6*nu.*up./h)*D, Z(:,1); ...
         -diag(hp)-diag(h)*D/2, -diag(up)/2-diag(u)*D, Z(:,1); ...
         dir'];

%{
# Boundary conditions

We have just the same boundary conditions as in [capillary_bridge_continuation.m](), with just one additional one, the fact that the velocity at inlet (left boundary) should equal to $U$. Since we have an equation that imposes the conservation of mass and the system is stationnary, we do not need to also ompose this at the outlet, this is automatically satisfied.
%}

      % Boundary conditions
      cloc=[1,n+1,2,n];
      f(cloc)=[ ...
          u(1)-U; ...
          h(1)-R; ...
          h(n)-R; ...
          v-pi*wx*h.^2];
      A(cloc,:)=[ ...
        I(1,:),Z(1,:), -1; ...
        Z(1,:),I(1,:), 0; ...
        Z(n,:),I(n,:), 0; ...
        Z(1,:), -2*pi*wx.*h', 0];


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
# Theory and Validation

The straight cylinder is always a stationnary solution of this syetm, this corresponds to a volume$$
v=\pi R^2L
$$
even when there is viscosity, since the surface is a free surface. When we increase progressively the velocity, the bridge tens to become unstable, why? Let us have a look at the waves that can propagate on the straight cylinder. One way to get them is from the dispertion relation of the straight bridge without velocity:
$$
c=\sqrt{\frac{\sigma}{\rho R}}\sqrt{\frac{(2\pi R/\lambda)^2-1}{2}}
$$ 
with $\lambda$ the wavelength and $c$ its celerity. You see that this is a dispersive system (the velocity of waves depends on the wavelength, and that this velocity is smaller for larger wavelengthes. So the slowest wave has its wavelength equal to the length of the capillary bridge $L$. 

Now assume that the velocity through the bridge is exactly equal to the velocity of this slowest wave, then the wave which propagate to the left is stationnary in the laboratory reference frame. This means that exactly for this velocity, there is a nontrivial steady state in addition to the straight cylinder. This is in fact the sign that the is a pitchfork bifurcation at this velocity, and since the bridge was stable for slower velocities, then it means that for velocity larger than this critical one, it becomes unstable. 

The Keller continuation is very good to find folds and no so good to find pichforks, so let's make the system a little unperfect to transform the pitchfork into a fold by having a volume very slightly smaller than that of the cylinder, and also a little viscosity to break the right/left symetry. This way, the continuation should go to the branch with a neck downstream.

On the figure, we show the $Uc$ for the bifurcation for wavelengthes equal to $L$ but also for 3/4, 1/2 and 1/3 (vertical red lines) because depending on the physical parameters you may go to an other branch. You see on the figure that the fold indeed is just below the slowest velocity.

One interesting thing to realize, is that when we increase progressively the length $L$ of the capillary bridge, the velocity at which it becomes unstable decreases until that point when it reaches zero, just when the length becomes $2 pi R$ the perimeter of the cylinder. This means that the Rayleigh-Plateau instability is just a particular case of th dynamic instability of the capillary Venturi(!).

%}

   % Showing the shape of the bridge
   subplot(1,2,1);
   plot(x,h,'b',x,u,'r');axis equal; axis([0,L,0,4*R]); 
   xlabel('x'); ylabel('r'); title('blue:shape red: velocity');
   
   subplot(1,2,2);
   plot(U,min(h),'b.'); 
   xlabel('U'); ylabel('min(h)'); title('Bifurcation branch');
   axis([0,5,0,1.2]);
   grid on; hold on
   
   % theory for validation
   lamb=L./[1,1.5,2,2.5,3,3.5];
   Uc=sqrt(sig/rho/R)*sqrt(((R*2*pi./lamb).^2-1)/2);
   for gre=1:length(Uc);plot([Uc(gre) Uc(gre)],[0,1.2],'r-');end
   drawnow
   
   % new direction
   dir=A\[Z(:,1);Z(:,1);1]; % the new direction
   dir=dir/norm(dir); % normalization
   
   ind=ind+1;
end

% print the figure
set(gcf,'paperpositionmode','auto')
print('-dpng','-r100','capillary_venturi_continuation.png')

%{

On the left panel, we show the shape of the bridge and the velocity through it. Of course, when the radius is small, the velocity must be large because of conservation of mass, this is the Venturi effect.

![The final figure](/sandbox/easystab/capillary_venturi_continuation.png)

# Exercices/Contributions

* Please try to find the other branch of the imperfect pitchfork (for this instead of plotting the branch by drawing *min(h)*, you should rather plot the value of $h$ at a chosen position, because the other branch has its neck going upward instead of downward)
* Please try to find the other imperfect pitchforks corresponding to the other symmetries of waves (for instance the one with wavelength $3L/4$

%}
