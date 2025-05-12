%{ 
# Simulation of Jet stream

Coded by [stab2014/Luis.m]() and [stab2014/PaulValcke.m](). This codes
shows the simulation of a flow jet through a pipe. We solved this problem
using the Unsteady Navier-Stokes equations. Velocity-pressure non-linear
coupling is solved by means of a Newton method, based on the Jacobian
matrix. In this code we use a first-order march in time for integrating
through time.

There exist no analytical solution of this flow, so a theoretical
validation cannot be performed.

# Preliminary thought about the flow

This flow is very interesting from the point of view of stability. A jet
stream flow through a pipe is naturally an unstable flow.
[Transition article from Wikipedia](http://en.wikipedia.org/wiki/Laminar-turbulent_transition)
shows the main idea, and a similar example: the Osborne Reynold's
experiment. At its initial stage, the flow remains stable, and the mixing
layer is horizontal. Then, the jet stream starts oscillating, until its
main structure is broken and several vortical structures arise.

Note that the Unsteady Navier-Stokes equations of this code has been
implemented using the following form:

$$ \partial_tu_i + u_j\partial_ju_i = -\partial_ip +
\frac{1}{Re}\partial^2_{jj}u_i $$

Therefore, the meaning of the Reynolds number will depend on what the user
specifies as dimensionless speed.

%}

clear all; figure('OuterPosition',[0, 0, 800, 400]); format compact

%{
# Parameters

The domain consists on a rectangular mesh, where you can adjust its
dimensions through $$Lx$$ and $$Ly$$. The jet inlet is placed in the
centerline of its length, and its diameter can be modified. Moreover, the
inlet flow imposed has the shape of a Poiseuille's flow.

Simulating high Reynolds numbers is always tricky, but this code implements
a simple way of increasing smoothly (actually linearly) the Reynolds number
from its inital state until a user-provided value. The user also specifies
the time rate of increase of the Reynolds number.
%}
%%%% parameters 
Re=1000; % Initial Reynolds Number.
ReMax = 6000; % Final Reynolds Number
ReAccelTime = 5; % Reynolds increase lapse-time.
Lx = 2; % Length of the pipe
Ly = 0.6; % Height of the pipe
DiamInlet = 0.025; % Diameter of jet stream inlet.
Uinlet = 1; % Flow velocity at inlet. For dimensionless calculation use 1.
Nx=200; % number of grid nodes in x
Ny=150; %number of grid nodes in y
pts=3; % number of points in finite difference stencils
alpha=1; % Under-relaxation. If the code diverges, reducing alpha may favour convergence.
maxIter = 200; % Newton's method maximum iterations
resSTOP = 1e-5; % Newton's residual
PlCase = 'Vorticity'; % User-specified: plots Vorticity or Velocity
dt = 0.1; % Time step
tmax = 30; % Time max


% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,pts);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,pts);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);
D.lap=D.yy+D.xx; % Laplacian

% Stores the indices of the inlet and its corresponding values
y_inlet_ind = find(Y(l.left)<=(Ly+DiamInlet)/2 &...
    Y(l.left)>=(Ly-DiamInlet)/2);
y_inlet = Y(y_inlet_ind);


%{
# Boundary conditions
%}

%%%% preparing boundary conditions
u=(1:NN)'; v=u+NN; p=v+NN;
II=speye(3*NN);

% Condition at the top
dir=[l.ctl;l.left;l.top;l.bot;l.cbl]; % where to put Dirichlet on u and v
dirInlet = y_inlet_ind;
loc=[u(dir); v(dir);
    u(l.right); 
    v(l.right)]; 
    
C=[II([u(dir);v(dir)],:); % ...     % Dirichlet on u,v 
    D.x(l.right,:), Z(l.right,:), Z(l.right,:); ...   % Neuman on u at outflow
    Z(l.right,:),  D.x(l.right,:),Z(l.right,:)]; % ...    % Neumann on v at outflow
%{
We chose an initial guess that statisfies the boundary conditions, this is good for Newton, and it makes it also very easy to impose the nonhomogeneous boundary conditions in the lop.
%}

% initial guess
U = zeros(NN,1);
U(dirInlet) = 1.5.*Uinlet.*(1-((2.*(y_inlet-min(y_inlet)))./DiamInlet - 1).^2);
V=zeros(NN,1);
P=ones(NN,1);
q0=[U(:);V(:);P(:)];

%{
Main loop resolution

Here we perform the march in time. For each time step, a Newton's method
iterative resolution is performed by means of the Jacobian of the
Navier-Stokes equations.
%}
% TIME MARCHING
% Newton iterations
disp('Newton loop')
ReDer = dt*(ReMax-Re)/ReAccelTime; % Calculates the increase rate of Reynolds
q=q0;
qNm1 = q0;
qM=q;
quit=0;count=0;time = 0;
for t = 0:dt:tmax
while ~quit     
 
    % the present solution and its derivatives
    UM = qNm1(u); % Previous time-step solution
    VM = qNm1(v);
    PM = qNm1(p);
    U=alpha.*q(u) + (1-alpha).*qM(u);
    V=alpha.*q(v) + (1-alpha).*qM(v);
    P=alpha.*q(p) + (1-alpha).*qM(p);
    Ux=D.x*U; Uy=D.y*U;
    Vx=D.x*V; Vy=D.y*V; 
    Px=D.x*P; Py=D.y*P;

%{
This is now the heart of the code: the expression of the nonlinear fonction that should become zero, and just after, the expression of its Jacobian, then the boundary conditions.
%}
    % nonlinear function
    f=[UM-U-(U.*Ux+V.*Uy+Px-(D.lap*U)/Re).*dt; ...
       VM-V-(U.*Vx+V.*Vy+Py-(D.lap*V)/Re).*dt; ...
      D.x*U+D.y*V];
    
    % Jacobian 
    A=[-spd(ones(NN,1))+dt.*(-(spd(Ux)+spd(U)*D.x + spd(V)*D.y )+(D.lap)/Re), -spd(Uy).*dt, (-D.x).*dt; ...
        -spd(Vx).*dt,  -spd(ones(NN,1))+dt.*(-(spd(Vy) + spd(V)*D.y + spd(U)*D.x )+(D.lap)/Re), (-D.y).*dt; ...
         D.x, D.y, Z];  
     
    % Boundary conditions 
    f(loc)=C*(q-q0);
    A(loc,:)=C;
   
    % convergence test  
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if (count>maxIter) || (res>1e5); NoConv=1;disp('no convergence');break;else NoConv=0; end
    if res<resSTOP; quit=1; disp('converged'); end
    
    % Newton step
    qM = q;
    q=qM-A\f;   
    count=count+1;
end
count = 0;
qNm1 = q;
quit = 0;
%{
# Plotting the results

For each time-step, a PNG file is saved of wether the Velocity or Vorticity
field (depends on what the user has specified). Then, an animated gif is
done.
%}
fprintf('time=%g\n',time);
    % Plot
    switch PlCase
        case 'Vorticity'
            Vorticity = -D.y*U + D.x*V;
            surf(X,Y,reshape(Vorticity,Ny,Nx)); view(2); shading interp; axis equal; colorbar; grid off;
            caxis([-10 10]);
        case 'Velocity'
            Module = sqrt((U.^2+V.^2));
            surf(X,Y,reshape(Module,Ny,Nx)); view(2); shading interp; axis equal; colorbar; grid off;
    end
    if NoConv, PlCase='***NO CONVERGENCE***'; end
    xlabel('x'); ylabel('y'); title(sprintf('Jet. %s field. t=%0.2fs ; Reynolds %0.0f',PlCase,time,Re));
    drawnow

    set(gcf,'paperpositionmode','auto')
    print('-dpng','-r90',sprintf('Jet_t%0.1fs.png',time));
    
time =time+dt;
if Re < ReMax, Re = Re + ReDer; end
if NoConv, break;end
end

% Creates animated gif
frame = 1;
for t = 0:dt:tmax
   im = imread(sprintf('Jet_t%0.1fs.png',t));
   [imind,cm] = rgb2ind(im,256);
   if frame ==1
      imwrite(imind,cm,'Jet.gif','gif','Loopcount',inf);
   else
       imwrite(imind,cm,'Jet.gif','gif','WriteMode','append',...
           'DelayTime',dt);
   end
   frame = frame+1;
end


%{
# Results

The following results are presented for the case Reynolds 6000, $$Lx= 2$$
and $$Ly = 0.6$$. The resulting velocity and vorticity fields are:

![Jet Stream - Velocity field](http://www.acro3d.com/public/joomla/images/Jet_Vel_Re6000.gif)

![Jet Stream - Vorticity field](http://www.acro3d.com/public/joomla/images/Jet_Vor_Re6000.gif)

There are several remarks that is worth making. First of all, after testing
many cases for different Reynolds and pipe diameters, we found that the
maximum size of the vortex is approximately equal to the diameter. Small
vortices tends to merge between them creating bigger ones. The biggest
vortices remains in place longer and are more stable structures than the
smaller ones, that tend to disappear or more often to merge.

%}
