%{ 
# Simulation of Unsteady Flow around a Cylinder

Coded by [Paul Valcke](./PaulValcke.m) and [Luis Bernardos](./Luis.m). This codes
shows the simulation of a flow around a cylinder. We solved this problem
using the Unsteady Navier-Stokes equations. Velocity-pressure non-linear
coupling is solved by means of a Newton method, based on the Jacobian
matrix. In this code we use a first-order march in time for integrating
through time.

There exist no analytical solution of this flow, so a theoretical
validation cannot be performed.

# Preliminary thoughts about the flow

The flow around a cylinder is a classic problem in fluid dynamics. There
exist no analytical solution for this unsteady flow. We intend to simulate
the famous Von Karman street, the vortex shedding downstream. However, we
did not succeed in simulating this instability, and therefore we encourage
future contributions. 

Indeed, phyisically the vorticity areas around the cylinder becomes unstable with time.
At first, a symmetrical quasi-steady flow configuration is stablished, and
then it will become unstable and asymmetrical. With this code, for some reason that the
authors didn't manage to figure out, when the instabilities occur and they
are amplified, the calculation diverges, or the vortex are dissipated. We
will talk about in more detail in the results section.

Note that the Unsteady Navier-Stokes equations of this code has been
implemented using the following form:

$$ \partial_tu_i + u_j\partial_ju_i = -\partial_ip +
\frac{1}{Re}\partial^2_{jj}u_i $$

%}

clear all; figure('OuterPosition',[0, 0, 800, 600]); format compact

%{
# Parameters

Here the user inputs all the user-specified parameters for calculations. It
can be stated to make Steady or Unsteady calculations, to plot Velocity
field (with vectors) or the Vorticity field.
%}
%%%% parameters 
Re=1000; % Starting reynolds number
ReMax = 6000; % Stablished reynolds number
ReAccelTime = 2; % Reynolds stablishing time 
Uinlet=1;
Rmax = 1; % Radius of the domain (outter side)
Rmin = 0.01; % Radius of the domain (inner side, the cylinder)
Ntheta=181; % number of grid nodes in theta
Nrho=120; %number of grid nodes in rho
pts=3; % number of points in finite difference stencils
alpha=1; % Under-relaxation
maxIter = 200;
resSTOP = 1e-5;
CalcMode = 'Unsteady'; % Steady or Unsteady
PlCase = 'Vorticity'; % Vorticity or Velocity
dt = 0.1; % Time step
tmax = 6; % Time of the simulation

Lx = 2*pi*Rmax;
Ly = Rmax-Rmin;

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('fd',0,1,Ntheta,pts);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,1,Nrho,pts);
[D,l,~,~,Z,I,NN]=dif2D(d,x,y);
%{
# Mesh
We use a cylindrical mesh. Therefore, we have to build it by means of a
$\rho$ and a $\theta$ vector. The spacing of these two vectors will let us to
shrink the grid where detailed calculations are required (and typically
high gradients exist). This is specially true at the boundary layer of the
cylinder, so that is why whe use the function [Lspacing](./Lspacing.m),
that imposes a trigonometric discretization in a given direction controlled
by a $lambda$ parameter. Then, we make use of [map2D](./map2D.m) function
for transforming the rectangular-created differentiation matrices into the
new deformed grid domain.

![Mesh](./Cylinder_Mesh.png)

%}
% Meshing
Tspac =  linspace(0,2*pi,Ntheta); % Spacing of theta, linear
Rhospac = Lspacing(Nrho,2,Rmax-Rmin); % Spacing of Rho, trigonometric
[Theta, Rho] = meshgrid(Tspac,Rmin + Rhospac);
[X, Y] = pol2cart(Theta,Rho); % We transfor polar to cartesian coordinates
D=map2D(X,Y,D); % We transform the differentiation matrices
D.lap=D.yy+D.xx; % Laplacian calculation


% We look for the corresponding indices of velocity_inlet and outflow:
velocity_inlet = l.top(round(0.25*numel(l.top)):round(0.75*numel(l.top)));
outflow = l.top(l.top<min(velocity_inlet) | l.top>max(velocity_inlet) );
% % Plots mesh
mesh(X,Y,X*0,'EdgeColor','k','DisplayName','Mesh');hold on; view(2); axis equal; grid off;
xlabel('x'); ylabel('y'); title(sprintf('Cylinder mesh. N_{\\theta}=%d  N_{\\rho}=%d',Ntheta,Nrho));
plot(X(l.bot),Y(l.bot),'ob','DisplayName','Cylinder-wall');
plot(X(l.left),Y(l.left),'or','DisplayName','Link');
plot(X(l.right),Y(l.right),'+r','DisplayName','Overlapping Link');
plot(X(velocity_inlet),Y(velocity_inlet),'og','DisplayName','Velocity-inlet')
plot(X(outflow),Y(outflow),'+m','DisplayName','Outflow')
hold off
legend('show','location','NorthWest')
set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','Cylinder_Mesh.png');

%{
# Boundary conditions
The boundary conditions imposed for this problem are the following: 

- Velocity Inlet: Dirichlet conditions, horizontal speed equal to Uinlet.
- Outflow: Neumann condition on velocity.
- Cylinder Wall: Dirichlet condition, no-slip (zero velocity).

Furthermore, we included a link condition at the overlapping mesh boundary.
This link condition is implemented in the Constraint C matrix, by means of
a continuity condition (using identity matrix) and continuity of the first
spacial derivatives, or tangency, using the differentiation matrices.
%}

%%%% preparing boundary conditions
u=(1:NN)'; v=u+NN; p=v+NN;
II=speye(3*NN);

%%%%%%%%%%%%%%
DD=blkdiag(D.y,D.y,D.y);
% Boundary Conditions
dir=[velocity_inlet;l.bot]; % where to put Dirichlet on u and v
lnk1=[u(l.left); v(l.left); p(l.left)]; % Link indices
lnk2=[u(l.right); v(l.right); p(l.right)]; % Overlapping link indices
p0loc=1;
loc=[u(dir); v(dir); p(p0loc); ... % Dirichlet
    lnk1;lnk2; % Link
    u(outflow);v(outflow)]; % Neumann
    
C=[II([u(dir);v(dir);p(p0loc)],:);% Dirichlet on u,v at walls 
    II(lnk1,:)-II(lnk2,:); % Continuity of the value
    DD(lnk1,:)-DD(lnk2,:); % Continuity of the normal derivative
    D.x(outflow,:), Z(outflow,:), Z(outflow,:); % Neumann on u2 at outflow
    Z(outflow,:), D.x(outflow,:), Z(outflow,:);  % Neumann on v2 at outflow
    ];


%{
We chose an initial guess that statisfies the boundary conditions, this is good for Newton, and it makes it also very easy to impose the nonhomogeneous boundary conditions in the lop.
%}

% initial guess
%U=Uinlet*(Rho-Rmin)/(Rmax-Rmin); % Alternatively, use this initial
%condition. Not really required.
U = zeros(NN,1);
U(velocity_inlet)=Uinlet; % Here we impose the velocity inlet speed.
V=0*X;
P=0*X;
q0=[U(:);V(:);P(:)];

%{
# Forced Perturbation

In an attempt to destabilize the trailing vortices, we implemented a simple
finite perturbation through modifying the imposed vertical velocity at the velocity
inlet boundary. Even making this, we did not manage to sufficiently
destabilize the vortices and simulate the Von Karman street. The user and
future contributor may find useful using this perturbation.


%}

time_pert_start = 6; % Starting time of the perturbation 
time_pert_end = 6.5; % End time of the perturbation.
Vpert = V;
Vpert(velocity_inlet)=0; % Modify this parameter for perturbating the velocity field from time_pert_start to time_pert_end
q0pert=[U(:);Vpert(:);P(:)];

switch CalcMode
    case 'Unsteady'
% MARCHE EN TEMPS
% Newton iterations
ReDer = dt*(ReMax-Re)/ReAccelTime;
disp('Newton loop')
q=q0;
qNm1 = q0;
qM=q;
quit=0;count=0;time = 0;
for t = 0:dt:tmax
while ~quit     
 
    % the present solution and its derivatives
    UM = qNm1(u); % Solution du pas précédant
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
    
    % Boundary conditions, including the perturbation condition
    if time >= time_pert_start && time <= time_pert_end
        f(loc)=C*(q-q0pert);
    else
        f(loc)=C*(q-q0);
    end
    A(loc,:)=C;
   
    % convergence test  
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if (count>maxIter) || (res>1e5); NoConv=1;disp('no convergence');break;else NoConv=0; end
    if res<resSTOP; quit=1; disp('converged'); end
    
    % Newton step
    qM = q;
    q=qM-A\f;% alternatively replace A\f by gmres(A,f,1000,1e-7); 
    count=count+1;
   
end
count = 0;
qNm1 = q;
quit = 0;

fprintf('time=%g\n',time);
    % Plot
    switch PlCase
        case 'Vorticity'
            Vorticity = -D.y*U + D.x*V;
            surf(X,Y,reshape(Vorticity,Nrho,Ntheta)); view(2); shading interp; axis equal; colorbar; grid off;
            caxis([-50 50]);
            xlow=-2*Rmin; ylow=-4*Rmin; xhigh=8*Rmin; yhigh=+4*Rmin; % Figure view point
        case 'Velocity'
            Module = sqrt((U.^2+V.^2));
            surf(X,Y,reshape(Module,Nrho,Ntheta),'facealpha',0.8); view(2); shading interp; axis equal; colorbar;
            hold on
            uu=reshape(U,Nrho,Ntheta); vv=reshape(V,Nrho,Ntheta);
            quiver(X(1:5:end,1:5:end),Y(1:5:end,1:5:end),uu(1:5:end,1:5:end),vv(1:5:end,1:5:end)); hold off
            xlow=-2*Rmin; ylow=-4*Rmin; xhigh=8*Rmin; yhigh=+4*Rmin; % Figure view point
    end
    xlim(gca, [xlow xhigh])
    ylim(gca, [ylow yhigh])
    axis(gca, [xlow xhigh ylow yhigh])
    if NoConv, PlCase='***NO CONVERGENCE***'; end
    xlabel('x'); ylabel('y'); title(sprintf('Cylinder. %s field. t=%0.2fs ; Reynolds %0.0f',PlCase,time,Re));
    drawnow

    set(gcf,'paperpositionmode','auto')
    print('-dpng','-r80',sprintf('Cylinder_t%0.1fs.png',time));
time =time+dt;
if Re < ReMax, Re = Re + ReDer; end
if NoConv, break;end
end

% Creates animated gif
FirstFrame = true;
gifname = sprintf('Cylinder_Re%d.gif',Re);
for t = 0:dt:tmax
   im = imread(sprintf('Cylinder_t%0.1fs.png',t));
   [imind,cm] = rgb2ind(im,256);
   if FirstFrame
      imwrite(imind,cm,gifname,'gif','Loopcount',inf);
      FirstFrame = false;
   else
       imwrite(imind,cm,gifname,'gif','WriteMode','append',...
           'DelayTime',dt);
   end
end
    case 'Steady'
% % STEADY
% Newton iterations
disp('Newton loop')
q=q0;
qM=q;
quit=0;count=0;
while ~quit     
 
    % the present solution and its derivatives
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
    f=[-(U.*Ux+V.*Uy)+(D.lap*U)/Re-Px; ...
       -(U.*Vx+V.*Vy)+(D.lap*V)/Re-Py; ...
      D.x*U+D.y*V];
    
    % Jacobian 
    A=[-( spd(Ux) + spd(U)*D.x + spd(V)*D.y )+(D.lap)/Re, -spd(Uy), -D.x; ...
        -spd(Vx),  -(spd(Vy) + spd(V)*D.y + spd(U)*D.x )+(D.lap)/Re, -D.y; ...
         D.x, D.y, Z];  
     
    % Boundary conditions 
    f(loc)=C*(q-q0);
    A(loc,:)=C;

    % plotting
    Module = sqrt((U.^2+V.^2));
    surf(X,Y,reshape(Module,Nrho,Ntheta)); view(2); shading interp; axis equal
%     quiver3(X,Y,ones(Ny,Nx).*max(Module),reshape(U,Ny,Nx),reshape(V,Ny,Nx),...
%         zeros(Ny,Nx),'k','AutoScale','off');
    xlabel('x'); ylabel('y'); title('Cylinder Steady'); grid off;hold off
    colorbar;
    drawnow
    
    % convergence test  
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if (count>maxIter) || (res>1e5); disp('no convergence');break; end
    if res<resSTOP; quit=1; disp('converged'); continue; end
    
    % Newton step
    qM = q;
    q=qM-A\f;   
    count=count+1;
end
    otherwise
        warning('Mode "%s" not recognized',CalcMode);
end

%{
# Results

As stated previously, we did not manage to simalate the Von Karman Street.
Here we present how the vorticity field is stablished and respects the
symmetry: 

![Symmetric vorticity stablishing](http://www.acro3d.com/public/joomla/images/cylinder_re6000.gif)

We tested many different Reynolds number and conditions. For example, the
following image shows a higher Reynolds number simulation (note: the label
says Reynolds 6000, but the geometry of the cylinder was different, so the
meaning of the Reynolds is not the same).

![Unstable vortices](http://www.acro3d.com/public/joomla/images/cylinder_re6000_vortexdissipation.gif)

The previous image shows how the vortices are destabilized and are
eventually dissipated. We think that the outflow boundary was too close,
and perhaps using a bigger mesh would be enough for simulating the Von
Karman street. Future contributors may try using a bigger mesh.

In the sake of variety, and in an attempt to make understand the future
contributor the problems that we encountered, we include also an image of a
simulation that diverges. We can clearly see that when the vortices start
to destibilize and potentially to "dettach", the code explodes. We did not
find the reason.

![Diverging simulation](http://www.acro3d.com/public/joomla/images/cylinder_re15000_div.gif)


%}