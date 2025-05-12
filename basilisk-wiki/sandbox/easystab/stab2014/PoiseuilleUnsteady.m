%{
# Poiseuille Unsteady-Navier-Stokes multi-domain test

Coded by [Paul Valcke](./PaulValcke.m) and [Luis Bernardos](./Luis.m). The
objective of this study is to test the Unsteady Navier-Stokes equations in
a multi domain. That is why we have chosen a quite simple flow: Poiseuille.
 The multi-domain constraint has no impact on the final solution, which is
 what we desired. Therefore, we can conclude that this multi-domain
 strategy is suitable for unsteady navier-stokes calculations.

We solve the
laminar unsteady Navier Stokes equations, by means of a Newton method using
the Jacobian of the equation. Hence, the solved equation (appart from
$\partial_iu_i=0$ of course) is:

$$\partial_t + u_j\partial_ju_i = -\partial_ip + \frac{1}{Re}\partial_{jj}^2u_i$$
%}


clear all; figure('OuterPosition',[100, 100, 600, 400]); format compact

%%%% parameters 
Re=1000; % reynolds number
Lx1 = 1; % Length of first domain
Lx2 = 1; % Length of second domain
Ly = 0.25;
DiamInlet = Ly; % Diameter of inlet
Uinlet = Re*1.79e-5/(1.225*DiamInlet);
Nx=30; % number of grid nodes in x
Ny=30; %number of grid nodes in y
pts=5; % number of points in finite difference stencils
alpha=1; % Under-relaxation
maxIter = 200;
resSTOP = 1e-5;
dt = 0.1;
tmax = 6;

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx1,Nx,pts);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,pts);
[D1,l1,X1,Y1,Z1,I1,NN1]=dif2D(d,x,y);

[d.x,d.xx,d.wx,x]=dif1D('fd',Lx1,Lx2,Nx,pts);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,pts);
[D2,l2,X2,Y2,Z2,I2,NN2]=dif2D(d,x,y);

D.x=blkdiag(D1.x,D2.x);
D.y=blkdiag(D1.y,D2.y);
D.xx=blkdiag(D1.xx,D2.xx);
D.yy=blkdiag(D1.yy,D2.yy);
Z=blkdiag(Z1,Z2);

NN=NN1+NN2;

X=[X1(:);X2(:)];
Y=[Y1(:);Y2(:)];

% Stores the indices of the inlet and its corresponding values
y_inlet_ind = find(Y(l1.left)<=(Ly+DiamInlet)/2 &...
    Y(l1.left)>=(Ly-DiamInlet)/2);
y_inlet = Y(y_inlet_ind);

D.lap=D.yy+D.xx;

%{
# Boundary conditions
We invoke the unknown $x_1$ for the first domain, and $x_2$ for the second
domain, where $x=u,v,p$.
%}

%%%% preparing boundary conditions
u1=(1:NN1)';  
u2=u1+NN2;
v1=u2+NN1;
v2=v1+NN2;
p1=v2+NN1;
p2=p1+NN2;
II1=speye(3*NN1);
II2=speye(3*NN2);
II=blkdiag(II1,II2);

% Dirichlet Conditions
dir1=[l1.ctl;l1.left;l1.top;l1.bot;l1.cbl]; % where to put Dirichlet on u1 and v1
lnk1 = l1.right; % Link
loc1=[u1(dir1); v1(dir1); % Dirichlet
    u1(lnk1); v1(lnk1); p1(lnk1)]; % Link

dir2=[l2.ctl;l2.top;l2.bot;l2.cbl]; % where to put Dirichlet on u2 and v2
lnk2 = l2.left;% Link
neu = l2.right; % Neumann
loc2=[u2(dir2); v2(dir2); % Dirichlet
      u2(lnk2); v2(lnk2); p2(lnk2); % Link
      u2(neu);  v2(neu)]; % Neumann

lnk = [lnk1;lnk2];
loc = [loc1;loc2];
dir = [dir1;dir2];

%{
# Multi-domain Link condition
The link condition between the two domains is performend in the Constraint
matrix $C$. Therefore, we make use of the Identity matrix, to impose that
the overlapping points must be equal. Similarly, we make use of the $D_x$
differentiation matrix in order to impose the continuity of the derivative
normal to the overlapping boundary: we call this "Link Tangency", and has
to be accomplished for the three unknowns (i.e, $u$, $v$, and $p$).
%}
C=[II([u1(dir1);v1(dir1);u2(dir2);v2(dir2)],:);% Dirichlet on u,v 
    II([u1(lnk1);v1(lnk1);p1(lnk1)],:)-II([u2(lnk2);v2(lnk2);p2(lnk2)],:); % Link Continuity
    D.x(lnk1,:)-D.x(lnk2+NN1,:), Z(lnk1,:), Z(lnk1,:); % Link Tangency on u
    Z(lnk1,:),D.x(lnk1,:)-D.x(lnk2+NN1,:), Z(lnk1,:); % Link Tangency on v
    Z(lnk1,:), Z(lnk1,:),D.x(lnk1,:)-D.x(lnk2+NN1,:); % Link Tangency on p
    D.x(neu+NN1,:), Z(neu+NN1,:), Z(neu+NN1,:); % Neumann on u2 at outflow
    Z(neu+NN1,:), D.x(neu+NN1,:), Z(neu+NN1,:)];  % Neumann on v2 at outflow

%{
We chose an initial guess that statisfies the boundary conditions, this is good for Newton, and it makes it also very easy to impose the nonhomogeneous boundary conditions in the lop.
%}

% initial guess
U1 = zeros(NN1,1);
U1(l1.left) = Uinlet; % We impose uniform speed. This way we will see the estabilish length of the flow.
V1=zeros(NN1,1);
P1=ones(NN1,1);

U2 = zeros(NN2,1);
V2=zeros(NN2,1);
P2=ones(NN2,1);

U = [U1;U2];
V = [V1;V2];
P = [P1;P2];
q0=[U1(:);U2(:);V1(:);V2(:);P1(:);P2(:)];

% Newton iterations
disp('Newton loop')
q=q0;
qNm1 = q0;
qM=q;
quit=0;count=0;time = 0;
for t = 0:dt:tmax
while ~quit     
 
    % the present solution and its derivatives
    UM = qNm1([u1;u2]); % Solution du pas précédant
    VM = qNm1([v1;v2]);
    PM = qNm1([p1;p2]);
    U=alpha.*q([u1;u2]) + (1-alpha).*qM([u1;u2]);
    V=alpha.*q([v1;v2]) + (1-alpha).*qM([v1;v2]);
    P=alpha.*q([p1;p2]) + (1-alpha).*qM([p1;p2]);
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
    if (count>maxIter) || (res>1e5); disp('no convergence');break; end
    if res<resSTOP; quit=1; disp('converged'); end
    
    % Newton step
    qM = q;
    q=qM-A\f;%gmres(A,f,1000,1e-7);%A\f;   
    count=count+1;
    
%     % NewtonPlot
%     Module = sqrt((U.^2+V.^2));
%     surf([X1,X2],[Y1,Y2],reshape(Module,Ny,2*Nx)); view(2); shading interp; hold on
%     xlabel('x'); ylabel('y');grid off;hold off
%     colorbar;
%     axis equal
%     drawnow
end
count = 0;
qNm1 = q;
quit = 0;
fprintf('t=%g\n',time);
    % plotting
    Module = sqrt((U.^2+V.^2));
    surf([X1,X2],[Y1,Y2],reshape(Module,Ny,2*Nx)); view(2); shading interp; hold on
    caxis([0 1.5*Uinlet]); % Color scale
    xlabel('x'); ylabel('y'); title(sprintf('Poiseuille t=%0.1fs ; Uinlet = %0.3fm/s; Reynolds %d',time,Uinlet,Re)); grid off;hold off
    colorbar;
    axis equal
    drawnow
    set(gcf,'paperpositionmode','auto')
    print('-dpng','-r80',sprintf('Poiseuille_t%0.1fs.png',time));
time =time+dt;
end

% Creates animated gif
FirstFrame = true;
gifname = sprintf('Poiseuille_Re%d.gif',Re);
for t = 0:dt:tmax
   im = imread(sprintf('Poiseuille_t%0.1fs.png',t));
   [imind,cm] = rgb2ind(im,256);
   if FirstFrame
      imwrite(imind,cm,gifname,'gif','Loopcount',inf);
      FirstFrame = false;
   else
       imwrite(imind,cm,gifname,'gif','WriteMode','append',...
           'DelayTime',dt);
   end
end




%{
# Results

![Poiseuille flow](./Poiseuille_Re1000.gif)

The flow successfully establishes and we can conclude that the multi-domain
calculation for 2D Unsteady-Navier-Stokes is feasible.
%}