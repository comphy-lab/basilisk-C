
clear all; clf; format compact;
warning('off')
%{ 

# Cartesian ,UNStationnary, Uncompressible, Navier-Stokes with Jacobian resolution

Coded by [stab2014/Luis.m]() and [stab2014/PaulValcke.m](), on the initial [NavierStokesStationnary]() code.

The structure of the code is the same as StationnaryNS. We changed the
expression of the function, the Jacobi and added an Euler March in Time
(first order).
If you want to understand the structure of the code, take a look at
StationnaryNS, then take a look at the comment of the function and the
jacobi

# Comments

We have no test for the moment, if you have any or experiment to compare, please contribute

You can use this program into find stationnary solutions, it will slowly
goes to the solution (if it exist). If you have a very complicate flow to
calculate, it might be a good way.

As the fluid is uncompressible, if it is laminar (then a stationnary
solution exist), a Poiseuille flow will be instantly created in right
conditions. A facingStep is a better test.
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHYSICAL PARAMETERS
Re=2000; % reynolds number

%GEOMETRIC PARAMETERS
Lx = 2;
Ly = 0.6;
DiamInlet = 0.025; % Diameter of inlet
tmax = 100; %WE ADD A MAXIMUM TIME

%NUMERICAL PARAMETERS=
Nx=200; % number of grid nodes in x
Ny=100; %number of grid nodes in y
dt = 0.1; %STEP BETWEEN TWO TIMES

maxIter = 200; % Maximum number of iteration for Newton method.
resSTOP = 1e-5; % Precision required to stop convergence

%FIXED PARAMETERS OR INTERMEDIATE FORMULAS
pts=5; % number of points in finite difference stencils
alpha=1; % No Under-relaxation
Uinlet = Re*1.79e-5/(1.225*DiamInlet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DIFFERENTIATION MATRIXES%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIATION 

[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,pts);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,pts);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% PLAN TRANSFORMATION
%You can transform X and Y as done in geometry
 D=map2D(X,Y,D);
 D.lap=D.xx+D.yy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BOUNDARIES GESTION%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% POSITIONS OF EVERY VECTORS
u=(1:NN)'; v=u+NN; p=v+NN;
II=speye(3*NN);

%%%% EXAMPLE OF DIRICHLET ENTRY : POISEUILLE FLOW ON ONE PART OF ONE BOUNDARY OF THE DOMAIN
%Here this over the height DiamInlet
y_inlet_ind = find(Y(l.left)<=(Ly+DiamInlet)/2 &...
Y(l.left)>=(Ly-DiamInlet)/2);
y_inlet = Y(y_inlet_ind);

% BOUNDARY POSITIONS
dir=[l.ctl;l.ctr;l.cbl;l.cbr;l.left;l.top;l.bot;]; % where to put Dirichlet on u and v
dirInlet = y_inlet_ind;

% Loc will be used to apply boundary condition to the function 
loc=[u(dir); v(dir); 
    u(l.right); ...
    v(l.right)];

%MATRIX OF CONSTRAINTS   
C=[II([u(dir);v(dir)],:); % ...     % Dirichlet on u,v 
    D.x(l.right,:), Z(l.right,:), Z(l.right,:); ...   % Neuman on u at outflow
    Z(l.right,:),  D.x(l.right,:),Z(l.right,:)]; % ...    % Neumann on v at outflow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INITIAL GUESS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = 0*ones(NN,1);
U(dirInlet) = 1.5.*Uinlet.*(1-((2.*(y_inlet-min(y_inlet)))./DiamInlet - 1).^2);

%Initial condition on V
V=zeros(NN,1);
V(y_inlet_ind) = 0;

%Initial condition on P
P=ones(NN,1);

%Initial Solution
sol0=[U(:);V(:);P(:)];


%%%%%%%%%%%%%%%%%RESOLUTION ALGORITHM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEWTON ALGORITHM
sol=sol0;
solM=sol;
solNm1 = sol0; % WE TAKE THE SOLUTION OF PREVIOUS EQUATION

quit=0;
count=0;
for t = 0:dt:tmax
    fprintf('t=%g\n',t);
while ~quit     
 
%%%% RELAXATION METHOD
% PREVIOUS TERMS OF THE RESOLUTION
    UM = solNm1(u);
    VM = solNm1(v);
    PM = solNm1(p);
    
% ACTUAL TERMES    
    U=alpha.*sol(u) + (1-alpha).*solM(u);
    V=alpha.*sol(v) + (1-alpha).*solM(v);
    P=alpha.*sol(p) + (1-alpha).*solM(p);
    
% DEFINITION OF THE DERIVATES    
    Ux=D.x*U; Uy=D.y*U;
    Vx=D.x*V; Vy=D.y*V; 
    Px=D.x*P; Py=D.y*P;

%DEFINITION OF THE FUNCTIONS YOU WANT TO RESOLVE
% Here we solve Navier stokes in 2D, (conservation of the impulsion on X,Y,
% and conservation of mass). You just replace spatial operator with their
% equivalent in matrixes, the temporal derivative is replaced by (U-UM)/dt.
% Instead of a division by dt, because we integrate the equation, we
% multiply everything by dt.
    f=[UM-U-(U.*Ux+V.*Uy+Px-(D.lap*U)/Re).*dt; ...
       VM-V-(U.*Vx+V.*Vy+Py-(D.lap*V)/Re).*dt; ...
      D.x*U+D.y*V];
    
%JACOBIAN OF THE EQUATIONS
% We take the previous Jacobi, we just multiply it with dt (like f). We add
% an Idendity because of the $\frac{d(U-Um)}{dU} (resp, $\frac{d(V-Vm)}{dV}) 
    A=[-spd(ones(NN,1))+dt.*(-(spd(Ux)+spd(U)*D.x + spd(V)*D.y )+(D.lap)/Re), -spd(Uy).*dt, (-D.x).*dt; ...
        -spd(Vx).*dt,  -spd(ones(NN,1))+dt.*(-(spd(Vy) + spd(V)*D.y + spd(U)*D.x )+(D.lap)/Re), (-D.y).*dt; ...
         D.x, D.y, Z];  
     
%BOUNDARY CONDITIONS
    f(loc)=C*(sol-sol0);
    A(loc,:)=C;
   
%CONVERGENCE TEST
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if (count>maxIter) || (res>1e5); disp('no convergence');break; end
    if res<resSTOP; quit=1; disp('converged'); end
    
%NEWTON STEP FOR RESOLUTION
    solM = sol;
    sol=solM-A\f;   
    count=count+1;
   
end


solNm1 = sol;
count = 0;
quit = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Place here figures you want at each step of the iterations

%Plot of the Module
    figure(1)
    Module = sqrt((U.^2+V.^2));
    surf(X,Y,reshape(Module,Ny,Nx)); view(2); shading interp; hold on
    axis equal
    quiver3(X,Y,ones(Ny,Nx).*max(Module),reshape(U,Ny,Nx),reshape(V,Ny,Nx),...
            zeros(Ny,Nx),'k');    
    xlabel('x'); ylabel('y'); title('Poiseuille Flow'); grid off;hold off
    colorbar;
    drawnow
    
    set(gcf,'paperpositionmode','auto')
    print('-dpng','-r80',sprintf('Jet_t%0.1fs.png',t));
end



%Place here figures you only want once at the end

% Creates animated gif
frame = 1;
for t = 3*dt:dt:tmax
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



