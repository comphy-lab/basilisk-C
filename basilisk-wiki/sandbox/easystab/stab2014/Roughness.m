%{

# UNFINISHED Roughness simulation page

Coded by [Paul Valcke](./PaulValcke.m) and [Luis Bernardos](./Luis.m). This
code has not been finished, it is an attempt to calculate the roughness
effect in a multi-domain mesh. The code programmatically creates several
meshes and imposes the link conditions between all domains also
programmatically. Through the parameter $N_{cav}$ The user can specify as
many cavities as desired, and the code will mesh and construct the C matrix
programmatically.

We did not have enough time for accomplishing this interesting job, and we
hope future contributors will like the idea and continue working in this
direction.

![The Mesh](./RoughnessMesh.png)
%}

clear all; figure('OuterPosition',[100, 100, 800, 400]); format compact

%%%% parameters 
Re=1000; % reynolds number
Lin = 1; % Length before roughness
DiamInlet = 1; % Diameter of main tube
Lcav = 0.2; % Length of roughness 
Dcav = 0.2; % Depth of roughness
Ncavs = 10; % Number of cavities defining roughness

Ny_main = 30; % Number of vertical divisions of main tube
Nx_inlet = 10; % Number of horizontal divisons of inlet
Ny_cav = 10; % Number of vertical divisons of cavity
Nx_cav = 10; % Number of horizontal divisons of cavity.

% Speed of inlet
Uinlet = Re*1.79e-5/(1.225*DiamInlet);
pts=5; % number of points in finite difference stencils
alpha=1; % Under-relaxation
maxIter = 200;
resSTOP = 1e-5;
dt = 0.1;
tmax = 6;

% differentiation
D = cell(3*Ncavs+1,1); l=D;X=D;Y=D;Z=D;I=D; NN = zeros(3*Ncavs+1,1);% Pre-allocation
[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lin,Nx_inlet,pts);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,DiamInlet,Ny_main,pts);
[D{1},l{1},X{1},Y{1},Z{1},~,NN(1)]=dif2D(d,x,y);

G = 0;
for i = 1:3:3*Ncavs
    G = G+1;
[d.x,d.xx,d.wx,x]=dif1D('fd',Lin+(G-1)*2*Lcav,Lcav,Nx_cav,pts);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,DiamInlet,Ny_main,pts);
[D{i+1},l{i+1},X{i+1},Y{i+1},Z{i+1},~,NN(i+1)]=dif2D(d,x,y);

[d.x,d.xx,d.wx,x]=dif1D('fd',Lin+(G-1)*2*Lcav,Lcav,Nx_cav,pts);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,-Dcav,Ny_cav,pts);
[D{i+2},l{i+2},X{i+2},Y{i+2},Z{i+2},~,NN(i+2)]=dif2D(d,x,y);

[d.x,d.xx,d.wx,x]=dif1D('fd',Lin+(G-1)*2*Lcav+Lcav,Lcav,Nx_cav,pts);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,DiamInlet,Ny_main,pts);
[D{i+3},l{i+3},X{i+3},Y{i+3},Z{i+3},~,NN(i+3)]=dif2D(d,x,y);
end

Dx=D{1}.x;
Dy=D{1}.y;
Dxx=D{1}.xx;
Dyy=D{1}.yy;
ZZ=Z{1};
for i = 2:3*Ncavs+1
Dx=blkdiag(Dx,D{i}.x);
Dy=blkdiag(Dy,D{i}.y);
Dxx=blkdiag(Dxx,D{i}.xx);
Dyy=blkdiag(Dyy,D{i}.yy);
ZZ=blkdiag(ZZ,Z{i});
end
Dlap=Dyy+Dxx;
NNsum=sum(NN);

% Plot mesh
for i = 1:3*Ncavs+1
    mesh(X{i},Y{i},0*X{i},'edgecolor','k'); view(2); grid off;hold on
end
xlabel('X');ylabel('Y'),title('Mesh of the roughness');hold off;axis equal;
drawnow

set(gcf,'paperpositionmode','auto')
print('-dpng','-r90','RoughnessMesh.png');
break

%{
# Boundary conditions
%}

%%%% preparing boundary conditions
u = cell(3*Ncavs+1,1); v=u; p=u; % pre-allocates
u{1}=(1:NN(1))';
for i = 2:3*Ncavs+1
    u{i}=u{i-1}+NN(i);
end
v{1}=u{end}+NN(1);
for i = 2:3*Ncavs+1
    v{i}=v{i-1}+NN(i);
end
p{1}=v{end}+NN(1);
for i = 2:3*Ncavs+1
    p{i}=p{i-1}+NN(i);
end
II=speye(3*NNsum);


% Dirichlet Conditions
dir = cell(3*Ncavs+1,1);lnkMM=dir;lnkM=dir;lnkP=dir;lnkPP=dir;loc=dir; % Preallocates
dir{1}=[l{1}.ctl;l{1}.left;l{1}.cbl;l{1}.top;l{1}.bot;l{1}.ctr;l{1}.cbr]; % where to put Dirichlet on u1 and v1
lnkMM{1}=[]; % Link with previous of the previous domain
lnkM{1}=[]; % Link with previous domain
lnkP{1}= l{1}.right; % Link with next domain
lnkPP{1}=[]; % Link with next to the next domain
loc{1}=[u{1}(dir{1}); v{1}(dir{1}); % Dirichlet
    u{1}(lnkP{1}); v{1}(lnkP{1}); p{1}(lnkP{1})]; % Link

G = 0;
for i = 2:3:3*Ncavs+1
    G = G+1;
    dir{i}=[l{i}.ctl;l{i}.top;l{i}.cbl;l{i}.cbr]; % where to put Dirichlet on u v
    lnkMM{i}=[];
    lnkM{i}=l{i}.left;  % Link with previous domain
    lnkP{i}=l{i}.bot; % Link with next domain
    lnkPP{i}=l{i}.right; % Link with next to the next domain
    loc{i}=[u{i}(dir{i}); v{i}(dir{i}); % Dirichlet
            u{i}(lnkM{i}); v{i}(lnkM{i}); p{i}(lnkM{i});
            u{i}(lnkP{i}); v{i}(lnkP{i}); p{i}(lnkP{i});
            u{i}(lnkPP{i}); v{i}(lnkPP{i}); p{i}(lnkPP{i})]; % Link
        
    dir{i+1}=[l{i+1}.ctl;l{i+1}.left;l{i+1}.cbl;l{i+1}.bot;l{i+1}.cbr;l{i+1}.right;l{i+1}.ctr]; % where to put Dirichlet on u v
    lnkMM{i+1}=[];
    lnkM{i+1}=l{i+1}.top; % Link
    lnkP{i+1}=[]; % Link with next domain
    lnkPP{i+1}=[]; % Link with next to the next domain    
    loc{i+1}=[u{i+1}(dir{i+1}); v{i+1}(dir{i+1}); % Dirichlet
            u{i+1}(lnkM{i+1}); v{i+1}(lnkM{i+1}); p{i+1}(lnkM{i+1})]; % Link
    
    dir{i+2}=[l{i+2}.ctl;l{i+2}.cbl;l{i+2}.bot;l{i+2}.cbr;l{i+2}.ctr]; % where to put Dirichlet on u v
    lnkMM{i+2}=l{i+2}.left;
    lnkM{i+2}=[]; % Link
    lnkP{i+2}=[];
    lnkPP{i+2}=[];
    loc{i+2}=[u{i+2}(dir{i+2}); v{i+2}(dir{i+2}); % Dirichlet
            u{i+2}(lnkMM{i+2}); v{i+2}(lnkMM{i+2}); p{i+2}(lnkMM{i+2})]; % Link
        
    if G == Ncavs
        neu=l{i+2}.right;
        loc{i+2}=[loc{i+2}; u{i+2}(neu); v{i+2}(neu)];
    else
        lnkP{i+2}=l{i+2}.right;
        loc{i+2}=[loc{i+2};u{i+2}(lnkP{i+2}); v{i+2}(lnkP{i+2}); p{i+2}(lnkP{i+2})];
    end
end

loc = cell2mat(loc(:));
% =================================
% CONSTRUCTION OF CONSTRAINT MATRIX
% =================================

C = II(u{1}(dir{1}),:);
% Dirichlet conditions:
for i = 2:3*Ncavs+1
    C = [C;II([u{i}(dir{i});v{i}(dir{i})],:)];
end
% Link Continuity
for i = 1:3*Ncavs
C=[C;II([u{i}(lnkP{i});v{i}(lnkP{i});p{i}(lnkP{i})],:)-II([u{i+1}(lnkM{i+1});v{i+1}(lnkM{i+1});p{i+1}(lnkM{i+1})],:);

% Link Tangency
     % on u
[Dx(lnkP{i}+sum(NN(1:i-1)),:)-Dx(lnkM{i+1}+sum(NN(1:i)),:),ZZ(lnkP{i}+sum(NN(1:i-1)),:),ZZ(lnkP{i}+sum(NN(1:i-1)),:)];
     % on v
[ZZ(lnkP{i}+sum(NN(1:i-1)),:),Dx(lnkP{i}+sum(NN(1:i-1)),:)-Dx(lnkM{i+1}+sum(NN(1:i)),:),ZZ(lnkP{i}+sum(NN(1:i-1)),:)];
     % on p
[ZZ(lnkP{i}+sum(NN(1:i-1)),:),ZZ(lnkP{i}+sum(NN(1:i-1)),:),Dx(lnkP{i}+sum(NN(1:i-1)),:)-Dx(lnkM{i+1}+sum(NN(1:i)),:)]];
if i<3*Ncavs
C=[C;II([u{i}(lnkPP{i});v{i}(lnkPP{i});p{i}(lnkPP{i})],:)-II([u{i+2}(lnkMM{i+2});v{i+2}(lnkMM{i+2});p{i+2}(lnkMM{i+2})],:); 
     [Dx(lnkPP{i}+sum(NN(1:i-1)),:)-Dx(lnkMM{i+2}+sum(NN(1:i)),:),ZZ(lnkPP{i}+sum(NN(1:i-1)),:),ZZ(lnkPP{i}+sum(NN(1:i-1)),:)];
     [ZZ(lnkPP{i}+sum(NN(1:i-1)),:),Dx(lnkPP{i}+sum(NN(1:i-1)),:)-Dx(lnkMM{i+2}+sum(NN(1:i)),:),ZZ(lnkPP{i}+sum(NN(1:i-1)),:)];
     [ZZ(lnkPP{i}+sum(NN(1:i-1)),:),ZZ(lnkPP{i}+sum(NN(1:i-1)),:),Dx(lnkPP{i}+sum(NN(1:i-1)),:)-Dx(lnkMM{i+2}+sum(NN(1:i)),:)]];
end

end

% Neumann
C=[C;[Dx(neu+sum(NN(1:3*Ncavs)),:), ZZ(neu+sum(NN(1:3*Ncavs)),:), ZZ(neu+sum(NN(1:3*Ncavs)),:)];
     [ZZ(neu+sum(NN(1:3*Ncavs)),:), Dx(neu+sum(NN(1:3*Ncavs)),:), ZZ(neu+sum(NN(1:3*Ncavs)),:)]];
 break
 % CUALQUIER PARECIDO CON LA REALIDAD ES PURA COINCIDENCIA...
% C=[II([u1(dir1);v1(dir1);u2(dir2);v2(dir2)],:);% Dirichlet on u,v 
%     II([u1(lnk1);v1(lnk1);p1(lnk1)],:)-II([u2(lnk2);v2(lnk2);p2(lnk2)],:); % Link Continuity
%     D.x(lnk1,:)-D.x(lnk2+NN1,:), Z(lnk1,:), Z(lnk1,:); % Link Tangency on u
%     Z(lnk1,:),D.x(lnk1,:)-D.x(lnk2+NN1,:), Z(lnk1,:); % Link Tangency on v
%     Z(lnk1,:), Z(lnk1,:),D.x(lnk1,:)-D.x(lnk2+NN1,:); % Link Tangency on p
%     D.x(neu+NN1,:), Z(neu+NN1,:), Z(neu+NN1,:); % Neumann on u2 at outflow
%     Z(neu+NN1,:), D.x(neu+NN1,:), Z(neu+NN1,:)];  % Neumann on v2 at outflow

%{
We chose an initial guess that statisfies the boundary conditions, this is good for Newton, and it makes it also very easy to impose the nonhomogeneous boundary conditions in the lop.
%}

% initial guess
U1 = zeros(NN1,1);
U1(l1.left) = Uinlet;%1.5.*Uinlet.*(1-((2.*(y_inlet-min(y_inlet)))./DiamInlet - 1).^2);
V1=zeros(NN1,1);
P1=ones(NN1,1);

U2 = zeros(NN2,1);
V2=zeros(NN2,1);
P2=ones(NN2,1);

U = [U1;U2];
V = [V1;V2];
P = [P1;P2];
sol0=[U1(:);U2(:);V1(:);V2(:);P1(:);P2(:)];

% Newton iterations
disp('Newton loop')
sol=sol0;
solNm1 = sol0;
solM=sol;
quit=0;count=0;time = 0;
for t = 0:dt:tmax
while ~quit     
 
    % the present solution and its derivatives
    UM = solNm1([u1;u2]); % Solution du pas précédant
    VM = solNm1([v1;v2]);
    PM = solNm1([p1;p2]);
    U=alpha.*sol([u1;u2]) + (1-alpha).*solM([u1;u2]);
    V=alpha.*sol([v1;v2]) + (1-alpha).*solM([v1;v2]);
    P=alpha.*sol([p1;p2]) + (1-alpha).*solM([p1;p2]);
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
    f(loc)=C*(sol-sol0);
    A(loc,:)=C;
   
    % convergence test  
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if (count>maxIter) || (res>1e5); disp('no convergence');break; end
    if res<resSTOP; quit=1; disp('converged'); end
    
    % Newton step
    solM = sol;
    sol=solM-A\f;%gmres(A,f,1000,1e-7);%A\f;   
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
solNm1 = sol;
quit = 0;
time =time+dt;
fprintf('t=%g\n',time);
    % plotting
    Module = sqrt((U.^2+V.^2));
    surf([X1,X2],[Y1,Y2],reshape(Module,Ny,2*Nx)); view(2); shading interp; hold on
    caxis([0 2*Uinlet]); % Color scale
    xlabel('x'); ylabel('y'); title(sprintf('Poiseuille t=%0.1fs ; Uinlet = %0.3fm/s; Reynolds %d',time,Uinlet,Re)); grid off;hold off
    colorbar;
    axis equal
    drawnow
    set(gcf,'paperpositionmode','auto')
    print('-dpng','-r80',sprintf('Poiseuille_t%0.1fs.png',time));
end

% Creates animated gif
FirstFrame = true;
gifname = sprintf('Poiseuille_Re%d.gif',Re);
for t = 3*dt:dt:tmax
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

%}
