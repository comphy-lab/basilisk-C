%{
# The Venturi flow

This is a flow in a cylindrical pipe with a localized neck. So we have the Navier-Stokes equations and we solve for a steady state with the proper boundary conditions. The neck is described using a mapping of the mesh just like we did in [diffmat_mapping.m](), and [peristalsis.m](). The new thing here is that it is done, like for [pipe_sym.m]() and [pipe.m](), in cylindrical coordinates.

This flow is interesting, becasue there is a quick accelleration of the fluid through the neck, and then detachment if the Reynolds number is large enough. This is strongly the case with the present parameters.

Dependency:

* [easypack.zip]()

%}
clear all; clf; format compact

%%%% parameters 
Re=1000; % reynolds number
Nx=101; % number of grid nodes in z
Ny=40; %number of grid nodes in r
Lx=30; % length in z of the domain [0,Lz]
pts=5; % number of points in finite difference stencils
amp=0.5; % radius at venturi
xpos=8; % position of the neck

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D('cheb',-1,2,Ny);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% mapping of the mesh
disp('Mapping of the mesh')
etay=1-(1-amp)*exp(-((x-xpos)/3).^2); 
Y=Y.*repmat(etay',Ny,1);  
D=map2D(X,Y,D);

% cylindrical Laplacian
ym1=spd(1./Y);   ym2=spd(1./Y.^2);
D.lap=D.yy+D.xx+ym1*D.y;

%{
# Boundary conditions

Since there are many boundary conditions, we prepare the vector *loc* and the constraint matrix *C* before the loop. Since the boundary conditions are linear, we can build them once for all before the Newton iterations. 

The boundary conditions are Dirichlet everywhere for $u$ (radial velocity) and $w$ (axial velocity), and homogeneous Neumann on $u,w$ out the outflow. For the pressure we impose its value 0 at the second gridpoint from the left of the top boundary, this is teold in *p0loc*. 

Then we need like we did for [peristalsis.m]() to impose some additional constraints on the pressure, these are imposed for the gridpoints in *neuploc*. how to impose these conditions is still a little mysterious to me, I do trial and errors until I get something that work. Please contribute if you have ideas more clear that this and write to me if you want to learn more.

A new thing in this code is that we use sparse matrices, this is very good to reduce the memory usage for our large 2D differentiation matrices, and reduces greatly the computation cost.
%}

%%%% preparing boundary conditions
NN=Nx*Ny;
u=(1:NN)'; w=u+NN; p=w+NN;
II=speye(3*NN);

neuploc=[l.ctl;l.ctr;l.ctr-Ny];  % where to impose the neumann condition on the pressure
p0loc=2*Ny; % where to impose zero pressure
dir=[l.ctl;l.ctr;l.left;l.top;l.bot]; % where to put Dirichley on u and w

loc=[u(dir); w(dir); p(p0loc); ...
    u(l.right); ...
    w(l.right); ...
    p(neuploc)];

C=[II([u(dir);w(dir);p(p0loc)],:); ...     % Dirichlet on u,w,and p
   D.x(l.right,:), Z(l.right,:), Z(l.right,:); ...   % Neuman on u at outflow
   Z(l.right,:),  D.x(l.right,:),Z(l.right,:); ...    % Neumann on w at outflow
   Z(neuploc,:), D.lap(neuploc,:)/Re, -D.x(neuploc,:)]; % neuman constraint on pressure
%{
We chose an initial guess that statisfies the boundary conditions, this is good for Newton, and it makes it also very easy to impose the nonhomogeneous boundary conditions in the lop.
%}
% initial guess
U=zeros(NN,1);
W=(1-(Y(:,1)*ones(1,Nx)).^2); % mean velocity 1 on the pipe of diameter 1,
P=-X/Re; P=P-P(p0loc); % pressure zero at p0loc
sol0=[U(:);W(:);P(:)];

% Newton iterations
disp('Newton loop')
sol=sol0;
quit=0;count=0;
while ~quit     
 
    % the present solution and its derivatives
    U=sol(u); W=sol(w); P=sol(p);
    Ux=D.x*U; Uy=D.y*U;
    Wx=D.x*W; Wy=D.y*W; 
    Px=D.x*P; Py=D.y*P;

%{
This is now the heart of the code: the expression of the nonlinear fonction that should become zero, and just after, the expression of its Jacobian, then the boundary conditions.
%}
    % nonlinear function
    f=[-(U.*Uy+W.*Ux)+(D.lap*U-ym2*U)/Re-Py; ...
       -(U.*Wy+W.*Wx)+D.lap*W/Re-Px; ...
      (ym1+D.y)*U+D.x*W];
    
    % Jacobian 
    A=[-(spd(W)*D.x+spd(U)*D.y+spd(Uy))+(D.lap-ym2)/Re, -spd(Ux), -D.y; ...
         -spd(Wy), -(spd(W)*D.x+spd(Wx)+spd(U)*D.y)+D.lap/Re, -D.x; ...
         ym1+D.y, D.x, Z];
     
    % Boundary conditions 
    f(loc)=C*(sol-sol0);
    A(loc,:)=C;
%{
Plotting
%}
    % plotting
    subplot(2,1,1);
    surf(X,Y,reshape(U-1,Ny,Nx)); view(2); shading interp; hold on
    
    sely=1:Ny; selx=1:6:Nx;
    ww=reshape(W,Ny,Nx); uu=reshape(U,Ny,Nx); 
    quiver(X(sely,selx),Y(sely,selx),ww(sely,selx),uu(sely,selx),'k');
    axis([0,Lx,-1,1]);
    xlabel('z'); ylabel('r'); title('radial velocity U'); grid off;hold off
    
    subplot(2,1,2);
    surf(X,Y,reshape(W,Ny,Nx)); view(2); shading interp; 
    xlabel('x'); ylabel('y'); title('axial velocity W'); grid off
    drawnow
    
    % convergence test  
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if count>50|res>1e5; disp('no convergence');break; end
    if res<1e-5; quit=1; disp('converged'); continue; end
    
    % Newton step
    sol=sol-A\f;   
    count=count+1;
end


set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','venturi.png')



%{

![The velocity field](venturi.png)


# Validation

Here we should compare these results to a Gerris simulation, or to the method by Pierre-Yves Lagree with the boundary layer equations.

# Exercices/Contributions

* Please add as well a mapping on the $z$ cordinate to be able to have better resolution close to the neck. For the moment, this resolution is the limiting factor for how steep the neck can be. you can do that just before the call to the function [mapping2D.m](), and your modifications will be take into account in the mapped differentiation matrices.
* Please simulate only the top half of the domain, with a symmetry boundary conditions along the axis, as is suggested also for [pipe_sym.m](). This is divide by two the ammount of gridpoints. You can also try to use the symmetry method of the differentiation matrices that is used in [pipe_sym.m]().
* Please compare these calculations with results from gerris with the same geometry----> [venturi_gerris.gfs]().
* Comparison between the Navier-Stokes equations and the Boundary-layer equations -------> [venturi_pyl.m]().

%}