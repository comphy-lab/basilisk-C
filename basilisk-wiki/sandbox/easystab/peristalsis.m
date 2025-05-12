%{
# Peristalsis

Peristalsis is the pumping mechanism in the guts: wall motion in the shape of a traveling wave can induce a pumping. If we look at this flow while moving with the frame of the traveling wave, the flow is steady. Thus we can use the Newton iteration just as we did for the meniscus to compute a steady solution of the Navier-Stokes equations. The motion of the wall is accounted with 2D differentiation matrices and a mapping just like we did before. We use Newton iterations to solve for the nonlinear equations. Here the peristalsis is done in a plane channel instead of a pipe (in 2D).

Dependency:

* [fdper.m]() for the finite difference periodic differentiation matrices
* [fddif.m]() for non-periodic
* [ufdwt.m]() for the stencils of finite difference 
* [spd.m]() to build sparse diagonal matrices 
%}


clear all; clf

% parameters
Re=1; % reynolds number
Nx=21; % gridpoints in x; must be odd
Ny=16; % gridpoints in y 
Lx=20; % domain size in x
Ly=2; % domain size in y
pts=5; % number of points in finite difference stencils
amp=0.5; % wall deformation amplitude
c=1; % wall wave speed

%{
# Periodicity in x

The solution that we compute is periodic in the x direction since we only solve for one wavelength but we assume that there are many of them. For this we build a differentiation matrix using finite difference. Since the domain and everything we will have to differentiate in x is periodic, we do not need to use decentered stencils close to the boundaries, we can use the value of the function "on the other side" of the domain. This is done in the function *fdper.m*.

For this code, we have chosen to use finite difference instead of Chebychev. Finite difference are good here because they make sparse differentiation matrices. This way the 2D differentiation matrices will be sparse, so they will be manipulated quite fast. But apart from this we could have used Chebychev in *y* and Fourrier in *x*.
%}

% 1D differentiation matrices
scale=Lx/2;
[x,dx] = fdper(Nx,1,pts); 
[x,dxx] = fdper(Nx,2,pts); 
dx=dx/scale; dxx=dxx/scale^2; x=x*scale; x=x-x(1);
ix=ones(1,Nx)*(x(2)-x(1)); 
  
scale=-Ly/2;
[y,dy] = fddif(Ny,1,pts); 
[y,dyy] = fddif(Ny,2,pts); 
dy=dy/scale; dyy=dyy/scale^2; y=y*scale;
iy=-([diff(y)',0]+[0,diff(y)'])/2; 

% 2D differentiation matrices
Dx=kron(dx,speye(Ny));
Dxx=kron(dxx,speye(Ny));
Dy=kron(speye(Nx),dy);
Dyy=kron(speye(Nx),dyy);
[X,Y]=meshgrid(x,y);

% vectors for coordinate selections 
NN=Nx*Ny;
dom=reshape(1:NN,Ny,Nx);
top=dom(1,:); top=top(:);
bot=dom(end,:); bot=bot(:);
u=dom(:);  v=u+NN;  p=v+NN; 

% domain mapping
eta=1+amp*cos(x*2*pi/Lx); 
etap=-amp*2*pi/Lx*sin(x*2*pi/Lx);
Xc=X; X=repmat(x',Ny,1); 
Yc=Y; Y=repmat(y,1,Nx).*repmat(eta',Ny,1);  

%{
# Periodicity for the mapping

For the mapping of the computational domain into the physical domain we need to compute the derivatives of the mapped coordinates with respect to the computational coordinates. here we have to be a little careful since everything is periodic in x, except x itself, which is periodic but with a linear increase. So when we compute the derivative of the mapping in x we do

>*xx=Dx*(X(:)-Xc(:))+1*

That is, we substract in x what is not periodic, and we add the known value of its derivative: 1. We do this as well for the second derivative.

Here we use *spd.m* a function to make a sparse diagonal matrix, the sparse equivalent of *diag*.
%}

% mapping coefficients   
xx=Dx*(X(:)-Xc(:))+1; xy=Dy*X(:);   
yx=Dx*Y(:); yy=Dy*Y(:);
yyy=Dyy*Y(:); xxx=Dxx*(X(:)-Xc(:)); yxy=Dy*yx;   
yxx=Dxx*Y(:); xyy=Dyy*X(:); xyx=Dx*xy;
jac=xx.*yy-xy.*yx;

% diff matrices
DMx=spd(yy./jac)*Dx-spd(yx./jac)*Dy;
DMy=-spd(xy./jac)*Dx+spd(xx./jac)*Dy;

DMxx=(spd((yy./jac).^2)*Dxx+spd((yx./jac).^2)*Dyy)+ ...
     spd((yy.^2.*yxx-2*yx.*yy.*yxy+yx.^2.*yyy)./jac.^3)*(spd(xy)*Dx-spd(xx)*Dy)+ ...
     spd((yy.^2.*xxx-2*yx.*yy.*xyx+yx.^2.*xyy)./jac.^3)*(spd(yx)*Dy-spd(yy)*Dx);

DMyy=(spd((xy./jac).^2)*Dxx+spd((xx./jac).^2)*Dyy)+ ...
     spd((xy.^2.*yxx-2*xx.*xy.*yxy+xx.^2.*yyy)./jac.^3)*(spd(xy)*Dx-spd(xx)*Dy)+ ...
     spd((xy.^2.*xxx-2*xx.*xy.*xyx+xx.^2.*xyy)./jac.^3)*(spd(yx)*Dy-spd(yy)*Dx);

DMxx2=-2*spd(yx.*yy./jac.^2)*(Dy*Dx); 
DMyy2=-2*spd(xx.*xy./jac.^2)*(Dy*Dx); 

Dx=DMx; Dxx=DMxx+DMxx2; 
Dy=DMy; Dyy=DMyy+DMyy2;
Dlap=Dxx+Dyy;


%{
# Integration operator
To compute the flux induced by the peristalsis, we wll need to integrate U on the domain, so we do not only need to differentiate, we also need to integrate. The 2D integration is done simply by combining the 1D integration vectors. Then we account for the mapping using the jacobian of the mapping transformation. This jacobian was also used in the computation of the mapped differentiation matrices in 2D.

%}

% integration operators
ww=iy(:)*ix(:)'; 
Dww=ww(:)'; 
Djac=jac(:);
Dww=Dww.*Djac'; 

% useful matrices
Z=spalloc(NN,NN,0); 
I=speye(NN); 
II=blkdiag(I,I,I); 
DDx=blkdiag(Dx,Dx,Dx); 

%{
# Initial guess

To initiate the Newton solver we use an initial guess where U is equal to minus the speed c of the wall wave. Indeed we move to the right with the wave, so the steady fluid in the laboratory frame seems like going to the left at velocity c. This U velocity thus satisfy the no-slip boundary conditions at the wall since the walls are fixed in the laboratory frame, it is just the wave that travels, not the particles at the walls. 

For V we build just a simple field that satisfies the boundary conditions. We have called *eta* the vector that describes the shape of the wall, so the velocity in the moving frame is related to the derivative in *x* of *eta* as written below.

For the pressure we use 0 as an initial guess.

This initial guess is built such as to satisfy the boundary conditions. This will make easier later on the impositino of the boundary conditions on our solutions, we will just impose that the boundary conditions of our solution must be the same as those of the initial guess, which we will store in the vector *base*. 
%}

% initial guess
U=-c+0*X;
V=-c*y*etap';
P=0*U;
base=[U(:);V(:);P(:)];
sol=base;

% newton iterations
quit=0; count=0;
while ~quit

    % base flow derivatives
    U=sol(u); Ux=Dx*U;  Uy=Dy*U;
    V=sol(v); Vx=Dx*V;  Vy=Dy*V;
    P=sol(p); Px=Dx*P;  Py=Dy*P;

%{
# Navier-Stokes

Here we have the steady Navier-Stokes equations, with the conservation of momentum for U and for V and the continuity equation. The solution has *f=0*. We as well buid the Jacobian of it.

%}
    % the nonlinear function
    f=[-(U.*Ux+V.*Uy)+Dlap*U/Re-Px; ...
       -(U.*Vx+V.*Vy)+Dlap*V/Re-Py; ...
       Ux+Vy];

    % jacobian 
    A=[-(spd(U)*Dx+spd(V)*Dy+spd(Ux))+Dlap/Re, -spd(Uy), -Dx; ...
       -spd(Vx), -(spd(U)*Dx+spd(V)*Dy+spd(Vy))+Dlap/Re, -Dy; ...
       Dx, Dy, Z];

%{
# Boundary conditions

We impose Dirichlet boundary conditions at all wall positions: the top and the bottom of u and v. We enforce that our solutino must have the same values of the boundary conditions as the initial guess. Here since we have a periodic domain in the *x* direction, there is no boundary conditions on the left and the right of the domain. 

We need as well some additional constraints for the pressure since it comes into the equations only through its x and y derivatives. The first thing is to impose its value somewhere in the domain, we do that at the first gridpoint (the top left point), where we say that the value of the pressure must be 0. 

We yet need to impose some additional constraint on the pressure, to remove the so-called "spurious modes". What is a spurrious mode? It is a non-zero distribution of pressure that has a zero x and y derivative inside the domain. Indeed the continuity equation is imposed at all gridpoints including the boundaries, except at the first gridpoint where we remove the continuity equation to impose the value of the pressure. But the pressure comes into the momentum equations only at the inside gridpoints (everywhere except at the boundaries) since at the boundaries we replace the momentum equations by the boundary conditions. Here we resolve all this by imposing at two boundary points of the continuity equation a projected version of the momentum conservation in V:

> *Dlap*v/Re-Dy*p*=0

%}
    % Boundary conditions
    loc=[u([top; bot]); v([top; bot]); p(1)]; 
    locp=[Ny;2*Ny]; 
    C=[II(loc,:); ... 
    Z(locp,:), Dlap(locp,:)/Re, -Dy(locp,:)]; 
    
    f([loc; p(locp)])=C*(sol-base);
    A([loc; p(locp)],:)=C;

    % plotting
    subplot(1,2,1);
    surf(X,Y,reshape(U,Ny,Nx)-100); shading interp; view(2); hold on; 
    quiver(X,Y,reshape(U,Ny,Nx)+c,reshape(V,Ny,Nx)); hold off;
    xlabel('x'); ylabel('y'); title('velocity field');
    drawnow 
  
    % convergence test
    res=norm(f);
    if res<1e-6; quit=1; continue; 
    elseif  res>1e8|count>20; 
    disp('% convergence failed'); quit=1; continue;
    end
  
    % Newton step
    sol=sol-A\f;
    count=count+1;
end

%{
# Validation

For the validation of this computation we need to compute the total flux induced by the wave. For this we add *c* the wave velocity to the U velocity field. We integrate this all the complete domain, we divide by the domain length in x to have the average velocity, and then we divide by 2 to have the flux on half of the domain, as was done in the reference article ???. To get the analytical formula for comparison, you take the Stokes equation, assume a long wavelength in x, which can give you an exact solution. Once integrated this gives a flux
$$
\phi=3a^2/(2+a^2)
$$ 
where *a* is the amplitude of the wall deformation which can vary from 0 (straight wall, no pumping) to 1 where the channel is completely ocluded.

For the extreme value of *a=1* it is very easy to get a solution, since there cannot be any flux through the contriction, that means that all the fluid is moving to the right at the velocity *c* of the traveling wave. The quantity of fluid carried by a wavelength is the wavelength Lx times the height Ly=2. Thus the maximum flux is just 1. 

We see on the graph two things: the long wave approximation tends nicely to one for large wall deformation amplitude, and the red dot showing our result is nicely on the blue curve. Indeed we have chosen a long wavelength with Lx=20. 
%}

% computing the flux and validation
flux=(Dww*(U+c)/Lx)/2;

ampvec=linspace(0,1,100);
fluxth=3*ampvec.^2./(2+ampvec.^2);
subplot(1,2,2);plot(ampvec,fluxth,'b',amp,flux,'ro');
xlabel('wall amplitude'); ylabel('flux'); title('validation');
legend('theory','computed')


%{
# The results
Here we put the figure

![Validation](/peristalsis.png)


# Exercices/contributions

* Please show that there is good agreement with the theory for many amplitudes of the wall deformation.
* Please investigate how the flux is not equal to the long wave approximation when the wavelength is short, and show how the one tends to the other when the wavelength is increased.
* Please investigate how the flow and the flux is affected by having a larger Reynolds number. 
* Please investigate the convergence of the solution when varying the grid resolution.
* Please also do the grid convergence test using different number of points in the finite difference stencil (here by default it is 5 in x and in y).
* Please replace finite difference in y by Chebychev and periodic finite difference in x by Fourier and see how the convergence is changed.
* Please compute the same thing only on the top half of the domain since everything is up/down symetric and check that it gives the same thing as the full domain computation (one easy way is to have a straight bottom wall with a symetry boundary condition.
* Please try a pumping with only one deformed wall and find a paper in the litterature that computes this for validation.
* Please change the wall wave shape to see what difference it makes. Here we have used a sinusoidal shape. You can try the short-step/long-straight-channels wave that is shown (with a long-wave theory) in ??? (my paper on peristalsis).
* Write a simple code, just to test the integration in a mapped 2D domain.
* Please try to put somewhere else the 3 constraints on pressure.
* Please add a pressure gradient along the channel to see the competition between Poiseuille flow induced by the pressure gradient and the pumping induced by the peristalsis (and find some litterature reference for validation). A simple case is: no peristalsis, only pressure gradient.
%}
