%{
# Advecting tracer particles

This code is just an example on how to have a velocity field moving along particles. A trouble with stability is that it is hard to realize when you get an eigenvalue, that this is a real flow with all its complexity and beauty. There are many modelisation steps to back-engineer. This is why I give a few examples of showing the actual velocity field of the eigenvectors, and it helps very much to let it advect particles, because we usually view the flow by the motion of particles instead of velocity vectors.
%}

clear all; clf; format compact

% parameters
n=50;       % gridpoints for the velocity field
np=20000;   % numer of tracer particles
Lx=2*pi;    % Domain size in x
Ly=2*pi;    % Domain size in y
dt=0.1;     % Time step
tmax=30;   % final time

% build the grid for the velocity field
x=linspace(0,Lx,n);
y=linspace(0,Ly,n);
[X,Y]=meshgrid(x,y);

%{
# The particles

We set many particles, located by their $x$ and $y$ coordinate. We chose this randomly. In addition to that, we will give a different color to the particles depending on wether they are on the right or on the left initially to help the eye see what happens and the global motion. For this selection, I use a logical array *selp*, which is true (value 1) if the particle is on the right and false (value 0) if it is on the left. We give the color blue to the "true" particles and the color red to the "false" ones. below in the plotting command.
%}

% Random initial positions of the particles
px=rand(np,1)*Lx;
py=rand(np,1)*Ly;

% to select the particles
selp=px>Lx/2;

% time loop
tvec=0:dt:tmax;
for ind=1:length(tvec);

%{
# The velocity field

We build an artificial velocity field based on the streamfunction
$$
\phi=\cos(x)\cos(y)
$$
so that we have 
$$
u=-\phi_y=\cos(x)\sin(y)
$$ 
and 
$$v=\phi_x=-\sin(x)\cos(y)$$
and to let it be a little more interesting, we have the velocity field change in time instead as
$$
\phi=\cos(x-t)\cos(y)
$$
This velocity field is now defined on the cartesian grid that we have built before the loop using the very usefull function *meshgrid*; we now need to interpolate the field on this grid to know what is the velocity vector at the position of each particle. We do this using *interp2* using the default linear interpolation scheme.
%}

    % the velocity field
    u=cos(X-tvec(ind)).*sin(Y);
    v=-sin(X-tvec(ind)).*cos(Y);

    % interpolate the velocity field
    vx=interp2(X,Y,u,px,py,'linear');
    vy=interp2(X,Y,v,px,py,'linear');

%{
# Periodicity

To make things simple, we have a periodic box and a periodic velocity field and the particles may cross the boundaries of the box, so we need to find a way to insert them back on the other side of the box when they are advected out. This is made simply by using the function *mod* for the output of the very simple explicit advection scheme where the particle is advanced using the actual velocity for the time step *dt*

>px=px+dt*vx;

%}
    % march one step
    px=mod(px+dt*vx,Lx);
    py=mod(py+dt*vy,Ly);

    % plotting
    plot(px(selp),py(selp),'b.',px(~selp),py(~selp),'r.')
    axis equal; axis([0,Lx,0,Ly]);
    xlabel('x'); ylabel('y'); title('advection of particles');
    drawnow
    %print('-dpng','-r100',['frames/particles_' num2str(ind) '.png'])
end


%{
![the last time of the simulation](particles.png)

The code when it runs displays the figure synchronously so you see the animation, but  please click on animation to see it saved for you: [particles.mov](http://www.lmm.jussieu.fr/~hoepffner/enseignement/particles.mov)

# Exercices/Contributions

* Please put more than two colors
* Please use another velocity field
* Please build a velocity field using point vortices of a perfect fluid that advect each other in time and advect the particles using this instationnary field. You can as well have sources and sink, but then it will be nice to introduce new particles in the sources to avoid area without particles arround them.


%}