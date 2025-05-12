%{
# Advection 2D : $f_t+Uf_x+Vf_y=0$
In this code, we do the marching in time of the advection equation $$f_{t}+Uf_{x}+Vf_{y}=0$$ with the boundary condition at the left by changing the code given [vibrating_string.m](). 

%}
clear all; clf

% parameters
Nx=30; % number of gridpoints in x
Ny=30; % number of gridpoints in y 
Lx=10; % domain length x
Ly=10; % domain length y
U=1;% x velocity
V=0;% y velocity
dt=0.01; % time step
tmax=10; % final time

x0=Lx/4; % x-coordinate at initial time
y0=0;% y-coordinate at initial time
l0=1; % length width

% differentiation matrices
scale=-2/Lx;
[x,DM] = chebdif(Nx,2); 
dx=DM(:,:,1)*scale;        
x=(x-1)/scale; 

scale=-2/Ly;
[y,DM] = chebdif(Ny,1); 
dy=DM(:,:,1)*scale;    
y=(y-1)/scale; 

[X,Y]=meshgrid(x,y);
% 2D differentiation matrices
Dx=kron(dx,eye(Ny));
Dy=kron(eye(Nx),dy);
[X,Y]=meshgrid(x,y);


I=eye(Nx*Ny);


%{
# System matrices

Here we build the matrices that give the discretization of the equations. The equations is
$$
f_{t}+Uf_{x}+Vf_{y}=0
$$ 
We can also write this equation in the following way
$$
f_{t}=-UD_xf-VD_yf
$$ 

We can put an identity matrix on the left of the equation, in this way, we can impose the
boundary condition more easliy.

So the equation becomes :
$$
If_{t}=-UD_xf-VD_yf
$$ 

Now we have a linear system with the form
$$
Ef_{t}=Af
$$ 
%}


% system matrices
E=I;
A=-U*Dx-V*Dy;



%{
# Boundary conditions
For 2D problems, by extracting the boundary cell locations (methode used in
[poisson2D.m](/sandbox/easystab/poisson2D.m)), we can set our boundary
conditions more efficiantly. Know that we have just one derivite on x in our equation, the
only condition we have to set is the left side of the array.
 

%}

% boundary conditions

% locations at the boundary
dom=reshape(1:Nx*Ny,Ny,Nx);
top=dom(1,2:end-1); top=top(:); 
bot=dom(end,2:end-1); bot=bot(:); 
left=dom(1:end,1); left=left(:); 
right=dom(1:end,end); right=right(:); 
II=eye(Nx*Ny); ZZ=zeros(Nx*Ny,Nx*Ny);
loc=[left];


A(loc,:)=II(loc,:);
E(loc,:)=0;

%{
# March in time 
Here we use the same methode with [vibrating_string.m](/sandbox/easystab/vibrating_string.m)
%}

% march in time matrix 
M=(E-A*dt/2)\(E+A*dt/2);

%{
# Initial condition
Here we need to describe the initial state of the system. We say here that the surface is initially deformed as a bell curve, and that the velocity is zero. From that state we then perform a loop to repeatedly advance the state of a short forward step in time. We use drawnow to show the evolution of the simulation as a movie when running the code. We store for validation the string position at the midle of the domain, to do this without worrying about the the grid points are, we interpolate 
f
 with the function interp1.
%}

% initial condition
q=Y*0+exp(-((X-x0)/l0).^2-((Y-y0)/l0).^2); 
q=q(:);

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(length(y),Nt);

for ind=1:Nt    
    q=M*q; % one step forward
    q_reshape=reshape(q,Ny,Nx);
    e(:,ind)=interp2(X,Y,q_reshape,Lx/2,y); % store center point
    
    % plotting
    figure (1)
    mesh(X,Y,q_reshape); 
   
    axis([0,Lx,0,Ly,0,1]);
    drawnow
end
legend('position'); title('Vibrating string')
xlabel('x'); ylabel('y');zlabel('f');

%{
# Validation
Now we show that the code does what we think it does, and also give
ourselves the means to tell how precise and robust it is.In order to verify
the solutions, we draw the position of the point in the middle of L at different times.
%}

% time evolution of central point
figure(2);
[TVEC,YT]=meshgrid(tvec,y);
mesh(TVEC,YT,e); 
axis([0,Lx,0,Ly,0,1]);

hold on

tt=linspace(0,tmax,500);
[TT,YTT]=meshgrid(tt,y);
etheo=YTT*0+exp(-((Lx/2-U*TT-x0)/l0).^2-((YTT-V*TT-y0)/l0).^2); 
surf(TT,YTT,etheo); 
axis([0,Lx,0,Ly,0,1]);
xlabel('temp/s');
zlabel('amplitude');
title('Valitation');
%{

<center>
![validation](/advection_2d.gif)
![validation](/validation.png)
</center>
%}

%{
Form the figure of validation, we find a superposition of the theory results and the numeral results
%}

%{
# ZHAO's contribution
\ See also [advection_1D.m](/sandbox/easystab/advection_2D.m)

\ Or link to my contribution page [zhao.m](/sandbox/easystab/stab2014/zhao.m)




# Notes

De la part de Jérôme

~~~matlab
domaine	        valeur	note
connectivité 	2	2
recyclage	2	2
graphiques	2	1
théories	4	0
Originalité	4	0
note /14	14	5
~~~

Ce que tu as fait est bien mais ça correspond seulement aux deux premières séances de TP. Dommage que tu as laissé tombé.


%}


