%{
#Vibrating string 2D

The idea is to use the [Vibrating_string](./easystab/vibrating_string.m) code which is in 1D and transfer it to two dimensions.

There are to possibilities to run the code. Either you want to validate it and you choose limit=2 
which will put boundary conditions to mimic the behaviour of the string in 1D and be sure we have 
the right evolution of the central point, either you chose limit=1 and you can observe the eigenmodes 
for a vibrating membrane and the march in time of a membrane where we apply a speed different from zero in 
one point of the membrane.
%}
clear all; clf

limit=1;  % limit=1 allows you to test the response of the membrane to a punctual shock
          % limit=2 validation of the case by comparison with the 1D behaviour

%%%% parameters and flags
c=1; % wave velocity
dt=0.5; % time step
tmax=30; % final time
Nx=20; % gridpoints in x 
Ny=20; % gridpoints in x  
Lx=2*pi; % domain size in x
Ly=2*pi; % domain size in y

%{
# 1D matrices
Here we just build the 1D differentiation matrice as we use to do.
 We build the dx and dy with different sizes and number of points.
%}

%1D differentiation matrices
scale=-2/Lx;
[x,DM] = chebdif(Nx,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 

scale=-2/Ly;
[y,DM] = chebdif(Ny,2); 
dy=DM(:,:,1)*scale;  
dyy=DM(:,:,2)*scale^2;    
y=(y-1)/scale; 

%{
# 2D matrices
We buid *Dx* and *Dy* out of *dx* and *dy* and identity matrices (of the right size). Remember that 
Octave is case-sensitive, that is `D`can be a different variable from `d`. We had already the 1D 
meshes `x` and `y`, and we now as well build mesh arrays using the function `meshgrid`. This function 
does something similar to `kron`and is very convenient for building arrays from functions of `x`and `y`, 
in a much easier way than one may think at first. Here too, `X`and `Y` are the 2D equivalents 
for `x`and `y`. 
%}

% 2D differentiation matrices
Dx=kron(dx,eye(Ny));
Dy=kron(eye(Nx),dy);
Dxx=kron(dxx,eye(Ny));
Dyy=kron(eye(Nx),dyy);

[X,Y]=meshgrid(x,y);

%{
# Test
To test our differentiation matrices, we will compare the analytical and numerical derivatives for 
trigonometric functions. You see here how convenient is `meshgrid`. Be aware that, when building `f` `fx` 
and `fy`, we use the `.*`multiplication which is the element by element multiplication, as opposed to 
the `*`multiplication which by default is the matrix multiplication.
%}

Z=zeros(Nx*Ny,Nx*Ny);
I=eye(Nx*Ny);
II=eye(2*Nx*Ny);
DDx=blkdiag(Dx,Dx);

%{
Here we remind you the equation we are studying

$$\partial^2_{ii}  \Phi_j  = \frac{1}{c^2} \partial_t^2 \Phi_j $$

%}

E=[I,Z; Z,I];
A=[Z,c^2*(Dxx+Dyy); I, Z];


% locations at the boundary
dom=reshape(1:Nx*Ny,Ny,Nx);
top=dom(1,1:end); top=top(:); 
bot=dom(end,1:end); bot=bot(:); 
left=dom(2:end-1,1); left=left(:); 
right=dom(2:end-1,end); right=right(:);

loc=[right;left;top;bot];
%{
These lines are meant to differentiate the two cases by changing the loc vector to apply the boundary conditions differently.
%}

if limit==1
    % We impose a Dirichlet condition on all borders
    C=II(loc,:);
end
if limit==2 
    % top & bot with dirichlet; right & left with Neuman (that's why it's DDX)
    C=[II([top;bot],:); DDx([right;left],:)]; 
end

A(loc,:)=C;
E(loc,:)=0;

% locations in the state
v=1:Nx*Ny;
f=Nx*Ny+1:2*Nx*Ny;

% march in time matrix 
M=(E-A*dt/2)\(E+A*dt/2);

%{
# Initial condition
We need to describe the state of the system at initial time. We say here that the string is initially 
deformed as a sinus (this satisfies the boundary conditions), and that the velocity is zero. 
From that state we then perform a loop to repeatedly advance the state of a short forward step in time. 
We use *drawnow* to show the evolution of the simulation as a movie when running the code. 
We store for validation the string position at the midle of the domain, to do this without 
worrying about the the grid points are, we interpolate $f$ with the function *interp1*. 
%}

% initial condition

if limit==1
q=[zeros(Nx*Ny,1); zeros(Nx*Ny,1)];   % Speed and location to 0 everywhere
q((Ny*Ny/4)+Nx/2)=3;                  % Dirac in SPEED in on point 
    
end
if limit==2
q=[zeros(Nx*Ny,1); sin(pi*Y(:)/Ly)]; %Position of the sinus
end

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
%e=zeros(Nt,1);

e(1:Nt)=0;

for ind=1:Nt    
   
    q=M*q; % one step forward
   
    e(ind)=q(Nx*Ny+(Nx*Ny/2)+(Nx/2));  % value of the speed at the center stored for the validation
                                     
    % plotting
    fX=reshape(Dx*f(:),Ny,Nx);
    figure(1)
    mesh(X, Y, reshape(q(f), Ny, Nx))       
    xlim([0,Lx])
    ylim([0,Ly])
    zlim([-1,1])
    drawnow
    pause(0.1)
end
%{
This is the result of the validation (limit=1)
%}
%{
![Validation](/sandbox/easystab/stab2014/vibrating_string_2D_validation3.png)
%}
%{
![Visualisation of the membrane for the verification](/sandbox/easystab/stab2014/vibrating_string_2D_validation2.png)
%}


%{
Here is the evolution of the membrane when hit punctualy (limit=2)
%}
%{
![March in time1](/sandbox/easystab/stab2014/vibrating_string_2D_marchintime1.png)
%}
%{
![March in time2](/sandbox/easystab/stab2014/vibrating_string_2D_marchintime2.png)
%}
%{
![March in time3](/sandbox/easystab/stab2014/vibrating_string_2D_marchintime3.png)
%}
%{
![March in time4](/sandbox/easystab/stab2014/vibrating_string_2D_marchintime4.png)
%}

%{
Now we show that the code does what we think it does, and also give 
ourselves the means to tell how precise and robust it is. 
The period of oscillation of this wave should be the time that it takes 
for a wave to travel at speed $c$ along its wavelength. Here we have 
chosen an initial condition whose wave length is twice the domain 
length $2L$ (this is "mode $1/2$"), thus the period of oscillations is 
T=2L/c
thus the central point of the string should evolve like
cos(2\pi t/T)
Here we build a second time vector, to plot the theoretical solution with
large resolution in time to have a smooth plot even when the time step of
the computation is large.
time evolution of central point
%}
if limit==2
Ttheo=2*Ly; 
tt=linspace(0,tmax,Nt);
etheo=cos(2*pi*tt/Ttheo); %theoretical position
valid(1:Nt)=0
valid(:)=e(:)-etheo(:)  %error between the 1D-theory and the 2D-code
figure(2)
plot(tt,valid)
xlabel('time'); ylabel('amplitude de l erreur');
title('erreur theorie1D/code 2D sur la position du point central')
end

% computing eigenmodes
[U,S]=eig(A,E);
s=diag(S);               %read the eigenvalues
[t,o]=sort(imag(s));    %ordonate the eigenvalues
s=s(o); U=U(:,o);        
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];  %delete the modes for which abs(s)>1000
rem=imag(s)<0; s(rem)=[]; U(:,rem)=[];  %delete the modes for which imag(s)<1000
rem=abs(s)<1e-4; s(rem)=[]; U(:,rem)=[];  %delete the modes for which abs(s)<1e-4

%Visualisation of the natural modes
figure(3)

for ind=1:9
subplot(3,3,ind);
surf(X,Y,abs(reshape(U(f,ind),Ny,Nx))); shading interp; view(2); axis tight
title(s(ind))
end
%{
Here you can visualise the nine first eigenmodes
%}

%{
![Visualisation of the Eigenmodes](/sandbox/easystab/stab2014/vibrating_string_2D_eigenmodes.png)
%}

%{
Coded by [Antoine and Aloïs](./alois&antoine)
%}