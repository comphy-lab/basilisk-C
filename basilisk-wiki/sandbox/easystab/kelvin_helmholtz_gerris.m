%{

# Kelvin-Helmholtz : validation with Gerris

This code computes the instability of the Kelvin-Helmholtz instability of a shear layer, just like [kelvin_helmholtz_hermite.m](), but here the idea is to generate the velocity field for a nonlinear simulation with [Gerris](http://gerris.dalembert.upmc.fr/main_page.html). 

The parameter file for Gerris is [kelvin_helmholtz.gfs](). We use a domain built out of two square boxes each of size *Lx*.

**Dependency:**

* [herdif.m]()
* [herroots.m]()
* [poldif.m]()
* [kelvin_helmholtz.gfs]() the parameter file for Gerris
* [v.gfv]() the view settings for gfsview2D
* You need to have Gerris and Gfsview installed on your system

**Files created:**

* *eevo* showing the evolution in time of the perturbation energy in Gerris
* *u.cgd* and *v.cgd* the initial condition for Gerris

%}


clear all; clf;
Lx=10;
alpha=2*pi/Lx;    % the wave number
Re=1000;    % the Reynolds number
N=100;      % the number of grid points
L=2*Lx;     % height of domain in gerris

%{
Here we do a scaling of the *y* axis and associated differentiation matrices such that the domain in this file corresponds to the height of the domain in Gerris.
%}

% Hermite differentiation matrices
scale=2*13.4064873381/L;
[y, DM] = herdif(N,2,1);
D=DM(:,:,1)*scale;
DD=DM(:,:,2)*scale^2;
y=y/scale;
Z=zeros(N,N); I=eye(N); 

% renaming the differentiation matrices
dy=D; dyy=DD;
dx=i*alpha*I; dxx=-alpha^2*I;
Delta=dxx+dyy;

% base flow
U=tanh(y); 
Uy=(1-tanh(y).^2);

% the matrices
S=-diag(U)*dx+Delta/Re;
A=[ ...
    S, -diag(Uy), -dx; ...
    Z, S, -dy; ...
    dx, dy, Z];
E=blkdiag(I,I,Z);

% computing eigenmodes 
disp('computing eigenmodes')
[UU,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); UU=UU(:,o);
rem=abs(s)>1000; s(rem)=[]; UU(:,rem)=[];

%{
# Initial condition for Gerris

From now on, we prepare the 2D field for the initial condition of Gerris. *Nx* is the nimber o grid points in *x*. And here we select the first eigenmode, which is the least stable one. 
%}

% showing the velocity field
Nx=100;
u=1:N; v=u+N; p=v+N; 
q=UU(:,1); 
Lx=2*pi/alpha;  x=linspace(-Lx/2,Lx/2,Nx);

%{
We need to translate vertically the axis for Gerris, because we use two boxes and the $y=0$ is in the center of box one which is the bottom one. This way, the shear layer is located just at the connection between the two boxes of Gerris.
%}
% translate y axis for gerris
y=y+Lx/2;

%{
Here we chose the relative amplitude of the eigenmode and the base flow (the base flow has max u velocity equal to 1). Chosing an amplitude too large means that the evolution will be nonlinear right from the beginning, which we don't want because we would like to compare the evolution of the energy in Gerris and see that as long as the perturbatino is small, it is the same as in this code. If we start on the other hand with an amplitude too small, grid resolution in Gerris and round-off error may disturb the evolution of our eigenmode.

To expand to physical space, we do like in [kelvin_helmholtz_hermite.m#velocity-field]().
%}
% scale mode amplitude 
q=0.01*q/max(abs(q(u))); 

% expand to physical space
qphys=2*real(q*exp(i*alpha*x));

% add the base flow to the perturbations
uu=qphys(u,:)+U*ones(1,Nx);
vv=qphys(v,:);
pp=qphys(p,:);

% show the velocity field
selx=1:5:Nx; sely=1:3:N;
quiver(x(selx),y(sely),uu(sely,selx),vv(sely,selx),'k'); hold on
surf(x,y,pp-10,'facealpha',0.5); shading interp;
axis equal;axis([x(1),x(end),y(1),y(end)]);
xlabel('x'); ylabel('y'); title('Velocity field of the initial condition');

%{
# Saving to disk
There are several ways to save fields for gerris, here we chose the simplest one, the [cgd format (cartesian grid data)](http://gerris.dalembert.upmc.fr/gfsfunction.html#Cartesian_Grid_Data_.28CGD.29_files)
%}

% save to cartesian grid data format for gerris (.cgd)
disp('saving initial condition')
% for u
fid=fopen('u.cgd','w');
fprintf(fid,'%s\n','2 x y');
fprintf(fid,'%u %u\n',Nx,N);
fprintf(fid,'%f ',x);fprintf(fid,'\n');
fprintf(fid,'%f ',y);fprintf(fid,'\n');
fprintf(fid,'%f\n',uu);  
fclose(fid);
% for v
fid=fopen('v.cgd','w');
fprintf(fid,'%s\n','2 x y');
fprintf(fid,'%u %u\n',Nx,N);
fprintf(fid,'%f ',x);fprintf(fid,'\n');
fprintf(fid,'%f ',y);fprintf(fid,'\n');
fprintf(fid,'%f\n',vv);  
fclose(fid);

%{
# Starting Gerris

The parameter file for the Gerris simulation is in [kelvin_helmholtz.gfs]().

We can start Gerris directly from here. The command to be executed in the shell is stored in the string *command*. The character "!" at the start of a command tells that the command should be executed in the shell. This command, stored in a string, is then evaluated using the function *eval*.

In Gerris it is possible to define some of the parameter values directly from the command line, like we do here to set the value of the *boxsize*. To set it automatically to the value of *Lx* in the code, we use *sprintf* which will insert in the string the value of *Lx* with two digits after the comma.

If starting Gerris like this does not work on your system, you can type this in the shell (and change the value *10* if you changed *Lx* in this code):

> gerris2D -m -Dboxsize=10. kelvin_helmholtz.gfs | gfsview2D v.gfv

The file [v.gfv]() parameterizes the view of the running dimulation by gfsview. If you don't have gfsview installed (this is more difficult to install than gerris itself), just remove the line "OutputSimulation" in *kelvin_helmholtz.gfs* and type in the shell

> gerris2D -m -Dboxsize=10. kelvin_helmholtz.gfs 

The Gerris simulation saves on the disc a file *eevo* which memorizes the time in the first column and the value of the energy of the perturbation as a second column.
%}

% start gerris simulation
disp('starting gerris')
command=['!gerris2D -m -Dboxsize=' sprintf('%3.2f',Lx) ' kelvin_helmholtz.gfs | gfsview2D v.gfv'];
disp(command);
eval(command)

%{

![A snapshot from the vorticity of the Gerris simulation](kelvin_helmholtz_snapshot.png)

# Validation

To validate these computations, we compare the evolution of the energy here and in Gerris. The eigenvalue of the mode we use here is stored in *s(1)*. Since the energy is the square of the mode amplitude, and the mode amplitudes grows in time like 
$$
\exp(s_rt)
$$
where $s_r$ is the real part of the eigenvalue $s$, then the energy grows like
$$
\exp(2s_rt)
$$
We scale this energy so that the linear one and the nonlinear one are initialy equal. 
%}

% validation
clf;
a=load('eevo');
t=a(:,1); e=a(:,2);
etheo=exp(2*real(s(1))*t); 
etheo=etheo/etheo(1)*e(1); 

semilogy(t,e,'b.-',t,etheo,'r');
xlabel('time'); ylabel('energy of perturbation')
legend('gerris','easystab');
title('Comparison Gerris/easystab')
grid on

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','kelvin_helmholtz.png')

%{
![Comparison of easystab and Gerris](kelvin_helmholtz_gerris.png)

# Exercices/Contributions

* Please change this code such that we can have several wavelengthes of the instability in the initial condition, to have a look at the phenomenon of vortex-pairing
* Please validate the code for several values of the viscosity and of the wavelength
* Please do these validation with Gerris for other instabilities, like [poiseuille_uvp.m]() or [couette_uvwp.m](), or [pipe_sym.m]() or even [rayleigh_benard.m]() for which you will need to initiate as well the temperature field. 


%}