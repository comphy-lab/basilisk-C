clear all; clf; format compact
%{ 
# GEOMETRIC TRANSFORMATION

USE IN GEOMETRY CHANGE

THIS PROGRAM SHOWS EXAMPLES OF GEOMETRY TRANSFORMATION OF DOMAIN

%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DIFFERENTIATION MATRIXES%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFERENTIATION 
% We create here our matrixes of differentiatons of our system, on a
% rectangle domain. 
%[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
%[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,5);
%[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

%{
# PLAN TRANSFORMATION
 We have different tools to change the form of our domain :
 If you want to apply a translation of every mesh of a line (on X, or Y)
 use  X=X+repmat(f(y)',Nx,1)'  to modify the global form. (on Y, use
 instead Y=Y+repmat(f(x)',Ny,1)'
 As example, look at sinusoidal pipe, or Cat tail for symmetrical result.
 If you want to increase or decrease the size of one part of the domain,
 use X=X.*+repmat(f(y)',Nx,1)' or Y=Y.*repmat(f(x)',Ny,1)'. For example,
 look at Ventury pipe

 You can of course do transformation with f(x,y). Use X and Y instead of x
 and y (they are matrixes of the coordinates instead of vecteur). As you
 change them, you may create X0=X and Y0=Y to keep the original form. You
 don't need any repmat as you work on objects with the same size.

 To create Cylinder coordinates, we need a minimum radius as it is hard do
 discretise r=0. It's better to use meshes as close a squares, so we need
 to have $$dr=r.d\theta$$. It might be usefull to use only one criteria of
 discretisation both for X and Y at the beginning.
 The trick is to consider X as the information of the radius, and Y as the
 information of the angle. 
 We first do a translation of the domain of Rmin. Then, we do a
 transformation to have a logarithmic dispersion of the lengths dr
 Then, we use $$ \theta= 2\pi Y/Ly $$, r=X, and we do our transformation
 at the same time : X=X.*cos(theta), Y=Y.*sin(theta).

 If you want to place another object than a Cylinder, you can place
 $$R(\theta)$$ instead of only R, so the position of the object in
 cylinder coordinates, as a vector or as a function. You may be interest
 in Legendre Polynoms to so such functions, see this interactive program I
 made here for this kind of problem : [URL].

 As you want to keep a cylinder domain at the exterior boundary, juste
 multiply your function by $$(Rmax-abs(X))/(Rmax-Rmin)$$ this way, the
 higher order polynoms are taken into account at R=rmin, and not at
 R=rmax.

 TAKING INTO ACCOUNT THE MODIFICATION
 We use Map2D to change our differentiations matrixes according to the
 transformation of X and Y

Here is some examples of forms you can modelise :
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
#Cat Tail
%}
% (will be used in a brusselator-like)
Lx=1;
Ly=1;
Nx=10;
Ny=10;
a=0.5;

[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

X=X;
Y=Y-Ly/2;
Y=Y.*(1+repmat(a*x',Ny,1));

D=map2D(X,Y,D);

subplot(3,2,1);
surf(X,Y,X*0); view(2);
xlabel('x'); ylabel('y'); title('Cat-tail');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
# SINUSOIDAL PIPE
%}
Lx=10;
Ly=1;
Nx=30;
Ny=10;
a=0.2;% magnitude of the sinus
N=4; %Number of periods

[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

X=X;
Y=Y-Ly/2+repmat(a*cos(N*pi*x/Lx)',Ny,1);

D=map2D(X,Y,D);

subplot(3,2,2);
surf(X,Y,X*0); view(2);
xlabel('x'); ylabel('y'); title('Sinusoidal pipe');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
#CREATING CYLINDER COORDINATES
%}
Lx = 1;
Ly = 10;
Nx=30; % number of grid nodes in x
Ny=40; %number of grid nodes in y

Rmin=0.3;
Angle=2*pi;

[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

X=-X-Rmin;
R0=X;
Shape=1;
X=R0.*cos(Angle*Y/Ly);
Y=-R0.*sin(Angle*Y/Ly);

D=map2D(X,Y,D);

subplot(3,2,3);
surf(X,Y,X*0); view(2);
xlabel('x'); ylabel('y'); title('Cylinder coordinates');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %{ 
 # AN OBJECT DEFINED BY A FUNCTION
 Please look at this program I made, wich is a good way to create form with Legendre polynoms : [external](http://www.openprocessing.org/sketch/192374)
 %}
 Lx = 1;
Ly = 10;
Nx=30; % number of grid nodes in x
Ny=40; %number of grid nodes in y

Rmin=1;
Rmax=2;
Angle=2*pi;

[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

X=-X-Rmin;
R0=X;
Shape=1+(Rmax-abs(X))/(Rmax-Rmin).*(0.1*cos(Angle*8*Y/Ly)+0.05*cos(Angle*16*Y/Ly));
X=X.*cos(Angle*Y/Ly).*Shape;
Y=-R0.*sin(Angle*Y/Ly).*Shape;

D=map2D(X,Y,D);

subplot(3,2,4);
surf(X,Y,X*0); view(2);
xlabel('x'); ylabel('y'); title('Cylinder around an object');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %{
 # RUGOUS CAVITY
 %}
Lx=10;
Ly=1;
Nx=30;
Ny=10;
a=0.1;% magnitude of the perturbation
N=4; %Number of periods

[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,5);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Ly,Ny,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

X=X;
Y=Y.*(1+repmat(a*cos(N*pi*x/Lx)'.^10,Ny,1));

D=map2D(X,Y,D);

subplot(3,2,5);
surf(X,Y,X*0); view(2);
xlabel('x'); ylabel('y'); title('Rugous cavity');

%{
# Venturi shape
%}

Lx = 30;
Ly=2;
Nx=50;
Ny=20;
Ymin=-Ly/2;
pts=5;
amp=0.8; %(diameter at the thickest part)
xpos=10; %(position of the thickest part)

[d.x,d.xx,d.wx,x]=dif1D('fd',0,Lx,Nx,pts);
[d.y,d.yy,d.wy,y]=dif1D('cheb',Ymin,Ly,Ny);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

etay=1-(1-amp)*exp(-((x-xpos)/3).^2); 
Y=Y.*repmat(etay',Ny,1);  

D=map2D(X,Y,D);

subplot(3,2,6);
surf(X,Y,X*0); view(2);
xlabel('x'); ylabel('y'); title('Venturi cavity');


        set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','geometry.png')

%{ 
#Figure 
![alt text](/geometry.png)
%}