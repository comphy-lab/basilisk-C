%{
# Advection diffusion global

Same as the [diffusion global](diffusion_global.m), here I tried to represent the advection diffusion without the march in time.
%}

clear all; clf

% parameters
Nx=41; % number of gridpoints
Nt=40; % gridpoints in time
Lx=10; % domain length in x
Lt=4; % duration in time
U=1; % advection velocity
mu=0.1;
xpos=2; % position of the initial condition

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,Lx,Nx);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Lt,Nt,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% rename y to time derivative
D.t=D.y;
D.tt=D.yy;
t=y;
T=Y;

l.left=[l.cbl;l.left;l.ctl];
l.right=[l.cbr;l.right;l.ctr];
l.start=l.bot;
l.end=l.top;

% system matrices
A=D.t+U*D.x-mu*D.xx;
b=zeros(NN,1);

% initial guess
f0=exp(-((X-xpos-U*Lt)/0.5).^2);

% boundary conditions
loc=[l.start,l.end];
C=I(loc,:);
 
b(loc)=C*f0(:); 
A(loc,:)=C;

f=A\b;
f=reshape(f,Nt,Nx);

% show evolution of f
subplot(2,1,1);
surf(X,T,f); shading interp; view(3)
xlabel('x'); ylabel('t');
title('time evolution of the string');

% compare final time with theory
subplot(2,1,2);
xx=linspace(0,Lx,100);
ftheo=exp(-((xx-(xpos+U*Lt))/0.5).^2);

plot(X(l.end),f(l.end),'b.',xx,ftheo,'r-');
legend('numerical','theory'); title('solution at final time')
xlabel('x'); ylabel('f(x,t=Lt)');

%{
![](/ad_dif_gl.png)


# Contributions
- Please try to remove the periodical boundaries
%}