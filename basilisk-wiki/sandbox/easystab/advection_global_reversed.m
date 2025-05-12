%{
# Advection equation without time marching (but reversed)

Here is like in [advection_global.m](), except that we impose the final condition instead of the initial condition.

We do as well the same test for the diffusion equation in [diffusion_global_reversed.m]().

%}

clear all; clf

% parameters
Nx=41; % number of gridpoints
Nt=40; % gridpoints in time
Lx=10; % domain length in x
Lt=4; % duration in time
U=1; % advection velocity
xpos=6; % position of the initial condition

%{
We build the differentiation matrices, just like we did in [venturi.m](). And we rename the $y$ direction by a $t$. This is just a renaming, there is no changing.
%}

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
A=D.t+U*D.x;
b=zeros(NN,1);

%{
Here this is the "final guess" instead of the initial guess
%}
% final guess
f0=exp(-((X-xpos)/0.5).^2);

%{
# Boundary conditions

So we impose the final condition to be a bump at position *xpos*, and we impose the right boundary condition, instead of the left one that we did in [advection_global.m]().
%}

% boundary conditions
loc=[l.end; l.right];
C=I(loc,:);
 
b(loc)=C*f0(:); 
A(loc,:)=C;

%{
# Solve system
%}

f=A\b;
f=reshape(f,Nt,Nx);

%{
# Validation

For the validation, we compare the initial time (the "bottom boundary") with the initial condition just translated of $-UL_t$, that is the motion of moving to the right at velocity $U$ during time $L_t$ (but backward in time...).
%}

% show evolution of f
subplot(2,1,1);
surf(X,T,f); shading interp; view(2)
xlabel('x'); ylabel('t');
title('time evolution of the string');

% compare initial time with theory
subplot(2,1,2);
xx=linspace(0,Lx,100);
ftheo=exp(-((xx-(xpos-U*Lt))/0.5).^2);

plot(X(l.start),f(l.start),'b.',xx,ftheo,'r-');
legend('numerical','theory'); title('solution at initial time')
xlabel('x'); ylabel('f(x,t=0)');

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','advection_global_reversed.png');

%{

![Validation](advection_global_reversed.png)

# Exercises/contributions

* Please add exercices and contributions

%}
