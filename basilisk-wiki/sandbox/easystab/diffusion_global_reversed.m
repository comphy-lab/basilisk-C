%{
# Diffusion equation without time marching (but reversed)

Here is like in [advection_global_reversed.m](), except that we do this for the diffusion equation (which should be much more ill-behaved, and more sensitive to errors).

%}

clear all; clf

% parameters
Nx=61; % number of gridpoints
Nt=60; % gridpoints in time
Lx=10; % domain length in x
Lt=4; % duration in time
xpos=5; % position of the final condition
mu=0.1; 

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
A=D.t-mu*D.xx;
b=zeros(NN,1);

%{
Here this is the "final guess".
%}

% final guess
f0=exp(-((X-xpos)/1).^2);

%{
# Boundary conditions

The "left", "right" and "end" boundaries should be like in the initial guess.
%}

% boundary conditions
loc=[l.end; l.right; l.left];
C=I(loc,:);
 
b(loc)=C*f0(:); 
A(loc,:)=C;

%{
# Solve system
%}

f=A\b;
f=reshape(f,Nt,Nx);

% show evolution of f
subplot(1,2,1);
surf(X,T,f); view(2)
xlabel('x'); ylabel('t');
title('time evolution of the string');

%{
# March in time of the initial condition

Here this is some kind of validation, we use a different method to go back from the initial condition to the final one, and compare to the global solution.

%}

% march in time the initial condition
I=eye(Nx);
E=I;
A=mu*d.xx; 

% boundary conditions
loc=[1;Nx];
E(loc,:)=0;
A(loc,:)=I(loc,:);

% march in time matrix 
dt=t(2)-t(1);
Mm=(E-A*dt/2);
M=Mm\(E+A*dt/2);

% initial condition
q=f(1,:)';

% marching loop
for ind=1:(Nt-1); q=M*q; end
    
% plotting
subplot(1,2,2);
plot(x,q,'b',X(l.end),f(l.end),'r.-');    
axis([0,Lx,0,1.2])    
legend('march in time','global'); title('diffusion reversed')
xlabel('x'); ylabel('final time');

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','diffusion_global_reversed.png');


%{

The subplot on the left shows the evolution in time of the solution. You see that this is ill-behaved because the solution becomes very large and oscillatory.

![Validation](diffusion_global_reversed.png)

# Exercises/contributions

* Please add exercices and contributions

%}
