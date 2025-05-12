%{
# MARCH IN TIME OF A 1D WAVE EQUATION

Based on the code [vibrating_string.m](/sandbox/easystab/vibrating_string.m) we do the march in time of a 1D wave equation. 

The equation we are studying is :
$$ U_{tt} + (1-r)U = U_{xx}
$$

We choose a solution in the form of :

$$ U(x,t) = Û e^{i \alpha x + \beta t} +cc
$$

We inject this solution in our equation and we can simplify it to get a second equation, which allows us to study the stability of the initial
equation: 

$$ \beta ^2 = r-1 - \alpha ^2
$$

This gives us three cases : 
$$ \beta ^2 < 0 if r < \alpha ^2 + 1
$$
$$ \beta ^2 > 0 if r > \alpha ^2 + 1
$$
$$ \beta ^2 = 0 if r = \alpha ^2 + 1
$$

In the third case, we get a neutral curve, as we can see below on the
graphic r=f(alpha).

So, the final system we must resolve is :

$$ Eq_{t}= Aq 
$$

$$
\begin{pmatrix}
I & Z \\
Z & I \\
\end{pmatrix}
\left(\begin{array}{l}
V_{t} \\
U_{t} \\
\end{array}
\right)
=
\begin{pmatrix}
Z & -(1-r)*I+D_{xx} \\
I & Z \\
\end{pmatrix}
\left(\begin{array}{l}
V \\
U \\
\end{array}\right)
$$

where $$ V = U_{t} $$
%}


%{
## Resolution
%}

clear all; clc; clf; %close all

% parameters
N=100;      % number of gridpoints
L=20;       % domain length
dt=0.1;     % time step
tmax=20;    % final time


%{
## Control parameter

You can change this parameter to see the stable character of the equation. You can refer to the stability curve at the end of this topic to choose the correct parameter.
%}
r=1.1;        % control parameter

%{
Here we use the differentaition matrices as seen in [diffmat.m](/sandbox/easystab/diffmat.m) to write the equation in he matrix form as shown earlier.
%}

% differentiation matrices
% scale=-2/L;
% [x,DM] = chebdif(N,2); 
% dx=DM(:,:,1)*scale; 
% dxx=DM(:,:,2)*scale^2;
% x=(x-1)/scale; 
Z=zeros(N,N); II=eye(2*N); I=eye(N);

[d.x,d.xx,d.wx,x]=dif1D('cheb',0,L,N);
dx=d.x;
dxx=d.xx;

% system matrices
E=II;
A=[Z, (r-1)*I+dxx; I, Z];

%{
We must choose boundary condtions that are in adequacy with our intial solution. Since we will be choosing a sinus form for the intial solution, we impose all the boundariy condtions to 0.
%}

% boundary conditions
loc=[1,N];
E(N+loc,:)=0;
A(N+loc,:)=II(N+loc,:);


% compute eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

mnum=2;
figure(1);
plot(x,U(1:N,mnum),'b',x,U(N+1:2*N,mnum),'r');
title('Position');
xlabel('x'); ylabel('Position / Vitesse');
legend('Vitesse','Position');


%{
## March in time

Here we use the same method as [vibrating_string_spring_dissipation.m]().
%}

% march in time matrix 
M=(E-A*dt/2)\(E+A*dt/2);

%{
## Initial condition

We choose an initial condition but you can generate a random function like in [diffusion_eigenmodes.m](/sandbox/easystab/diffusion_eigenmodes.m).
%}

% initial condition
q=[0*sin(2*pi*x/L); sin(pi*x/L)]; 
ql=q;

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);


%{
## Plotting
%}

figure(2)
filename = 'vibrating_string_zones_de_stabilites_r0.gif';

for ind=1:Nt    

     q=M*q; 
     
    % plotting
    plot(x,q(N+1:N*2));
    axis([0,L,-10,10]); grid on
    xlabel('x'); ylabel('U');
    title('Position')
    drawnow 
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if ind == 1;
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
    
end

%graphic r=f( \alpha )
figure(3)
alpha=[0:0.3:N];
r_vecteur=alpha.^2 + 1;
plot(alpha(1:15),r_vecteur(1:15),'b-*');
title('Zones de stabilites')
xlabel('alpha'); ylabel('r parametre de controle');
text(0.5, 16, '\color{green} Zone instable','FontSize',18); 
text(3,3,'\color{red} Zone stable','FontSize',18); 
text(0.5,6,'\color{blue} Courbe neutre \rightarrow','FontSize',18);
print('-dpng','-r80','vibrating_string_zones_de_stabilites');
set(gcf,'paperpositionmode','auto');

%{
## Figures 
We present you here the stability curve and the displacement of our system generated for r = 0, then r = 2. Let us note that for r = 0, the system is stable and for r > 1, it diverges.

![Instability curve](vibrating_string_zones_de_stabilites.png)

You can see here that above the neutral curve, in other words, when r > 1 the system is unstable and stable when r < 1. When r = 1, the system is both stable and unstable.

For r = 0 : the solution converges 

![r = 0 (stable)](vibrating_string_zones_de_stabilites_r0.png)

The system is stable. If you chose a higher tmax (we chose 20), you will see that the oscillations diminish. The system comes back to an equilibrium.

For r = 2 : the solution diverges

![r = 2 (unstable)](vibrating_string_zones_de_stabilites_r2.png)

Confirming what we see in the neutral curve graph, the system diverges when r > 1.

%}