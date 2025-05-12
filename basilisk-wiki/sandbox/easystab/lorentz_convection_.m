%{

# Numerical integration of the Lorenz system and plotting of the
corresponding convection patterns


The Lorenz system is a classical dimension-3 dynamical 
system modelling thermal convection in a rectangular cell.

The system is defined as follows :

$$
\frac{ dx}{dt} = P (y-x)
$$
$$
\frac{ dy}{dt} = r x - y - x z;
$$
$$
\frac{ dz}{dt} = x y - bz.
$$

- r is the control parameter and represents the temperature difference between upper and lower surfaces (Rayleigh number)

- P reprents the Prandtl number (classical value 10)

- b is a geometrical parameter related to the aspect ratio of the cell (classical value 8/3)

This program will compute one trajectory of the system and generate a movie
of the correspoding convection pattern.

What you can modify when playing with the program is :

- control parameter r (line 47) 

- initial conditions X0 : change definition in lines 55

%}


%{ 
## Definition of the parameters and initialization of the figures
%}


% control parameter
r= 15; 

% Physical parameters
P = 10;
b = 8/3;

%% Initial conditions 

X0 = [1;1;1];


% Numerical parameters

T0 = 0; % initial time
Tstep = .25; % time-step 
Nt = 100; % max-number of integrations
Tpause = .5; %time (in seconds between two plots) 

% Initialisation of the figures
close all;
figure(1);
xlabel('x'), ylabel('y'), zlabel('z');
title('Trajectory in phase-space for two initial conditions'); 
grid;
view(20,10); % set angle for 3D view

figure(20);
xlabel('x');ylabel('y');
title('Structure of the convection cells');
global XX YY;
Xarray = linspace(0,4,60);
Yarray = linspace(0,1,15);
[XX,YY] = meshgrid(Xarray,Yarray);


disp('Position the figures and press enter to launch');
pause;

%{ 
## Performing the temporal integration 
%}

for N = 1:Nt
    Dt = [T0 : Tstep/100 : T0+Tstep] ;  % bornes d'integration
    [t, X] = ode45(@(t,X)(f_lorenz(X,r,P,b)), Dt, X0); % integration for "blue"
    figure(1); hold on;% 
    plot3(X(:,1), X(:,2), X( :,3),'b-');    

    X0 = X(end,:); %Xb0 = Xb(end,:); 
    T0 = t(end);% initial condition for next run

    PlotConvection(X(end,1),X(end,2),X(end,3));
    pause(Tpause);
end

%{ 
## Results 
%}

    figure(1);
    set(gcf,'paperpositionmode','auto');
    print('-dpng','-r80','Lorenz1.png');
     figure(20);
    set(gcf,'paperpositionmode','auto');
    print('-dpng','-r80','LorenzConv.png');
    


%{      
![**Figure 1 :** Trajectories of the Lorenz system in phase space with two close initial conditions](/sandbox/easystab/david/Lorenz.png)
![**Figure 20 :** [convection patterns](/sandbox/easystab/david/Lorenz2.png)
 %}

   
%{
## Definition of the function governing the dynamical system 
%}
    
    
function f = f_lorenz(X,r,P,b)
% cette fonction definit le systeme de Lorenz

x = X(1); y = X(2); z = X(3);


f = zeros(3,1); % vecteur colonne
f(1) = P*(y-x);
f(2) = r*x - y - x*z;
f(3) = x*y - b*z;

end
 

function PlotConvection(X,Y,Z)

global XX YY;
ampX = .02;ampY = 0.02;ampZ = .005;

figure(20);
UU = -X*ampX*(sin(XX*pi).*cos(YY*pi));
VV = X*ampX*(cos(XX*pi).*sin(YY*pi));
TT = (1-YY)-ampZ*Z*sin(2*pi*YY)+ampY*Y*(cos(XX*pi).*sin(YY*pi));
contourf(XX,YY,TT);hold on;
quiver(XX,YY,UU,VV,ampX*30,'k');
axis equal;

end