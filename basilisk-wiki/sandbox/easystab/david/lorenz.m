%{

# Numerical integration of the Lorenz system


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

This program will draw two trajectories with slightly different initial conditions.

What you can modify when playing with the program is :

- control parameter r (line 47) 

- initial conditions X0 and Xb0 : change definition in lines 55-56 or uncomment lines 61-63

%}

function [] = lorenz()

%{ 
## Definition of the parameters and initialization of the figures
%}


% control parameter
r= 24; 

% Physical parameters
P = 10;
b = 8/3;

%% Initial conditions 
%% X0 is initial condition of "blue" trajectory and Xb is initial condition of "red" trajectory

% First idea : two points close to [1,1,1]
X0 = [1;1;1];
Xb0 = [1+1e-5;1;1] ; 

%% Second idea : two points close to the equilibrium point
%% the system has an equilibrium point [Xf,Xf,Xf^2/b] ; with Xf = sqrt((r-1)*b) ;
%% we will use initial condition close to this point
% Xf = sqrt((r-1)*b); 
% X0 = [Xf+1e-2;Xf;Xf^2/b+1e-3] 
% Xb0 = [Xf;Xf+1e-2;Xf^2/b-1e-3] 

% Numerical parameters

T0 = 0; % initial time
Tstep = 1; % time-step 
Nt = 100; % max-number of integrations
Tpause = .3; %time (in seconds between two plots) 

% Initialisation of the figures

figure(1);
xlabel('x'), ylabel('y'), zlabel('z');
title('Trajectory in phase-space for two initial conditions'); 
grid;
view(20,10); % set angle for 3D view

figure(2);
subplot(2,1,1);
xlabel('t');ylabel('x (-) , y(--) , z (..)');
title('Time-history of the solutions');
subplot(2,1,2);
title('Distance between the two trajectories');
xlabel('t');ylabel('|X(blue)-X(red)|');

disp('Position the figures and press enter to launch');
pause;

%{ 
## Performing the temporal integration 
%}

for N = 1:Nt
    Dt = [T0 : Tstep/100 : T0+Tstep] ;  % bornes d'integration
    [t, X] = ode45(@(t,X)(f_lorenz(X,r,P,b)), Dt, X0); % integration for "blue"
    [tb, Xb] = ode45(@(t,X)(f_lorenz(X,r,P,b)), Dt, Xb0); % integration for "red"
    figure(1); hold on;% 
    plot3(X(:,1), X(:,2), X( :,3),'b-');
    plot3(Xb(:,1), Xb(:,2), Xb(:,3),'r-');
    

    X0 = X(end,:); Xb0 = Xb(end,:); T0 = t(end);% initial condition for next run
    figure(2);
    subplot(2,1,1);hold on;
    plot(t,X(:,1),'b-',t,X(:,2),'b--',t,X(:,3),'b:', tb,Xb(:,1),'r-',tb,Xb(:,2),'r--',tb,Xb(:,3),'r:');
    
    Xbb = interp1(tb,Xb,t);
    diff = sqrt( (Xbb(:,1)-X(:,1)).^2 + (Xbb(:,2)-X(:,2)).^2 + (Xbb(:,3)-X(:,3)).^2); 
    subplot(2,1,2);hold on;
    semilogy(t,diff,'k');
    
    pause(Tpause);
end

%{ 
## Results 
%}

    figure(1);
    set(gcf,'paperpositionmode','auto');
    print('-dpng','-r80','Lorenz1.png');
     figure(2);
    set(gcf,'paperpositionmode','auto');
    print('-dpng','-r80','Lorenz2.png');
    
    end

%{      
![**Figure 1 :** Trajectories of the Lorenz system in phase space with two close initial conditions](/sandbox/easystab/david/Lorenz1.png)
![**Figure 2 :** Time-history of the two trajecrories and relative distance](/sandbox/easystab/david/Lorenz2.png)
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
 
