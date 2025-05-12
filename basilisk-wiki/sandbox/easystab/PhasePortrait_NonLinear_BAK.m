%{ 
# Phase portrait of a NON-LINEAR dynamical system of dimension 2.

This program draws a phase portrait for a dynamical system written under
the form 

$$
\frac{d}{dt} X(t)  = F(X(t)). 
$$

With 
$$
X = \left[ \begin{array}{c} x_1 \\ x_2 \end{array} \right] 
$$

and $F$ is a function defined at the end of the program.
 
The drawing of the phase portrait is done using the same representation
as for the linear case treated by the program [PhasePortrait_Linear.m]().

## Example

The following figure was generated for a dynamical system corresponding to a simple pendulum. Use the program interactively to consider other cases, change parameters and plot more trajectories !

![**Figure 1 :** Phase portrait of a 2x2 nonsystem. This example corresponds to the simple pendulum. 
Two trajectrories, both ending on the stable point, are drawn. Use the program interactively to draw more trajectories 
and/or study another model systems !
](PhasePortrait_NonLinear/PhasePortraitNonLinear.svg)

## **How to use this program :**

- Select your case by setting the variable 'typeproblem' at line 38.
- Change the parameters in corresponding definition section, to be foundin subfunction F_PortraitdePhase defined at the end of this file (in the corresponding block of the "switch/case" selection) 
- Run the program with Octave or Matlab
- Click on the figure to plot trajectories.

The program may be customized by adding cases with new values of 'typeproblem', just add the definition in a new switch/case statement at the end !

# Main program:

%}
        
function [] = PhasePortrait_NonLinear()	

% parametres :
global typeproblem xmin xmax ymin ymax Tmax Tmaxbak;
 
typeproblem = 'Exam2022';  

% Select the problem to consider among the following cases :
% 'Pendulum', 'RotatingPendulum', 'InvertedPendulum','SaddleNode', 'Transcritical', 
% 'Pitchfork', 'Brusselator', 'VanDerPol', 'LotkaVolterra', 'BuffaloWolf', 'Trefethen', 
% 'Exam2021'; 
% 'Custom' [ , ... ]
%typeproblem = 'LotkaVolterra';

%plotEnergy=0; % use 1 to plot the energy in figure 2  
 
% Dimensions of the figure. We give default values here. These can also be adjusted in the definition of the function. 
    close all;
    xmin = -2; xmax = 2;
    ymin = -2; ymax = 2;
    Z = F_PortraitdePhase([0.;0.]); % (trick) call to the function to get the values defined there if not using default..
    eps = 1e-6;
% Pour le trace de trajectoires en cliquant sur la figure : nombre de trajectoires, Tmax
    Ntraj = 10;
    Tmax = 50; % max time for integratation of trajectory
    Tmaxbak = 5; % max time for integratation of trajectory (backwrards)
    Tmaxplot = 5; % max time for plot in figure 2.
    ODEoptions = odeset('Refine',4,'RelTol',1e-7);
% Fin des parametres


%{ 
## Illustration of the flow with vectors

This is done again with the *quiver* command.

%}
    xP = linspace(xmin,xmax,21);
    yP = linspace(ymin,ymax,21);
	[xG, yG] = meshgrid(xP,yP);
    	for i=1:length(xP)
        for j= 1:length(yP)
            F = F_PortraitdePhase([xP(i),yP(j)]);
            ux(j,i) = F(1);
            uy(j,i) = F(2);
        end
    end
    
   figure(1); 

   hold on
   quiver(xG, yG, ux, uy, 'Color', 'b');
    title({['Phase portrait for ',typeproblem ],'Left-click to draw a trajectory ; right-click to stop'},'FontSize',14);
    xlabel('x_1');ylabel('x_2');
    axis([xmin xmax ymin ymax]);
    hold on;
    plot([xmin xmax],[0 0],':k',[0,0],[ymin,ymax],':k');

%    if(plotEnergy==1)
%    figure(1);subplot(4,2,7:8);hold on;
    %title('Amplitude vs. time');
%    xlabel('t');ylabel('(x_1^2+x_2^2)^{1/2}');
%    end   
    
pause(1);

%{
## Drawing of a few trajectories selected by clicking on the figure :
%}  
 INTERACTIVE = isempty(getenv('OCTAVE_AUTORUN'));
   % to detect if the code is running on the server or interactively
   if INTERACTIVE
     % this part is for using the program in interactive mode, 
     % several trajectories can be drawn by clicking on the figure
   figure(1); 
   [xp,yp,button] = ginput(1)
   while(button==1);
    disp('click on figure to launch a trajectory...')
    xinit = [xp ; yp];
    dFdx1 = (F_PortraitdePhase([xp+eps,yp])-F_PortraitdePhase([xp-eps,yp]))/(2*eps);
    dFdx2 = (F_PortraitdePhase([xp,yp+eps])-F_PortraitdePhase([xp,yp-eps]))/(2*eps);
    disp([' Drawing trajectory starting from [x1,x2] = [ ',num2str(xp), ' , ' num2str(yp), ' ] ']);
    disp([' Gradient matrix at starting point :   [ [ ' ,num2str(dFdx1(1)), ' , ',num2str(dFdx2(1)), ' ] ;']);
    disp(['                                     [ ' ,num2str(dFdx1(2)), ' , ',num2str(dFdx2(2)), ' ] ]']);
 
    
   [t,xtraj] = ode45(@(t,X)F_PortraitdePhase(X),[0, Tmax],xinit,ODEoptions);
   [tbak,xtrajback] = ode45(@(t,X)F_PortraitdePhase(X),[0,-Tmaxbak],xinit,ODEoptions);
    plot(xtraj(:,1),xtraj(:,2),'r',xtraj(1,1),xtraj(1,2),'ro');
    hold on;
    plot(xtrajback(:,1),xtrajback(:,2),'m:');
    axis([xmin xmax ymin ymax]);
    xp=xtraj(end,1); yp=xtraj(end,2);
    dFdx1 = (F_PortraitdePhase([xp+eps,yp])-F_PortraitdePhase([xp-eps,yp]))/(2*eps);
    dFdx2 = (F_PortraitdePhase([xp,yp+eps])-F_PortraitdePhase([xp,yp-eps]))/(2*eps);
    disp([' This trajectory ends at point [x1,x2] = [ ',num2str(xp), ' , ' num2str(yp), ' ] ']);
    disp([' Gradient matrix at this point :         [ [ ' ,num2str(dFdx1(1)), ' , ',num2str(dFdx2(1)), ' ] ;']);
    disp(['                                         [ ' ,num2str(dFdx1(2)), ' , ',num2str(dFdx2(2)), ' ] ]']);
    disp(' ');
 %    if(plotEnergy==1)
       %figure(2);
 %      subplot(4,2,7:8);
 %      plot(t,sqrt(xtraj(:,1).^2+xtraj(:,2).^2));
       %figure(1);
%     end
    [xp,yp,button] = ginput(1);
   end   
    figure(1);
    xlabel('x_1');  ylabel('x_2');title('Phase portrait of a nonlinear 2x2 system');
   
  else
    % this part is for using the program in non-interactive mode, 
    % it is execcuted when running automatically on the server.
    % a single trajectory will be drawn and a figure will be generated for the website
    xinit = [0;2];
    [t,xtraj] = ode45(@(t,X)F_PortraitdePhase(X),[0, Tmax],xinit,ODEoptions);
    [tbak,xtrajback] = ode45(@(t,X)F_PortraitdePhase(X),[0,-Tmaxbak],xinit,ODEoptions);
    plot(xtraj(:,1),xtraj(:,2),'r',xtrajback(:,1),xtrajback(:,2),'m',xtraj(1,1),xtraj(1,2),'ro'); 
    xinit = [0;.2];
    [t,xtraj] = ode45(@(t,X)F_PortraitdePhase(X),[0, Tmax],xinit,ODEoptions);
    [tbak,xtrajback] = ode45(@(t,X)F_PortraitdePhase(X),[0,-Tmaxbak],xinit,ODEoptions);
    plot(xtraj(:,1),xtraj(:,2),'r',xtrajback(:,1),xtrajback(:,2),'m',xtraj(1,1),xtraj(1,2),'ro'); 
  end
  set(gcf,'paperpositionmode','auto');
  print('-dpng','-r80','PhasePortraitNonLinear.png');
  
end



%{
# Definition of the function F :
%} 

function F = F_PortraitdePhase(X)

% Fonction pour tracer le portrait de phase d'un systeme 2-2 de la forme 
% d X /d T = F(X). 
%
global typeproblem xmin xmax ymin ymax;
x1 = X(1); x2 = X(2);



switch (typeproblem)
    
%{ 

### Case "Pendulum" :
        
Equation for the motion of a damped pendulum.        

$$
 m L^2 \ddot \theta = - \mu_f \dot \theta - m g \sin \theta
$$

Considering this problem as a dynamical system of order 2 for 
$X = [x_1 ; x_2] = [\theta, \dot \theta]$ leads to : 

$$\frac{ d}{dt} [x_1 ; x_2]  = [ x_2 ; - \mu x_2 - \omega_0^2 \sin(x_1) ]$$

with $\omega_0^2 = \sqrt{\frac{g}{L}}$ ; $\mu = \frac{\mu_f}{M L^2}$.   
    
**Exercice :** 
Observe the structure of the phase portrait for various
values of the damping parameter mu. What do we observe for mu=0 ?
           
%}
    case('Pendulum')
        
        mu = .1;
        omega0 = 1; % omega_0 = sqrt(g/L)
        F = [x2;-omega0^2*sin(x1)-mu*x2];
        % adjust the figure axes 
        xmin = -3*pi ; xmax = 3*pi; ymin = -2*pi; ymax = 2*pi;
        
%{ 

### Case "RotatingPendulum" :
        
Equation for the motion of a pendulum with imposed precession 

$$
 m L^2 \ddot \theta = - \mu_f \dot \theta - m g L \sin \theta + m L^2 \Omega^2 \sin \theta \cos \theta 
$$

Considering this problem as a dynamical system of order 2 for 
$X = [x_1 ; x_2] = [\theta, \dot \theta]$ leads to : 

$$\frac{ d}{dt} [x_1 ; x_2]  = [ x_2 ; - \mu x_2 - \omega_0^2 \sin x_1 (1 - R \cos x_1 )  ]$$

with $\omega_0^2 = \sqrt{\frac{g}{L}}$ ; $\mu = \frac{\mu_f}{M L^2}$ ; $R = \frac{\Omega^2}{\omega^2}$.   
        
**Exercice :** 
Observe the structure of the phase portrait for various values of the rotation parameter R. 
Draw the corresponding bifurcation diagram.
           
%}
    case('RotatingPendulum')
        
        omega0 = 1; mu = .5;
        R = 1.1;
        F = [x2;-omega0^2*sin(x1)*(1-R*cos(x1))-mu*x2];
        % adjust the figure axes 
        xmin = -pi/2 ; xmax = pi/2; ymin = -.5; ymax =.5;
       
        %{ 

### Case "InvertedPendulum" :
        
Equation for the motion of an inverted pendulum with a spring (offset angle alpha) 

$$
 m L^2 \ddot \theta = - \mu_f \dot \theta + m g L \sin \theta - K (\theta-\alpha)  
$$

Considering this problem as a dynamical system of order 2 for 
$X = [x_1 ; x_2] = [\theta, \dot \theta]$ leads to : 

$$\frac{ d}{dt} [x_1 ; x_2]  = [ x_2 ; - \mu x_2 + \omega_0^2 \sin x_1 - k ( x_1- \alpha)    ]$$

with $\omega_0^2 = \sqrt{\frac{g}{L}}$ ; $\mu = \frac{\mu_f}{m L^2}$ ; $k = \frac{K}{m L^2}$.   
        
**Exercice :** 
Observe the structure of the phase portrait for various values of the offset angle alpha. 
Draw the corresponding bifurcation diagram.
           
%}
    case('InvertedPendulum')
        
        omega0 = 1; mu = .5;k = .9;
        alpha = 0.04;
        
        F = [x2;+omega0^2*sin(x1)-k*(x1-alpha)-mu*x2];
        % adjust the figure axes 
        xmin = -pi/2 ; xmax = pi/2; ymin = -.5; ymax =.5;
       
  

%{ 
### Case 'VanDerPol' :

The VanDerPol oscillator is a well-known model of self-sustained
oscillator. It is defined as follows :

$$ m \ddot x +\omega_0^2 x   = (r - \delta x^2) \dot x$$
        

**Exercice :** 
Observe the structure of the phase portrait for various choices of the parameters. 
Reconstruct the corresponding bifurcation diagram in terms of the bifurcation parameter r.     
        
%}
    case('VanDerPol');    
        omega0 = 1;
        r = 0.1;
        delta = 1;
        F = [x2 ; -omega0^2*x1+(r-delta*x1^2)*x2];
        % adjust the figure axes 
        xmin = -3 ; xmax = 3; ymin = -3; ymax = 3;       


        
%{ 
### Case 'Brusselator' :

The 'Brusselator' (more precisely the homogeneous brusselator) is a simple model 
for a chemical reaction displaying unsteadiness. It is defined as follows :
        
$$ 
\dot x_1 = -(\beta + 1) x_1 + x_1^2 x_2 + \alpha  ;
$$
$$ 
\dot x_2 = \beta x_1 - x_1^2 x_2 ;
$$
        
          
      
**Exercice :** 
Observe the structure of the phase portrait for various choices of the parameters. 
Construct the bifurcation diagram as function of the bifurcation parameter beta.
        
%}
   case('Brusselator')
       
       alpha = 1; 
       beta = 3; 
       F = [-(beta+1)*x1+x1^2*x2+alpha ; beta*x1-x1^2*x2] ;
        
        
        % adjust the figure axes 
        xmin = 0 ; xmax = 6; ymin = 0; ymax = 6;   


%{ 
### Case 'LotkaVolterra' :

The Lotka-Volterra is a simple system representing the competition of two
animal species (for instance hares and lynxes). It is defined as follows :
        
$$ 
\dot x_1 = x_1( \alpha - \beta x_2) ;
$$
$$ 
\dot x_2 = x_2( -\gamma + \delta x_1) ;
$$
        
          
      
**Exercice :** 
Observe the structure of the phase portrait for various choices of the parameters. 
Is this problem conservative or not ? 
        
%}
   case('LotkaVolterra')
       % we take the same values as wikipedia...
       alpha = 2/3; 
       beta = 4/3; 
       gamma = 1; 
       delta = 1;
       F = [(alpha-beta*x2)*x1 ; (delta*x1-gamma)*x2];
        % adjust the figure axes 
        xmin = 0 ; xmax = 2.5; ymin = 0; ymax = 2.5;   

        %{ 
### Case 'BuffaloWolf' :

This problem is a 
        
$$ 
\dot x_1 = r x_1 - A x_1 (x_2+ E x_2^2) - B x_1^2;
$$ 
$$ 
\dot x_2 = -C x_2 + D x_1 (x_2+ E x_2^2) ) ;
$$
        
          
      
**Exercice :** 
Observe the structure of the phase portrait for various choices of the parameters. 
Reconstruct the corresponding bifurcation diagram as the parameter r is
varied
        
%}
    case('BuffaloWolf')
    r = 2.3;
    A = .3;
    B = 0.1;
    C = 1;
    D = .2;
    E = 0.15;
    F = [(r*x1 -A*x1*(x2+E*x2^2) - B*x1^2) ; (-C*x2 + D*x1*(x2+E*x2^2)) ];

    % adjust the figure axes 
 xmin = -1 ; xmax = 20; ymin = -1; ymax = 10;   

        %{ 
### Case 'Trefethen' :

This is a simple model which is asymptotically stable but may
 lead to sucritical transition with very small initial conditions thanks to
 transient growth.
 
 
        
 
 
$$ 
\dot x_1 = -(1/R) x_1  + x_2 -\sqrt{x_1^2+x_2^2}  x_2;
$$
$$ 
\dot x_2 = -(2/R) x_2 + \sqrt{x_1^2+x_2^2} x_1;
$$
        
          
      
**Exercice :** 
observe the behaviour of the system for very small values of the initial
 condition (try with various values of $R$?)
        
%}
    case('Trefethen')
    R = 10;
    NL = 1;% select 0 or 1
    normX = sqrt(x1^2+x2^2)  ; 
    F1 = -(1/R)*x1 + x2 - NL*normX*x2; 
    F2 = -(2/R)*x2 + NL*normX*x1;
    F = [F1;F2];
 % adjust the figure axes 
 xmin = -.05 ; xmax = .05; ymin = -.05; ymax = .05;   
 
        %{ 
### Case 'SH_5.1' :

This case is a variant of Trefethen's original model,
and is taken from Schmid & Henningson (Exercice 5.1, p. 525)
 
%}
    case('SH_5.1')
    R = 10;
    NL = 0;% select 0 or 1
    normX = sqrt(x1^2+x2^2)  ; 
    F1 = -(1/R)*x1 + NL*normX*x2; 
    F2 = -(2/R)*x2 + x1 - NL*normX*x1;
    F = [F1;F2];
 % adjust the figure axes 
 xmin = -1 ; xmax = 1; ymin = -1; ymax = 1;   
       
 %{ 
### Case 'SaddleNode' :

We investigate a saddle-node bifurcation in dimension 2 by considering the
following model:

$$ m \ddot x + \mu \dot x  = r - x^2 $$
        
If the friction is dominant over the inertia this equation reduces to the classical 
one-dimensional saddle-node bifurcation.
        
**Exercice :** 
Observe the structure of the phase portrait for various
choices of the parameters. Reconstruct the corresponding bifurcation diagram.     
        
%}
   case('SaddleNode')
         mass = .5; friction = 1; gravity = 1; 
         r=  -.1; 
         F = [x2 ; -friction*x2/mass + (r-x1^2)*gravity];
        % adjust the figure axes 
        xmin = -3 ; xmax = 3; ymin = -1; ymax = 1;
 %{ 
### Case 'Pitchfork' :

We investigate a transcritical bifurcation in dimension 2 by considering the
following model:

$$ m \ddot x + \mu \dot x  = r x  - \delta x^3 $$
        
If the friction is dominant over the inertia this equation reduces to the classical 
one-dimensional pitchfork bifurcation.
        
**Exercice :** 
Observe the structure of the phase portrait for various choices of the parameters. 
Reconstruct the corresponding bifurcation diagram.     
        
%}
   case('Pitchfork')
         mass = .5; friction = 2; 
         r= .5; delta = 1;
         F = [x2 ; -friction*x2/mass+(r*x1-delta*x1^3)/mass];
         % adjust the figure axes 
         xmin = -1 ; xmax = 1; ymin = -1; ymax = 1;         

%{ 
### Case 'TransCritical' :

We investigate a transcritical bifurcation in dimension 2 by considering the
following model:

$$ m \ddot x + \mu \dot x  = r x  - x^2 $$
        
If the friction is dominant over the inertia this equation reduces to the classical 
one-dimensional transcritical bifurcation.
        
**Exercice :** 
Observe the structure of the phase portrait for various choices of the parameters. 
Reconstruct the corresponding bifurcation diagram.     
        
%}

   case('TransCritical')
         mass = .5; friction = 1; 
         r= .5; 
         F = [x2 ; -friction*x2/mass+(r*x1-x1^2)/mass];
         % adjust the figure axes 
         xmin = -1 ; xmax = 1; ymin = -.5; ymax = .5;       
         
%{
### Case 'Exo3.1' :
%}

   case('Exo3.1')
         F = [x1*(1-x1) ; x2*(2-4*x1)];
         % adjust the figure axes 
         xmin = -.5 ; xmax = 1.5; ymin = -.5; ymax = 1.5;         
           
         
%{
### Case 'Exam2021' :
%}
   case('Exam2021')
         r = -.9; A = 100; B = 1;
         F = [x1*(r+1-(x1-1)^2)+A*x2 ; -B*x2];
         % adjust the figure axes 
         xmin = -.5 ; xmax = 1.6; ymin = -.01; ymax = .02;         

 case('Exam2022')
         r = -.9; A = 100; B = 1;
         F = [(x1-r)*(r+1-(x1-1)^2)+A*x2 ; -B*x2];
         % adjust the figure axes 
         xmin = -1.5 ; xmax = 1.6; ymin = -.01; ymax = .02; 


   case('Custom')
        % add your own case here !
  
   otherwise
       error(['Error in PhasePortrait_NonLinear : unknown type ', typeproblem]);
 
end%swith

end%function
