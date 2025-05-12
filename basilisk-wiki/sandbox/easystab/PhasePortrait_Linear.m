%{  

*This document belongs to the [easystab](http://basilisk.fr/sandbox/easystab/README) project, check the main page of the project to understand the general philosophy of the project.*

*The present program was specifically designed as a support for the course ["Introductions to hydrodynamical instabilities"](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) for the "DET" Master's cursus in Toulouse. For the underlying theory please see [Lecture notes for chapters 1-2](http://basilisk.fr/sandbox/easystab/LectureNotes_DynamicalSystems.md)*

*To use this program: click "raw page source" in the left column, copy-paste in the Matlab/Octave editor, and run.*

# Phase portrait of a LINEAR dynamical system of dimension 2

This program investigates mathematically and numerically the behaviour of 
a linear dynamical system with two degrees of freedom, written under the following form :


$$
d/dt X = A X
$$

With 
$$
X = \left[ \begin{array}{c} x_1 \\ x_2 \end{array} \right] 
$$

and

$$
A = \left[ \begin{array}{cc} A_{11} & A_{12} \\ A_{21} & A_{22} 
 \end{array} \right]
$$

## Example

The following figure was generated using sample values for coefficients of matrix A. Use the program interactively to change the values and plot more trajectories !

![Phase portrait of the linear 2x2 system for A = [[ -.5, 1],[0,-.7]] ; sample trajectory ](PhasePortrait_Linear/Figure1.png)

## **How to use this program :**

- Run the program with Octave or Matlab
- In the "Command Window", enter the coefficients of the matrix A as requested,
- Click on the figure to plot trajectories.

# Main program :
%}

function [] = PhasePortrait_Linear()
clear all; close all; 


INTERACTIVE = isempty(getenv('OCTAVE_AUTORUN'));	% to detect if the code is running on the server or interactively

if INTERACTIVE
  % When running the code interactively : select the parameters by hand
  fprintf('\n Draw the phase portrait for a linear dynamical system of dimension 2 defined as : '); 
  fprintf('\n      f_1 = A11 x1 + A12 x2 ')
  fprintf('\n      f_2 = A21 x1 + A22 x2 \n')
	A11 = input(' A11 = ');
	A12 = input(' A12 = ');
  A21 = input(' A21 = ');
	A22 = input(' A22 = ');
else 
  % If running the code automatically on server ; taking default values
  A11 = -.5;
  A12 = 1;
  A21 = -1;
  A22 = -.7;
end

  
%{

## Plot of the "Flow" in phase space

For a 2d system, the phase space corresponds to the $(x_1,x_2)$-plane. The
flow in phase-space is visulasized by a vector field which is plotted using
the Matlab/Octave command **quiver**.

The flow will be plotted in the figure with blue arrows 
 
%}

% Parameters for plots    
% range for figure 1
    xmin = -2; xmax = 2;
    ymin = -2; ymax = 2;
	
   [xG, yG] = meshgrid(linspace(xmin, xmax, 21), linspace(ymin, ymax, 21));
	ux =  A11*xG+A12*yG;
	uy =  A21*xG+A22*yG;
	figure(1);
     title({'Phase portrait of the linear system dX/dt = A X',...
     ['A11 = ',num2str(A11),' ; A12 = ',num2str(A12),' ; A21 = ',num2str(A21),' ; A22 = ',num2str(A22)],...
     ' Left-click to draw a trajectory ; right-click to stop'},'FontSize',14)
    hold on
	quiver(xG, yG, ux, uy, 'Color', 'b');
    axis([xmin xmax ymin ymax]);
    
    
 %{ 

   ## Eigenvalues/Eigenmodes of the matrix A:
    
   The eigenvalues/eigenmodes are computed with the **eig** command. 
    
 %}
    
 A = [[A11, A12];[A21,A22]];
    
 disp('Matrix A :');
 A
 
 disp([' Det(A) = ',num2str(A11*A22-A12*A21)]);
 disp([' Tr(A)  = ',num2str(A11+A22)]);
 disp('');
    
   [em,D] = eig(A);
   lambda = diag(D); % eigenvalues
    l1 = D(1,1); l2 = D(2,2); 
   
   
   
    % warning : matlab orders the eigenvalues by increasing magnitude, not decreasing real part as in the theory !

   e1 = em(:,1);
   e2 =  em(:,2);
   disp(['Eigenvalues : lambda1 = ',num2str(l1),' ; lambda2 = ',num2str(l2)]) ;
   disp(' ');
   disp('Eigenmodes :');
   e1
   e2
    
 %{
 ## Classification of the equilibrium point (0,0) as function of its eigenvalues :
    
A number of cases of subcases have to be considered, depending if the eigenvalues are real or complex, different or identical, and if the real parts are positive and negative, etc...
    
In the figure, red lines will display the unstable eigenspaces, magenta lines the stable eigenspaces,
green lines the neutral eigenspaces, orange lines the algebraically unstable directions.

In the case where eigenvalues are complex conjugate the directions will be
plotted in dash-dot style.
    
%}
    
    
 if (imag(l1)==0)
        %{ 
        ### Two real eigenvalues
        %}
        
   if(l1<0)&&(l2<0)&&(l1~=l2)
        disp('Two real, negative, distinct eigenvalues : this is a stable node');
           plot([-e1(1),e1(1)],[-e1(2),e1(2)],'m--');
           plot([-e2(1),e2(1)],[-e2(2),e2(2)],'m--');
       elseif(l1<0)&&(l2<0)&&(l1==l2)&&(abs(A12)+abs(A21)>0)
           disp('Two real, negative, equal eigenvalues, nondiagonal matrix : this is an improper stable node');
           plot([-e1(1),e1(1)],[-e1(2),e1(2)],'m--');
       elseif(l1<0)&&(l2<0)&&(l1==l2)&&(abs(A12)+abs(A21)==0)
           disp('Two real, negative, equal eigenvalues, diagonal matrix : this is a stable star');
           plot([-1,1],[0,0],'m--');plot([0,0],[-1,1],'m--');
           plot([-sqrt(.5),sqrt(.5)],[-sqrt(.5),sqrt(.5)],'m--');plot([-sqrt(.5),sqrt(.5)],[sqrt(.5),-sqrt(.5)],'m--');
       elseif(l1>0)&&(l2>0)&&(l1~=l2)
           disp('Two real, positive, distinct eigenvalues : this is an unstable node');
           plot([-e1(1),e1(1)],[-e1(2),e1(2)],'r--');
           plot([-e2(1),e2(1)],[-e2(2),e2(2)],'r--');
       elseif(l1>0)&&(l2>0)&&(l1==l2)&&(abs(A12)+abs(A21)>0)
           disp('Two real, positive, equal eigenvalues, nondiagonal matrix : this is an improper unstable node');
           plot([-e1(1),e1(1)],[-e1(2),e1(2)],'r--');
       elseif(l1>0)&&(l2>0)&&(l1==l2)&&(abs(A12)+abs(A21)==0)
           disp('Two real, positive, equal eigenvalues, diagonal matrix : this is an unstable star');
           plot([-1,1],[0,0],'r--');plot([0,0],[-1,1],'r--');
           plot([-sqrt(.5),sqrt(.5)],[-sqrt(.5),sqrt(.5)],'r--');plot([-sqrt(.5),sqrt(.5)],[sqrt(.5),-sqrt(.5)],'r--');
       elseif(l1*l2<0)
           disp('Two real eigenvalues with opposite sign : this is a saddle');
           if(l1<0)
           plot([-e1(1),e1(1)],[-e1(2),e1(2)],'m--');
           plot([-e2(1),e2(1)],[-e2(2),e2(2)],'r--');
           else
           plot([-e1(1),e1(1)],[-e1(2),e1(2)],'r--');
           plot([-e2(1),e2(1)],[-e2(2),e2(2)],'m--');
           end
       elseif(l1==0)&&(l2>0)
           disp('one null eigenvalue and one positive eigenvalue : this is an unstable degenerated node')
           plot([-e1(1),e1(1)],[-e1(2),e1(2)],'g--');
           plot([-e2(1),e2(1)],[-e2(2),e2(2)],'r--');
        elseif(l1==0)&&(l2<0)
           disp('one null eigenvalue and one negative eigenvalue : this is a nonhyperbolic point of codimension 1')
           plot([-e1(1),e1(1)],[-e1(2),e1(2)],'g--');
           plot([-e2(1),e2(1)],[-e2(2),e2(2)],'m--');
        elseif(l1==0)&&(l2==0)&&(abs(A12)+abs(A21)==0)
           disp('two null eigenvalues, diagonal matrix : this is a nonhyperbolic point of codimension 2')
           plot([-1,1],[0,0],'g--');plot([0,0],[-1,1],'g--');
           plot([-sqrt(.5),sqrt(.5)],[-sqrt(.5),sqrt(.5)],'g--');plot([-sqrt(.5),sqrt(.5)],[sqrt(.5),-sqrt(.5)],'g--');
         elseif(l1==0)&&(l2==0)&&(abs(A12)+abs(A21)>0)
           disp('two null eigenvalues, non diagonal matrix : this is an algebraically unstable point')
           plot([-e1(1),e1(1)],[-e1(2),e1(2)],'g--');
           plot([e1(2),-e1(2)],[-e1(1),e1(1)],'o--');
       end
    else
        %{
        ### Two complex eigenvalues
        %}
         er = real(e1);ei = imag(e1);
         
   if(A11==A22)&&(A21==-A12) % this is a circular focus/center
      if(real(l1)>0)
         disp('two complex conjugate eigenvalues with positive real part, all directions equivalent : this is an unstable circular focus')
           plot([-1,1],[0,0],'r-.');plot([0,0],[-1,1],'r-.');
           plot([-sqrt(.5),sqrt(.5)],[-sqrt(.5),sqrt(.5)],'r-.');plot([-sqrt(.5),sqrt(.5)],[sqrt(.5),-sqrt(.5)],'r-.');
            
   elseif(real(l1)<0)
       disp('two complex conjugate eigenvalues with negative real part, all directions equivalent : this is a stable circular focus')
             plot([-1,1],[0,0],'m-.');plot([0,0],[-1,1],'m.');
             plot([-sqrt(.5),sqrt(.5)],[-sqrt(.5),sqrt(.5)],'m-.');plot([-sqrt(.5),sqrt(.5)],[sqrt(.5),-sqrt(.5)],'m-.');
        elseif(real(l1)==0)
             disp('two complex conjugate eigenvalues with null real part :  all directions equivalent : this is a circular centre')
             plot([-1,1],[0,0],'g-.');plot([0,0],[-1,1],'g.-');
             plot([-sqrt(.5),sqrt(.5)],[-sqrt(.5),sqrt(.5)],'g-.');plot([-sqrt(.5),sqrt(.5)],[sqrt(.5),-sqrt(.5)],'g-.');
          end
             
         
   else % elliptical focus/center
             
   phi= -atan2(2*er'*ei,(er'*er-ei'*ei))/2;
   ec = er*cos(phi)-ei*sin(phi);
   es = -er*sin(phi)-ei*cos(phi);
   disp('main directions for an elliptical focus/center:')
         ec
         es
        if(real(l1)>0)
             disp('two complex conjugate eigenvalues with positive real part : this is an unstable elliptical focus')
             plot([-ec(1),ec(1)],[-ec(2),ec(2)],'r-.');
             plot([-es(1),es(1)],[-es(2),es(2)],'r-..');  
        elseif(real(l1)<0)
             disp('two complex conjugate eigenvalues with negative real part : this is an stable elliptical focus')
            plot([-ec(1),ec(1)],[-ec(2),ec(2)],'m-.');
             plot([-es(1),es(1)],[-es(2),es(2)],'m-..');  
        elseif(real(l1)==0)
             disp('two complex conjugate eigenvalues with null real part : this is an elliptical centre')
              plot([-ec(1),ec(1)],[-ec(2),ec(2)],'g-.');
             plot([-es(1),es(1)],[-es(2),es(2)],'g-..');  
        end
    end  
    end
       
    
%{ 
## Drawing of a few trajectories to complete the phase portrait
 
   Click on the figure to draw trajectories : the part of the
   trajectory for t>0 will be plotted in red and the part for t<0 in magenta.
    
%}
          
 Tmax = 20./max([abs(l1),abs(l2)]); % max time for integratation of trajectory
 Tmaxplot = Tmax/2; % max time for plot in figure 2.
 if INTERACTIVE
   % interactive mode : plotting trajectory by clicking on figure
    [xp,yp,button] = ginput(1);
   while(button==1);
    xinit = [xp;yp];
    
   [t,xtraj] = ode45(@(t,x)([A11*x(1)+A12*x(2); A21*x(1)+A22*x(2)]),linspace(0, Tmax,500),xinit);
    [tb,xtrajback] = ode45(@(t,x)(-[A11*x(1)+A12*x(2); A21*x(1)+A22*x(2)]),linspace(0, Tmax,500),xinit);
    figure(1);
    plot(xtraj(:,1),xtraj(:,2),'r',xtrajback(:,1),xtrajback(:,2),'m',xtraj(1,1),xtraj(1,2),'ro');
    %figure(2); plot(t,sqrt(xtraj(:,1).^2+xtraj(:,2).^2));hold on;
    figure(1);
    [xp,yp,button] = ginput(1);
   end
   
 else % non-interactive mode : draw one sample trajectory 
    xinit = [1;1];
    [t,xtraj] = ode45(@(t,x)([A11*x(1)+A12*x(2); A21*x(1)+A22*x(2)]),linspace(0, Tmax,500),xinit);
    [tb,xtrajback] = ode45(@(t,x)(-[A11*x(1)+A12*x(2); A21*x(1)+A22*x(2)]),linspace(0, Tmax,500),xinit);
    figure(1);
    plot(xtraj(:,1),xtraj(:,2),'r',xtrajback(:,1),xtrajback(:,2),'m',xtraj(1,1),xtraj(1,2),'ro');
    hold on;
    xinit = [-1;1];
    [t,xtraj] = ode45(@(t,x)([A11*x(1)+A12*x(2); A21*x(1)+A22*x(2)]),linspace(0, Tmax,500),xinit);
    [tb,xtrajback] = ode45(@(t,x)(-[A11*x(1)+A12*x(2); A21*x(1)+A22*x(2)]),linspace(0, Tmax,500),xinit);
    plot(xtraj(:,1),xtraj(:,2),'r',xtrajback(:,1),xtrajback(:,2),'m',xtraj(1,1),xtraj(1,2),'ro');
    title({'Phase portrait of the linear system dX/dt = A X',...
        ['A11 = ',num2str(A11),' ; A12 = ',num2str(A12),' ; A21 = ',num2str(A21),' ; A22 = ',num2str(A22)],...
     ' Sample trajectories'},'FontSize',13)
 end
   
set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','Figure1.png');
 
  
end  

     
%{ 

Exercices :
    
- Please play with the code and try to observe all particular cases 

    
%}
