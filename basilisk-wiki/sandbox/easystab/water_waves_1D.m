%{
# WATER WAVES IN A 1D SWIMMING POOL

This code is based upon [vibrating_string.m](/sandbox/easystab/vibrating_string.m).

We change the boundary conditions so that our system is a model of the water waves in a 1D swimming pool (we use homogenous Neumann conditions instead of Dirichlet conditions).


# March in time of the vibrating string
%}

clear all; clf

% parameters
N=50; % number of gridpoints
L=2*pi; % domain length
c=1; % wave velocity
dt=0.01; % time step
tmax=20; % final time

% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 
II=eye(2*N);

%{ 
# System matrices 
%}
% system matrices
E=[I,Z; Z,I];
A=[Z,c^2*dxx;I,Z]; 

% locations in the state
v=1:N;
f=N+1:2*N;

%{
# Boundary conditions
The original system has one variable and two time derivatives, so we should apply two boundary conditions. We will simply say that the ends of the state are not attached but the first space derivative is zero. This is a homogeneous Neumann condition applied on the first and last point of the state vector for the state position.
%}

% boundary conditions
loc=[f(1),f(N)];
C=[Z(1,:),dx(1,:); Z(N,:), dx(N,:)]; % neumann condition
E(loc,:)=0;
A(loc,:)=C;

%{ 
# March in time
%}
% march in time matrix 
M=(E-A*dt/2)\(E+A*dt/2);

%{
# Initial condition
We need to describe the state of the system at initial time. We say here that the string is initially deformed as a cosinus (this satisfies the boundary conditions), and that the velocity is zero. From that state we then perform a loop to repeatedly advance the state of a short forward step in time. We use *drawnow* to show the evolution of the simulation as a movie when running the code. We store for validation the string position at the midle of the domain, to do this without worrying about the the grid points are, we interpolate $f$ with the function *interp1*. 
%}

% initial condition
q=[zeros(N,1); cos(2*pi*x/L)]; % we want to see an entire wavelength hence we choose a period of 2*pi

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

%{
## Plotting

We want to see the evolution of the eigenmodes. So, we plot the eigenmodes calculated earlier in a time loop using "drawnow" function. Furthermore, we save the resulting graph as an animated image. To do that, we used the code found in [http://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab](http://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab). If you do not want to do this, put on comment the lines following "drawnow" and "filename = 'brusselator_eigenmode_lambda4.gif'". 
%}

figure(1)
filename = 'water_wave_1D_model_position_veolocity.gif';
for ind=1:Nt    
    q=M*q; % one step forward
    e(ind)=interp1(x,q(f),L/2); % store center point
    
    % plotting
    plot(x,q(f),'b',x,q(v),'r--'); grid on; axis([0,L,-2,2]); 
    legend('position','velocity'); title('Water Wave');
    xlabel('x'); ylabel('f');
    drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if ind == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
end

%{
# Validation
%}

% time evolution of central point
figure(2)
Ttheo=L/c; % we change the period so that the theory superimposes on our numerical solution 
tt=linspace(0,tmax,500);
etheo=-cos(2*pi*tt/Ttheo); % there is a '-' in front of the cosinus every 4*pi
plot(tvec,e,'b.-',tt,etheo,'r-');
title('Validation');
legend('numerical','theory');
xlabel('time'); ylabel('f(L/2)');
print('-dpng','-r80','water_wave_1D_model_theory_numerical');
set(gcf,'paperpositionmode','auto');

%{
# Figures 

![Evolution of the position and the velocity](/sandbox/easystab/stab2014/water_wave_1D_model_position_veolocity.gif)

As you can see, when the wave is at it's maximum height, the veolicity diminishes, pulling the wave downwards (actually creating the wave).  

![Comparision theory / numerical solution](/sandbox/easystab/stab2014/water_wave_1D_model_theory_numerical.png)

When the number of grid points increases, the difference between the analitical and the numerical solution tends towards 0.
%}