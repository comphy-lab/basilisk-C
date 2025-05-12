%{
# The Falkner Skan boundary layer

In this code, we solve the Falkner-Skan equation to get the boundary layer velocity profile with favourable or adverse perssure gradient like we did in [falkner_skan.m]().

Here we do the Keller arclength continuation to follow the branch when the pressure gradient parameter $\beta$ is progressively reduced (like we did in [meniscus_overturn_continuation.m]()). There is a fold on this nonlinear branch (there are two profile solutions for $\beta$ negative).

**Dependency**

* [chebdif.m]()
* [falkner_skan_continuation_ref.png]() a figure for the validation

%}

clear all; clf;
%setenv("GNUTERM","X11")

% parameters
L=25; % box height
N=150; % number of gridpoints
beta=1; % initial pressure gradient parameter
delta=0.01; % continuation length

% differentiation matrices
scale=-2/L;
[y,DM] = chebdif(N,3); 
dy=DM(:,:,1)*scale;    
dyy=DM(:,:,2)*scale^2;    
dyyy=DM(:,:,3)*scale^3;
y=(y-1)/scale;  	
I=eye(N); Z=zeros(N,N);
w=([diff(y)',0]+[0,diff(y)'])/2; % integration weights

%{
Use the Blasius solution from Pohlhausen order 4 polynomial as initial guess.
%}

% initial guess from order 4 Polhausen profile
d99=6*sqrt(35/37);
eta=y/d99; deta=dy*d99; % rescaled vertical coordinate
u0=1-(1-eta).^3.*(1+eta); u0(eta>1)=1;
A=deta; A(1,:)=I(1,:); u0(1)=0; % set up integration 
g0=A\u0; % compute integral

%initial guess
sol=[g0;beta];
dir=[Z(:,1);-1];   % initial direction

% Continuation loop
ind=1;quitcon=0;
bbeta=[]; ggyy=[];
while ~quitcon 
    
   solprev=sol;
   sol=sol+dir*delta; % extrapolation of solution
   
   % Newton iterations
   quit=0;count=0;
   while ~quit
    
      % the present solution
      g=sol(1:N); gy=dy*g; gyy=dyy*g; gyyy=dyyy*g;
      beta=sol(end);

      % nonlinear function
      f=[gyyy+g.*gyy+beta*(1-gy.^2); ...
          dir'*(sol-solprev)-delta];
    
      % analytical jacobian
      A=[dyyy+diag(g)*dyy+diag(gyy)-2*beta*diag(gy)*dy, 1-gy.^2; ...
         dir'];
     
      % Boundary conditions
      loc=[1,2,N];
      C=[I(1,:),0; dy(1,:),0; dy(N,:),0];
      f(loc)=C*sol-[0;0;1];
      A(loc,:)=C;
    
      % convergence test
      res=norm(f);
      disp([num2str(count) '  ' num2str(res)]);
      if count>50|res>1e5; disp('no convergence'); break; end
      if res<1e-5; quit=1; disp('converged'); continue; end
    
      % Newton step
      sol=sol-A\f;
      count=count+1;
   end
   
   % showing the results
   subplot(1,2,1)
   plot(gy,y,'b-',u0,y,'r--'); xlabel('u'); ylabel('y'); 
   title('velocity profile');
   axis([-0.5,1,0,L]); grid on; 
   
   subplot(1,2,2)
   bbeta=[bbeta,beta];
   ggyy=[ggyy,gyy(1)];
   plot(bbeta,ggyy,'b.-'); grid on;
   xlabel('beta'); ylabel('wall shear'); title('bifurcation branch');
   drawnow
   
   % new direction
   dir=A\[Z(:,1);1]; % the new direction
   dir=dir/norm(dir); % normalization
   
   % stop loop above given d1
   d1=w*(1-gy);
   if d1>15; quitcon=1; end
   
   % adjust continuation step
   if count==1; delta=delta*1.5; end
   if count>5; delta=0.5*delta; end
   
   ind=ind+1;
end

% loading the scanned figure for comparison and plotting
subplot(1,2,2);
a=imread('falkner_skan_continuation_ref.png'); ss=size(a);
xx=linspace(-0.4,2,ss(2));
yy=linspace(2.8,-0.4,ss(1));
image(xx,yy,a); axis xy; hold on

plot(bbeta,ggyy,'r-','linewidth',3); grid on;
xlabel('beta'); ylabel('wall shear'); 
title('validation');


set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','falkner_skan_continuation.png');



%{
# Result and Plot

The plot corresponds to the result from Stewartson 1964. There are several curves, that correspond to the parameter S in the book describing the compressibility (see exercices below). Our case is the incompresible one with $S=1$. Our curve is in red.

![The velocity profile](falkner_skan_continuation.png)

# Exercices/Contributions

* Please compute the Pohlhausen approximation of the $f''_0(\beta)$ curve 
* Please check that  $f''_0 \simeq 1.544(−\beta)^{3/4}$ for small $\beta$  (Brown & Stewartson 66)
* Please check that $f'''+1-f'^2=0,$ and $f''_0=1.15,$ and $\int_0^\infty (1-f')d\eta=0.779$ 
for large $\beta$ (convergent case)   
* Please solve the compressible problem $f'''+ff''+\beta (S-f'^2) =0$ and $S''+f S'=0$ (Stewartson 64) 
 
# Links
* the Blasius solution [blasius.m]() $\beta=0$
* the Falkner Skan solution [falkner_skan.m]() for one value of $\beta<0$ (with separation)
* the Hiemenz solution [hiemenz.m]() for $\beta=1$
* the Falkner Skan solution [falkner_skan_continuation.m]() for all the $\beta$ 
  
# Bibliography
* see [blasius.m]()
* Brown & Stewartson "On The Reversed Flow Solutions Of The Falkner-Skan Equation" Mathematica 1966,
* Stewartson "Compressible Boundary Layer" Oxford 1964
 
%}