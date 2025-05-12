%{
# Gravity free surface waves
  We showed in [the general code](http://basilisk.fr/sandbox/easystab/free_surface_gravity.m#validation) that the wave velocity at the surface of the sea is: 
  
$$
c=\sqrt(g \tanh(\alpha L)/{\alpha})
$$

%}

clear all; clf;

alphamin=5;
alphamax=50;
j=0;        % counter's initialization 

for alphax=alphamin:alphamax
alpha=alphax/1000;    % wavenumber in x
n=10;      % number of gridpoints
j=j+1;      % counter
L=1;        % Fluid height in y
rho=1;      % fluid density
mu=0.0001;    % fuid viscosity
g=1;        % gravity

% differentiation matrices
scale=-2/L;
[y,DM] = chebdif(n,2); 
D=DM(:,:,1)*scale;    
DD=DM(:,:,2)*scale^2;    
y=(y-1)/scale; 
I=eye(n); Z=zeros(n,n);

% renaming the matrices
dy=D; dyy=DD;
dx=i*alpha*I; dxx=-alpha^2*I;
Delta=dxx+dyy;

% System matrices
A=[mu*Delta, Z, -dx, Z(:,1); ...
   Z, mu*Delta, -dy, Z(:,1); ...
   dx, dy, Z, Z(:,1); ...
   Z(1,:),I(n,:),Z(1,:),0];

E=blkdiag(rho*I,rho*I,Z,1);

% boundary conditions
loc=[1,n,n+1,2*n];
C=[I(1,:),Z(1,:),Z(1,:),0; ... 
   Z(1,:),I(1,:),Z(1,:),0; ...
   Z(1,:),Z(1,:),-I(n,:),rho*g; ...
   dy(n,:),dx(n,:),Z(1,:),0]; 

E(loc,:)=0;  
A(loc,:)=C;

% compute eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];


%{
# Validation for long wavelength

If we have long wavelength compared to the fluid depth, then we have:

$$
\alpha L << 1 
$$

so applying the Taylor expansion:

$$
tanh(\alpha L)=\alpha L
$$

so the wave velocity tends to :

$$
c = \sqrt(g L)
$$

This shows that the velocity of the waves does not depend on the wavelength when this one become large: the system is not dispersive.

%}


%validation
ctheo_longwavelength(j)=1;  %approached velocity
cnum=abs(imag(s(1)))/alpha; %numerical velocity
cnumvec(j)=cnum;            %vector of numerical velocities
alphavec(j)=alpha;         %vector of wavenumbers

end

ctheo_general=sqrt(g*tanh(alphavec*L)./alphavec); %theorical velocity

subplot(1,2,1)
plot(alphavec,ctheo_longwavelength,'g-',alphavec,cnumvec,'r.');
xlabel('alpha');ylabel('wave velocity'); title('validation')
legend('c_{theoretical approched}','c_{num}')
grid on

subplot(1,2,2)
plot(1./alphavec,(ctheo_longwavelength-ctheo_general),'r-');
xlabel('alpha');ylabel('différence between general & approached solution'); title('validation')
grid on

set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','free_surface_gravity-long_wave_length.png');

set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','free_surface_gravity.png');


%}

%{ 
Here are the graphs for the wave velocity as a fonction of the wave number. You can see :

- in the first one the graph of the numerical and the approached theoretically approached solution

- in the second one, the difference between the exact and the approached solution

![wave velocity as a function of long wavenumbers $\alpha$](/free_surface_gravity-long_wavelength.png)