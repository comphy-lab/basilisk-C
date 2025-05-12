%{
# Gravity free surface waves
  We show in [the general code](http://basilisk.fr/sandbox/easystab/free_surface_gravity.m#validation) that the wave velocity at the surface of the sea is: 
  
$$
c=\sqrt(g \tanh(\alpha L)/{\alpha})
$$

%}

clear all; clf;
alphamin=5;
alphamax=10;
for alpha=alphamin:alphamax  ;  % wavenumber in x
n=100;      % number of gridpoints
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
# Validation for short wavelength

If we have short wavelength compared to the fluid depth, then we have:

$$
tanh(\alpha L)=1
$$

so the wave velocity tends to :

$$
c = \sqrt({g}/{\alpha})
$$


thus the group velocity tends to:

$$
c_g=\sqrt{g \alpha}
$$

This shows that the group velocity of the waves depends on the wavenumber when this one become large: the system is dispersive.

%}


% validation
alphavec=linspace(alphamin,alphamax,100);
ctheo_general=sqrt(g*tanh(alphavec*L)./alphavec);
ctheo_shortwavelength=sqrt(g./alphavec);

cnum=abs(imag(s(1)))/alpha;

cnumvec(alpha+1-alphamin)=cnum;
alphavec2(alpha+1-alphamin)=alpha;
end

subplot(1,2,1)
plot(alphavec,ctheo_shortwavelength,'g-',alphavec2,cnumvec,'b.');
xlabel('alpha');ylabel('wave velocity'); title('validation')
legend('c_{theoretical approched}','c_{num}')
grid on

subplot(1,2,2)
plot(alphavec,(ctheo_shortwavelength-ctheo_general),'r-');
xlabel('alpha');ylabel('différence between general & approached solution'); title('validation')
grid on

set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','free_surface_gravity-short_wave_length.png');

%{ 
Here are the graphs for the wave velocity as a fonction of the wave number. You can see :

- in the first one the graph of the numerical and the approached theoretically approached solution

- in the second one, the difference between the exact and the approached solution

![wave velocity as a function of long wavenumbers $\alpha$](/free_surface_gravity-short_wave_length.png)