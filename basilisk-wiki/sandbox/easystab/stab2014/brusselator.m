%{
# THE BRUSSELATOR
%}

%{
U and V are the concentrations of two chemical componants that varies with time and in the direction x. $\lambda$ is a chemical paramter we chose, which will have an impact on the stability of the system. You can change the value of $\lambda$. Here is the equation system that describes the brusselator :

$$
U_{t} = 1 - (\lambda + 1)U + 2U_{xx} + U²V          (i)
$$
$$
V_{t} = \lambda *U + V_{xx} + U²V          (ii)
$$
%}

%{
First, we calculate the linear and stationary form of the system. If the system is stationary, then 
$$ U_{t} = V_{t} = 0 $$ 
Then from (i) we have,
$$ V = \lambda /U $$ 
which we replace in (ii), which gives leaves us with : $$ 1 = (\lambda + 1)U - \lambda *U  
$$
We then disrupt the system by posing :

$$
U = U_b + u
$$
$$
V = V_b + v
$$
%}

%{
We now have the linear form of the brusselator (by negelecting terms containing U², product UV, etc.) : 

$$
U_t = 2U_{xx} + (\lambda - 1)U + V
$$
$$
V_t = - \lambda *U + V_{xx} - V
$$
%}

%{
THIS PROJECT CONTAINS 2 CODES : EIGENMODES & MARCH IN TIME. DO NOT FORGET TO COMMENT ONE CODE BEFORE RUNNING THE ONE YOU WANT.
%}

%{
# Eigenmodes of the brusselator (linear form)

Let us use [diffusion_eigenmodes.m](/sandbox/easystab/diffusion_eigenmodes.m) to calculate the eigenvectors and eigenmodes of our system. By using the differentiation matrices as we learnt from [diffmat.m](/sandbox/easystab/diffmat.m) and [diffmat_2D.m](/sandbox/easystab/diffmat_2D.m), we can turn the linear system into a matrix system :

$$ Eq_{t}= Aq
$$

$$
\begin{pmatrix}
I & Z \\
Z & I \\
\end{pmatrix}
\left(\begin{array}{l}
U_{t} \\
V_{t} \\
\end{array}
\right)
=
\begin{pmatrix}
( \lambda -1)I+2D_{xx} & I \\
- \lambda *I & D_{xx}-I \\
\end{pmatrix}
\left(\begin{array}{l}
U \\
V \\
\end{array}\right)
$$
%}

clear all; clf

% parameters
N=100; % number of gridpoints
L=1; % domain length

lambda = 0; % chemical parameter

% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 
II=eye(2*N);

% system matrices
E=II;
A(1:N,1:N)=(lambda-1)*I+2*dxx;
A(N+1:2*N,1:N)=-lambda*I;
A(1:N,N+1:2*N)=I;
A(N+1:2*N,N+1:2*N)=dxx-I;

%{
## Boundary Conditions
For now, we will use homogenous Dirichlet conditions as the boundary conditions. We will later try replacing them by periodic boundary condtions.
%}

% boundary conditions
loc=[1,N];
E(loc,:)=0; 
E(N+loc,:)=0;
A(loc,:)=II(loc,:);
A(N+loc,:)=II(N+loc,:);
         
%{
We compute the eigenmodes using the function *eig*. We then sort the modes according to decaying real part of the eigenvalue. With this choice, the first eigenvalue will be the one with the largest real part. We then remove the eigenmodes for which the eigenvalue is larger than 1000. We do this because since we have the matrix $E$ to impose the boundary conditions, the system is algebro-differential, with the concequence that there will be some infinite eigenvalues corresponding to the fact that the constraints are imposed infinitely fast (their dynamics is infinitely rapid).

If you want to see the eigenvalues and the eigenvectors, plot the real and imaginary parts components of s and U, s, being the eigenvalue and U being the eigenvector.
%}

% computing eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); 
s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

%{
We generate a random function for the initial condition via the function "randn".
%}

% initial condition

% show the evolution of an initial condition
n=5; % number of eigenmodes in the initial condition
a=randn(n,1); % the weights

%{
## Plotting

We want to see the evolution of the eigenmodes. So, we plot the eigenmodes calculated earlier in a time loop using "drawnow" function. Furthermore, we save the resulting graph as an animated image. To do that, we used the code found in [http://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab](http://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab). If you do not want to do this, put on comment the lines following "drawnow" and "filename = 'brusselator_eigenmode_lambda4.gif'". 
%}

% time loop
figure(1)
filename = 'brusselator_eigenmode_lambda4.gif';

for t=1:1:60;
    
    % the present state
    q=U(1:N,1:n)*(a.*exp(s(1:n)*t));
    plot(x,q);
    hold on
    
    % draw each of the components
    for gre=1:n
       plot(x,U(1:N,gre)*a(gre)*exp(s(gre)*t),'r--');
    end
    hold off
    
    axis([0,L,-2,2]); grid on
    xlabel('domain'); ylabel('Evolution of the Eigenmodes'); title(t); 
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    if t == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
end

%{
## Figures 
We present you here the first 5 eignemodes of the brusselator generated for $\lambda$ = 2, then $\lambda$ = 4. Let us note that for $\lambda$ = 2, the system is stable and for $\lambda$ > 2, it diverges.

![$\lambda$ = 2 (stable)](brusselator_eigenmode_lambda2.gif)
![$\lambda$ = 4 (unstable)](brusselator_eigenmode_lambda4.gif)


HERE ENDS THE RESOLUTION OF THE BRUSSELATOR WITH EIGENMODES. THE FOLLOWING PART IS THE RESOLUTION WITH THE MARCH OF TIME. IT IS A SECONDE CODE.

%}

%{
# March in Time of the brusselator (nonlinear form)

We use [advection1D.m](/sandbox/easystab/advection1D.m) to do the march in time of the brusselator. We keep the the linear system Eq_{t}= Aq  previously described and add the two nonlinear parts which (stored in a vector) are :

$$ 
nl_{1} = \lambda *U² + 2U*V + U²*V
$$
$$
nl_{2} = - nl_{1}
$$
Let us remind ourselves that these nonlinear parts were neglected in the linear form (wfet we disrupted the system with new forms for U and V).

The system is now :

$$ Eq_{t}= Aq + nl 
$$
$$ where nl is the vector containing the non linear parts
$$
First, let us deal with the linear part of the system, which is exactly the same as seen in the "March in Time". What differs is the resolution.
%}

%{
## Resolution
%}

clear all; clc; clf; close all

% parameters
N=100; % number of gridpoints
L=20; % domain length
dt=0.1; % time step
tmax=50; % final time
x0=L/8; %x-coordinate at initial time
l0=0.5; length width

lambda=2; % chemical parameter

% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale; 
dxx=DM(:,:,2)*scale^2;
x=(x-1)/scale; 
Z=zeros(N,N); II=eye(2*N); I=eye(N);

% system matrices
E=II;
A(1:N,1:N)=(lambda-1)*I+2*dxx;
A(N+1:2*N,1:N)=-lambda*I;
A(1:N,N+1:2*N)=I;
A(N+1:2*N,N+1:2*N)=dxx-I;

% boundary conditions
loc=[1,N];
E(loc,:)=0; 
E(N+loc,:)=0;
A(loc,:)=II(loc,:);
A(N+loc,:)=II(N+loc,:);

%{
## March in time

Here we use the same method as [vibrating_string_spring_dissipation.m]().
%}

% march in time matrix 
Mm=(E-A*dt/2);
M=Mm\(E+A*dt/2);

%{
## Initial condition

We chose an initial condition rather than generating a random function like earlier but you can do either.
%}

% initial condition
q=[2*sin(2*pi*x/L); 3*sin(8*pi*x/L)]; 
ql=q;

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);
nl=zeros(2*N,1);

%{
## Adding the nonlinear parts to the calculus and Plotting

If you follow the same methode as seen in [vibrating_string_spring_dissipation.m](), the system becomes :
$$ q(t+dt)=\frac{\frac{(E-A*dt/2)q}{E+A*dt/2}}{nl*dt} 
$$

As for the plotting, we use the method previously used.
%}

figure(1)

for ind=1:Nt    
    
    % Nonlinear terms
     nl(1:N)=2*q(1:N).*q(N+1:2*N) + lambda*q(1:N).*q(1:N) + q(1:N).*q(1:N).*q(N+1:2*N);
     nl(N+1:2*N)=-(2*q(1:N).*q(N+1:2*N) + lambda*q(1:N).*q(1:N) + q(1:N).*q(1:N).*q(N+1:2*N));
     nl(loc)=0; % homogenous boundary conditions of the non linear term
     nl(loc+N)=0;
     
     ql=M*ql;
     q=M*q+Mm\nl*dt; % one step forward
     q=M*q; 
     
    % plotting
    plot(x,q(1:N),'b',x,q(N+1:N*2),'b.-',x,ql(1:N),'r',x,ql(N+1:N*2),'r.-'); 
    %plot(x,q(1:N),'b',x,q(N+1:N*2));
    axis([0,L,-2,2]); grid on
    legend('Linear part U','Linear part V','Nonlinear part U','Nonlinear part V'); title('March in Time of the Brusselator')
    xlabel('x'); ylabel('q');
    drawnow 
    
end

print('-dpng','-r80','brusselator_march_in_time_lambda2');
set(gcf,'paperpositionmode','auto');

%{
## Figures 
We present you here the position of the brusselator for $\lambda$ = 2, then $\lambda$ = 4. Let us note that for $\lambda$ = 2, the system is stable and for $\lambda$ > 2, it diverges.

For $\lambda$ = 2 : the solution converges. The system comes back to its equilibrium. Run the code with a higher time if you do not see this.

![$\lambda$ = 2 (stable)](brusselator_march_in_time_lambda2.png)

For $\lambda$ = 4 : the solution diverges.

![$\lambda$ = 4 (unstable)](brusselator_march_in_time_lambda4.png)

%}

%}