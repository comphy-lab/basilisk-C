
%{
# Mathematical analysis and Numerical resolution of the "Ginsburg-Landau" Equation 
# Part I : linear study


This program is a pedagogical introduction to the numerical resolution of
partial-difference equations systems displaying instability phnomena, 
starting from the eigenmode analysis of the linearised system to predict the stability, 
and ending up to the illustration of the dynamics in the full nonlinear case.

## Definition of the problem :

The Ginzburg-Landau equation is a simple model used for various processes in physics,
chemistry or biology. We define this equation as follows;

$$
\frac{\partial \phi}{\partial t} = 
\sigma \phi 
+ \kappa \frac{\partial^2 \phi}{\partial x^2}
-\phi^{N_{nl}}
$$

Here :

- $\phi(x,t)$ is the unknown function defined on an interval $x \in
\Omega$, which can be finite ($\Omega= [a,b]$) or infinite.

- $\sigma$ is the local growth rate. The latter may be constant ($\sigma = \sigma_0$) or a function of space with the form $\sigma = \sigma(x) = \sigma_0-\sigma_2 x^2$. 

- $\kappa$ is a diffusion coefficient. 

- $N_{nl}$ is the order of the nonlinearity (2 or 3).

- Boundary conditions on a,b may be Dirichlet ($\psi=0$) or Neumann
($\partial \psi / \partial x = 0$).

A possible interpretation of this equation, in biology, is as follows.
Consider the evolution of a population of organisms (bacteria for instance)
evolving in a one-dimensional medium (a tube for instance). $\Psi(x,t)$ 
will correspond to the local density of organisms. $\sigma(x)$ is the
local concentration of nutriments (assumed stationnary). The diffusion term
represents the motion of the organisms modelled as "random walkers", and the 
nonlinear term represent the competition between the organisms. In biology,
$\phi(x,t)$ is positive and nonlinearities can be modelled with $N_{nl}=2$.

In the next chapters, the Ginsburg-Landau will also be used to model flow
instabilities, and in this context it is more relevant to take $N_{nl}=3$.



In this program, sigma(x) is assumed as a polynomial of order 2 :
$$
\sigma(x) = \sigma_0 + \sigma_1 x + \sigma_2 x^2
$$

This program (and the next one allows to play with the parameters ;
the main parameters to chose are the coeffic

%}

clear all; clf
close all;

% physical parameters

% Here is where to specify the parameters when playing with the program !

xmin = 0 % location of left boundary
xmax = pi  % location of right boundary
kappa=.1 % diffusion coefficient
sigma0 = .3 %  parameter sigma_0
sigma1 = 0; % parameter sigma_1 
sigma2 = 0 %  parameter sigma_2
NL = 3 % order of nonlinearities ; should be 2 or 3
%amp =1e-4 % Initial amplitude
%typecondinit=2; % 1 for single-mode initial condition, 2 for random initial condition.




% numerical parameters for space discretisation
N=100; % number of gridpoints
bctype = 'Dirichlet'; %boundary conditions may be Dirichlet (phi =0) or Neumann (dphi/dx = 0) 
                      %(this choice applies at both boundaries)

% Show the domain and plot sigma(x)

x = linspace(xmin,xmax,N);
sigma = sigma0+sigma1*x+sigma2*x.^2; 


figure(1);
plot(x,sigma,'k');hold on;
xlabel('x');ylabel('\sigma(x)');title('Domain, growth rate and BC');
xlim([xmin-.5,xmax+.5]);%
ylim([min(sigma)-.5,max(sigma)+.5]);
plot([xmin xmin],[min(sigma)-.5,max(sigma)+.5],'r','Linewidth',2)
plot([xmax xmax],[min(sigma)-.5,max(sigma)+.5],'r','Linewidth',2)
plot([xmin xmax],[0,0],'k:')
text(xmin-.1,0,bctype(1))
text(xmax+.1,0,bctype(1))
pause(1);
                      
                      
%{ 
## Study of the linearised system.

### Mathematical analysis of the problem

The problem has an trivial solution $\phi(x,t) = 0$. We first consider
the linear stability of this solution by investigating the behaviour
of small-amplitude perturbations of this solution (i.e. $\phi \ll 1$).
This means that we can drop the nonlinear term and consider the linear
problem

$$
\frac{\partial \phi}{\partial t} = 
\sigma \phi 
+ \kappa \frac{\partial^2 \phi}{\partial x^2}
$$

We look for particular solutions in the form of eigenmodes, with the
following ansatz
$$
\phi(x,t) = \hat{\phi}(x) exp (\lambda t)
$$

Inserting in the starting equations, we are lead to the eigenvalue problem

$$
\lambda \hat{\phi} = \sigma(x) \hat{\phi} + \kappa \frac{\partial^2 \hat \phi}{\partial x^2}
$$



This problem is of Sturm-Liouville type, so it possesses a discrete number
of eigenvalues/eigenmodes $[\lambda_n, \hat{\phi}_n]$.


The problem has analytical solutions in a few cases (see last section of the program)
but the objective of this program is to do a numerical resolution.



%}

%{
### Numerical resolution of the eigenmode problem

We will transform the linear problem into a generalized eigenvalue problem
into the form

$$
\lambda E X = A X
$$ 

The spatial discretisation is done using the tools provided by the [easystab](http://basilisk.fr/sandbox/easystab/README) project
( see for instance the programs [diffmat.m]() and [diffmat_dif1D]() as a starting point to understand the way to use this set of programs.)

Here we use the function [dif1D]() to construct the following objects :

- The grid $x$ (a column-vector containing the gridpoints $x_{1..N}$),

- The first-order and second-order differentiation matrices $dx$ and $dxx$ (square matrices),

- The "weight" $w$ (a line-vector used to cumpute integrals over the domain ; not used in the present program but used in the nonlinear part;
                    see [integration.m]() for more details and examples)

%}

discretization = 'fds'; % fds for simple finite differences ; you may try other cases
[dx,dxx,wx,x]=dif1D(discretization,xmin,xmax-xmin,N);
Z=zeros(N,N); I=eye(N); II=eye(2*N);

%{
We now use these "basic bricks" to build the matrices A and E corresponding to the discretized Ginsburg-Landau equation :
%}

E=I;
sigma = sigma0+sigma1*x+sigma2*x.^2; 
A=kappa*dxx+diag(sigma); 

% then construct boundary conditions
if(bctype =='Dirichlet')
    E(1,:)=0; E(N,:)=0; 
    A(1,:)=I(1,:); A(N,:)=I(N,:);
elseif(bctype=='Neumann')
    E(1,:)=0; E(N,:)=0; 
    A(1,:)=dx(1,:); A(N,:)=dx(N,:);
end
    
%{
We compute the eigenmodes using the function *eig*. We then sort the modes according to decaying real part of the eigenvalue. 
With this choice, the first eigenvalue will be the one with the largest real part. We then remove the eigenmodes for which the 
eigenvalue is larger than 1000. We do this because since we have the matrix $E$ to impose the boundary conditions, 
the system is algebraic differential system, with the consequence that there will be some unphysical infinite eigenvalues.
%}

% computing eigenmodes and a few manipulations to sort them in proper order
[U,S]=eig(A,E); 
s=diag(S);  [t,o]=sort(-real(s)); 
s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
for ind = 1:4
  if(max(U(:,ind))<0.95)
    U(:,ind) = -U(:,ind); % to make sure that the first eigenmode is positive 
  end 
end

disp('computed leading eigenvalues:');
s(1:4)


% show the eigenvalues
figure(2);
subplot(2,1,1);
co='rbmck';%color codes
for ind=1:4
   plot(real(s(ind)),imag(s(ind)),[co(ind),'*']);hold on;
end
xlim([min(s(5)-1,0),max(1,max(real(s)+1))]);% adjust the range to show the first five
xlabel('Re(\lambda)');ylabel('Imag(\lambda)');title('Spectrum');
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4');
grid on; 

% show the eigenvectors
subplot(2,1,2);
for ind=1:4
   plot(x,real(U(:,ind)),co(ind))   %x,imag(U(:,ind)),[co(ind) '--']);
   hold on
end
xlabel('x');  ylabel('eigenvectors');title('eigenvectors');
legend('\psi_1(x)','\psi_2(x)','\psi_3(x)','\psi_4(x)');
grid on; 
set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','GinsburgLandau_eigenmodes.png');




%{ 
![**Figure 1 :** The eigenvalues and eigenvectors](/sandbox/easystab/david/GinsburgLandau_eigenmodes.png)


### Comparison with theoretical results 

Analytical solutions are available for a few particular cases. 
We discuss only three of them :

- If the growth rate is uniform and Dirichlet conditions are used
at both boundaries ($x =0 ; L$)
Then the eigenvalues/eigenmodes are as follows :
$$
[\lambda_n;\hat \phi_n(x)] = [ \sigma_0 - \kappa n^2 \pi^2/L^2 ;  
\sin (n \pi x / L ) ]  \quad (\mbox{ with } n=1,2,...)
$$

- If the growth rate is uniform and Neumann conditions are used at both boundaries 
($x = 0 \pm L$) 
Then the eigenvalues/eigenmodes are as follows :
$$
[\lambda_n;\hat \phi_n(x)] = [ \sigma_0 - \kappa n^2\pi^2/L^2 ;  \cos ( n \pi x / L ) ] 
 \quad (\mbox{ with } n=0,1,...)
$$

- If the domain is infinite ( $x \in [-\infty,\infty]$ and the growth rate
is a second-order polynomial ($\sigma = \sigma_0 - \sigma_2 x^2$ with
$\sigma_2>0$) , then the eigenvalues/eigenmodes are as follows :
$$
[\lambda_n;\hat \phi_n(x)] = [ \sigma_0 - (2n+1) \kappa /l_c^2;  exp(- x^2/2l_c^2)H_n(x) ] 
 \quad (\mbox{ with } n=0,1,...)
$$
where $l_c = (\kappa/\sigma_2)^{1/4}$ and $H_n(x)$ is the Hermite
polynomial of order $n$ (Note that the mathematical analysis of this case is similar 
to that of the quantum harmonic oscillator equation ; 
many online resources are available on this problem).


# Exercises/contributions

Run the program with the following sets of parameters (lines 69 to 75)

## First case : 

   sigma =  sigma_0 (constant), 
   
   kappa = 0.1, 
   
   [xmin,xmax] = [0,pi], Dirichlet conditions

  - Compute the eigenvalues using the program for different values of sigma_0.
  
  - Compare with the theoretical results
 
  - What is the minimum threshold sigma_0 = sigma0_c for instability ? what can we expect to happend in the nonlinear 
  regime for such parameters ?
  
## Second case : 

   sigma = sigma_0 + sigma_2 x^2 with sigma _2 = -.25  (sigma_1 = 0) 
   
   [xmin, xmax] = [-5,5], Dirichlet conditions
   
   kappa = 0.1

 - Compute the eigenvalues using the program for different values of sigma_0 starting from 1
  
 - Compare with the theoretical results (NB the theory assumes [xmin,xmax] = [-infty,infty], is this relevant here ?)
       
 - What is the minimum thresholdsigma_0 = sigma0_c for instability ? what can we expect to happend in the nonlinear 
  regime for such parameters ?
  

## Third case : 

  sigma = sigma_0 + sigma_1 x + sigma_2 x^2 with sigma1 = 0.5 and sigma _2 = 2 ; 
  
  [xmin, xmax] = [-2,2], Dirichlet conditions
  
   kappa = 0.1
   

 - Compute the eigenvalues using the program for different values of sigma_0 starting from -4
 
 - What is the minimum threshold $\sigma_0 = {\sigma_0}_{c,1}$ for instability ? 
 
 - What is the threshold $\sigma_0 =  {\sigma_0}_{c,2}$ for the onset of a second unstable mode ?
  
 - What can be expected to happpen if $\sigma_0 \ge {\sigma_0}_{c,1}$ ? and for $\sigma_0 \ge {\sigma_0}_{c,2}$ ?



%}