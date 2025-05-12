%{
# Numerical resolution of the "Ginzburg-Landau" Equation 


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

- $\sigma$ is the local growth rate. The latter may be constant ($\sigma = \sigma_0$) or a function of space with the form $/sigma = \sigma(x) = \sigma_0-\sigma_2 x^2$. 

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

In the next chapters, the Ginzburg-Landau will also be used to model flow
instabilities, and in this context it is more relevant to take $N_{nl}=3$.


%}

clear all; clf
close all;

% physical parameters
L=1 % domain length
x0 = 0 % location of left boundary
kappa=.1 % diffusion coefficient
sigma0 = 3.5 %maximum amplification rate (CONTROL PARAMETER)
sigma2 = 0 % variation of sigma with x
NL = 2 % order of nonlinearities ; should be 2 or 3
amp =1e-4 % Initial amplitude
typecondinit=2; % 1 for single-mode initial condition, 2 for random initial condition.

% numerical parameters for space discretisation
N=100; % number of gridpoints
bctype = 'Dirichlet'; %boundary conditions may be Dirichlet or Neumann
Tstep = .5;%time between two plots in figure 2
Tmax = 8 % maximum time for plots



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

- The "weight" $w$ (a line-vector used to cumpute integrals over the domain ; see [integration.m]() for more details and examples)

%}

discretization = 'fds'; % fds for simple finite differences ; you may try other cases
[dx,dxx,wx,x]=dif1D(discretization,x0,L,N);
Z=zeros(N,N); I=eye(N); II=eye(2*N);

%{
We now use these "basic bricks" to build the matrices A and E corresponding to the discretized Ginsburg-Landau equation :
%}

E=I;
sigma = sigma0-sigma2*x.^2; 
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

% computing eigenmodes
[U,S]=eig(A,E); 
s=diag(S);  [t,o]=sort(-real(s)); 
s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

% show the eigenvalues
figure(1)
subplot(2,1,1);
plot(real(s),imag(s),'b.')
xlim([min(s(5)-1,0),max(1,max(real(s)+1))]);% adjust the range to show the first five
xlabel('real part');ylabel('imaginary part');title('eigenvalues');
grid on; 

% show the eigenvectors
subplot(2,1,2);
co='rbmck';
for ind=1:4
plot(x,real(U(:,ind)),co(ind),x,imag(U(:,ind)),[co(ind) '--']);
hold on
end
xlabel('x');  ylabel('eigenvectors');title('eigenvectors');
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

In this program we compare the numerical results with the theory,
restricting to the most-amplified eigenvalue $\lambda_1$.

%}



disp('Leading eigenvalue :') ; s(1)

if(sigma2==0)
    if(bctype =='Dirichlet')
    disp('Theoretical value (bounded, Dirichlet) : ') ; sigma0-kappa*(pi/L)^2
    disp('Threshold : ') ; 
    sigmas = kappa*(pi/L)^2
    elseif(bctype =='Neumann')
    disp('Theoretical value (bounded, Neumann : ') ; sigma0
    disp('Threshold : ') ; 
    sigmas = 0
    end
else
    disp('Theoretical value (unbounded, inhogeneous) ') ; sigma0-sqrt(kappa*sigma2)
    disp('Threshold : ') ; 
    sigmas = sqrt(kappa*sigma2)
end



%{

### Linear initial value problem

We consider an initial condition of amplitude "amp" and structure
corresponding either to the leading mode or randomly selected.

%}



if(typecondinit==1)
        Psi = amp*real(U(:,1)); 
elseif(typecondinit==2)
        Psi =amp*rand(1,N); Psi(1)=0;Psi(N)=0;
        Psi = Psi';% must be a column vector
end

%{

The general solution of the problem has the following form :
$$
\phi(x,t) = \sum_{n=1}^\infty c_n exp(\lambda_n t) \hat \phi_n(x)
$$
where the coefficients $c_n$ are obtained by projecting the initial condition
upon the corredonding modes :

$$
c_n = \frac{\int_\Omega \phi(x,0) \hat \phi_n(x) dx}{\int_\Omega \hat \phi_n(x)^2 dx}
$$

We compute the $c_n$ coefficients for the first four modes, and we plot
the amplitudes $A_n(t) = c_n exp (\lambda_n t)$ as function of time.

%}

for n=1:4
    c(n) = (wx*(Psi.*U(:,n)))/(wx*U(:,n).^2);
end


figure(2);subplot(2,1,1);
plot(x,Psi);hold on; title('Initial condition');xlabel('x');ylabel('\psi(x,0)');
figure(2);
subplot(2,1,2);
ttab = linspace(0,Tmax,100);
for n=1:4
    semilogy(ttab,abs(c(n)*exp(s(n)*ttab)),[co(n),'--']); hold on;
end
ylim([1e-8,1]);
title('Amplitudes A_n(t) ; linear (dashed curves)');xlabel('t');ylabel('A_n(t)');
print('-dpng','-r80','GinsburgLandau_Amplitudes_Linear.png');

disp('Program paused after linear calculations ; press enter to continue');
pause;


%{

![**Figure 2:** Initial condition and time-evolution of the amplitudes 
of the first four modes according to linear theory](/sandbox/easystab/david/GinsburgLandau_Amplitudes_Linear.png)



## Nonlinear dynamics : 

### Single-mode approximation

In the nonlinear range, the solution can still be projected onto the linear
eigenmodes, but the time-evolution will not be exponential. So we take the following
ansatz :


$$
\phi(x,t) = \sum_{n=1}^\infty A_n(t) \hat{\phi}_n(x) 
$$

Thanks to the orthogonality property of the eigenmodes, we can obtain amplitude equations for each of the $A_n(t)$ by projecting
upon the corresponding mode :

$$
\frac{d A_n}{dt} = \lambda_n A_n - 
\frac{\sum_\Omega \left(\sum_1^\infty A_n(t) \hat \phi_n(x) \right)^{N_{nl}} \phi_n(x) dx}{\int_\Omega \phi_n(x)^2 dx} 
$$


We note $r = \lambda_1 = \sigma_0 - \sigma_s$ the "control parameter", with $\sigma_s$ the linear threshold.
If we assume that the leading mode remains dominant (A_1(t) \gg (A_2(t),A_3(t),...) for all times,
we can drop the higher modes. The previous equation thus leads to 

$$
\frac{\partial A_1}{\partial t} = r A_1 - \beta A_1^{N_{nl}}
$$

with 
$$
\beta = \frac{\int (\hat \phi_1)^{N_{nl}+1} dx}{\int (\hat \phi_1)^2 dx}
$$

We recognize the classical bifurcation equation, whose solution is the
following :

- For $N_{nl}=2$ (transcritical bifurcation for $r>0$)

$$
A_1(t) = \frac{r/\beta}{1+g_1 exp (-r t)} \quad \mbox{ with } g_1
= \frac{r}{\beta A_{1,0}} - 1
$$


- For $N_{nl} = 3$ (supercritical pitchfork bifurcation for $r>0$)
$$
A_1(t) = \frac{\sqrt{r/\beta}}{\sqrt{1+g_1 exp (-2 r t)}} \quad \mbox{ with } g_1
= \frac{r}{\beta A_{1,0}^2} - 1
$$
%}

r = s(1);
ttabsinglemode = linspace(0,Tmax,100);
if(NL==3)
    beta = wx*U(:,1).^4/(wx*U(:,1).^2);
    Ainf = sqrt(r/beta);
    g1 = r/(beta*c(1)^2)-1;
    A1singlemode = sqrt(r/beta)./sqrt(1+g1*exp(-2*r*ttabsinglemode)); 
elseif(NL==2)
    beta = wx*U(:,1).^3/(wx*U(:,1).^2);
    Ainf = r/beta;
    g1 = r/(beta*c(1))-1;
    A1singlemode = (r/beta)./(1+g1*exp(-r*ttabsinglemode)); 
end

disp('One-mode approximation predictions :')
disp(['   Linear amplification rate of dominant mode : ',num2str(r)]);
disp(['   Nonlinear coefficient beta = ',num2str(beta)]);
disp(['   Saturation amplitude = ',num2str(Ainf)]);
disp(['   Time scale to reach saturation : ',num2str(log(abs(Ainf/c(1)))/r)]);
disp(' ');



figure(2);subplot(2,1,2);%hold on;
semilogy(ttabsinglemode,abs(A1singlemode),'k+');hold on;
ylim([amp^2,5*abs(Ainf)]);
title('Amplitudes A_n(t) ; linear (dashed), one-mode approximation (symbols)')
disp('Press Enter to launch time-integration');
pause;

%{

### Nonlinear dynamics : time integration

In the last part of this program we illustrate qualitatively the effect of nonlinearities
by performing time-integration (direct numerical simulation) of the system.

Since the focus is to illustrate qualitatively the dynamics, we use a very
simple integration method, namely forward Euler. This imposes a small
time step to ensure the stability of the numerical scheme, namely :
%}

deltax = L/(N-1);
dt = deltax.^2/(5*kappa);

figure(2);subplot(2,1,1); 
title('Solution at several instants')
figure(2);subplot(2,1,2);
ylim([amp^2,5*abs(Ainf)]);
title('Amplitudes A_n(t) ; linear (dashed), one-mode approx. (symbols), numerical solution (full lines)')

%{

 The time-stepping loop is done as follows (and the amplitudes of the first four modes will be stored in the matrix Atab)

%}


for it = 1:(Tmax/dt);
    Psi = E*Psi+dt*(A*Psi-Psi.^NL);
    Psi(1)=0;Psi(N)=0;
    ttab(it) = it*dt;
    
    for ind = 1:4
        Atab(ind,it) = wx*(Psi.*U(:,ind))/(wx*U(:,ind).^2);
    end
    
    if(mod(it,round(Tstep/dt))==0) 
        figure(2);subplot(2,1,1);
        plot(x,Psi);
        figure(2); subplot(2,1,2);
        for ind = 1:4
             semilogy(ttab(1:it),abs(Atab(ind,1:it)),co(ind),'LineWidth',2)
        end
        pause(0.1);%to allow refreshing of figures
    end

end


print('-dpng','-r80','GinsburgLandau_Amplitudes_NonLinear.png');
disp('Program Ginsburg_Landau ended');


%{

![**Figure 3 :** Solution $\Phi(x,t)$ at several instants, and time-evolution
of the amplitudes of the first four modes ](/sandbox/easystab/david/GinsburgLandau_Amplitudes_NonLinear.png)


# Exercises/contributions

- Please play with the parameters and observe the dynamics 

- Please verify the theoretical solutions given above 

- Please do the same for other model 1D equations (swift-Ohenberg for instance)
