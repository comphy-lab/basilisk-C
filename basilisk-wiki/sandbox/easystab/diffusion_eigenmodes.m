%{


*This document belongs to the [easystab](http://basilisk.fr/sandbox/easystab/README) project, check the main page of the project to understand the general philosophy of the project.*

*The present program was specifically designed as a support for the course ["Introductions to hydrodynamical instabilities"](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) for the "DET" Master's cursus in Toulouse. For the underlying theory please see [Lecture notes for chapter 3](http://basilisk.fr/sandbox/easystab/NumericalMethodsForEigenvalueProblems.md)*

*To use this program: click "raw page source" in the left column, copy-paste in the Matlab/Octave editor, and run.*

# Eigenmodes of the diffusion equation




This is a pedagogical introduction for the computation of the eigenmodes of a 1D system. 
We then use these eigenmodes to see the evolution in time of an initial condition. 
%}

clear all; clf
% parameters
N=20; % number of gridpoints
L=pi; % domain length
mu=1; % diffusion coefficient


%{
# Theory

In this program we consider the 1D diffusion equation :

$$
\frac{\partial T}{\partial t}  = \mu \frac{\partial^2 T}{\partial x^2} 
\qquad \textrm{ for } \quad x \in [0,L]
$$
with boundary conditions :
$$
T_{x=0}  = T_{x=L} = 0
$$

We consider solutions in eigenmode form :
$$
T = \hat{T} e^{\lambda t}
$$

Eigenpairs $(\lambda,  \hat{T})$ are solutions of the following problem :

$$
\lambda \hat{ T } = \mu \frac{\partial^2 \hat{T}}{\partial x^2} 
$$

The analytical solution is as follows :
$$
\lambda_n = - \frac{\mu n^2 \pi^2}{L^2} ; \quad \hat{T}_n = \sin ( n \pi x /L ) 
$$

%}

stheory = -mu*pi^2/L^2*[1:1:N].^2;



%{
# Numerical resolution

## Discretization
 
 We use the function *dif1D* from the easystab project to construct the grid 
 *x* and the differentiation matrices *dx* and *dxx* (the 'weight' variable *wx* is not used here) 


%}  

%%
discretization = 'fds'; 
% try either 'fds' (finite differences), 'fd' (fd with sparse matrix storage) or  'cheb' (chebyshev)
[dx,dxx,wx,x]=dif1D(discretization,0,L,N); 

% NB to see the structure of the matrix try this:
%figure ; spy(dxx) 


%{
## Construction of the matrices 
%}
%%
Z=zeros(N,N); I=eye(N); 
B=I;
A=mu*dxx; 

%% boundary conditions
loc=[1,N];
B(loc,:)=0; 
A(loc,:)=-I(loc,:);
         
%{

## Resolution of the eigenvalue problem

We compute the eigenmodes using the function *eig*. We then sort the modes according to decaying real part of the eigenvalue. With this choice, the first eigenvalue will be the one with the largest real part. We then remove the eigenmodes for which the eigenvalue is larger than 1000. We do this because since we have the matrix $B$ to impose the boundary conditions, the system is algebro-differential, withthe concequence that there will be some infinite eigenvalues corresponding to the fact that the constraints are imposed infinitely fast (their dynamics is infinitely rapid).
%}

%% computing eigenmodes with direct method (matrix diagnonalization) 
[U,S]=eig(A,B);

% NB to compute only a selection of eigenmodes using Arnoldi method try this:
% [U,S]=eigs(A,B,10,'LR')


%% sort the eigenmodes
s=diag(S);  [t,o]=sort(-real(s)); 
s=s(o); U=U(:,o);


%% show the eigenvalues/eigenmodes
figure(1);hold off;
plot(real(stheory),imag(stheory),'ro');xlim([0,L]);
hold on;
plot(real(s),imag(s),'bx');
xlim([-25,1])
xlabel('real part');ylabel('imaginary part');title('eigenvalues');
legend('theory','computed')
grid on; 
saveas(gcf,'diffusion_spectrum','svg');

%{
![**Figure :** Spectrum ](diffusion_eigenmodes/diffusion_spectrum.svg)
%}


%% show the eigenvectors
figure(2); hold off;
co='brkmc';
for ind=1:3
plot(x,real(U(:,ind)),'+-',x,imag(U(:,ind)),[co(ind) '--']);
hold on
end
xlabel('x');  ylabel('T(x)');title('eigenvectors');
grid on; 
saveas(gcf,'diffusion_eigenmodes1','svg');
%{
![**Figure :** Structure of the three first eigenvectors](diffusion_eigenmodes/diffusion_eigenmodes1.svg)
%}

INTERACTIVE = isempty(getenv('OCTAVE_AUTORUN'));
if (INTERACTIVE)
  disp('click any key to continue...');
  pause;
end

%{

# Time evolution

Here we use the eigenmodes to show the time evolution of an initial condition. If the initial condition of our system is one of its eigenmodes, we know exactly all the time evolution, it will simply be the eigenvector multiplied by the exponential of the associated eigenvalue time the time. Going further one step, this means that if the initial condition is a combination of the eigenvectors, we get the evolution as sum sum of the eigenvectors weighted by their individual exponential time dependency.

With the initial condition
$$
T(x,0)=\sum_i \alpha_i \hat{T}_i(x)
$$
with $\alpha_i$ the amplitude of each eigenmode in the initial condition. The evolution in time is thus
$$
T(x,t)=\sum_i \alpha_i \exp(\lambda_i t) \hat{T}_i(x)
$$

%}

%% show the evolution of an initial condition
%figure(2)


%{
Here we simply build a random initial condition by combining the *n* least stable eigenvectors
%}
n=10; % number of eigenmodes in the initial condition
a=randn(n,1);a(1)=1; % the weights


figure(3);
% time loop
for t=linspace(0,1,51)
  
   % the present state
   q=U(:,1:n)*(a.*exp(s(1:n)*t));
   plot(x,q);
   grid on; xlim([0,L]);hold on;    
   % draw each of the components
   for gre=1:n
      plot(x,U(:,gre)*a(gre)*exp(s(gre)*t),'r--');
   end
   hold off
   title(['t = ',num2str(t)])
   legend('T(x,t)','eigenmode components')
   drawnow;
   if INTERACTIVE&&(t==0)
        disp('click any key to continue...');
        pause;
        disp('launching time evolution...');
   end
   pause(0.1)
   if (t==0.0)
         set(gcf,'paperpositionmode','auto');      
         print('-dsvg','-r80','diffusion_eigenmodes2.svg');
   end
   if (t==0.25)

          set(gcf,'paperpositionmode','auto');      
          print('-dsvg','-r80','diffusion_eigenmodes3.svg');
   end    
end





%{
![Initial condition and projection on eigenmodes](diffusion_eigenmodes/diffusion_eigenmodes2.svg)


Here we see in blue the state at the time of the snapshot, and in red the eigenvectors that we have used when building the initial condition, each multiplied with $exp(s_i t)$ where $s_i$ is the associated eigenvalue. This shows that quickly, the evolution of the state of the system tends to the evolution of the least stable eigenmode.

![A snapshot of the time evolution for t=0.25](diffusion_eigenmodes/diffusion_eigenmodes3.svg)

# Exercises/contributions

* Please do the same for other model 1D equations (advection, advection-diffusion...)
* Please do this for an unstable system (look for instance at the Ginsburg-landau example)


%}
