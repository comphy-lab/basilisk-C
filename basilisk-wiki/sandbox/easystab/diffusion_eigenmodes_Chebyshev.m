%{
# Eigenmodes of the diffusion equation

This is a pedagogical introduction for the computation of the eigenmodes of a 1D system.

This program is adapted from diffusion_eigenmodes.m but uses 
"Chebychev pseudo-spectral" discretization

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
discretization = 'cheb'; % (Chebyshev discretization)
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
saveas(gcf,'diffusion_spectrum_cheb','svg');

%{
![**Figure :** Spectrum ](diffusion_eigenmodes_Chebyshev/diffusion_spectrum_cheb.svg)
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
saveas(gcf,'diffusion_eigenmodes1_cheb','svg');
%{
![**Figure :** Structure of the three first eigenvectors](diffusion_eigenmodes_Chebyshev/diffusion_eigenmodes1_cheb.svg)
%}
