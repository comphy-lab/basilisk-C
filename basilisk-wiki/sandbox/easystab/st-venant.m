%{

* Overall description

%}

clear all; clf;

% Parameters
N=100;       % Number of gridpoints
b=0;         % fluid viscosity
alpha=0.5 ;  % Wavenumber in x
g=9.81;      % Gravity 
H=100;       % Depth



% Differentiation matrices

 I=eye(N);
 II=blkdiag(I,I);
 Z=zeros(N);
%{
# System matrices

* Equations

In the following of easystab, we do this in all the codes, we do not build the system matrices with $i\alpha I$, instead, we 
create a new differentiation matrix in $x$: dx=i*alpha*I to simplify the notation. We do the same for the second derivative in $x$: dxx=-alpha^2*I.


%}

% Naming the matrices
dx=i*alpha*I; dxx=-alpha^2*I;

% Location vectors
u=1:N; eta=u+N;
% System matrices

A=[-b*I,-g*dx; ...
   -H*dx,Z];
E=II;
   

%{
# Boundary conditions
We consider our domain as infinite, we do not need boundary conditions.   
   
%}

%{
# Eigenmodes

We compute the eigenvalues and the corresponding eigenvectors of the matrix system $(E,A)$. Above we said quickly that our linear systems are in the shape
$$
q_t=Aq
$$
thus the solution is an exponential like
$$
q(t)=q_0\exp(st)
$$
but in fact the truth is a little more subtle than this. This exponential solution is
true for a scalar system (when $q$ is not a vector, but is just a scalar value). For matrix systems there are many particular solutions that behave like this exponential, these are the eigenvectors. In general for our purpose there are as many eigenvectors as degrees of freedom in $q$. So there are $N$ different values of $i$
for which
$$
s_iEq_i=Aq_i
$$
where $s_i$ is a (possibly complex) eigenvalue and $q_i$ is the associated eigenvector. In most cases that we are interested here, the eigenvectors form a *basis* of the space of $q$. Thus for whatever initial condition $q_0$ you chose, you can write it as a sum of the eigenvectors
$$
q_0=\sum_{i=1}^N a_i q_i
$$
so if you know the eigenvectors, you can know the evolution of any initial condition to your system. The evolution is then
$$
q(t)=\sum_{i=1}^N a_i q_i\exp(s_i t) 
$$
where we have used the fact that we know exactly how each eigenvector will evolve in time using the exponential with the associated eigenvalue $s_i$.

For most of the cases that we are interested, the question is not 

> "how will an initial condition evolve?", 

but is rather 

> "is it possible that the system will explode?"

Usually there is just one, or just a few of unstable eigenmodes (an eigenvector with associated eigenvalue with positive real part such that it grows in time). So, if there is one such eigenmode, however small the initial condition amplitude $a_i$ on this mode, it will grow exponentially in time and the system will most likely explode. This is why we are usually only interested in just the few "least stable" eigenmodes, the ones that have the largest growth rate (of the least decay rate). This is why below for the validation, we just show the four first eigenvalues.

This approach is the "classical" approach to stability, and we will see later that 
there are some "modern" cases where this approach is too simple and must be revised, see [transient_growth.m]().

%}

% Compute eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

%{
# Validation

Comparison between theoretical solution of the dispersion relation and
numerical ones.
%}

p=[1 b g*H*alpha^2];

r=roots(p);

hold on
plot(real(r(1)),imag(r(1)),'g*',real(r(2)),imag(r(2)),'r*');
plot(real(s),imag(s),'b*','markersize',50);
% plot(real(s3),imag(s3),'y.',real(s4),imag(s4),'m.',real(s5),imag(s5),'c*');
grid on
hold off


