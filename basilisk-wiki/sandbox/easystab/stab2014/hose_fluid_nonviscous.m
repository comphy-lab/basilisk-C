%{
# Hose wherein flows a non viscous liquid

Here we study here the “Hose wherein flows a non viscous liquid”, there may be a instability which the tube is deformed in the form of a wave.../hose_fluid_nonviscous.m. For a description, please see the EX5. [http://www.lmm.jussieu.fr/~hoepffner/enseignement/stab2011_td2.pdf](http://www.lmm.jussieu.fr/~hoepffner/enseignement/stab2011_td2.pdf) and [http://www.lmm.jussieu.fr/~hoepffner/enseignement/stab2012_td4correc.pdf](http://www.lmm.jussieu.fr/~hoepffner/enseignement/stab2012_td4correc.pdf)

f is the position of points of hose that varies with time and in the direction $x$. ksi and beta are two parameters which will have the impact on the stability of the system. ksi corresponds to the tube voltage and beta corresponds to the speed of the fluid in the pipe. normally the tube with a high voltage is stability and the fluid with a high speed is instability. We do this ex like the code [wave_like.m](../wave_like.m). Here is the equations:
$$
\begin{array}{l}
fxxxx+(1-ksi)fxx+2sqrt(\beta )fxt+ftt=0
\end{array}
$$
To run this code you may need Chebychev differentiation matrix [here](../chebdif.m), and Fourier differentiation matrix [here](../fourdif.m).
%}


clear all; clf

% parameters
N=50; % number of gridpoints
L=4*pi; % domain length
dt=0.05; % time step
tmax=20; % final time
ksi=0.9; % parameter corresponds to the tube voltage
beta=0; % parameter corresponds to the speed of the fluid

% differentiation matrices
[dx,dxx,d.wx,x]=dif1D('fou',0,L,N,5);
dxxxx=dxx*dxx;

Z=zeros(N,N); I=eye(N); 
II=eye(2*N);
ddx=blkdiag(dx,dx);

%{
# System matrices

Here we build the matrices that give the discretization of the equations. The equations is
$$
f_{tt}=c^2f_{xx}
$$ 
This is the classical vibrating string equation, also known as D'Alembert's equation. This is also known as the *wave equation*, here just in 1D. The vibration of the string will depend upon how heavy it is (the mass per unit length), and how tensed it is (the forced applied onto the string by the two attachment points). The heavier it is, the slower will th wave be, and the more tensed, the faster, so these two parameters compete for the wave speed *c*. Here we have written the  equations already with the parameterization that let *c* appear.

This is an equation with two time derivatives and to apply always the same methods, we will transform it into an equation with a single derivative, simply by augmenting the state vector: instead of describing the state of the system by only the position of the string, we will store the position of the string and also the velocity of the string (for every gridpoint). The state vector is thus
$$
q=\left(
\begin{array}{c}
v \\
f \\
\end{array}\right)
$$
once this done, we need to tell the system that $v$ is indeed the velocity, that is
$$
f_t=v
$$ 
and then the systems dynamics in terms of both $v$ and $f$
$$
v_t=c^2f_{xx}
$$
we thus have the array representation for the system
$$
\left(
\begin{array}{cc}
I & 0 \\
0 & I \\
\end{array}\right)
\left(\begin{array}{c}
v_t \\
f_t \\
\end{array}\right)=
\left(
\begin{array}{cc}
0 & c^2D_{xx} \\
I & 0 \\
\end{array}\right)
\left(\begin{array}{c}
v \\
f \\
\end{array}\right)
$$
Here we have put explicitely a large identity matrix on the left because this will be useful to impose the boundary conditions. We thus have a linear system of the form
$$
Eq_t=Aq
$$
%}

% system matrices
E=[-I,Z; Z,I];
A=[2*sqrt(beta)*dx,dxxxx+(1-ksi)*dxx; I, Z]; 

% locations in the state
v=1:N;
f=N+1:2*N;

% computing eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); 
s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

% show the eigenvalues
plot(real(s),imag(s),'b.')
xlim([-1,1])
ylim([-10,10])
xlabel('real part');ylabel('imaginary part');title('eigenvalues');
grid on;hold on

%{
# Boundary conditions
The original system has one variables and four space derivatives, so we should apply 4 boundary conditions. We will simply say that the position of the string at its both ends should be zero: two attachment points. This is a homogeneous Dirichlet condition applied on the first and last point of the state vector for the state position. And the Newmann condition applied on the first and last point too.
For this problem, we consider our domain as infinity, so do not need boundary conditions.
%}
% % boundary conditions
% dir=[f(1),f(N)];
% neu=[v(1),v(N)];
% loc=[dir; neu];
% E(loc,:)=0; 
% A(loc,:)=[II(dir,:); ddx(neu,:)];

%{
# March in time

To write a relation between the state at a given time and the state a little later, we will do a numerical approximation of the time derivative. We use the Crank-Nicolson scheme, which is quite simple and robust. We just integrate in time the evolution equation from time $t$ to time $t+dt$
$$
\int_t^{t+dt} E q_t dt=\int_t^{t+dt} Aq dt
$$ 
the integral on the left hand side is just $E(q(t+dt)-q(t))$ since $E$ does not change in time in this example, and the right hand side we approximate with the trapezoidal rule for integration
$$
E(q(t+dt)-q(t))\approx A dt (q(t+dt)+q(t))/2
$$ 
since now we want to express the new state $q(t+dt)$ as a function of the old state $q(t)$, we rearrange
$$
(E-Adt/2)q(t+dt)\approx (E+A dt/2) q(t)
$$ 
and we build the march-in-time matrix by inverting the left hand side matrix to transform this implicit system into an explicit system
$$
q(t+dt)\approx (E-Adt/2)^{-1}(E+A dt/2) q(t)
$$ 
This matrix inverse is well-defined since we have already imposed the propoer boundary conditions. 
%}
% march in time matrix 
M=(E-A*dt/2)\(E+A*dt/2);

%{
# Initial condition
We need to describe the state of the system at initial time. We say here that the string is initially deformed as a cosnus and two ends stay zero, the velocity is zero. From that state we then perform a loop to repeatedly advance the state of a short forward step in time. We use drawnow to show the evolution of the simulation as a movie when running the code.
%}
% initial condition
q=[zeros(N,1); (1-cos(2*pi*x/L))/2]; 

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

for ind=1:Nt    
    q=M*q; % one step forward
    figure(3)
    % plotting
    plot(x,q(f),'b',x,q(v),'r--');    
    axis([0,L,-2,2])    
    drawnow
end
legend('position','velocity'); title('Vibrating string')
xlabel('x'); ylabel('f');

%{
# Validation

To validate this computation, we compare it to the eigenvalues we get by assuming a domain invariant in the *x* direction, solving mode by mode.We do the wavelike assumption
$$
f(x,t)=\hat{f}\exp(i\alpha x+st)+C.C.
$$
So the system of equations becomes an eigenvalue problem for eigenvalue $s$, the eigenvalues are the roots of the caracteristic polynomial of this system.
And then we can get the dispersion relation.
$$
s^2+s(2i\alpha*sqrt(\beta))+\alpha^2(\alpha^2-(1-ksi))=0
$$
we compute using the function *roots*. After we plot it and compare between the theoretical solutions and the numerical ones. 


%}
% validation
for n=[1,2,3,4,5];
    alpha=n*2*pi/L; % wavenumber
    stheo=roots([1,2*alpha*sqrt(beta)*i,alpha^2*(alpha^2-(1-ksi))]);
    figure(1)
    plot(real(stheo),imag(stheo),'ro');
end

figure(2)
co='brk';
numvec=[1:3];
for ind=1:length(numvec);
    num=numvec(ind);
    plot(x,imag(U(f,num)),co(ind));
    hold on
end
legend('mode 1','mode 2','mode 3');
xlabel('x');  ylabel('eigenvectors');title('eigenvectors');
grid on; 

%{
Here is the figure at time *tmax*:
<center>
![](/sandbox/easystab/stab2014/hose_fluid_nonviscous.png)
</center>
<center>
![Eigenvalues](/sandbox/easystab/stab2014/hose_fluid_nonviscous_eigenvalues.png)
![Eigenvectors](/sandbox/easystab/stab2014/hose_fluid_nonviscous_eigenvectors.png)
</center>
$$
Eigenvalues and eigenvectors
$$
%}
