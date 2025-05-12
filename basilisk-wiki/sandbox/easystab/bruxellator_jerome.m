%{

We study here the "Bruxellator" in a 1D version, as an exercice to compare the eigenmodes on a finite periodic domain with the theory of normal modes. For a description, please see [http://www.lmm.jussieu.fr/~hoepffner/enseignement/stab2011_td2.pdf](http://www.lmm.jussieu.fr/~hoepffner/enseignement/stab2011_td2.pdf) and [http://www.lmm.jussieu.fr/~hoepffner/enseignement/stab2012_td4correc.pdf](http://www.lmm.jussieu.fr/~hoepffner/enseignement/stab2012_td4correc.pdf)

U and V are the concentrations of two chemicals that varies with time and in the direction $x$. $\lambda$ is a chemical parameter we chose, which will have an impact on the stability of the system. Here is the nonlinear set of equations:
$$
\begin{array}{l}
U_{t} = 1 - (\lambda + 1)U + 2U_{xx} + U²V   (i)\\
V_{t} = \lambda*U + V_{xx} + U²V          (ii)
\end{array}
$$
%}

clear all; clf; format long

% parameters
N=100; % number of gridpoints
L=10; % domain length
lambda = 3; % chemical parameter

% differentiation matrices
[d.x,d.xx,d.wx,x]=dif1D('fou',0,L,N);
Z=zeros(N,N); I=eye(N); II=eye(2*N);

%{

# Linearization

First, we calculate the steady and uniform state of the system. If the system is stationary, then 
$U_t = V_t = 0$. Then from (i) we have, $V = \lambda U$ which we replace in (ii), which leaves us with : $1 = (\lambda + 1)U - \lambda U $. We then perturb the system as follows :
$$
U = U_b + u , V = V_b + v
$$
We now have the linearized equations: 
$$
\begin{array}{l}
U_t = 2U_{xx} + (\lambda - 1)U + V\\
V_t = -\lambda U + V_{xx} - V
\end{array}
$$
which can be rewritten in the matrix form
$$ 
Eq_{t}= Aq
$$
with
$$
\begin{pmatrix}
I & Z \\
Z & I \\
\end{pmatrix}
\left(\begin{array}{l}
u_{t} \\
v_{t} \\
\end{array}
\right)
=
\begin{pmatrix}
(\lambda-1) I+2D_{xx} & I \\
-\lambda I & D_{xx}-I \\
\end{pmatrix}
\left(\begin{array}{l}
u \\
v \\
\end{array}\right)
$$

%}
% system matrices
E=II;
A=[(lambda-1)*I+2*d.xx, I; ...
    -lambda*I, d.xx-I];

%{
# Boundary conditions

We have used Fourier discretization in $x$, so it means that we already assume that the solution is periodic. This is useful since the theory that we want to check assumes a periodic solution. We nevertheless impose something because we do not want to allow a constant solution to be an eigenmode of the system. So we use the integration matrix *d.wx* (see [integration_2D.m]() for more information) to impose that the constraint that the integral of *v* should be zero.

%}

% boundary conditions
u=1:N; v=N+1:2*N;
loc=[1];
C=[Z(1,:),d.wx];
E(loc,:)=0; A(loc,:)=C;

% computing eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); 
s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

% show the eigenvalues
subplot(1,2,1);
plot(real(s),imag(s),'b.')
xlim([-20,1])
xlabel('real part');ylabel('imaginary part');title('eigenvalues');
grid on; hold on

%{
# Validation

To validate this computation, we compare it to the eigenvalues we get by assuming a domain invariant in the *x* direction, solving mode by mode, like we did in [wave_like.m](). We do the wavelike assumption
$$
u(x,t)=\hat{u}\exp(i\alpha x+st)+C.C.
$$
so the system of equations becomes an eigenvalue problem for eigenvalue $s$
$$
s\left(\begin{array}{l}
\hat{u} \\
\hat{v} \\
\end{array}\right)
=
\begin{pmatrix}
\lambda-1-2\alpha^2 & 1 \\
-\lambda  & \alpha^2-1 \\
\end{pmatrix}
\left(\begin{array}{l}
\hat{u} \\
\hat{v} \\
\end{array}\right)
$$
The eigenvalues are the roots of the caracteristic polynomial of this system
$$
s^2+s(2+3\alpha^2-\lambda)+1+\alpha^2(3-\lambda)+2\alpha^4=0
$$
which we compute using the function *roots*. We plot on top of eachother the eigenvalues from the 1D domain and the theoretical eigenvalues. 

For this we need to account for all the periodicities that are alowed by the grid in 1D. Since we have a domain with periodic boundary conditions, we allow waves of wavelength $L, L/2, L/3, L/4 \dots$, that is wavenumbers $2\pi/L, 4\pi/L, 6\pi/L, \dots$. 

We then show as well the eigenvectors, and we see that they indeed correspond to the periodicities we are talking about.

%}
% validation
for n=[1,2,3,4,5];
    k=n*2*pi/L; % wavenumber
    stheo=roots([1,2+3*k^2-lambda,1+k^2*(3-lambda)+2*k^4]);
    plot(real(stheo),imag(stheo),'ro');
end

% show the eigenvectors
subplot(1,2,2);
co='brkmcrbbrkmcrbbrkmcrb';
numvec=[1:10];
for ind=1:length(numvec);
    num=numvec(ind);
    plot(x,real(U(u,num)),co(ind));
    hold on
end
legend('mode 1','mode 2','mode 3');
xlabel('x');  ylabel('eigenvectors');title('eigenvectors');
grid on; 

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','bruxellator_jerome.png');

%{


![Eigenvalues and eigenvectors](bruxellator_jerome.png)

# Random initial condition

Now we give a feeling for the time evolution of the system by showing the evolution in time of random initial conditions built using the *n* first eigenmodes of the system. Remeber that the initial condition is real, so we should give equal random weights on the eigenmodes that are complex conjugate to each other. For this, we simply remove the eigenmodes whose eigenvalue has a negative imaginary part, then we take the real part of the weighted sum of the other eigenmodes. 
%}


% remove eigenmodes with negative imaginary part
rem=imag(s)<-1e-5; s(rem)=[];U(:,rem)=[];

% loop to repeat
for tre=1:10

    % to build the initial condition
    n=10; % number of eigenmodes in the initial condition
    a=randn(n,1)+i*randn(n,1); % the weights
    tmax=5;
    
    % to set the y scale
    if max(real(s(1:n)))>0 % if unstable
      ymax=max(real(U(u,1:n)*(a.*exp(s(1:n)*tmax))));
    else % if stable
      ymax=max(real(U(u,1:n)*(a.*exp(s(1:n)*0))));
    end
    
% time loop
for t=linspace(0,tmax,100);
    
    % the present state
    q=real(U(u,1:n)*(a.*exp(s(1:n)*t)));
    plot(x,q);
    hold on
    
    % draw each of the eigenmode component
    for gre=1:n
       plot(x,real(U(u,gre)*a(gre)*exp(s(gre)*t)),'r--');
    end
    hold off
    axis([0,L,-ymax,ymax]); grid on
    title(t)
    drawnow;
end
end






