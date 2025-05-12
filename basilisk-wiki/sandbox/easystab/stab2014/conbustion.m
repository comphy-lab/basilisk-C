%{
# Eigenmodes of the diffusion equation

This is a second pedagogical introduction for the computation of the eigenmodes of a 1D system. We then use these eigenmodes to see the evolution in time of an initial condition. 


U is the density of the flame that varies with time and in the direction x. mu is the diffusion parameter we chose, which will have an impact on the stability of the system.Tau parameter is positive and  parameterizes the violence of the combustion reaction . Here is the equation system that describes the conbustion systeme :
$$
\begin{array}{l}
U_{t} =\mu U_{xx} -  (1 / \tau) (1-U)U   \\
\end{array}
$$
%}
%{
# Linearization

First, we calculate the linear and stationary form of the system. If the system is stationare,then we find two solution: $$U_{1}=0 \;\; and\;\; U_{2}=1$$





We then disrupt the system by posing :
$$
U = U_b + u \;\;\;\;\;\;\;\;  with\;\; b=1\; or\; 2
$$



We now have the linear form of the conbustion systeme : 
$$
\begin{array}{l}
U_{t} =\mu U_{xx} -  (U/\tau) \;\;\;\; (i) \\
\end{array}
$$
$$
\begin{array}{l}
U_{t} =\mu U_{xx} +  (U/\tau) \;\;\;\; (ii) \\
\end{array}
$$
We choose the steady solution (first)


%}

clear all; clf
% parameters
N=100; % number of gridpoints
L=10; % domain length
mu=0.1 ;%0.1; % diffusion coefficient
for tau=[0.0001 0.001 0.01 0.1  50 120 400  700 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 20000 25050 100000 ]%;10000;%0.1;
% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 
II=eye(2*N);
%{

Let us use [diffusion_eigenmodes.m](/sandbox/easystab/diffusion_eigenmodes.m) to calculate the eigenvectors and eigenmodes of our system. By using the differentiation matrices as we learnt from [diffmat.m](/sandbox/easystab/diffmat.m) and [diffmat_2D.m](/sandbox/easystab/diffmat_2D.m), we can turn the linear system into a matrix system :

$$ 
Eq_{t}= Aq 
$$
$$ 
with\;\;E=I\;\;and\;\;A= \mu dxx-1/ \tau
$$

%}

% system matrices
E=I;
A=mu*dxx-1/tau   ; 
%{
# Boundary conditions

We choose to keep the same Dirichlet Boundary conditions. 

%}
% boundary conditions
loc=[1,N];
E(loc,:)=0; 
A(loc,:)=I(loc,:);
         
%{
Please take a look at [diffusion_eigenmodes.m](/sandbox/easystab/diffusion_eigenmodes.m) if you want to understand the section computing eigenmodes

%}
% computing eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); 
s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
    
%}

%{
# Validation

To validate our code we did a loop that varies tau from 0.0001 to 100000. And we plot the eigenvectors to fallow their evolution. So  the first value of tau represents our problem, then we can see eigenvectors looks like more and more to  eigenvectors of Diffusion problem. Indeed for $$ tau=100000$$ we have the same eigenvectors.
%}
if tau==0.0001
    % show the eigenvalues
figure(1)
subplot(2,1,1);
plot(real(s),imag(s),'b.')
xlim([-5,1])
xlabel('real part');ylabel('imaginary part');title('eigenvalues');
grid on; 

% show the eigenvectors
subplot(2,1,2);
co='brkmc';
for ind=1:3
plot(x,real(U(:,ind)),co(ind),x,imag(U(:,ind)),[co(ind) '--']);
hold on
end
xlabel('x');  ylabel('eigenvectors');title(['eigenvectors tau=',num2str(tau)])  % title('eigenvectors' num2str(tau));
grid on; 

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','diffusion_eigenmodes1.png');

else
    figure(3)



% show the eigenvectors

co='brkmc';
for ind=1:3
plot(x,real(U(:,ind)),co(ind),x,imag(U(:,ind)),[co(ind) '--']);
hold on
end
xlabel('x');  ylabel('eigenvectors');title(['eigenvectors tau=',num2str(tau)])  % title('eigenvectors' num2str(tau));
grid on; 

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','diffusion_eigenmodes1.png');
hold off
%{


# Time evolution
 
Here we use the eigenmodes to show the time evolution of an initial condition. If the initial condition of our system is one of its eigenmodes, we know exactly all the time evolution, it will simply be the eigenvector multiplied by the exponential of the associated eigenvalue time the time. Going further one step, this means that if the initial condition is a combination of the eigenvectors, we get the evolution as sum sum of the eigenvectors weighted by their individual exponential time dependency.

With the initial condition
$$
q(0)=\sum_i \alpha_i q_i
$$
with $\alpha_i$ the amplitude of each eigenmode in the initial condition. The evolution in time is thus
$$
q(t)=\sum_i \alpha_i \exp(s_i t)q_i
$$

%}
end

%break
% show the evolution of an initial condition
n=10; % number of eigenmodes in the initial condition

%{
Here we simply build a random initial condition by combining the *n* least stable eigenvectors
%}
a=randn(n,1); % the weights

% time loop
if tau==0.0001
    figure(5)
for t=linspace(0,40,200);
    
    % the present state
    q=U(:,1:n)*(a.*exp(s(1:n)*t));
    plot(x,q);
    hold on
    
    % draw each of the components
    for gre=1:n
       plot(x,U(:,gre)*a(gre)*exp(s(gre)*t),'r--');
    end
    hold off
    axis([0,L,-2,2]); grid on
    title(t)
    drawnow;
    
    xlabel('x');  ylabel('eigenvectors');title('Marche en temps de la conbustion pour tau=0.0001')
     
end
end
end

%{
# Result

![Above, we can see the last picture produced by the code. We can notice it is the same that the diffusion code ](/combustion_last_tau.png)

![](/combustion1st_tau.png)

![combustion picture of march in time](/combustion_marcheintime.png)

%}