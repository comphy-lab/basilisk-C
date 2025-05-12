%{
# Eigenmodes of the unsteady reaction  diffusion equation

We use our problem  but now we take this equation $$ U_{t} =\mu U_{xx} +  (U/\tau) $$ . This équation is conditionally steady indeed we can get positiv eigenvalue. So let see what happens with the evolution in time when a eigenvalue is positive .
%}

clear all; clf

% parameters
N=100; % number of gridpoints
L=10; % domain length
mu=0.1 ; % diffusion coefficient
 tau=[ 120] 
% differentiation matrices
scale=-2/L;
[x,DM] = chebdif(N,2); 
dx=DM(:,:,1)*scale;    
dxx=DM(:,:,2)*scale^2;    
x=(x-1)/scale; 
Z=zeros(N,N); I=eye(N); 
II=eye(2*N);

% system matrices
E=I;
A=mu*dxx  +1/tau ; 

% boundary conditions
loc=[1,N];
E(loc,:)=0; 
A(loc,:)=I(loc,:);
         

% computing eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); 
s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
control=find(s>0);  
 if norm(control)>0 % if a eigenvalue is positiv the system is unsteady
     disp(tau)
     disp('The problem is unsteady')
 else
    disp('The problem is steady') 
 end
    
%{
# Eigenmodes

For a linear scalar system like
$$
u_t=au
$$
with $u(t)$ a scalar varying in time and $a$ a constant, the explicit evolution of this system is
$$
u(t)=u(0)\exp(at)
$$
where $u(0)$ is the initial value of $u$.

For our systems it is the same, except that they are matrix-systems and the varuable is a vector
$$
q_t=Aq
$$
but this can be put into the same form by diagonalizing the system. To diagonalize means simply change the reference to write the system in the basis of its eigenvectors. In this very special basis, the matrix $A$ is diagonal. This means that in this basis, all the degrees of freedom are independent from each-other, and behave like the scalar example above.


%}




% set(gcf,'paperpositionmode','auto');
% print('-dpng','-r80','diffusion_eigenmodes1.png');


   



%{
![The eigenvalues and eigenvectors](diffusion_eigenmodes1.png)

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


%break
% show the evolution of an initial condition
n=10; % number of eigenmodes in the initial condition

%{
Here we simply build a random initial condition by combining the *n* least stable eigenvectors
%}
a=randn(n,1); % the weights

% time loop
%if tau==0.0001
    figure(2)
for t=linspace(0,20,300);
    
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
end




%{
Here we see in blue the state at the time of the snapshot, and in red the eigenvectors that we have used when building the initial condition, each multiplied with $exp(s_i t)$ where $s_i$ is the associated eigenvalue. This shows that quickly, the evolution of the state of the system tends to the evolution of the least stable eigenmode.

![A snapshot of the time evolution](/sandbox/easystab/stab2014/combustion_instable.png)




%}