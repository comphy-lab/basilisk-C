%{
# The brusselator - ongoing project

U and V are the concentrations of two chemical componants that varies with time and in the direction x. lambda is a chemical paramter we chose, which will have an impact on the stability of the system. Here is the equation system that describes the brusselator :

$$
U_{t} = 1 - (\lambda + 1)U + 2U_{xx} + U²V          (i)
$$
$$
V_{t} = \lambda*U + V_{xx} + U²V          (ii)
$$

First, we calculate the linear and stationary form of the system. If the system is stationary, then $$ U_{t} = V_{t} = 0 $$ . Then from (i) we have, 
$$ 
V = \lambda \ U 
$$ 
which we replace in (ii), which gives leaves us with : 
$$ 
1 = (\lambda + 1)U - \lambda*U 
$$ 
We then disrupt the system by posing :
$$
U = U_b + u
$$
$$
V = V_b + v
$$

We now have the linear form of the brusselator : 
$$
U_t = 2U_{xx} + (\lambda - 1)U + V
$$
$$
V_t = -\lambda*U + V_{xx} - V
$$

%}

clear all; clf

% parameters
N=100; % number of gridpoints
L=4*pi; % domain length
lambda = 2; % chemical parameter

% differentiation matrices
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,L,N);
Z=zeros(N,N); I=eye(N); II=eye(2*N);

% system matrices
E=II;
A=[(lambda-1)*I+2*d.xx, I; ...
    -lambda*I, d.xx-I];

% boundary conditions
u=1:N; v=N+1:2*N;
loc=[u(1),u(N),v(1),v(N)];
E(loc,:)=0;
A(loc,:)=II(loc,:);

%{
# Eigenmodes of the brusselator (lineaer form)

Let us use [diffusion_eigenmodes.m](/sandbox/easystab/diffusion_eigenmodes.m) to calculate the eigenvectors and eigenmodes of our system. By using the differentiation matrices as we learnt from [diffmat.m](/sandbox/easystab/diffmat.m) and [diffmat_2D.m](/sandbox/easystab/diffmat_2D.m), we can turn the linear system into a matrix system :

$$ 
Eq_{t}= Aq
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
(\lambda-1)I+2D_{xx} & I \\
-\lambda I & D_{xx}-I \\
\end{pmatrix}
\left(\begin{array}{l}
U \\
V \\
\end{array}\right)
$$
%}

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

% validation
for n=[1,2,3,4,5];
    k=n*pi/L;
    stheo=roots([1,2+3*k^2-lambda,1+k^2*(3-lambda)+2*k^4]);
    plot(real(stheo),imag(stheo),'ro');
end

% show the eigenvectors
subplot(1,2,2);
co='brkmcrb';
numvec=[1,3,5,7];
for ind=1:length(numvec);
    num=numvec(ind);
    plot(x,real(U(1:N,num)),co(ind));
    hold on
end
legend('mode 1','mode 2','mode 3');
xlabel('x');  ylabel('eigenvectors');title('eigenvectors');
grid on; 

set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','bruxellator_jerome.png');

%{

![](brusselator_jerome.png)

%}