clear all; 
clf

%{ 
# Study of the stability of a simple non constant differential equation   

It's well known in first years of university that if ytou take the classical differential equation :
$$
\frac{1}{\omega_0^2} \frac{\delta^2 U}{\delta t^2}+ \tau \frac{\delta U}{\delta t}+ U =
f(t)
$$

if the coefficient of the left side are all from the same sign, then the
system is stable. If they are not, he system is unstable.

I started wondering, what happen if the sign of the first derivate was
sometime positive and sometimes negative ?
As eigeinmode is a very powerful tool to analyse stability, and a
non-linear differential equation keep the same eigeinmode during his
evolution, I tried to analyse this system with eigeinmodes to determine his
stability. 

The simplest function for me was to multiply the damp factor by a sinus of
the time, with a certain frequency. After few computations, I realised that
there was a range of value where the system was unstable. 

I plotted the eigeinvalues in function of the damping frequency, I tried to
understand the dynamics of the system.

Actually, I cannot understand it, we see that some eigeinmodes do
bifurcations, must I have no ways to understand properly the diagram
traced.
Moreover, the system is 0D (one point with an evolution in time). I have in
fact no idea of what does the Eigenvectors represent, and since we already
compute the evolution in time, I also have no idea what are the
interpretation of eigenvalues...

Some theory would be usefull to develop, I'm actually looking at this kind
of system and I will update this programm as soon as I know more about

#How the code work ?

This is a mix between  and   

%}

%%%%%%%%%%%%%%%%%%%%%%PARAMETERS OF THE MODEL%%%%%%%%%%%%%%%%%%%
% numerical parameters
Tmax=10; % domain length
N=100; % number of points
t=linspace(0,Tmax,N)'; %DON'T CHANGE THIS

% physical parameters
omega0=1;
tau=1;

%vector of the excitation, choose an expression for b
b=zeros(N,1);  %No sollicitation

amax=5;
%Nonconstant vectors
for a=0.0:0.01:amax
mattau=tau*diag(sin(a*t));
matomega0=omega0*diag(1);

%Boundary condition
type1= 1 ; %1 dirichlet, 2 Neumann, 3 f''
place1= 0 ; % in the physical dimension (we calculate the indice after)
val1=1;
type2= 2 ; % 1 dirichlet, 2 Neumann, 3 f''
place2= 0 ; % in the physical dimension (we calculate the indice after)
val2=0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=t(2)-t(1); % the grid size
indice1=floor(place1/h)+1;
indice2=floor(place2/h)+1;

    % first derivative
    D=zeros(N,N);
    D(1,1:3)=[-3/2, 2, -1/2]/h;
    for ind=2:N-1
        D(ind,ind-1:ind+1)=[-1/2, 0, 1/2]/h;
    end
    D(end,end-2:end)=[1/2, -2, 3/2]/h;

    % second derivative
    DD=zeros(N,N);
    DD(1,1:3)=[1, -2, 1]/h^2;
    for ind=2:N-1
        DD(ind,ind-1:ind+1)=[1, -2, 1]/h^2;
    end
    DD(end,end-2:end)=[1, -2, 1]/h^2;

%%%%%%%%%%%%%%%%%%%%%%%Creation of the matrix%%%%%%%%%%%%%%%%%%%%%%%%%%
I=eye(N);
A=I+mattau*D+(matomega0^(-2))*DD;


%%%%%%%%%%%%%%%%%%%%%%%%%%initial conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if type1==1, A(1,:)=I(indice1,:) ; end
 if type1==2, A(1,:)=D(indice1,:) ; end
 if type1==3, A(1,:)=DD(indice1,:); end
 
 if type2==1, A(N,:)=I(indice2,:) ; end
 if type2==2, A(N,:)=D(indice2,:) ; end
 if type2==3, A(N,:)=DD(indice2,:); end
 
 %{Attribution of the values}%
 b([1,N])=[val1,val2];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%Resolution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U=A\b;

%%%%%%%%%%%%%%%%%%%%%%%%%%EIGENVALUES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V,S]=eig(A,eye(N)); %We compute all the eigenvectors of the system
s=diag(S);          %We extract the eigenvalues 
rem=abs(s)>2;       %We remove too big or too weak eigeinvectors
s(rem)=[];
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TRACE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(2,1,1);
if max(U./t)>1.1, plot(a,real(s),'r.'); end; %We plot unstable modes in another color
if max(U)<1.1, plot(a,real(s),'b.'); end;
xlabel('a');ylabel('values');title('eigenvalues');
grid on; 
hold on 

subplot(2,1,2);
if mod(a,0.05)==0,
if max(U)>1.1, plot(t,U,'r'); end;
if max(U)<1.1, plot(t,U,'b'); end;
axis([0,Tmax,-5,5])
hold on
end
end
hold off
  set(gcf,'paperpositionmode','auto')
    print('-dpng','-r90',sprintf('unstable.png'));

%{
# Figures !

Here is the figure the program draw, in blue stable values of a and in red unstable values. Honestly, I don't understand my results.
![alt text](/unstable.png)

%}
