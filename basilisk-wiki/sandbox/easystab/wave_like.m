%{
# The assumption of a wavelike perturbation

This code is a pedagogical one used to show how to trasform a 3D problem into
a linear system in 1D by assuming that the perturbation has a wavelike structure in all the "homogeneous directions".

When we do a slideshow to introduce the instabilities in fluids, we usually show waves. Why is it so that fluid flow should destabilize 
as a wave? It is not a general fact, it is because we have historically studied the stability of fluid flow with what we call "homogeneous" directions.
A direction of space is called homogeneous if the system is invariant by translation in this direction. For instance if you have a fluid
at rest between two infinite horizontal flat plates then it is clear that if you move of a given distance horizontally, ou will see nothing changing: this is invariant by translation in the horizontal direction.
On the other hand if you move vertically, you will find yourself closer to a plate and farther from the othe one. So, this system is not invariant by translatino in the vertical direction.

An other example is the flat Poiseuille flow, the flow of a viscous fluid between two infinite parrallel flat plates, induced by a pressure gradiant in one of the horizontal direction.
The velocity profile is a parabola in the vertical direction. This fow is as well invariant by translation in the two horizontal direction. Now if for instance you consider a square duct, 
then only the axis of the duct is a homogeneous direction.

Since we solve linear equations, and since everything is constant in the homogeneous directions, we can do Fourier transforms in these directions, and we can study the dynamics for each of the wavelength of the Fourier
transform. If we do this for every frequency, then we now everything about the flow. This means that it is enough to consider one by one 
each wavelength. This is why we say that we do the "assumption of wavelike perturbation".

Assume then that we have a flow in 2D with velocity and pressure 
$$
u(x,y,t), v(x,y,t), p(x,y,t).
$$ 
Each of these unkown is a function of three 
coordinates. Now saying that we study waves means that we are interested in the behaviour of the flow in a particular form
$$
u(x,y,t)=\tilde{u}(y,t)\cos(\alpha x+\phi(y,t))
$$
where we have kept in $\tilde{u}$ the dependency with $y$ since $y$ is not a homogeneous direction we cannot assume a wavelike behaviour in general in this direction.
And we have for now kept the dependency with time $t$. In this expression, $\alpha$ is called the wavenumber, it is equal to $2\pi$ divided by the wavelength, 
and $\phi(y,t)$ is a phase shift, and here we have assumed that this phase shift can vary with the vertical $y$ direction and time. By doing this transformation we have greatly
simplified our problem because we have removed one dimension of the unknown.

This formulation with a cosinus is nice, but we can yet do something even simpler and compact using the complex exponential
$$
u(x,y,t)=\tilde{u}(y,t)\exp(i\alpha x +st)+Conjugate(\tilde{u}(y,t)\exp(i\alpha x +st))
$$
where $Conjugate()$ denotes the complex conjugate. By using this notation, we have allowed are unknown $\tilde{u}$ to be a complex function,
and this way we have absorbed the phase shift $\phi$ into it. Indeed, since $\tilde{u}$ is complex, its complex phase can change along $y$,
and since we multiply it to a complex exponential, we will have the phase of the product change in $y$ thus removing the need of an additional variable $\phi$.
We have added the complex conjugate to the complex expression of the wave since the wave is a real thing thus its description should be real. For the $x$ dependency with have $i\alpha$ with $i$ the complex unity and $\alpha$ a real number. This way there is just oscillations in $x$, and no exponential growth, since
$$
\exp(a+ib)=\exp(a)(\cos(b)+i\sin(b)).
$$

Usually we use the notation
$$
u(x,y,t)=\tilde{u}(y,t)\exp(i\alpha x)+C.C.
$$
instead where $C.C.$ means "Complex Conjugate". We always add this complex conjugate to make the wave real. 
A complex number plus its complex conjugates is equal to twice the real part of the number.

Now we see how to write derivatives with this notation. Since nothing depends on $y$ in the exponential, then the $y$ derivative of $u$ 
is simply the $y$ derivative of $\tilde{u}$
$$
u_y=\tilde{u}_y\exp(i\alpha x)+C.C.
$$
and for the $x$ derivative is even simpler because we have assumed a wave in $x$, thus we know in advance the dependency of $u$ with $x$
$$
u_x=i\alpha\tilde{u}\exp(i\alpha x)+C.C.
$$
thus differentiating in $x$ just ammounts to multiplying with the wavenumber $i\alpha$.

Now let's have a look at the time dependency in $u$. Since we will always look at linear systems we will manipulate dynamics looking like
$$
q_t=Aq
$$
where we equate the first derivative in time of our unknown $q$ with a matrix $A$ times the unknown. 
This first order linear system have solutions with an exponential dependency in time in the form
$$
q(t)=q_0\exp(st)+C.C.
$$
where $q_0$ is the initial condition of the evolution and $s$ is a (possibly complex) eigenvalue of the matrix $A$. 
This way we see that it will be in many cases possible to assume as well this simple exponential dependency of our state with time. Here we
have again added the complex conjugate $C.C.$ to allow for a complex initial condition and a complex eigenvalue.

Thus we now can simplify even more the assumed shape of our unknown with an exponential behaviour in $x$ and also in $t$
$$
u(x,y,t)=\hat{u}(y)\exp(i\alpha x+st)+C.C.
$$
And if we wonder how will the time derivative look like, we see that we will have
$$
u_t=s\hat{u}(y)\exp(\alpha x+st)+C.C.
$$
since here is no $t$ dependency in $\hat{u}$, the derivative in time just ammounts to multiplying by $s$. 

In most of the cases that we will consider in easystab, the wavenumber of the wave will be chosen, and $s$ will be an unknown of the system. Thus most of our work will be to do the computations in order for a given $\alpha$ to find out what is the value of $s$. This choice is called the "temporal stability". The other choice where you choose for yourself the temporal properties $s$ and get the spatial properties $\alpha$ is called the "spatial stability". 

If $s$ has a negative real part, then the exponential $\exp(st)$ with be decreasing in amplitude in time, this is a **stable** wave, and if $s$ has a positive 
real part, then $\exp(st)$ wil be growing in amplitude in time, this is an **unstable** wave.


%}

clear all; clf;

% parameters
N=50;      % number of gridpoints
L=1;        % Fluid height in y
rho=1;      % fluid density
mu=1;    % fuid viscosity
alpha=4 ;    % wavenumber in x
g=1;        % gravity

% 1D differentiation matrices
[D,DD,wy,y]=dif1D('cheb',0,L,N,3);
I=eye(N); Z=zeros(N,N);

%{
# System matrices

In the following, we consider a fluid at rest between two infinite horizontal plates. The base flow is thus 0, and if we linearize the Navier-Stokes equations and continuity about this zero base flow we get the system
$$
\begin{array}{l}
\rho u_t=-p_x+\mu\Delta u\\
\rho v_t=-p_y+\mu\Delta v\\
u_x+v_y=0.\\
\end{array}
$$
We have seen that differrentiation with respect to $x$ ammounts to multiplication by $i\alpha$, and differentiation 
with respect to $t$ ammounts to multiplication by $s$, thus we have 
$$
\begin{array}{l}
s\rho \hat{u}=-i\alpha \hat{p}+\mu (-\alpha^2 \hat{u}+\hat{u}_{yy})\\
s\rho \hat{v}=-\hat{p}_y+\mu(-\alpha^2 \hat{v}+\hat{v}_{yy} \\
i\alpha\hat{u}+\hat{v}_y=0\\
\end{array}
$$
And you will now ask: 

> "but what happened with the C.C.?"

Well, the good news is that you are lucky, you don't have to worry about the C.C.

If we now represent the system in matrix form $sEq=Aq$ we have
$$
q=\left(\begin{array}{c}
\hat{u} \\ \hat{v}\\ \hat{p}\\ 
\end{array}\right)
, \quad
E=\left(\begin{array}{cccc}\rho I&0&0\\0&\rho I&0\\0&0&0\end{array}\right)
, \quad
A=\left(\begin{array}{ccc}
\mu(-\alpha^2I+\partial_{yy})&0&-i\alpha I\\
0&\mu(-\alpha^2I+\partial_{yy})&-\partial_y\\
i\alpha I &\partial_y&0
\end{array}\right)
$$

In the following of easystab, we do this in all the codes, we do not build the system matrices with $i\alpha I$, instead, we 
create a new differentiation matrix in $x$: dx=i*alpha*I to simplify the notation. We do the same for the second derivative in $x$: dxx=-alpha^2*I.


%}

% renaming the matrices
dy=D; dyy=DD;
dx=i*alpha*I; dxx=-alpha^2*I;
Delta=dxx+dyy;

% location vectors
u=1:N; v=u+N; p=v+N;

% System matrices
A=[mu*Delta, Z, -dx; ...
   Z, mu*Delta, -dy; ...
   dx, dy, Z];

E=[rho*I, Z, Z; ...
   Z, rho*I, Z; ...
   Z, Z, Z];

%{
# Boundary conditions

We simply have no-penetration at both walls so $v$ vanishes at $y=0$ and $y=L$, and for $u$ with have that its $y$ derivative must be zero (homogeneous Neumann boundary condition). We could have chosen just no-slip, but slip is good because there is a simple analytical solution for the validation.  

The constraint matrix $C$ for the boudanry conditions is thus
$$
C=\left(\begin{array}{ccc}
\partial_y|_0&0&0\\
\partial_y|_L&0&0\\
0&I|_0&0\\
0&I|_L&0\\
\end{array}\right)
$$
such that 
$$
Cq=0
$$
is the compact expression of our four boundary conditions. 
   
   
%}

% boundary conditions
loc=[u(1),u(N),v(1),v(N)];
II=eye(3*N); DD=blkdiag(dy,dy,dy);
C=[DD([u(1),u(N)],:); II([v(1),v(N)],:)];

E(loc,:)=0;  
A(loc,:)=C;

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

% compute eigenmodes
[U,S]=eig(A,E);
s=diag(S);  [t,o]=sort(-real(s)); s=s(o); U=U(:,o);
rem=abs(s)>1000; s(rem)=[]; U(:,rem)=[];

%{
# Validation

We said early in this file that $y$ is not a homogeneous direction so that it is not possible to assume a wavelike behaviour in this direction because the walls make the system not invariant by translation in $y$. But in fact here there is no base flow and by chance the solution we get if we do a wavelike assumption in $y$ can satisfy the boundary conditions if we chose properly the wavelength in $y$. This is for this reason that we have chosen a slip at the walls instead of a no-slip. It is probably also not too difficult to work out the analytical solution with no-slip, but the derivation below is quite interesting on its own, to open the view on wavelike assumptions.
This is what we do here to get some model for the viscous attenuation
of oscillations of this flow.

We write 
$$
u(x,y,t)=\hat{u}\exp(i\alpha x+i\beta y+st)+C.C.
$$
where now $\hat{u}$ is just a complex scalar value. We will have to chose the wavenumber $\beta$ 
in $y$ such that we have half a wavelength (the wavelength is twice the distance between the walls), 
or one full wavelength, of three fourth of a wavelength and so on, so that the 
wave can satisfy the no-penetration at the walls (so that the walls are at nodes of the wave in $y$, that is, $v$
is zero at $y=0$ and $y=L$). 

If we insert this assumption into the equations, we get rid of all the derivatives
$$
\begin{array}{l}
s\rho \hat{u}=-i\alpha \hat{p}-\mu(\alpha^2+\beta^2) \hat{u}\\
s\rho \hat{v}=-i\beta \hat{p}-\mu(\alpha^2+\beta^2) \hat{v}\\
i\alpha\hat{u}+i\beta \hat{v}=0\\
\end{array}
$$
And we call $k^2=\alpha^2+\beta^2$. We remove the pressure by multiplying the $u$ equation with $\beta$
and the $v$ line with $\alpha$, and we substract the $v$ one to the $u$ one. We then use the continuity equation to have the 
relation between $u$ and $v$
$$
v=-\frac{\alpha}{\beta}u
$$
so that at the end we get
$$
s=-\frac{\mu}{\rho} k^2
$$
which is the theoretical result that we compare below to our calculation. We run a loop on the increasing
values of the $y$ wavenumber $\beta$ and compute the decay rate for an array of values of the $x$ wavenumber $\alpha$. This equation tells us first that all modes are stable. The eigenvalue $s$ is purely real, and is negative, so all modes decay in time. If the dynamic viscosity $\mu/\rho$ is larger, the modes decay faster. The wavenumber also play a role in the decay rate of the modes, when $\alpha$ of $\beta$ is large, it means that the wavelength is small, and it is clear now that short wavelenght are more touched by the damping effect of viscosity so the decay rate will be fast for short wavelengthes.

it is interesting to see that the two boundary conditions at the walls induce a *quantization* of the possible wavelengthes in the $y$ direction, so we consider some given quantized values of $\beta$. This quantization corresponds exactly to the quantization of the eigenmodes of the numerical system. On te graph we show the theory and the numerics for the four least stable modes.


We see that the fit is nice.
%}

% validation
subplot(2,2,1);
alphavec=linspace(0,10,100);
betavec=2*pi./[2*L,L,2/3*L,1/2*L];
for ind=1:length(betavec)
stheo=-mu*(alphavec.^2+betavec(ind).^2);
plot(alphavec,stheo,'r-'); hold on
end

plot(alpha,s(1:4),'b.');
xlabel('alpha');ylabel('exponential growth rate'); title('validation')
grid on

%{
# Velocity field of the eigenvectors

At the start of this code above, we spent a long time to introduce the idea of wavelike perturbation, in order to reduce the size and complexity for our calculations. Now that all the calculations are done and that we have validated them, it is time to do the work of simplification backward and recover the real velocity field corresponding to the least stable modes. The eigenvectors are the vectors
$$
q_i=\left(
\begin{array}{c}
\hat{u}(y) \\ \hat{v}(y) \\ \hat{p}(y)
\end{array}\right)
$$
so to get the velocity field, we just to multiply this with the exponential in $x$, plus the complex conjugate. For the theory it is better to see it as C.C. and for the computations it is more straightforward to take twice the real part. For other examples of this, please see [free_surface_gravity_particles.m](). 

%}

% show the velocity field of the first three eigenvectors
Lx=2*pi/alpha;
x=linspace(0,Lx,20);
for ind=1:3 
    subplot(2,2,1+ind);
    q=U(:,ind);
    qphys=2*real(q*exp(i*alpha*x));
    quiver(x,y(1:2:N),qphys(u(1:2:N),:),qphys(v(1:2:N),:));
    axis equal; axis([0,Lx,0,L]);xlabel('x');ylabel('y'); 
    title(['mode' num2str(ind)]);
end

set(gcf,'paperpositionmode','auto');
print('-dpng','-r100','wave_like.png');

%{

![The figure](wave_like.png)

On the figure, we show the velocity field for the three least stable eigenmodes. All off course have the same wavelength $2\pi/\alpha$ in $x$, but we see that as we go to more decaying modes, they have more and more oscillations in the wall-normal direction. The first mode corresponds to two counter-rotating vortices occupying the whole height between the walls. The second mode is four counter-rotating vortices, since now from $y=0$ to $L$, space is made for two vortices, and so on for the following modes.


# Exercices/Contributions

* Please do the same thing as in here for the brusselator:
$$
\begin{array}{l}
u_t=1-(\lambda+1)u+2u_{xx}+u^2v\\
v_t=\lambda u+v_{xx}-u^2v\\
\end{array}
$$
where $u$ and $v$ are two variables depending on $x$ and $t$. Start by finding the uniform steady state (constant in time and space). Then linearize about this steady state, then compute the eigenvalues. Compare with the theoretical analysis where you assume an infinite domain. Take periodic boundary conditions to be able to compare your code and you theory


* Please do the same as above, but with the reaction-diffusion equation
$$
u_t=\mu u_{xx}+\frac{1}{\tau}u(1-u)
$$
where $\mu$ and $\tau$ are two positive constants.

* Please do the same for the excited wave equation
$$
u_{tt}+(1-r)u=u_{xx}
$$
with $r$ a parameter.

* Please do the same analysis for a 1D model of a pipe conveying fluid -->[hose_fluid_nonviscous.m](/sandbox/easystab/stab2014/hose_fluid_nonviscous.m)
$$
f_{xxxx}+(1-\xi)f_{xx}+2\sqrt{\beta}f_{xt}+f_{tt}=0
$$
* Please do the same analysis for the Saint-Venant equation -->[st-venant.m]()
$$
\begin{array}{l}
u_t=-g\eta_x-bu\\
\eta_t=-(Hu)_x\\
\end{array}
$$
where $H$ is the depth of the sea and $\eta$ is the perturbation to the sea depth and $u$ is the associated fluid velocity.

* Please do the samae analysis for the Saint-Venant equation with Coriolis force
$$
\begin{array}{l}
u_t-fv=-g\eta_x-bu,\\
v_t+fu=-g\eta_y-bv,\\
\eta_t=-Hu_x-Hv_y\\
\end{array}
$$




%}
