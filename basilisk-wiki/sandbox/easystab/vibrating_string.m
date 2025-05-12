%{
# March in time of the vibrating string

In this code, we simulate the evolution in time of an initial condition of a model of a guitar string, simply an elastic string tensed under two fixed attachement points. What is interesting with this system from the technical point of view is that the dynamic equation has a second derivative in time and we show here how to transform this into a larger system with just a single time derivative.

We also show hee a very simple way to do marching in time, by simply building a matrix *M* which if we do a matrix-vector multiplication with a given state of the string, we get the state of the string a little time later (equal to the chosen time step).

%}

clear all; clf

% parameters
N=50; % number of gridpoints
L=2*pi; % domain length
c=1; % wave velocity
dt=0.5; % time step
tmax=20; % final time

% 1D differentiation matrices
[dx,dxx,wx,x]=dif1D('cheb',0,L,N,3);
Z=zeros(N,N); I=eye(N); II=eye(2*N);

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
E=[I,Z; Z,I];
A=[Z,c^2*dxx; I, Z]; 

% locations in the state
v=1:N;
f=N+1:2*N;

%{
# Boundary conditions
The original system has one variable and two time derivatives, so we should apply two boundary conditions. We will simply say that the position of the string at its both ends should be zero: two attachment points. This is a homogeneous Dirichlet condition applied on the first and last point of the state vector for the state position.
%}

% boundary conditions
loc=[f(1),f(N)];
E(loc,:)=0; 
A(loc,:)=II(loc,:);
         
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
We need to describe the state of the system at initial time. We say here that the string is initially deformed as a sinus (this satisfies the boundary conditions), and that the velocity is zero. From that state we then perform a loop to repeatedly advance the state of a short forward step in time. We use *drawnow* to show the evolution of the simulation as a movie when running the code. We store for validation the string position at the midle of the domain, to do this without worrying about the the grid points are, we interpolate $f$ with the function *interp1*. 
%}

% initial condition
q=[zeros(N,1); sin(pi*x/L)]; 

% marching loop
tvec=dt:dt:tmax; 
Nt=length(tvec);
e=zeros(Nt,1);

for ind=1:Nt    
    q=M*q; % one step forward
    e(ind)=interp1(x,q(f),L/2); % store center point
    
    % plotting
    subplot(1,2,1);
    plot(x,q(f),'b',x,q(v),'r--');    
    axis([0,L,-2,2])    
    drawnow
end
legend('position','velocity'); title('Vibrating string')
xlabel('x'); ylabel('f');

%{
# Validation
Now we show that the code does what we think it does, and also give ourselves the means to tell how precise and robust it is. The period of oscillation of this wave should be the time that it takes for a wave to travel at speed $c$ along its wavelength. Here we have chosen an initial condition whose wave length is twice the domain length $2L$ (this is "mode $1/2$"), thus the period of oscillations is 
$$
T=2L/c
$$
thus the central point of the string should evolve like
$$
cos(2\pi t/T)
$$
Here we build a second time vector, to plot the theoretical solution with large resolution in time to have a smooth plot even when the time step of the computation is large.
%}

% time evolution of central point
subplot(1,2,2);
Ttheo=2*L/c;
tt=linspace(0,tmax,500);
etheo=cos(2*pi*tt/Ttheo);
plot(tvec,e,'b.-',tt,etheo,'r-');
legend('numerical','theory');
xlabel('time'); ylabel('f(L/2)');


%{
Here is the figure at time *tmax*:
![Validation](vibrating_string.png)

# Exercises/contributions

* Please use as initial condition the mode 1 (one wavelength in one domain length) and adapt the theoretical solution. -->[vibrating string mode1](stab2014/vibrating_string_mode1.m)
* Please do a convergence study of the time step --> [vibrating string time step convergence](stab2014/vibrating_string_time_convergence)
* Please replace the Crank-Nicolson marching scheme by forward Euler and also Backward Euler and compare. (Forward Euler ammounts to approximating the time integral of the right hand side by multiplying $Aq(t)$ by the time step, this gives an explicit scheme; and backward Euler is approximating it with $Aq(t+dt)dt$, this gives an implicit scheme (for which we need to invert a matrix).
--> [vibrating_string_forward_backward_euler.m](stab2014/vibrating_string_forward_backward_euler.m)
* Please write a time marching scheme that does not need to replace the second time derivative with a single time derivative.
* Please use Runge-Kutta for the time marching. --> [vibrating_string_RungeKutta.m](stab2014/vibrating_string_RungeKutta.m)
* Please do the time marching with the same scheme as here, but without computing and storing the inverse of a matrix.
* Please add some dissipation in the system (for instance the friction with the air) and show how this affects the time evolution.
* Please change the boundary conditions so that our system is a model of the water waves in a 1D swimming pool (the ends of the state are not attached but the first space derivative is zero) [water_waves_1D.m]()
* Please change the system so that one of the ends of the string is not fixed but is attached to a spring (something like the board of a guittar). You can also put some dissipation for this spring. 
-->[vibrating_string_spring_dissipation.m](stab2014/vibrating_string_spring_dissipation.m)
* Please change the system to model the Larsen effect leading to the Jimmy Hendrix instability (the guitar vibration is fed back through a strong amplifier and in turn re-exites the mechanical vibration).
* Please change the boundary condition to model a string that is forced through one of the attachment points. For instance start with a string at rest and then let one of the attachment points move up and down and induce waves on the string. Now this is now longer a homogenous boundary conditions at the ends.-->[vibrating string forced.m]()
* Please change the system to model a string whose mass density is not everywhere the same (the density is a non-constant parameter).
* Please change the system to model two tensed strings connected by a heavy knot (here two difficulties: model the fact that the knot has a weight, and connect two strings). --> [2 strings one mass](stab2014/two_vibrating_strings_one_mass.m)
* Please modify the model to account for a string whose properties can change in time (like the wave speed $c$ because someone puts more tension while the string is vibrating (this is for instance the tremolo effect for a guitar). For this you may need to account in the time marching scheme for the fact that the matrices are not stationnary.
* Please change the code to march in time the wave equation in 2D
$$
f_{tt}=c^2 (f_{xx}+f_{yy})
$$
(first on a square then on a disc).
* Please add some suggestions on new exercices/contributions for this code.
* Please do the marching in time of the advection equation $f_t+Uf_x=0$ with just one boundary condition at the left end (if U is positive)
  [advection_1D.m]()
* Please do the marching in time of the diffusion equation $f_t=\mu f_{xx}$ --> [stab2014/Diffusion.m]()
* Please do the marching in time of the advection-diffusion equation $f_t+Uf_x=\mu f_{xx}$ --> [stab2014/adv-diff.m]()
 
**For all these, please show that the code does what you expect: validate!**
%}
