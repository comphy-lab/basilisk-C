%{
# The shape of a meniscus

This code is a simple first test to computing the shape of a liquid meniscus for a bath of liquid in contact with a solid vertical wall. The liquid is in contact with the wall with a contact angle that we call here $theta$. We can build an equation to describe how the fluid surface goes to this contact angle at the wall to a flat interface far from the wall. The set-up is described on the figure.


% here the figure

This system has two physical effects: first the gravity induces inside of the liquid a hydrostatic pressure distribution: pressure increases linearly with depth below the interface. The second effect is due to the surface tension at the liquid-air interface which creates a pressure jump through a curved interface, described by the Laplace law:
$$
\Delta p=\sigma/R
$$
where $R$ is the radius of curvature of the interface and $\sigma$ is the surface tension. The pressure is lower on the concave side of the surface. To write a differential equation for the interface, we will simply say that the pressure just below the curved interface must be equal to the hydrostatic pressure. For this we will need to describe the radius of curvature of the interface based on the function $h(x)$ that describe where the interface is located. The curvature is related to the second derivative $h''$, but with a correction term if the interface is not flat:
$$
1/R=\frac{h''}{(1+h'^2)^{3/2}}.
$$ 
We can now write the equation for $h$:
$$
\sigma \frac{h''}{(1+h'^2)^{3/2}}=\rho g h.
$$
This is the nonlinear equation that we will need to solve numerically to get the shape of the meniscus. But before to solve this nonlinear equation, let's start with a simpler case, for instance assuming that the contact angle is close to $\pi/2$, this is the low-slope approximation for which the curvature can be approximated by
$$
1/R\approx h''
$$ 
since now $h'^2$ is small. We get the linear equation
$$
\sigma h''=\rho g h.
$$
The boundary conditions are that $h$ must be flat far from the wall, let's say $h(L)=0$ where $L$ is the domain length, and that the slope must obey the contact angle at the wall. The solutino is then
$$
h(x)=h_0 \exp(x/\sqrt{\sigma/\rho g})
$$ 
where $h_0$ is the meniscus height. Since the contact angle is $\theta$, we must have the relation from the slope at the wall
$$
h_0=\sqrt{2/(\rho g)\sigma(1-\sin(\theta)}
$$
his gives us a weak-slope solution that we can use for two things: first as an initial guess for the Newton iterations, and then to compare to our nonlinear solutino at the end of the convergence.
%}

clear all; clf;

% parameters
L=10;            % domain length
N=100;           % number of grid points
rho=1;           % density
sig=1;           % surface tension
g=1;             % gravity
theta=pi/2*0.2;  % contact angle

% Differentiation matrices and integration
[D,DD,wx,x]=dif1D('cheb',0,L,N,3);
Z=zeros(N,N); I=eye(N); 

%initial guess
h0=sqrt(sig/(rho*g))/tan(theta);
initguess=h0*exp(-x/sqrt(sig/(rho*g)));
q=initguess;

% Newton iterations
quit=0;count=0;
while ~quit

%{
Now we are inside the loop of the Newton iterations. The first thing to do is to extract from the vector q the different variables that we need for computing the nonlinear function. here it is quite simple because we just have one single variable, $h$ that describes the height of the interface, but in other codes, we may have more variables, like for instance when we compute the meniscus with allowing overturning of the interface, where we use a parametric descriptino of the interface $x(s), y(s), \theta(s)$ where $s$ is arclength along the interface, starting from the wall. 

here we as well compute and store the derivatives of $h$ for later use in the nonlinear function $f$.
%}

    % the present solution and its derivatives
    h=q; hp=D*h;  hpp=DD*h; a=1+hp.^2;
    
    % nonlinear function
    f=rho*g*h-sig*hpp.*a.^-1.5; 

%{
# Computing the Jacobian
For the Newton iterations, we need to build the jacobian of the function $f$, it is a square matrix stored in the variable A. The easiest way to obtain the expression for A is to do a perturbation of $f$
$$
f(h+\delta)\approx f(h)+A\delta
$$
where $\delta$ is a function just like $h$ but with a low amplitude such that we will be able to neglect the nonlinear terms and A is a linear operator: the Jacobian (the matrix equivalent of the derivative for a scalar function). here we have
$$
f= \rho g h -\sigma\frac{h''}{(1+h'^2)^{3/2}}
$$ 
thus
$$
f(h+\delta)=\rho g (h+\delta) -\sigma\frac{(h+\delta)''}{(1+(h+\delta)'^2)^{3/2}}
$$ 
and we use Taylor expansion and neglect nonlinear terms to get the expression
$$
f(h+\delta)=f(h)+\rho g \delta-\sigma(-\frac{3}{2}\frac{h'h''}{(1+h'^2)^{5/2}}\delta'+\frac{1}{(1+h'^2)^{3/2}}\delta'')
$$
thus we have the Jacobian
$$
A=\rho g I-\sigma(-\frac{3h'h''}{(1+h'^2)^{5/2}}D+\frac{1}{(1+h'^2)^{3/2}}D^2)
$$
where we have denote $I$ the identity and $D$ and $D^2$ the first and second derivative opertors. 
%}

    % analytical jacobian
    A=rho*g*I-sig*(diag(-3*hp.*hpp.*a.^-2.5)*D+diag(a.^-1.5)*DD);

%{
# Boundary conditions
Now the boundary conditions. Since $f$ is a second order function, we need to impose two boundary conditions. We replace the first and the last equations of $f$ by equations that impose the boundary conditions
$$
h'(x=0)=-\frac{1}{\tan(\theta}
$$
and $h(x=L)=0$. Since we have replaced these two equations in $f$, we must as well replace the first and last equations in $A$ with the jacobian of these equations. Since we have here linear noundary conditions, their Jacobian is easy to find.
%}

    % Boundary conditions
    loc=[1 N];
    C=[D(1,:); I(N,:)];
    f(loc)=C*q-[-1/tan(theta); 0];
    A(loc,:)=C;
    
%{
Here we plot the current solution of the problem, and we compare it to two things: the solution to the linear system and the meniscus height for the nonlinear system.
%}

    % convergence test
    res=norm(f);
    plot(x,h,'b-',x,initguess,'r--',0,sqrt(2*sig*(1-sin(theta))/(rho*g)),'r.',0,h(1),'bo'); 
    legend('h','linear solution','meniscus height');xlim([-0.01,L]);drawnow;
    disp([num2str(count) '  ' num2str(res)]);
    if count>50|res>1e5; disp('no convergence'); break; end
    if res<1e-5; quit=1; disp('converged'); continue; end

%{
# The Newton step
And now is the core of the Newton iteration, where we find the new guess for sol (that is: $f$), by solving a linear systme involving the jacobian. The idea is the following: we are over-optimistic to believe that the linear approximation of $f$ is good, so that the true solution of the system should be $h+\delta$ where the jump to the solution $\delta$ is yet unknown and satisfy
$$
f(h+\delta)=f(h)+A\delta=0
$$
thus we should have
$$
\delta=-A^{-1} f(h)
$$
and thus the new optimistic solution should be $h+\delta$. off course, it is not good to compute the inverse of A and then multiply with $h$, instead we do not inverse $A$ but solve the linear system using the backslash operator.
%}

    % Newton step
    q=q-A\f;
    count=count+1;
end

%{
# The figure
Here we see the result of the computation

![Shape of the nonlinear meniscus](/meniscus.png)

# Exercices/Contributions

* Please draw a graph that shows how the height of the meniscus depends upon the contact angle
* Please change the boundary condition from the contact angle to the meniscus height (and thus let the contact angle adjust accordingly)
* Please modify this code to compute the shape of a drop on a flat substrate
* Please modify this code to compute the shape of a drop on an (slightly) inclined substrate
* Please modify this code to compute the capillary rise between two vertical walls
* Please modify this code to compute the capillary rise inside a capillary tube (you will need to account also for the second component of the curvature)


%}
