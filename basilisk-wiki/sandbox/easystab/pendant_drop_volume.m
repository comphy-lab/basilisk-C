%{
# The shape of a pendant drop and the critical volume for dripping

The aim of this code is to compute the shape of a pendant drop situated
at the lower end of a vertical pipe, assumed plane. With the use of a 
continuation procedure, the bifurcation diagram for increasing volume 
of the droplet is also comptued.

Similarly as in the problem of the shape of the [meniscus.m](http://basilisk.fr/sandbox/easystab/meniscus.m), 
the gravity induces a hydrostatic pressure distribution and the surface 
tension imposes a pressure jump at the interface because of its 
curvature. Furthermore, a pressure corresponds to the volume of the droplet.

The equation for the height of the interface $h$ reads:
$$
\sigma \frac{h^{''}}{(1+h^{'^2})^{3/2}}=\rho g h - p \, .
$$

For this problem so as for others where overturning phenomena are 
present, it is convenient to consider a parametric description with the
interface being located at $x(s)$ and $y(s)$ in a cartesian coordinate 
system with $dx^2+dy^2=l^2ds^2$ and $s\in[0,1]$, $l$ being the length of 
the interface. The governing equation becomes then:
$$
\sigma \frac{x^{''}y^{'}-y^{''}x^{'}}{(x^{'^2}+y^{'^2})^{3/2}}=\rho g y - p \, .
$$
An additional equation is the conservation of the volume:
$$
\int_0^1 y x^' ds = V
$$
and because of the metric we also have:
$$
x^{'^2}+y^{'^2}=l^2
$$
%}

clear all; clf;

% parameters
L=1;            % domain length
N=200;          % number of points
rho=1;          % density
sig=1;          % surface tension
g=-1;           % gravity
hi = 0.1;       % thickness of uniform film
V=hi*L;         % volume of uniform film
delta=0.2;      % continuation length

% differentiation and integration
scale=-2/L;
[s,DM] = chebdif(N,2);
D=DM(:,:,1)*scale;
DD=DM(:,:,2)*scale^2;
s=(s-1)/scale;
Z=zeros(N,N); I=eye(N);
INT=([diff(s)',0]+[0,diff(s)'])/2;

%{
# Boundary conditions
        
Since we have a second order differential equation and two first order equations we
need four boundary conditions. 
We impose Dirichlet boundary conditions, representing the fixed
positions of the droplet at the pipe wall. 
Note that one condition will be used as additional equation in the system.
%}

% Boundary conditions
x1=0; xn=L;
y1=0; yn=0;

%{
As a initial guess for the Newton iteration it is sufficient to consider a uniform film of 
thickness $h_i=V/L$ with $V$ being the volume and $L$ the width of the 
pipe since we start with a small volume. The initial pressure is $p_i = V/L \rho g$.
%}

%initial guess
Pinit=V/L*rho*g;
initguess_x=s*L;
initguess_y=hi*ones(N,1);
linit=sum(sqrt(diff(initguess_x).^2+diff(initguess_y).^2));

%{
The solution vector contains the coordiantes of the interface, the
pressure, the length of the interface and the volume. The volume of the droplet is the
homotopy parameter, which is varied to obtain the bifurcation diagram.
Initially, the direction for the continuation is just along the volume.
%}

sol=[initguess_x;initguess_y;Pinit;linit;V];
dir=[zeros(N,1);zeros(N,1);0;0;1];   % initial direction

%{
In the next, the first outer loop is for the continuation over the volume $V$. 
The inner loop is instead needed to find the steady-state solution
with a Newton method.
%}

disp('Continuation loop')
ind=1;quitcon=0;
while ~quitcon
    solprev=sol;
    
    %{
    # Keller's pseudo-arclength continuation
    
    To construct the bifurcation diagram we employ the property
    that a steady state cannot appear or disappear.
    There are different procedures for the continuation. The Keller's
    pseudo-arclength continuation is used in this code. This method allows 
    the direction of the variation of the control
    parameter to change sign and permits therefore to continuate the
    bifurcation diagram after the fold.
    Once the non-linear solution for a given parameter set is found, a guess 
    for the new one is made and a Netwon iteration is performed to find it. 
    The future value for the volume is not imposed by found again by the
    Newton method and is indeed the last component of the residual function
    $f$. The new solution (V,H) has to satisfy:
    $$
    \mathbf{t} \cdot ((V,H)_{new}^T-(V,H)_{old}^T) = \delta \, ,
    $$
    where $\mathbf{t}$ is the direction of the continuation and
    $\delta$ the continuation length.
    For more informations see the [Lecture Notes on Numerical Analysis of Nonlinear Equations, by Eusebius Doedel.](http://indy.cs.concordia.ca/auto/notes.pdf).
    %}
    
    sol=sol+dir*delta; % new prediction of solution
    
    % Newton iterations
    disp('Newton loop')
    quit=0;count=0;
    while ~quit

        %{
        # Newton loop

        Now we are in the Newtoon loop to find the shape of the droplet with a specific volume (still unknow). See
        [meniscus.m](http://basilisk.fr/sandbox/easystab/meniscus.m) for
        more details.
        The main difference is that here the function $f$ of which the root
        has to be found is is a vector quantity whose components are the 
        governing equation
        $$
        \rho g y - p - \sigma
        \frac{x^{''}y^{'}-y^{''}x^{'}}{(x^{'^2}+y^{'^2})^{3/2}} \, ,
        $$
        the metric equation
        $$
        \quad x^{'^2}+y^{'^2}-l^2 \, ,
        $$
        the constraint for the volume
        $$
        \int_0^1 y x^' ds - V \, ,
        $$
        the equation for one boundary 
        $$
        x(0)=0 \, ,
        $$
        condition and the additional equation
        for the continuation loop:
        $$
        \mathbf{t} \cdot ((V,H)_{new}^T-(V,H)_{old}^T) - \delta \, .
        $$
        It is therefore a function of the interface coordinates $x$, $y$, 
        the pressure $p$ and also of the length of the former, $l$, via the 
        curvilinear parameter.

        *Notation*: the underscript `p` indicate a derivative w.r.t. `s`
        %}
        
        % the present solution and its derivatives
        x=sol(1:N); y=sol(N+1:2*N); P=sol(2*N+1); l=sol(2*N+2); V=sol(end); xp=D*x;  xpp=DD*x; yp=D*y;  ypp=DD*y; a=xp.^2+yp.^2;

        % nonlinear function
        f=[-rho*g*y.*a.^1.5-sig.*(xpp.*yp-ypp.*xp)+P*a.^1.5; ...
        a-l^2; ...
        INT*(y.*xp)-V; ...
        x(1)-x1; ...
        dir'*(sol-solprev)-delta];

        %{
        # Computing the Jacobian
        
        Since $f$ is here a vector quantity, the Jacobian is a bloc matrix, 
        whose blocs
        are the discretized Jacobian matrices of the different components
        w.r.t. the variables. They are found by perturbing all variables in the
        considered equation with small parameters and neglecting the nonlinear 
        terms, in the same way as what described in [meniscus.m](http://basilisk.fr/sandbox/easystab/meniscus.m).
        %}

        % analytical jacobian
        A=[ 3*P*diag(a.^0.5.*xp)*D-3*rho*g*diag(y.*xp.*a.^0.5)*D-sig*diag(yp)*DD+sig*diag(ypp)*D,...
        3*P*diag(a.^0.5.*yp)*D-3*rho*g*diag(y.*yp.*a.^0.5)*D+sig*diag(xp)*DD-sig*diag(xpp)*D-rho*g*diag(a.^1.5),...
        a.^1.5,zeros(N,1),zeros(N,1); ...
        2*diag(xp)*D,2*diag(yp)*D,zeros(N,1),-2*l*ones(N,1),zeros(N,1); ...
        y'.*I(N,:)-y'.*I(1,:)-INT.*yp',INT.*xp',0,0,-1;...
        I(1,:),zeros(1,N),0,0,0; ...
        dir'];

        %{
        The Jacobian has to be modified accordingly to the boundary conditions as well.
        Remember that one of them was already imposed as a
        constraint in the construction of the Jacobian as an additional equation.
        The remaining three Dirichlet boundary conditions are:
        $$
        x(1)=x_n \quad , \quad y(0) = y_1 \quad , \quad y(1) = y_n
        $$
        %}

        % Boundary conditions
        loc = [1 N N+1];
        f(loc)=[x(N)-xn; y(1)-y1; y(N)-yn];
        A(loc,:)=[I(N,:),zeros(1,N),0,0,0; zeros(1,N),I(1,:),0,0,0; zeros(1,N),I(N,:),0,0,0];

        %{
        If the residual is smaller than the tolerance, the shape of the droplet
        for this volume is found. If this is not the case, an additional
        Newton iteration has to be performed.
        %}

        % convergence test
        res=norm(f);
        %plot(x,y,'b-',initguess_x,initguess_y,'r--'); drawnow;
        disp([num2str(count) '  ' num2str(res)]);
        if count>50|res>1e5; disp('no convergence'); endLoop =1; quitcon=1; break; end
        if res<1e-5; quit=1; disp('converged'); continue; end
        
        %{
        # The Newton step
        
        The core of the Netwon iteration is the computation of the new guess
        for the solution.
        %}
        
        % Newton step
        sol=sol-A\f;
        count=count+1;
    end
    
    %{
    # The new direction for the continuation
    
    The new direction for the continuation is found by:
    %}
    
    % New direction
    dir=A\[zeros(N,1);zeros(N,1);0;0;1]; % new direction
    dir=dir/max(abs(dir)); % normalization
    ind=ind+1;
    
    %{
    # Bifurcation diagram and droplet shape
    
    The shape of the pendant drop for a given volume is finally found and its 
    height can be plotted as a function of its volume to construct the bifurcation diagram.
    A saddle node bifurcation is clearly visible from the diagram. After
    the fold, the branch is unstable and the volume has to be reduced to
    find a steady-state solution. Note that the stability of the solution cannot be
    determined with this method. Here the stability character of the branches is
    guessed from physical arguments.
    %}
    
    x=sol(1:N);
    y=sol(N+1:2*N);
    V=sol(end);
    H = max(y);
    
    figure(1)
    plot(x,-y,'b-',[0 L],[0 0],'-k'); axis equal; %axis([-L/2 3/2*L 0 4*L]);
    
    figure(2)
    plot(V,H,'bo')
    hold on
    grid on
    xlabel('V')
    ylabel('H')
    drawnow;
    
end

%{
# The figure

Here we see the result of the computation.

![Bifurcation diagram and the different shapes of droplet](/pendant_drop_volume.png)

%}

%{
# Exercices/Contributions

Modify the previous code to consider the problem of a pendant droplet connected
to a uniform film (zero contact angle) or to a horizontal substrate ( given contact angle depending on materials).
This is for example the situation of the dripping from
the ceiling. 
Start by implementing Neumann boundary conditions to
take into account the inclination of the interface at its extremities.
%}

