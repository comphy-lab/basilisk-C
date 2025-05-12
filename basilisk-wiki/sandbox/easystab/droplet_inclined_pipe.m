%{
# A droplet coming out a pipe under gravity effects

The aim of this code is to determine the shape of a droplet coming out a
vertical or sligthly inclined pipe. We use Newton method in order to solve the nonlinear
governing equations and Keller arc continuation method on the droplet volume in order to follow
solution branches in the bifurcation plane.

The position of the interface is found solving the stress equilibrium at
the interface. Similarly at the problem
[pendant_drop_volume.m](http://basilisk.fr/sandbox/easystab/pendant_drop_volume.m), the stresses acting on the interface are given
by the surface tension forces and the gravity which introduce an
hydrostatic pressure inside the droplet. Thus we write the stress
equilibrium as 
$$
\rho g h - p = \sigma \kappa \, .
$$
where $\sigma$ is the surface tension and $\kappa$ is the interface
curvature (here we consider only the curvature in the $(x,y)$ plane).

Writing explicitly the curvature we obtain
$$
\rho g y - p = \sigma \frac{x^{''}y^{'}-y^{''}x^{'}}{(x^{'^2}+y^{'^2})^{3/2}} \, ,
$$
where $x$ $y$ are the coordinates of the interface expressed in a parametric form. 
We decide to write the curvature in a parametric form in order to have a
well posed problem also in the case where the interface is overturning 
(two $y$ coordiantes of the interface corresponding to a single $x$ coordinate).

Furthermore we impose conservation of volume
$$
\int_0^1 y x^' ds = V
$$ 
and we define the metrics of the parametric description of the interface
$$
x^{'^2}+y^{'^2}=l^2
$$
where $l$ is the length of the interface.
%}

clear all; clf;

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

    % parameters
    V = 2;                          % starting volume
    L = 1;                          % domain length
    n = 100;                        % number of grid points
    rho = 1;                        % density
    sigma = 1;                      % surface tension
    g = 1;                          % gravity
    theta = -0.0001;               % pipe inclination
    delta = 0.2;                    % continuation step

    % Differentiation and integration
    [D,DD,wt,t]=dif1D('cheb',0,1,n,3);
    Z=zeros(n,n); I=eye(n); 
    
    %{
    # Boundary conditions
    
    Since we have a system containing a second and a first order differential
    equation we need three boundary conditions. Imposing the $x$ and $y$
    coordinate for the starting and ending point of the interface we get
    four conditions. We use three of them as boundary conditions and we can
    still adopt the last one as fourth governing equation which close the system of
    four unknowns $x$, $y$, $p$, $l$.
    %}

    % Boundary conditions
    x1 = -L*cos(theta); xn = L*cos(theta);
    y1 = L*sin(theta); yn = -L*sin(theta);

    %{
    Assuming a small volume and thus small interface deformations we can
    retrieve an intial guess by linearization which reads $x = L(2t-1)$,
    $y = \frac{p}{2\sigma}(L^2-x^2)$, $p = \frac{3\sigma V}{2L^3}$ and $l =
    \int_0^1 \sqrt{x'^2+y'^2}dt$ where $0<t<1$ by definition.
    %}
    
    %initial guess
    p = 3*sigma*V/L^3/2;
    x = 2*t*L-L;
    y = p/2/sigma*(L^2-x.^2);
    l = sum(sqrt(diff(x).^2+diff(y).^2));
    
    
    %{
    The solution vector contains the coordinates of the interface, the
    pressure, the length of the interface and the volume. The volume of the droplet is the
    homotopy parameter, which is varied in order to find different branches of the solution.
    %}
    
    sol = [x; y; p; l; V];
    dir=[zeros(2*n+2,1); 1];   % initial direction
    
    %{
    In the following, the first outer loop is for the continuation over the volume $V$. 
    The inner loop is instead needed to find the steady-state solution
    with a Newton method.
    %}
    
    %initialization
    f = zeros(numel(sol),1);

    disp('Continuation loop')
    ind=1;quitcon=0;
    while ~quitcon

        solprev = sol;
        
        %{
        # Keller's pseudo-arclength continuation
        
        In order to follow different branches of the solution and turn
        around bifurcations we implement Keller's pseudo-arclength
        continuation. For a more detailed explanation of the method look at [Lecture Notes on Numerical Analysis of Nonlinear Equations, by Eusebius Doedel.](http://indy.cs.concordia.ca/auto/notes.pdf).
       
        %}
        
        sol = sol+dir*delta; % new prediction of solution

        % Newton iterations
        disp('Newton loop')
        quit=0;count=0;
        while ~quit
            
        %{
        # Newton loop

        Now we are in the Newtoon loop to find the shape of the droplet with a specific volume. See
        [meniscus.m](http://basilisk.fr/sandbox/easystab/meniscus.m) for
        more details.
        The implementation follows the one adopted in [pendant_drop_volume.m](http://basilisk.fr/sandbox/easystab/pendant_drop_volume.m)
        where we find same governing equation and same type of
        boundary conditions.
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
        the equation for the remaning boundary condition
        $$
        x(0) - x_1 \, ,
        $$
        and the additional equation for the continuation loop:
        $$
        \mathbf{t} \cdot ((V,H)_{new}^T-(V,H)_{old}^T) - \delta \, .
        $$
        The system is therefore a function of the interface coordinates $x$, $y$, 
        the pressure $p$ and also of the length of the former, $l$, via the 
        curvilinear parameter.

        *Notation*: the underscript `p` indicate a derivative w.r.t. `s`
        %}

        % the present solution and its derivatives
        hx = sol(1:n); hpx = D*hx; hppx = DD*hx; hy = sol(n+1:2*n); hpy = D*hy; hppy = DD*hy;
        a = hpx.^2+hpy.^2;  % this is comfortable because ita appears very often

        %current values of pressure, interface length and volume
        p = sol(end-2);
        l = sol(end-1);
        V = sol(end);

        % NONLINEAR FUNCTIONS
        %continuity of normal stress
        f(1:n) = p*a.^1.5-rho*g*hy.*a.^1.5-(hppx.*hpy-hpx.*hppy)*sigma;
        %metrics
        f(n+1:2*n) = hpx.^2+hpy.^2-l^2;
        %volume conservation CAREFUL WITH FLIPPED DROPLET, YOU NEED TO FLIP AROUND THE ORIGIN!!!
        f(2*n+1) = wt*(hy.*hpx)-V;
        %BC
        f(2*n+2) = hx(1)-x1;
        %continuation
        f(2*n+3) = dir'*(sol-solprev)-delta;

        %{
        # Computing the Jacobian
        
        Since $f$ is here a vector quantity, the Jacobian is a block matrix, 
        whose blocks
        are the discretized Jacobian matrices of the different components
        w.r.t. the variables. They are found by perturbing all variables in the
        considered equation with small parameters and neglecting the nonlinear 
        terms, in the same way as what described in [meniscus.m](http://basilisk.fr/sandbox/easystab/meniscus.m).
        %}
        
        % analytical jacobian
        A = [diag(3*p*a.^0.5.*hpx)*D-diag(sigma*hpy)*DD+diag(sigma*hppy)*D-diag(3*rho*g*hy.*hpx.*a.^0.5)*D ...
            diag(3*p*a.^0.5.*hpy)*D-diag(sigma*hppx)*D+diag(sigma*hpx)*DD-diag(a.^1.5*rho*g)-diag(3*rho*g*hy.*hpy.*a.^0.5)*D...
            a.^1.5 zeros(n,2); ... %from continuity of normal stress
            diag(2*hpx)*D diag(2*hpy)*D zeros(n,1) -2*l*ones(n,1) zeros(n,1); ...   %from metrics
            [-wt.*hpy' wt.*hpx' zeros(1,2) -1] + [-hy(1) zeros(1,n-2) hy(n) zeros(1,n+3)];...   %from volume conservation
            1 zeros(1,2*n+2); ...   %from BOUNDARY CONDITION that I use as equation
            dir']; %from continuation

        %{
        The Jacobian has to be modified accordingly to the boundary conditions.
        Remember that one of them was already imposed in the construction of the Jacobian as an additional equation.
        The remaining three Dirichlet boundary conditions are:
        $$
        x(1)=xn \quad , \quad y(1) = y1 \quad , \quad y(1) = yn
        $$
        %}
            
        % Boundary conditions
        f([n n+1 2*n]) = [hx(n)-xn; hy(1)-y1; hy(n)-yn];
        A([n n+1 2*n],1:end) = [[I(n,:) zeros(1,n+2) 0];...
            [zeros(1,n) I(1,:) zeros(1,2) 0]; [zeros(1,n) I(n,:) zeros(1,2) 0]];

        %{
        If the residual is smaller than the tolerance, the shape of the droplet
        for this volume is found. If this is not the case, an additional
        Newton iteration has to be performed.
        %}
            
        %% convergence test
        res=norm(f);
        subplot(2,1,1)
        plot(hx,hy,'b-',[x1 xn],[y1 yn],'k',[x1-2*sin(theta) x1],[y1-2*abs(cos(theta)) y1],...
            'k',[xn-2*sin(theta) xn],[y1-2*abs(cos(theta)) yn],'k'); axis equal; drawnow;
        disp([num2str(count) '  ' num2str(res)]);
        if count>100||res>1e5; disp('no convergence'); delta = delta/2; break; end
        if res<1e-5; quit=1; disp('converged'); continue; end

        %{
        # The Newton step
        
        Computation of the solution, it is obtained by solving a linear system.
        %}
        
        % Newton step
        sol = sol-A\f;

        count=count+1;
        end

        %this is to have the well defined value of the parameter that I
        %monitor
        xprime = D(1,:)*sol(1:n);
        yprime = D(1,:)*sol(n+1:2*n);
        if xprime>0&&yprime>0
            alpha = atan(yprime/xprime);
        elseif xprime>0&&yprime<0
            alpha = atan(yprime/xprime)+2*pi;
        elseif xprime<0&&yprime>0
            alpha = atan(yprime/xprime)+pi;
        elseif xprime<0&&yprime<0
            alpha = atan(yprime/xprime)+pi;
        end
        
        %{
        # The new direction for the continuation

        The new direction for the continuation is found by:
        %}

        % New direction
        dir=A\[zeros(2*n+2,1);1]; % nouvelle direction
        dir=dir/norm(dir); % normalization
        ind = ind+1;
       
        %{
        # Bifurcation diagram and droplet shape

        In our problem we choose the droplet volume $V$ as homotopy parameter
        and we monitor the contact angle $\alpha$ of the droplet when
        touching the pipe, this gives an indication on the bending of the interface.
        The forces acting on the system are surface 
        tension and gravity, the first one trying to keep the droplet all 
        together against the action of its own weigth. When the system is 
        perfectly symmetric (not inclined pipe) the droplet is distribuited evenly respect to
        the pipe axis, in this case raising the droplet volume the interface will start
        to bend in order to better contrast the action of the weigth with surface tension forces, till
        reaching a point where a further bending is not enough anymore and
        a volume reduction is required leading to a bifurcation in our 
        diagram. On the other hand when the system is asymmetric (even
        very slightly inclined pipe) the droplet mass is not distribuited evenly but
        it will be more on one side, from this physical
        observation we can figure out that for the same volume the
        required bending of the interface will be much higher that in the
        case before, thus the bifurcation will be reached for a smaller volume. The
        symmetry of this second phenomenon around the angle $\theta$ will
        give rise to a symmetric bifurcation depending on this parameter
        (pitch-fork bifurcation). In the figure below we plot the contact angle $\alpha$ of the droplet as a function
        of the volume as homotopy parameter, finding the bifurcations previously described. A saddle node
        bifuraction appears around $V=5.2$ if the pipe is not inclined ($\theta=0$) and a
        pitch-fork bifuraction around $V=4.1$ when the pipe is sligthly inclined ($\theta=\pm 10^{-5}$).
        %}

       subplot(2,1,2)
       hold on
       plot(V,alpha,'ob')
       hold off
       grid on
       title('bifurcation diagram')
       xlabel('volume')
       ylabel('\alpha')

       %break if it diverges
       if res>1e5
           break
       end
       
    end
    
%{
# The figure

Result of the computation.

![Bifurcation diagram and the two possible droplet evolution](/droplet_inclined_pipe.png)

%}

%{
# Exercices/Contributions

Play with the code trying different values of the inclination angle
$\theta$, try to guess how it would look like the bifurcation diagram if we
determine $\alpha$ as a function of $\theta$. In order to study more in detail the influence of the angle $\theta$ on
this physical system one could use it as the homotopy parameter keeping the volume constant. Modify
the code in this way and than try to use both the angle $\theta$ and
volume as continuation parameters. This will allow for a deeper
understanding of the influnce of these parameters on the change in
toplogy of the bifurcation diagram.
%}
