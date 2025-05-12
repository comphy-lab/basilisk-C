%{
# Burgers equation without time marching

Just like in [advection_global.m](), but here for a nonlinear system, the Burgers equation. 
Since the system is now nonlinear, we use the Newton iterations to converge progressively to the solution.

We do as well the time marching for a comparison.

%}

clear all; clf

% parameters
Nx=61; % number of gridpoints
Nt=100; % gridpoints in time
Lx=10; % domain length in x
Lt=20; % duration in time
xpos=Lx/4; % position of the initial condition
mu=0.1; % diffusion

% differentiation 
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,Lx,Nx);
[d.y,d.yy,d.wy,y]=dif1D('fd',0,Lt,Nt,5);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

% rename y to time derivative
D.t=D.y;
D.tt=D.yy;
t=y;
T=Y;

l.left=[l.cbl;l.left;l.ctl];
l.right=[l.cbr;l.right;l.ctr];
l.start=l.bot;
l.end=l.top;


%initial guess
u0=exp(-((X-xpos)/1).^2);
sol=u0(:);

% Newton iterations
quit=0;count=0;
while ~quit


    % the present solution and its derivatives
    u=sol; ux=D.x*u;  uxx=D.xx*u; ut=D.t*u;
    
%{
# The function and its Jacobian

The Burgers equation is
$$
u_t+uu_x=\mu u_{xx}
$$
so the function that we want to become zero is
$$
f=u_t+uu_x-\mu u_{xx}
$$
and the jacobian is
$$
A=\partial_t+u\partial_x+u_xI-\mu\partial_{xx}
$$
we should not be afraid to have time derivatives in the Jacobian, since they do not have formally any difference with the spatial derivatives.
%}
    % nonlinear function
    f=ut+u.*ux-mu*uxx; 

    % analytical jacobian
    A=D.t+(spd(u)*D.x+spd(ux))-mu*D.xx;

    % boundary conditions
    loc=[l.start; l.left; l.right];
    C=I(loc,:); 
    f(loc)=C*(u-u0(:)); 
    A(loc,:)=C;
    
    % showing present solution
    subplot(1,2,1);
    surf(X,T,reshape(u,Nt,Nx)); view(2); 
    xlabel('x'); ylabel('t');title('Global solution');
    drawnow
      
    % convergence test
    res=norm(f); disp([num2str(count),' ' num2str(res)])
    if count>50|res>1e5; disp('no convergence'); break; end
    if res<1e-5; quit=1; disp('converged'); continue; end

    % Newton step
    sol=sol-A\f;
    count=count+1;
end


%{
# Time marching

Here we compare with the solution that we get using time marching. For this we treat the diffusion term by Crank-Nicolsin, just like in [vibrating_string.m](), and for the nonlinear terms, we use forward Euler. Since forward Euler is explicit, this avoids to solve a nonlinear system at each time step.

%}
% comparison with time marching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% system matrices
I=eye(Nx);
E=I;
A=mu*d.xx; 

% boundary conditions
loc=[1;Nx];
E(loc,:)=0;
A(loc,:)=I(loc,:);

%{
Just like in [vibrating_string.m]() we use the Crank-Nicolson scheme for the linear terms, this gives
$$
M_m u(t+dt)=M_p u(t)-uu_x dt
$$
where the integral of the nonlinear term from $t$ to $t+dt$ was approximated by $uu_xdt$ like in the forward Euler scheme, and with the matrices
$$
M_m=(E-Adt/2), M_p=(E+Adt/2)
$$
%}
% march in time matrix 
dt=t(2)-t(1);
Mm=(E-A*dt/2);
M=Mm\(E+A*dt/2);

% initial condition
q=u0(1,:)';

% marching loop
for ind=1:(Nt-1)    
    nl=-q.*(d.x*q);
%{
Here we should not forget that we have to impose the boundary conditions. Here the nonlinear terms act like a nonhomogeneous forcing. We should remember not to put a forcing on the boundary conditions. This would be good for nonhomogeneous boundary conditions, but this is not what we want to do here. So we put some zeros in the lines of the constraints (stored in *loc*). 
%}
    nl(loc)=0;
    
    q=M*q+Mm\nl*dt; % one step forward
end
    
% plotting
subplot(1,2,2);
plot(x,q,'b',X(l.end),u(l.end),'r.-');    
axis([0,Lx,0,0.5])    
legend('march in time','global'); title('Vibrating string')
xlabel('x'); ylabel('final time');


set(gcf,'paperpositionmode','auto');
print('-dpng','-r80','burgers_global.png');



%{
On the left subplot, you see the golbal solutino obtained using Newton iterations. And on the right you have the final state, comparing the two methods.

![Validation](burgers_global.png)

# Exercises/contributions

* Please do the time marching of the Burgers equations using as well Crank_nicolson for the nonlinear term, and so solving at each time step a nonlinear equation using the Newton iterations.

%}
