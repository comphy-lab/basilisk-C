%{
# Testing the domain geometry derivative

This is like [domain_derivative_1D.m]() but here in 2D. The only change is that instead of imposing
$$
h_{xx}+2=0
$$
in the domain, we impose
$$
h_{xx}+h_{yy}+2=0
$$
and we impose as a Dirichlet boundary condition on the left and right sides of the domain, the value of the theoretical solution of the problem
$$
h=y(1-y)
$$
%}

disp('%%%%%%%')
clear all; clf;

% parameters
Nx=20;           % number of grid points in x
Ny=20;            % number of grid points in y
p=-1;         % desired slope at top boundary
L=0.9;      % the length in y of the computational domain

% differentiation
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,1,Nx);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,L,Ny);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

ZZ=blkdiag(Z,zeros(Nx,Nx));
II=blkdiag(I,eye(Nx,Nx));

l.h=(1:NN)';
l.eta=NN+[1:Nx]';
l.top=[l.ctl; l.top; l.ctr];
l.bot=[l.cbl; l.bot; l.cbr];

%initial guess
eta0=zeros(Nx,1);
sol0=[Y(:).*(1-Y(:));eta0];
sol=sol0;

% Newton iterations
quit=0;count=0;
while ~quit

    % the present solution and its derivatives
    h=sol(l.h); hy=D.y*h;  hxx=D.xx*h; hyy=D.yy*h; 
    eta=sol(l.eta);
    
    % nonlinear function
    f=[hxx+hyy+2; hy(l.top)+eta.*hyy(l.top)-p*ones(Nx,1)]; 

    % analytical jacobian
    A=[D.xx+D.yy, ZZ(l.h,l.eta); ...
       D.y(l.top,:)+diag(eta)*D.yy(l.top,:), diag(hyy(l.top))];

    % Boundary conditions
    loc=[l.left; l.right; l.bot; l.top];
    C=II([l.bot;l.right;l.left],:);
    
    f(loc)=[C*(sol-sol0); ... % the linear boundary conditions
            h(l.top)+diag(eta)*hy(l.top)]; % the nonlinear boundary conditions
    A(loc,:)=[C; ...
              II(l.top,:)+diag(eta)*D.y(l.top,:)*II(l.h,:)+diag(hy(l.top))*II(l.eta,:)];
%{
# Linear extrapolation

I do here the linear extrapolation in one single line of code by using array operations. This is a little tricky if you are not used to it. If you prefer, you can replace this by a loop of x, and at each loop do the linear extrapolation like we did  in 1D in [domain_derivative_1D.m]().

%}
    % Show present solution
    ne=10;
    [XX,YY]=meshgrid(x,linspace(L,L+max(eta),ne)); % grid for linear extrapolation
    hh=(h(l.top)*ones(1,ne))'+(hy(l.top)*ones(1,ne))'.*(YY-L); % linear extrapolation
    mesh(X,Y,reshape(h,Ny,Nx),'edgecolor','b'); hold on
    mesh(XX,YY,hh,'edgecolor','r'); 
    plot3(x,L+eta,x*0,'k','linewidth',2); hold off % draw boundary
    legend('h','linear extrapolation','estimated boundary position','location','north')
    view([-101,16])
    drawnow;
    
    % convergence test
    res=norm(f,inf);
    disp([num2str(count) '  ' num2str(res)]);
    if count>50|res>1e5; disp('no convergence'); break; end
    if res<1e-9; quit=1; disp('converged'); continue; end

    % Newton step
    sol=sol-A\f;
    count=count+1;
end

set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','domain_derivative_2D.png')

%{
# The figure

![](domain_derivative_2D.png)

# Exercices/Contributions

* Please 
* Please 

%}