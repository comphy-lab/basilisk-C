%{
# Testing the domain geometry derivative

This is just like in [domain_derivative_1D_adapt.m]() but here in 2D.

%}

disp('%%%%%%%')
clear all; clf;

% parameters
Nx=10;           % number of grid points in x
Ny=20;            % number of grid points in y
p=-1;         % desired slope at top boundary
L=0.95;      % the length in y of the computational domain

% differentiation
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,1,Nx);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,L,Ny);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);
D0=D; Y0=Y;

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

%{
# Domain adaptation

We do here just like we did in 1D, making a loop running over $x$. 
%}
    % addapt the domain
    H=reshape(h,Ny,Nx); Hy=reshape(hy,Ny,Nx);
    mesh(X,Y,H,'edgecolor','b');hold on
    
    for gre=1:Nx
        ll=Y(Ny,gre);
        y=Y(:,gre);
        Y(:,gre)=Y(:,gre)*(1+eta(gre)/ll); % stretched grid
        yy=linspace(1.0001*ll,2*ll,100)'; % the grid outside the domain
        hh=H(Ny,gre)+(yy-ll)*Hy(Ny,gre); % linear extrapolation outside of the domain
        H(:,gre)=interp1([y; yy],[H(:,gre); hh],Y(:,gre),'splines'); % interpolation to the new grid
    end
    D=map2D(X,Y,D0); % the new differentiation matrices
    sol0=[Y(:).*(1-Y(:));eta0]; % update sol0 for the boundary conditions
    sol(l.h)=H(:); % update h
    sol(l.eta)=0; % set eta to zero

    mesh(X,Y,H,'edgecolor','r'); hold off; 
    xlabel('x');ylabel('y'); zlabel('h'); title('before and after interpolation')
    drawnow
    
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
print('-dpng','-r80','domain_derivative_2D_adapt.png')

%{
# Validation

For the validation, we compare the solution to the analytical solution
$$
h=y(1-y)
$$
%}
% validation
htheo=Y.*(1-Y);
err=norm(h-htheo(:),2)


%{
# The figure

![](domain_derivative_2D_adapt.png)

# Exercices/Contributions

* Please 
* Please 

%}