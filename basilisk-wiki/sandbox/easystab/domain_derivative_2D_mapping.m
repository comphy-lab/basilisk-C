%{
# Testing the domain geometry derivative

This is just like in [domain_derivative_1D_mapping.m]() but here in 2D.
The subtleties lie in the expression for the Jacobian, especially the part concerning perturbations of the stretching function $\eta$. To learn about this, please check out [stretching_formula.m]() and [jacobian_formula.m]().

%}

disp('%%%%%%%')
clear all; clf;

% parameters
Nx=10;           % number of grid points in x
Ny=20;            % number of grid points in y
p=-1;         % desired slope at top boundary
L=1;      % the length in y of the computational domain

% differentiation
[d.x,d.xx,d.wx,x]=dif1D('cheb',0,1,Nx);
[d.y,d.yy,d.wy,y]=dif1D('cheb',0,L,Ny);
[D,l,X,Y,Z,I,NN]=dif2D(d,x,y);

ZZ=blkdiag(Z,zeros(Nx,Nx));
II=blkdiag(I,eye(Nx,Nx));
T=kron(speye(Nx,Nx),ones(Ny,1));
yy=Y(:);

l.phi=(1:NN)';
l.e=NN+[1:Nx]';
l.top=[l.ctl; l.top; l.ctr];
l.bot=[l.cbl; l.bot; l.cbr];

%initial guess
phi0=Y.*(1-Y);
e0=1.1+0.5*sin(x*pi);
q0=[phi0(:);e0];
q=q0;

% Newton iterations
quit=0;count=0;
while ~quit
    
    % the present solution and its computational-domain derivatives
    phi=q(l.phi); phiy=D.y*phi; phixy=D.x*phiy;  phixx=D.xx*phi; phiyy=D.yy*phi;
    e=q(l.e); ex=d.x*e; exx=d.xx*e; E=T*e; Ex=T*ex; Exx=T*exx;
    
    % Physical space differentation matrices
    Dp.lap=D.xx ...
        +spd(-2*yy.*E.^-1.*Ex)*(D.x*D.y) ...
        +spd(yy.^2.*E.^-2.*Ex.^2+E.^-2)*D.yy ...
        +spd(yy.*(2*E.^-2.*Ex.^2-E.^-1.*Exx))*D.y;
    Dp.x=D.x+spd(-yy.*E.^-1.*Ex)*D.y;
    Dp.y=spd(E.^-1)*D.y;
   
    Dp.lapphie=-2*spd(phixy.*yy)*T*(diag(e.^-1)*d.x-diag(e.^-2.*ex)) ...
        +spd(phiyy)*2*(spd(yy.^2)*T*(diag(e.^-2.*ex)*d.x-diag(e.^-3.*ex.^2))-T*diag(e.^-3)) ...
        +spd(phiy.*yy)*T*(diag(4*e.^-2.*ex)*d.x-diag(e.^-1)*d.xx+diag(-4*e.^-3.*ex.^2+e.^-2.*exx));
    Dp.xphie=spd(yy.*phiy)*T*(diag(-e.^-1)*d.x+diag(e.^-2.*ex));
    Dp.yphie=spd(phiy)*T*diag(-e.^-2);
    
    % Compute the physical space derivatives
    phiplap=Dp.lap*phi;
    phipx=Dp.x*phi;
    phipy=Dp.y*phi;
    
    % nonlinear function
    f=[phiplap+2; phipy(l.top)-p*ones(Nx,1)]; 

    % analytical jacobian
    A=[Dp.lap, Dp.lapphie; ...
       Dp.y(l.top,:), Dp.yphie(l.top,:)];

    % Boundary conditions
    loc=[l.left; l.right; l.bot; l.top];
    C=II(loc,:);
    f(loc)=C*(q-q0);
    A(loc,:)=C; 
             
    % Show present solution
    Yp=reshape(E.*Y(:),Ny,Nx);
    mesh(X,Yp,reshape(phi,Ny,Nx)); hold off; 
    xlabel('x');ylabel('y'); zlabel('phi'); 
    axis([0 1 0 1.5 0 0.3]);
    view(-120,30);
    drawnow
    
    % convergence test
    res=norm(f,inf);
    disp([num2str(count) '  ' num2str(res)]);
    if count>50|res>1e5; disp('no convergence'); break; end
    if res<1e-9; quit=1; disp('converged'); continue; end

    % Newton step
    q=q-A\f;
    count=count+1;
end

set(gcf,'paperpositionmode','auto')
print('-dpng','-r80','domain_derivative_2D_mapping.png')

%{
# Validation

For the validation, we compare the solution to the analytical solution
$$
h=y(1-y)
$$
%}
% validation
phitheo=Y.*(1-Y);
err=norm(phi-phitheo(:))


%{
# The figure

This is the final figure produced by the Newton iterations. Once converged. Please do the computation yourself to see the initial guess and the process of the convergence.

![](domain_derivative_2D_mapping.png)

# Exercices/Contributions

* Please 
* Please 

%}
