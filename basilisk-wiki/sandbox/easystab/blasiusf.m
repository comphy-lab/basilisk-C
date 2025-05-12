%{

Simply the code for the Blasius boundary layer profile [blasius.m](), made into a function.

Input:

* L is the domain height (typically 10)
* N is the numbre of grid points (Typically 100)

Output

* y are the grid poins
* u the velocity profile on this grid

%}

function [y,u]=blasiusf(L,N);

% differentiation matrices
scale=-2/L;
[y,DM] = chebdif(N,3); 
dy=DM(:,:,1)*scale;    
dyy=DM(:,:,2)*scale^2;    
dyyy=DM(:,:,3)*scale^3;
y=(y-1)/scale;  	
I=eye(N); Z=zeros(N,N);

% initial guess from order 4 Polhausen profile
scalepol=6*sqrt(35/37);
eta=y/scalepol; deta=dy*scalepol; % rescaled vertical coordinate
u0=1-(1-eta).^3.*(1+eta);
u0(eta>1)=1;

A=deta; A(1,:)=I(1,:); u0(1)=0; % set up integration 
g0=A\u0; % compute integral

sol=g0;

% Newton iterations
quit=0;count=0;
while ~quit
    
    % the present solution and its derivatives
    g=sol; gy=dy*g; gyy=dyy*g; gyyy=dyyy*g;
    
    % nonlinear function
    f=2*gyyy+g.*gyy;
    
    % analytical jacobian
    A=2*dyyy+diag(g)*dyy+diag(gyy)*I;

    % Boundary conditions
    loc=[1,2,N];
    C=[I(1,:); dy(1,:); dy(N,:)];
    f(loc)=C*g-[0;0;1];
    A([1,2,N],:)=C;
    
    % convergence test
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if count>50|res>1e5; disp('Blasius: no convergence'); break; end
    if res<1e-5; quit=1; disp('Blasius: converged'); continue; end
    
    % Newton step
    sol=sol-A\f;
    count=count+1;
end

% the velocity profile
u=dy*g;
