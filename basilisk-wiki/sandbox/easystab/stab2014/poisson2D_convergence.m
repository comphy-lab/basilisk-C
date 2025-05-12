%{
# Convergence study of the Poisson 2D equation

Through this study, we are going to analyze the convergence of the poisson
equation using two different methods of resolution: chebichev's
differentiation matrices and finite difference differentiation matrices.


%}

clear all; clf


%{
In the first stage, we declare all the requried parameters of this study.
Note that this is a convergence analysis, so we need to make a loop through
the number of grid points. However, as it is a 2D problem, our grid extends
throughout both X and Y. We intend to make a single-variable analysis, so
we impose a given aspect ratio between X and Y.

The variable Nx, which will govern the loop, is a vector ranging all the
cases to be stored.
%}

%%%% parameters and flags 
Nx=5:5:40; % gridpoints in x (admits a vector for convergence study)
nNx = numel(Nx); % Number of cases for convergence study
Aspect_Ratio = 1; % Aspect ratio of the grid
Ny=round(Nx.*Aspect_Ratio); % gridpoints in y
Lx=1; % domain size in x
Ly=1; % domain size in y

% Comparation vectors of the order of finite difference methods
h1x = Lx./(Nx-1);
h1y = Ly./(Ny-1);

%{
The variables err and err_FD are used to store the residual of the norm of
the error vector for each analysis case of the main loop.
%}
err = nan(1,nNx); % Preallocates err vectors
err_FD = err;
for n = 1:nNx

% Grid for finite difference
x_FD = (0:h1x(n):h1x(n)*(Nx(n)-1));
y_FD = (0:h1y(n):h1y(n)*(Ny(n)-1));
[X_FD,Y_FD]=meshgrid(x_FD,y_FD);

%{
This code makes use of the functionality of dif1D and dif2D for creating
the matrices of Finite Difference. Please note that the finite differences
scheme used is $O(h^2)$.
%}

% differentiation 
[d.x,d.xx,d.wx]=dif1D('fd',0,Lx,Nx(n),5);
[d.y,d.yy,d.wy]=dif1D('fd',0,Ly,Ny(n),5);
D_FD=dif2D(d);


%{
In this part we create the 1D chebichev and then 2D chebichev
differentiation matrices. It is very important to note that using chebichev
differentiation matrices involves changing the grid discretization. This
has been taken into consideration, and the updated values of $x$ and $y$
corresponds to the new chebichev's discretization.
%}

%1D differentiation matrices Cheb
scale=-2/Lx;
[x,DM] = chebdif(Nx(n),2);
dx=DM(:,:,1)*scale;
dxx=DM(:,:,2)*scale^2;
x=(x-1)/scale;

scale=-2/Ly;
[y,DM] = chebdif(Ny(n),2); 
dy=DM(:,:,1)*scale;    
dyy=DM(:,:,2)*scale^2; 
y=(y-1)/scale; 

% 2D differentiation matrices Cheb
Dx=kron(dx,eye(Ny(n)));
Dxx=kron(dxx,eye(Ny(n)));
Dy=kron(eye(Nx(n)),dy);
Dyy=kron(eye(Nx(n)),dyy);
[X,Y]=meshgrid(x,y);

%{
We impose the same boundary conditions for both finite difference and
chebichev methods. We create two main matrices for solving the problem, one
corresponding to the finite difference and the other for chebichev.
%}
% locations at the boundary
dom=reshape(1:Nx(n)*Ny(n),Ny(n),Nx(n));
top=dom(1,1:end); top=top(:); 
bot=dom(end,1:end); bot=bot(:); 
left=dom(2:end-1,1); left=left(:); 
right=dom(2:end-1,end); right=right(:); 

% System matrix
A=Dxx+Dyy; % Chebichev method
A_FD=D_FD.xx+D_FD.yy; % Finite Difference method


% Forcing
k=2; l=1;
b=-pi^2*(k^2+l^2)*sin(pi*k*X).*sin(pi*l*Y);
b=b(:);
b_FD=-pi^2*(k^2+l^2)*sin(pi*k*X_FD).*sin(pi*l*Y_FD);
b_FD=b_FD(:);

% boundary conditions
II=eye(Nx(n)*Ny(n)); ZZ=zeros(Nx(n)*Ny(n),Nx(n)*Ny(n));
loc=[top; bot; left; right];
A(loc,:)=II(loc,:);
A_FD(loc,:)=II(loc,:);
b(loc)=0;
b_FD(loc)=0;

%{
Then we solve the linear systems, we calculate the norm of the residuals
and store them for plotting.
%}

% solving the linear system
f=A\b;
f_FD=A_FD\b_FD;
solexact=sin(pi*k*X).*sin(pi*l*Y);
solexact_FD=sin(pi*k*X_FD).*sin(pi*l*Y_FD);
err(n) = norm(solexact(:) - f,2);
err_FD(n) = norm(solexact_FD(:) - f_FD,2);

end

%{
For checking purposes, we plot the last case of our main loop. Therefore we
can see the huge difference in error, but also the zones where the biggest
error occurs.
%}


% plotting the result
subplot(2,2,1);
mesh(X,Y,reshape(f,Ny(n),Nx(n)));
xlabel('x'); ylabel('y'); zlabel('f')
title(sprintf('Poisson problem\nChebichev; Ntot = %d',Nx(n)*Ny(n)));
subplot(2,2,2);
mesh(X,Y,reshape(f,Ny(n),Nx(n))-solexact);
xlabel('x'); ylabel('y'); zlabel('f')
title('error - Chebichev');

subplot(2,2,3);
mesh(X_FD,Y_FD,reshape(f_FD,Ny(n),Nx(n)));
xlabel('x'); ylabel('y'); zlabel('f')
title(sprintf('Poisson problem\nFinite Difference; Ntot = %d',Nx(n)*Ny(n)));
subplot(2,2,4);
mesh(X_FD,Y_FD,reshape(f_FD,Ny(n),Nx(n))-solexact_FD);
xlabel('x'); ylabel('y'); zlabel('f')
title('error - Finite Difference');

namefig1 = 'Poisson_Conv_err';
set(gcf,'paperpositionmode','auto');
print('-dpng','-r100',namefig1);

%{
![Poisson 2D for 40 grid elements for X and 40 grid elements for Y](./Poisson_Conv_err.png)
%}

% Plots the convergence study
if nNx > 1
    ErrFig = figure('Name','Convergence Study');
    plot(Nx.*Ny,err,'r.-');
    hold on
    plot(Nx.*Ny,err_FD,'g.-');
    hold off
    xlabel('Number of points');
    ylabel('Norm of the error');
    title('Convergence study');
    set(gca,'YScale','log','XScale','log');
    legend('Chebichev','Finite Difference O(h^2)','Location','East');
end

namefig2 = 'Poisson_Conv_err_evol';
set(gcf,'paperpositionmode','auto');
print('-dpng','-r100',namefig2);

%{
Additionnally, we represent the evolution of the error as a function of
the number of gridpoints. We can see that the chebichev method reaches a
minimum, and then the error starts growing. Therefore we can conclude that
the machine tolerance has been reached, and further refinement of the grid
results in higher errors. 
% 2D differentiation matrices Finite Difference
Dxx_FD=kron(dxx_FD,eye(Ny(n)));
Dyy_FD=kron(eye(Nx(n)),dyy_FD);
%}

%{
![Poisson 2D evolution of the error comparison](./Poisson_Conv_err_evol.png)
%}

