%{
# The Venturi flow

Here is a code based on [venturi.m]() where we can switch between Navier-Stokes and the Boundary layer equations using the parameter *ns* (1 for Navier-Stokes and 0 for boundary layer)

NS
$$w \frac{\partial w}{\partial z} + u \frac{\partial w}{\partial r} = 
-\frac{\partial p}{\partial z} + \frac{1}{Re} (\frac{\partial^2 w}{\partial z^2} +
\frac{  1}{  r}\frac{\partial }{\partial r}r\frac{\partial w}{\partial r})
$$
$$w \frac{\partial u}{\partial z} + u \frac{\partial u}{\partial r} = 
-\frac{\partial p}{\partial r} + \frac{1}{Re} (\frac{\partial^2 u}{\partial z^2} +
\frac{  1}{  r}\frac{\partial }{\partial r}r\frac{\partial u}{\partial r})
$$
$$\frac{\partial w}{\partial z}   +
\frac{  1}{  r}\frac{\partial }{\partial r}ru=0$$

RNSP (Prandtl equations in the whole tube)
$$w \frac{\partial w}{\partial z} + u \frac{\partial w}{\partial r} = 
-\frac{\partial p}{\partial z} + \frac{1}{Re} ( 
\frac{  1}{  r}\frac{\partial }{\partial r}r\frac{\partial w}{\partial r})
$$
$$0 = 
-\frac{\partial p}{\partial r}
$$
$$\frac{\partial w}{\partial z}   +
\frac{  1}{  r}\frac{\partial }{\partial r}ru=0$$
these are just the same equation whith just $p(z)$ and negligible $\frac{\partial^2 w}{\partial z^2}$, the $1/Re$ can be then absorbed in the longitudinal scale of $z$.

And we do the comparison of the wall-shear stress or the two cases.

Dependency:

* [chebdif.m]()
* [fddif.m]()
* [mapping2D.m]()

%}

clear all; 

% loop on ns
for ns=[1,0]
format compact

%%%% parameters 
Re=1000; % reynolds number
Nz=151; % number of grid nodes in z
Nr=40; %number of grid nodes in r
Lz=30; % length in z of the domain [0,Lz]
pts=5; % number of points in finite difference stencils
amp=0.5; % radius at venturi
zpos=8; % position of the neck

% Chebychev in r
scale=1;
[r,DM] = chebdif(Nr,2);
dr=DM(:,:,1)/scale;
drr=DM(:,:,2)/scale^2;
r=r*scale;
ir=-([diff(r)',0]+[0,diff(r)'])/2; 

% finite difference in z
scale=Lz/2;
[z,dz] = fddif(Nz,1,pts);
[z,dzz] = fddif(Nz,2,pts);
dz=dz/scale;
dzz=dzz/(scale^2);
z=(z+1)*scale;
iz=([diff(z)',0]+[0,diff(z)'])/2; 

% 2D differentiation and integration
Dr=kron(speye(Nz),dr);
Drr=kron(speye(Nz),drr);
Dz=kron(dz,speye(Nr));
Dzz=kron(dzz,speye(Nr));
ww=ir(:)*iz(:)'; 
wz=ones(Nr,1)*iz(:)'; 
wr=ir(:)*ones(1,Nz); 

% mapping of the mesh
disp('Mapping of the mesh')
[Z,R]=meshgrid(z,r);
etar=1-(1-amp)*exp(-((z-zpos)/3).^2); 
R=R.*repmat(etar',Nr,1);  
[Dz,Dr,Dzz,Drr,ww,wz,wr]=mapping2D(Z,R,Dz,Dr,Dzz,Drr,ww,wz,wr);

% cylindrical Laplacian
rm1=spd(1./R);   rm2=spd(1./R.^2);
lap=Drr+ns*(Dzz)+rm1*Dr;
   
% vectors for selecting coordinates 
NN=Nr*Nz;
dom=reshape(1:NN,Nr,Nz);
out=dom(2:Nr-1,end); out=out(:);
in=dom(2:Nr-1,1); in=in(:);
top=dom(1,2:Nz-1); top=top(:);
bot=dom(end,2:Nz-1); bot=bot(:);
cor=dom(1,[1 end]); cor=cor(:);
u=(1:NN)'; w=u+NN; p=w+NN;

%%%% preparing boundary conditions
Ze=spalloc(NN,NN,0); I=speye(NN); II=speye(3*NN);

neuploc=[cor;cor(2)-Nr];  % where to impose the neumann condition on the pressure
p0loc=Nr+1; % where to impose zero pressure
dir=[cor;in;top;bot]; % where to put Dirichley on u and w

loc=[u(dir); w(dir); p(p0loc); ...
    u(out); ...
    w(out); ...
    p(neuploc)];

C=[II([u(dir);w(dir);p(p0loc)],:); ...     % Dirichlet on u,w,and p
   Dz(out,:), Ze(out,:), Ze(out,:); ...   % Neuman on u at outflow
   Ze(out,:), Dz(out,:),Ze(out,:); ...    % Neumann on w at outflow
   Ze(neuploc,:), lap(neuploc,:)/Re, -Dz(neuploc,:)]; % neuman constraint on pressure

% initial guess
U=zeros(NN,1);
W=(1-(R(:,1)*ones(1,Nz)).^2); % mean velocity 1 on the pipe of diameter 1,
P=-Z/Re; P=P-P(p0loc); % pressure zero at p0loc
base=[U(:);W(:);P(:)];

% Newton iterations
disp('Newton loop')
sol=base;
quit=0;count=0;
while ~quit     
 
    % the present solution and its derivatives
    U=sol(u); W=sol(w); P=sol(p);
    Uz=Dz*U; Ur=Dr*U;
    Wz=Dz*W; Wr=Dr*W; 
    Pz=Dz*P; Pr=Dr*P;
    
    
    % the nonlinear function
    f=[ns*(-(U.*Ur+W.*Uz)+(lap*U-rm2*U)/Re)-Pr; ...
       -(U.*Wr+W.*Wz)+lap*W/Re-Pz; ...
      (rm1+Dr)*U+Dz*W];
    
    % Jacobian 
    A=[ns*(-(spd(W)*Dz+spd(U)*Dr+spd(Ur))+(lap-rm2)/Re), ns*(-spd(Uz)), -Dr; ...
         -spd(Wr), -(spd(W)*Dz+spd(Wz)+spd(U)*Dr)+lap/Re, -Dz; ...
         rm1+Dr, Dz, Ze];
     
    % Boundary conditions 
    f(loc)=C*(sol-base);
    A(loc,:)=C;

    % plotting
    subplot(3,1,1);
    surf(Z,R,reshape(U-1,Nr,Nz)); view(2); shading interp; hold on
    
    selr=1:Nr; selz=1:6:Nz;
    ww=reshape(W,Nr,Nz); uu=reshape(U,Nr,Nz); 
    quiver(Z(selr,selz),R(selr,selz),ww(selr,selz),uu(selr,selz),'k');
    axis([0,Lz,-1,1]);
    xlabel('z'); ylabel('r'); title('radial velocity U'); grid off;hold off
    
    subplot(3,1,2);
    surf(Z,R,reshape(W,Nr,Nz)); view(2); shading interp; 
    xlabel('z'); ylabel('r'); title('axial velocity W'); grid off
    selr=1:3:Nr; selz=1:3:Nz;
    drawnow
    
    % convergence test  
    res=norm(f);
    disp([num2str(count) '  ' num2str(res)]);
    if count>50|res>1e5; disp('no convergence');break; end
    if res<1e-5; quit=1; disp('converged'); continue; end
    
    % Newton step
    sol=sol-A\f;   
    count=count+1;
end


% plot the wall shear stress
subplot(3,1,3);
switch ns; case 1; co='b'; case 0; co='r'; end
plot(Z(top),Dr(bot,:)*W,co); hold on
xlabel('z'); ylabel('Wr'); title('wall shear stress');
grid on
end

%{
![velocity field and wall shear stress](venturi_pyl.png)

notice that the wall shear stresses are superposed

## Bibliography

Franz Chouly, P.-Y. Lagrée (2012): 
"[Comparison of computations of asymptotic flow models in a constricted channel.](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/choulylagree12.pdf)",
Applied Mathematical Modelling 36 (2012), pp. 6061-6071

%}
