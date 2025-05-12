%{
# Mapping 2D made into a function

Please see [diffmat_mapping.m]() to learn how the mapping of meshes is done. This allows to solve equations on meshes which are not rectangular. The use of the present function is examplified in [diffmat_mapping_map2D.m]().

%}
function D=map2D(X,Y,D);
% a function the does a numerical mapping on a 2D domain
% D.x, D.y, D.xx and D.yy are the differentiation matrices on the original
% cartesian grid
% and X, Y are the mapped domain (the domain for which you whish to have
% the differentiation matrices given as output

Nx=size(X,2);
Ny=size(X,1);

% Derivatives of the mapping
xx=D.x*X(:); xy=D.y*X(:);   yx=D.x*Y(:);  yy=D.y*Y(:);
yyy=D.yy*Y(:);   xxx=D.xx*X(:);   yxy=D.y*yx;  yxx=D.xx*Y(:);   xyy=D.yy*X(:);   xyx=D.x*xy;
jac=xx.*yy-xy.*yx;

% diff matrices
DMx=spd(yy./jac)*D.x-spd(yx./jac)*D.y;
DMy=-spd(xy./jac)*D.x+spd(xx./jac)*D.y;

DMxx=spd((yy./jac).^2)*D.xx ...
     +spd((yx./jac).^2)*D.yy ...
     +spd((yy.^2.*yxx-2*yx.*yy.*yxy+yx.^2.*yyy)./jac.^3)*(spd(xy)*D.x-spd(xx)*D.y) ...
     +spd((yy.^2.*xxx-2*yx.*yy.*xyx+yx.^2.*xyy)./jac.^3)*(spd(yx)*D.y-spd(yy)*D.x);
DMxx2=-2*spd(yx.*yy./jac.^2)*(D.y*D.x);

DMyy=spd((xy./jac).^2)*D.xx ...
     +spd((xx./jac).^2)*D.yy ...
     +spd((xy.^2.*yxx-2*xx.*xy.*yxy+xx.^2.*yyy)./jac.^3)*(spd(xy)*D.x-spd(xx)*D.y) ...
     +spd((xy.^2.*xxx-2*xx.*xy.*xyx+xx.^2.*xyy)./jac.^3)*(spd(yx)*D.y-spd(yy)*D.x);
DMyy2=-2*spd(xx.*xy./jac.^2)*(D.y*D.x);

% the differentiation matrices
D.x=DMx; D.xx=DMxx+DMxx2; 
D.y=DMy; D.yy=DMyy+DMyy2;

% integration operators
D.w=D.w(:).*jac; % double integration 
D.wx=D.wx(:).*xx; % integration in x
D.wy=D.wy(:).*yy; % integration in y

D.w=reshape(D.w,Ny,Nx);
D.wx=reshape(D.wx,Ny,Nx);
D.wy=reshape(D.wy,Ny,Nx);
