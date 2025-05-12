%{
# A function of the mapping of 2D differentiation matrices

Here is done inside a function basically the same thing as was done in the code of [diffmat_mapping.m](), [peristalsis.m]().

%}

function [Dx,Dy,Dxx,Dyy,ww,wx,wy]=mapping2D(X,Y,Dx,Dy,Dxx,Dyy,ww,wx,wy);
% a function the does a numerical mapping on a 2D domain
% Dx, Dy, Dxx and Dyy are the differentiation matrices on the original
% cartesian grid
% and X, Y are the mapped domain (the domain for which you whish to have
% the differentiation matrices given as output

% Derivatives of the mapping
xx=Dx*X(:); xy=Dy*X(:);   yx=Dx*Y(:);  yy=Dy*Y(:);
yyy=Dyy*Y(:);   xxx=Dxx*X(:);   yxy=Dy*yx;  yxx=Dxx*Y(:);   xyy=Dyy*X(:);   xyx=Dx*xy;
jac=xx.*yy-xy.*yx;

% diff matrices
DMx=spd(yy./jac)*Dx-spd(yx./jac)*Dy;
DMy=-spd(xy./jac)*Dx+spd(xx./jac)*Dy;

DMxx=spd((yy./jac).^2)*Dxx ...
     +spd((yx./jac).^2)*Dyy ...
     +spd((yy.^2.*yxx-2*yx.*yy.*yxy+yx.^2.*yyy)./jac.^3)*(spd(xy)*Dx-spd(xx)*Dy) ...
     +spd((yy.^2.*xxx-2*yx.*yy.*xyx+yx.^2.*xyy)./jac.^3)*(spd(yx)*Dy-spd(yy)*Dx);
DMxx2=-2*spd(yx.*yy./jac.^2)*(Dy*Dx);

DMyy=spd((xy./jac).^2)*Dxx ...
     +spd((xx./jac).^2)*Dyy ...
     +spd((xy.^2.*yxx-2*xx.*xy.*yxy+xx.^2.*yyy)./jac.^3)*(spd(xy)*Dx-spd(xx)*Dy) ...
     +spd((xy.^2.*xxx-2*xx.*xy.*xyx+xx.^2.*xyy)./jac.^3)*(spd(yx)*Dy-spd(yy)*Dx);
DMyy2=-2*spd(xx.*xy./jac.^2)*(Dy*Dx);

% the differentiation matrices
Dx=DMx; Dxx=DMxx+DMxx2; 
Dy=DMy; Dyy=DMyy+DMyy2;

% integration operators
ww=reshape(ww(:).*jac,size(X)); 
wx=reshape(wx(:).*xx,size(X)); 
wy=reshape(wy(:).*yy,size(X)); 


