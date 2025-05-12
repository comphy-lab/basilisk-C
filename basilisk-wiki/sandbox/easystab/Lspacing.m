%{

A function used by [poisson_2D_hollow_cylinder.m]()

%}

function [ x, x_linear ] = Lspacing( Nx, lambda, Lx )
%% LSPACING function
%   INPUTS:
%   Nx: Number of divisions.
%   lambda: Kind of division.
%   Lx: Total Length of the spaced item.
%   OUTPUTS:
%   x: Distance positioning of the control point. Vector.
%   x_linear: Regularly (linear) spaced control points. Vector.

Xtr = (0:(1/Nx):1)';
x_linear = zeros((Nx-1),1);
for ii = 1:Nx
    x_linear(ii,1) = (Xtr(ii+1)+Xtr(ii))/2;
end

if lambda <= 3 && lambda >= 2
   L3 = x_linear;
   L2 = 1 + sin((pi/2).*(-1+x_linear));
   A = 3 - lambda;
   B = lambda - 2;
   x = B.*L3 + A.*L2;
elseif lambda <= 2 && lambda >= 1
   L2 = 1 + sin((pi/2).*(-1+x_linear));
   L1 = 0.5.*(1+cos(pi.*(1+x_linear)));
   A = 2 - lambda;
   B = lambda - 1;
   x = B.*L2 + A.*L1;
elseif lambda <= 1 && lambda >= 0
   L1 = 0.5.*(1+cos(pi.*(1+x_linear)));
   L0 = x_linear;
   A = 1 - lambda;
   B = lambda;
   x = B.*L1 + A.*L0;
elseif lambda <= 0 && lambda >= -1
   L0 = x_linear;
   L1 = 0.5.*(1+cos(pi.*(1+x_linear)));
   A = - lambda;
   B = lambda - (-1);
   x = B.*L0 + A.*L1;
elseif lambda <= -1 && lambda >= -2
   L1 = 0.5.*(1+cos(pi.*(1+x_linear)));
   L2 = sin((pi/2).*x_linear);
   A = (-1) - lambda;
   B = lambda - (-2);
   x = B.*L1 + A.*L2;
elseif lambda <= -2 && lambda >= -3
   L2 = sin((pi/2).*x_linear);
   L3 = x_linear;
   A = (-2) - lambda;
   B = lambda - (-3);
   x = B.*L2 + A.*L3;
else
    error('lambda = %0.2f out of limits. Lambda should be between [3,-3].',lambda);
end

x = x.*Lx;
x_linear = x_linear.*Lx;
end
