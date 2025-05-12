%{
# LSPACING function
 
Coded by [Paul Valcke](./PaulValcke.m) and [Luis Bernardos](./Luis.m). This
function is used for controlling the discretization of a mesh in a given
direction. It is especially suited for shrinking the mesh near the walls.
It accepts a lambda $\lambda\in[-3,3]$ parameter that interpolates the
linearly-spaced $x$ vector with a sinus or cosinus spacing.
%}
function [ x, x_linear ] = Lspacing( Nx, lambda, Lx )
%{
INPUTS:

- Nx: Number of divisions.
- lambda: Kind of division.
- Lx: Total Length of the spaced item.

OUTPUTS:

- x: Distance positioning of the control point. Vector.
- x_linear: Regularly (linear) spaced control points. Vector.
%}

Xtr = (0:(1/Nx):1)'; % Transition vector
x_linear = zeros((Nx-1),1);
for ii = 1:Nx
    x_linear(ii,1) = (Xtr(ii+1)+Xtr(ii))/2;
end

%{
# Interpolation between linear and sinus($\lambda\in[3,2]$)

The output spacing ranges from linear to a sinus spacing:

$$ x = 1 + sin(\frac{\pi}{2}(-1+x_l)) $$
%}
if lambda <= 3 && lambda >= 2 % interpolation between
   L3 = x_linear;
   L2 = 1 + sin((pi/2).*(-1+x_linear));
   A = 3 - lambda;
   B = lambda - 2;
   x = B.*L3 + A.*L2;
%{
# Interpolation between sinus and cosinus($\lambda\in[2,1]$)

The output spacing ranges from sinus spacing:

$$ x = 1 + sin(\frac{\pi}{2}(-1+x_l)) $$
   
To a cosinus spacing:
   
$$ x = \frac{1}{2}(1+cos(\pi(1+x_l))) $$
%}   
elseif lambda <= 2 && lambda >= 1
   L2 = 1 + sin((pi/2).*(-1+x_linear));
   L1 = 0.5.*(1+cos(pi.*(1+x_linear)));
   A = 2 - lambda;
   B = lambda - 1;
   x = B.*L2 + A.*L1;
%{
# Interpolation between cosinus and linear($\lambda\in[1,0]$)

The output spacing ranges from cosinus $x =
\frac{1}{2}(1+cos(\pi(1+x_l)))$ to linear sinus spacing.
%}    
elseif lambda <= 1 && lambda >= 0
   L1 = 0.5.*(1+cos(pi.*(1+x_linear)));
   L0 = x_linear;
   A = 1 - lambda;
   B = lambda;
   x = B.*L1 + A.*L0;
%{
# Interpolation between linear and cosinus($\lambda\in[0,-1]$)

The output spacing ranges from linear to cosinus $x =
\frac{1}{2}(1+cos(\pi(1+x_l)))$.
%} 
elseif lambda <= 0 && lambda >= -1
   L0 = x_linear;
   L1 = 0.5.*(1+cos(pi.*(1+x_linear)));
   A = - lambda;
   B = lambda - (-1);
   x = B.*L0 + A.*L1;
%{
# Interpolation between cosinus and sinus ($\lambda\in[-1,-2]$)

The output spacing ranges from cosinus $x =
\frac{1}{2}(1+cos(\pi(1+x_l)))$ to sinus $x=sin(\frac{\pi}{2}x_l)$ .
%} 
elseif lambda <= -1 && lambda >= -2
   L1 = 0.5.*(1+cos(pi.*(1+x_linear)));
   L2 = sin((pi/2).*x_linear);
   A = (-1) - lambda;
   B = lambda - (-2);
   x = B.*L1 + A.*L2;
%{
# Interpolation between sinus and linear ($\lambda\in[-2,-3]$)

The output spacing ranges from sinus $x=sin(\frac{\pi}{2}x_l)$ to linear .
%} 
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