%{
The Clenshaw-Curtis integration weights for the chebychev grid. This is used in [dif1D.m](). See [integration.m]() to see them tested in 1D and [integration_2D.m]() in 2D.
%}

function IW=clencurt(N)
%
% Computes the integration weigths for pseudo-chebychev on domain [-1 1]
%
% INPUTS:
% N  : the number of points 
%
% OUTPUT:
% IW : vector of the integration weigths.

nW=0:1:N-1;
jW=0:1:N-1;

bW=ones(1,N); bW(1)=0.5; bW(N)=0.5;
cW=2*bW;
bW=bW/(N-1);

S=cos(nW(3:N)'*jW*(pi/(N-1)));
IW=bW.*[(2+(cW(3:N).*((1+(-1).^nW(3:N))./(1-nW(3:N).^2)))*S)];
