%{
#Function  dif1D

%}

function [d,dd,w,s]=dif1D(type,s0,L,N,opt);

%{ 

This function creates a mesh s, differentiation matrices of first and second
order d and dd, as well as the weight w. 
 
 On output: s is a column-vector containing the mesh points, starting from the point at the left.
            d and dd are the derivation matrices (quare matrices, in full or sparse storing format)
            w is a line-vector containing the weights. 
            
Last parameter 'opt' is optional. For Finite-Difference it is interpreted as
the number of points of the stencil (default :3) and for
Stretched Chebyschev it is interpreted as the parameter b
controlling the stretching.

%}

switch type

%{
## Finite differences (optimized case)

For a domain [s0,s0+L]

Last parameter 'pts' is optional, if absent it is set to 3 (relevant value in almost all cases)
%}
  
   case 'fd' 
       
  if(nargin==4) 
    pts = 3;
  else
   pts = opt;
  end
  scale=L/2;
  [s,d] = fddif(N,1,pts); 
  [s,dd] = fddif(N,2,pts); 
  s=s*scale; 
  d=d/scale; dd=dd/scale^2; s=s-s(1)+s0;
  w=([diff(s'),0]+[0,diff(s')])/2;    
  
  
  %{ 
  ## Finite differences (basic implementation, equivalent to the previous)
  %}   
    
 case 'fds' 
    

  dx = L/(N-1); % grid spacing
  s = [  s0 : L/(N-1) : s0+L ]';  % creates an equispaced mesh with N point on the interval
  for i=2:N-1 % centered three-point formulas 
    d(i,i) = 0;
    d(i,i+1) =  1/(2*dx); 
    d(i,i-1) = -1/(2*dx);
    dd(i,i) = -2/dx^2;
    dd(i,i+1) =  1/(dx^2); 
    dd(i,i-1) =  1/(dx^2);        
  end
 % uncentred formulas for first and last grid points 
   d(1,1) = -3/(2*dx); d(1,2) = 4/(2*dx) ; d(1,3) = -1/(2*dx);
   d(N,N-2) =  1/(2*dx); d(N,N-1) = -4/(2*dx); d(N,N)   =  3/(2*dx);
   dd(1,1) = 2/dx^2; dd(1,2) = -5/dx^2; dd(1,3) = 4/dx^2; dd(1,4) = -1/dx^2; 
   dd(N,N) = 2/dx^2; dd(N,N-1) = -5/dx^2; dd(N,N-2) = 4/dx^2; dd(N,N-3) = -1/dx^2; 
  % "weight" to compute integrals using the trapeze rule
  for i=2:N-1
    w(i) = dx;
  end
  w(1) = dx/2; w(N) = dx/2;
  
%{ 
## Finite differences on a periodic domain 
 
 Domain is defined as [s0,s0+L]
%}
  
 case 'fp' % periodic finite differences 
  scale=L/2;
  [s,d] = fdper(N,1,opt); 
  [s,dd] = fdper(N,2,opt);
  s=s*scale; 
  d=d/scale; dd=dd/scale^2; s=s-s(1)+s0;
  w=ones(1,N)*(s(2)-s(1)); 

%{ 
## Chebyshev collocation 

on domain [s0,s0+L]
%}  

 case 'cheb' % chebychev 
  scale=-L/2;
  [s,DM] = chebdif(N,2);
  d=DM(:,:,1);  
  dd=DM(:,:,2);
  s=s*scale; 
  d=d/scale; dd=dd/scale^2; s=s-s(1)+s0;
  w=L*clencurt(N)/2; 
  
  
%{ 
## Fourier collocation 
   
   (requires a periodic domain [s0,s0+L])
  
%} 
  case 'fou' % Fourier 
  scale=L/(2*pi); 
  [s, d] = fourdif(N, 1);
  [s, dd] = fourdif(N, 2);
  s=s*scale; 
  w=ones(1,N)*(s(2)-s(1));
  d=d/scale; dd=dd/scale^2; s=s-s(1)+s0;

%{ 
## Hermite collocation on infinite domain.
  
  Domain is centred on s0, L controls the density 
  
  (increase L to cluster the points towards s0) 

%} 

  case 'her' % hermite 
  scale=1;
  [s,DM] = herdif(N,2,L);
  d=DM(:,:,1);  
  dd=DM(:,:,2);
  s=s*scale; 
  d=d/scale; dd=dd/scale^2; s=s+s0;
  w=zeros(1,N); % weight is not implemented yet !!!
  
%{
## Stretched Chebyshev ; infinite ; exponential
  
  * s0 is the center of the domain.
  * L controls the clustering of points around the origin : about half points will be localised in [s0-L,s0+L]
  * optional parameter b (default value 0.999) controls the stretching (b should be smaller than 1).
%}
  
  case 'chebInfExp' % chebychev with coordinate stretching to treat large domains
      
 
   if(nargin==4) 
      b = 0.999;
   else
      b = opt;
   end
   
   [dxi,ddxi,wxi,xi]=dif1D('cheb',-b,2*b,N); % Standart Chebyshev for xi \in [-b,b]
   
  s = L*atanh(xi)+s0; % coordinate stretching (exponential)
  G(:, 1) =                          - (xi.^2 - 1) /L;      % G1 = ds / dxi
  G(:, 2) =                   2 *xi .* (xi.^2 - 1) /L^2;    % G2 = d^2 s /d xi^2
  G(:, 3) =    - 2 * (3 * xi.^2 - 1) .* (xi.^2 - 1) / L^3;
  G(:, 4) = 8 * xi .* (3 * xi.^2 - 2) .* (xi.^2 - 1) / L^4;
  
  d = G(:,1).*dxi;
  dd = G(:,1).^2.*ddxi+G(:,2).*dxi;
  w = wxi./G(:,1)';
  
%{
## Stretched Chebyshev ; infinite ; algebraic
  
  * s0 is the center of the domain.
  * L controls the clustering of points around the origin : about half points will be localised in [s0-L,s0+L]
  * optional parameter b (default value 0.999) controls the stretching (b should be smaller than 1).
  %}
  
  case 'chebInfAlg' % chebychev with coordinate stretching to treat large domains
      
  if(nargin==4) 
      b = 0.999;
   else
      b = opt;
  end
  
  [dxi,ddxi,wxi,xi]=dif1D('cheb',-b,2*b,N); % Standart Chebyshev for xi \in [-b,b]
  
  s = L * xi ./ sqrt(1 - xi.^2); 
  K = sqrt(s.^2 + L^2);
  G(:, 1) =                                    L^2 ./ K.^3;
  G(:, 2) =                          - 3 * s * L^2 ./ K.^5;
  G(:, 3) =             3 * (4 * s.^2 - L^2) * L^2 ./ K.^7;
  G(:, 4) = - 15 * s .* (4 * s.^2 - 3 * L^2) * L^2 ./ K.^9;
 
  d = G(:,1).*dxi;
  dd = G(:,1).^2.*ddxi+G(:,2).*dxi;
  w = wxi./G(:,1)';
  
 %{
## Stretched Chebyshev ; semi-infinite ; exponential
  
  * s0 is the left bound of the domain.
  * L controls the clustering of points around the origin : about half points will be localised in [s0,s0+L]
  * optional parameter b (default value 0.999) controls the stretching (b should be smaller than 1).
  %}
  
 case 'chebSemiInfExp' % chebychev with coordinate stretching to treat large domains
 
 if(nargin==4) 
      b = 0.999;
   else
      b = opt;
 end
     
  [dxi,ddxi,wxi,xi]=dif1D('cheb',0,b,N); % Standart Chebyshev for xi \in [0,b]
  
  s = L*atanh(xi)+s0; % coordinate stretching (exponential)
  G(:, 1) =                          - (xi.^2 - 1) /L;      % G1 = ds / dxi
  G(:, 2) =                   2 * xi .* (xi.^2 - 1) /L^2;    % G2 = d^2 s /d xi^2
  d = G(:,1).*dxi;
  dd = G(:,1).^2.*ddxi+G(:,2).*dxi;
  w = wxi./G(:,1)';

  %{
## Stretched Chebyshev ; semi-infinite ; algebraic
  
  * s0 is the left bound of the domain.
  * L controls the clustering of points around the origin : about half points will be localised in [s0,s0+L]
  * optional parameter b (default value 0.999) controls the stretching (b should be smaller than 1).
  %}
  
  case 'chebSemiInfAlg' % chebychev with coordinate stretching to treat large domains
      
  if(nargin==4)  
      b = 0.999;  
  else
      b = opt; 
  end
   
  [dxi,ddxi,wxi,xi]=dif1D('cheb',0,b,N); % Standart Chebyshev for xi \in [0,1]

  s = L * xi ./ sqrt(1 - xi.^2); 
  K = sqrt(s.^2 + L^2);
  G(:, 1) =                                    L^2 ./ K.^3;
  G(:, 2) =                          - 3 * s .* L^2 ./ K.^5;
 
  d = G(:,1).*dxi;
  dd = G(:,1).^2.*ddxi+G(:,2).*dxi;
  w = wxi./G(:,1)';
  

end
