% test de la matrice de chebychev



% parameters 
L=4*pi;             % domain length 
N=50;           % number of grid points


% differentiation and integration 
scale=-2/L; [x,DM] = chebdif(N,2);
D=DM(:,:,1)*scale; DD=DM(:,:,2)*scale^2; x=(x-1)/scale;

% build a function 
f=cos(x);
fx=D*f;
plot(x,fx,'b',x,-sin(x),'k') 