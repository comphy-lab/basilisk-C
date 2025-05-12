function [x,D]=fdper(N,order,pts);
% build equispaced grid on [-1,1], and 
% pts points finite diference matrix for N mesh points 
% for periodic domains

% $$$ N=10
% $$$ order=1
% $$$ pts=9


x=linspace(-1,1,N+1)';
x=x(1:end-1);
h=x(2)-x(1);

% subroutine for finite difference weights
W=ufdwt(h,pts,order);
t=(pts+1)/2;


R=zeros(1,N); R([end-t+2:end,1:t])=W(t,:);
C=zeros(1,N); C([t:-1:1,end:-1:end-t+2])=W(t,:);
D=toeplitz(C,R);

