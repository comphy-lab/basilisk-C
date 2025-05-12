function [x,D]=fddif(N,order,pts);
% build equispaced grid on [-1,1], and 
% five points finite diference matrix for N mesh points 

%pts=5;
x=linspace(-1,1,N)';
h=x(2)-x(1);

% subroutine for finite difference weights
W=ufdwt(h,pts,order);
t=(pts+1)/2;

%central difference in the middle
D=spdiags(ones(N,1)*W(t,:),-t+1:t-1,N,N);

for indd=1:t-1
  D(indd,1:pts)=W(indd,:);   % at the left boundary
  D(N-indd+1,end-pts+1:end)=W(end-indd+1,:); % at the right boundary
end
