%{

%}

function a=spd(b);
% makes a sparse diagonal matrix from the columns of b
n=numel(b);
a=spdiags(b(:),0,n,n);
