a=0.1
n=256
x = linspace(-0.5,0.5,n);
y = linspace(-0.5,0.5,n);
[X,Y] = meshgrid(x,y);
f=((X.^2+Y.^2).^2 + 4*a*X.*(X.^2+Y.^2) - 4*a^2*Y.^2<=0);
output_matrix(f,n,x,y,'test_input_matrix.bin')
