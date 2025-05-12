%{
# Presentation of the differentiation matrix
In this program, we compare the second derivative computed by
multiplication twice by the first derivative matrix or once by the second
derivative matrix. We use the program blabla for the differentiation
matrix.

%}

clear all; clf

% parameters
L=2*pi; % domain length
N=15; % number of points

% the grid
x=linspace(0,L,N)';
h=x(2)-x(1); % the grid size

% first derivative
D=zeros(N,N);
D(1,1:3)=[-3/2, 2, -1/2]/h;
for ind=2:N-1
    D(ind,ind-1:ind+1)=[-1/2, 0, 1/2]/h;
end
D(end,end-2:end)=[1/2, -2, 3/2]/h;

% second derivative
DD=zeros(N,N);
DD(1,1:3)=[1, -2, 1]/h^2;
for ind=2:N-1
    DD(ind,ind-1:ind+1)=[1, -2, 1]/h^2;
end
DD(end,end-2:end)=[1, -2, 1]/h^2;


%{

Now that the differentiation matrices are built and stored in D and DD.In
the program below we compute the second derivative of cos(x) for the
comparaison.

%}

% test of the derivatives
f=cos(x); %initial function
fpp=-sin(x)'; %second derivative

plot(x,DD*cos(x),'m.--',x,D*D*cos(x),'b',x,-cos(x),'g.-');
legend('DD*cos(x)','D*D*cos(x)','real')
print('-djpeg','compare.jpg'); % save the figure
xlabel('x');ylabel('y');

%{

# test of the computation
And here is the figure that is produced by the code:

![Comparison of the numerical and exact derivative for the first and second derivative of a cosinus](/diffmat.png)

matrix or once by the second derivative matrix. 


%}