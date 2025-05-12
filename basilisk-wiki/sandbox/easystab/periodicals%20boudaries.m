%{
# The differentiation matrix using periodic boundaries

This code realizes a suggestion of [diffmat.m]() of building differentiation matrices with finite difference for a periodic domain, thus avoiding the use of decentered stencils.

%}

clear all; clf

% parameters
L=2*pi; % domain length
N=15; % number of points

%{
To make it work with periodicals boundaries, it is necessary to adapt the grid from $x=linspace(0,L,N)'$ to $x=linspace(0,L,N+1)';$ and precise that the end of the first grid is equal to the next beginning like that : $x=x(1:end-1);$
%}
% the grid
x=linspace(0,L,N+1)'; x=x(1:end-1);
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
#Principle of periodicals boundaries
Before the grid was represented by a single line such as :
$$ \left(f_1,f_2,\dots,f_{n-1},f_n \right)$$
Now the grid is an infinite loop. Which means that the first item of the second line is the last item of the first one and so on. After the creation of the differents line, you concatenate them into one infinite line. It is represented like that :
$$ \left(f_1,f_2,\dots,f_{n-1},f_n \right) \left(f_n,f_1,f_2,\dots,f_{n-1},f_n \right) \dots  \left(f_n,f_1,f_2,\dots,f_{n-1},f_n \right)$$
With the concatenation :
$$  \left(f_1,f_2,\dots,f_{n-1},f_n,f_1,f_2,\dots,f_{n-1},f_n,f_1,f_2,\dots,f_{n-1},f_n \right) $$
With this type of boundaries, two of the scheme are changed. Thus : 
$$
f_{x,1}=\frac{f_2-f_N}{2h} $$

$$ f_{x,N}=\frac{f_1-f_{N-1}}{2h} 
$$
That implies, that the matrix of differantiation is now : 
$$\left(
\begin{array}{l}
f_{x,1} \\
f_{x,2} \\
\vdots \\
f_{x,N-1} \\
f_{x,N} 
\end{array}
\right)
=
\frac{1}{2\Delta x}\begin{pmatrix}
0 & 1 & 0 & \dots & 0 & -1 \\
 -1 & 0 & 1 & & &  \\
& \ddots & \ddots & \ddots & &\\
& & \ddots & \ddots & \ddots \\
& & & -1 & 0 & 1 \\
1 & 0 & \dots & 0 & -1 & 0 
\end{pmatrix}
\left(\begin{array}{l}
f_1 \\
f_2 \\
\vdots \\
f_{N-1} \\
f_N 
\end{array}\right)
$$
%}

% periodicals boundaries
A=zeros(N,N);
for ind=2:N-1
    A(ind,ind-1:ind+1)=[-1/2, 0, 1/2]/h;
end
A(1,[2,N])=[1/2,-1/2]/h; 
A(N,[1,N-1])=[1/2,-1/2]/h;

%{
Now, let this what difference does that make, with the usual matrix. 
%}

% convergency test 
plot(x,A*cos(x)+sin(x),'r.--',x,D*cos(x)+sin(x),'g.--');
legend('periodicals boundaries ','no periodicals boundaries ');
xlabel('x'); ylabel('error'); title('Error for the two methods')
print('-dpng','-r100','diffmaterr.png'); % save the figure

%{

# Test of the computation
And here is the figure that is produced by the code:

![Difference between periodicals and no periodicals boundaries](diffmaterr.png)

The accuracy of the computation has been increased 

Return to the differentiation matrix page [diffmat.m]()

%}