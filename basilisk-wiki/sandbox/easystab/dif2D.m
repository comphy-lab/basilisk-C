%{
In this code, we build the 2D differentiation matrices of a rectangular domain. 

To learn how to use this function, please see [diffmat_2D_dif2D.m](). 

%}

function [D,l,X,Y,Z,I]=dif2D(d,x,y)

Nx=length(x);
Ny=length(y);
NN=Nx*Ny;

%{
# 2D differentiation matrices
To build the 2D differentiation matrices from the 1D ones, we use the *kron* operator. Please see [diffmat_2D.m#d-matrices-1]().
%}

% differentiation
D.x=kron(d.x,speye(Ny));
D.xx=kron(d.xx,speye(Ny));
D.y=kron(speye(Nx),d.y);
D.yy=kron(speye(Nx),d.yy);

%{
# Integration weights in 2D
To learn about the integration weights, please see [integration_2D.m]().
%}

% integration
D.w=d.wy(:)*d.wx(:)'; 
D.wx=ones(Ny,1)*d.wx(:)'; 
D.wy=d.wy(:)*ones(1,Nx); 

D.w=reshape(D.w,Ny,Nx);
D.wx=reshape(D.wx,Ny,Nx);
D.wy=reshape(D.wy,Ny,Nx);

%{
# Location vectors
To impose the boundary conditions we need to know where in the vector representation of the 2D variable are the different boundaries and corners. For this we build here a few location vectors stored as the fields of the structure *l*. To learn about the use of location vectors, please see [pedagogy#location-vectors-and-matrices](). To learn how the array variable is transformed into a vector, see [diffmat_2D.m](); there you will as well understand the structure of the mesh in 2D.

Here we discuss the different fields of *l* and how they are built. You see below the structure of the mesh:

![The different boundaries of the 2D mesh](locationfig.png)

We start by building the long vector of all the indices *1:NN*, where *NN=Nx*Ny* is the total number of degrees of freedom of the 2D mesh. This vector we transform into the shape of the mesh with the function *reshape*. We store this in the array *dom* (standing for "domain"). Once this done, we have the link between the positions in the array and the locations in the vector. Since the indices in the array increase with position, we have that:

* the first line of the array *dom* corresponds to the bottom boundary cells, 
* the last line of the array *dom* corresponds to the top boundary cells,
* the first column of the array *dom* corresponds to the left boundary cells,
* the last column of the array *dom* corresponds to the right boundary cells.
* The corners of the array *dom* correspond to the corner cells.

This means that when you display the array as Octave/Matlab does on screen, you see your mesh upside-down. To avoid this we could have chosen to have the $y$ position of cells decreasing with index instead of having them increase with index, but I prefered to keep the same convention for $x$ and $y$. 

In this generic function, I chose to store the side boundaries without the corners. Depending on your physical system, you may chose to include the corners in the location vectors for instance in the right and left boundaries, and keep the top and bottom boundaries without corners, see for instance [???]()

%}

% locations
dom=reshape(1:NN,Ny,Nx);
right=dom(2:Ny-1,end);    l.right=right(:);
left=dom(2:Ny-1,1);       l.left=left(:);
top=dom(end,2:Nx-1);      l.top=top(:);
bot=dom(1,2:Nx-1);    l.bot=bot(:);
cor=dom([1 end],[1 end]);     l.cor=cor(:);

l.cbl=cor(1); % bottom left corner
l.ctl=cor(2); % top left corner
l.cbr=cor(3); % bottom right corner
l.ctr=cor(4); % top right corner

%{
To learn about meshgrid, please see [diffmat_2D.m#meshgrid]()
%}

[X,Y]=meshgrid(x,y);

%{
It is usefull when building the system matrices, to have a zero matrix and an identity matrix. Here we build them with the size of the total number of degrees of freedom *NN*. These are sparse matrices, meaning that in memory is only kept the nonzero elements of the matrices and their location in the matrices. In 2D and 3D it is good to use sparse matrices: use less memory and faster computations because we do not need to multiply many zeros. For the issues on computational efficiency, please see [pedagogy#computational-efficiency]().
%}

Z=spalloc(NN,NN,0); I=speye(NN); 
