/**
# Miscalleneous functions

## Norm calculation 

*/
#ifdef VOF
coord normal (Point point, scalar c) {
/**
A function to rescale normals so that they are unit vectors w.r.t. the 2-norm
(by default, the 1-norm is adopted for efficiency purposes).*/
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}


/**
A function to compute 2-norm normal in every cell. */

void compute_normal (scalar f, vector normal_vector) {
  foreach() {
    coord n = normal (point, f);
    foreach_dimension() 
      normal_vector.x[] = n.x;
  }
  boundary((scalar*){normal_vector});
}
#endif

/**
## Norm-2

$|\cdot|_2$ norm
*/

double norm2 (coord u) {
  double sum = 0;
  foreach_dimension()
    sum += sq(u.x);
  return sqrt(sum);

}

/**
## Probe
Returns the value contained in a cell
*/

double capteur(Point point, scalar s){
  return s[];
}


/**
#Interpolations functions

* 1D linear
* 2D bilinear with ENO correction
* 2D bicubic

*/

/**
## Interpolation stencil

Used to choose which stencil will be used for interpolation.
*/
void InterpStencil (coord p, int Stencil[]){
  if(p.x<0){Stencil[0] = -1;} 
  else{Stencil[0] = 0;} 

  if(p.y<0) {Stencil[1] = -1;}
  else {Stencil[1] = 0;}
}


/**
## 1D - linear interpolation
*/
double linearInterpolation(double x1, double f_x1, double x2, double f_x2, 
  double x)
{
  double result = (x - x1)/(x2-x1)*f_x2  + (x2-x)/(x2-x1)*f_x1;
  return result;
}

/**
##Bilinear with quadratic correction
*/
#if dimension == 2
double mybilin (Point point, scalar s, int Arr[], coord p_interp,
  double coeff[]){

  double dx = p_interp.x - Arr[0];
  double dy = p_interp.y - Arr[1];
  double dxdy = dx*dy;
  // f11
  coeff[0] = 1. - dx - dy + dxdy;
  // f21
  coeff[1] = dx - dxdy;
  // f12 
  coeff[2] = dy - dxdy;
  // f22
  coeff[3] = dxdy;

//quadratic correction. from Min and Gibou 2006.
  
  double phixx0 = s[-2,0] - 2*s[-1,0] + s[];
  dx = phixx0*(dx)*(1.-dx)/2.;

  phixx0 = s[0,-1] - 2*s[] + s[0,1];
  dy= phixx0*(dy)*(1.-dy)/2.;
  // double phixx0 = s[Arr[0]-1,Arr[1]] - 2*s[Arr[0],Arr[1]] 
    // + s[Arr[0]+1,Arr[1]];
  // double phixx1 = s[Arr[0],Arr[1]] - 2*s[Arr[0]+1,Arr[1]] 
    // + s[Arr[0]+2,Arr[1]];
  // double phixx2 = s[Arr[0]-1,Arr[1]+1] - 2*s[Arr[0],Arr[1]+1] 
    // + s[Arr[0]+1,Arr[1]+1];
  // double phixx3 = s[Arr[0],Arr[1]+1] - 2*s[Arr[0]+1,Arr[1]+1] 
    // + s[Arr[0]+2,Arr[1]+1];
  // double phixval[4] = {phixx0,phixx1,phixx2,phixx3};  
  // fprintf(stderr, "%g %g %g %g\n", phixx0, phixx1, phixx2, phixx3);
  // dx = phixval[minlist(phixval)]*(dx)*(1.-dx)/2.;
  // fprintf(stderr, "%g\n", phixval[minlist(phixval)]);


  // phixx0 = s[Arr[0],Arr[1]-1] - 2*s[Arr[0],Arr[1]] 
  //   + s[Arr[0],Arr[1]+1];
  // phixx1 = s[Arr[0],Arr[1]] - 2*s[Arr[0],Arr[1]+1] 
  //   + s[Arr[0],Arr[1]+2];
  // phixx2 = s[Arr[0]+1,Arr[1]-1] - 2*s[Arr[0]+1,Arr[1]] 
  //   + s[Arr[0]+1,Arr[1]+1];
  // phixx3 = s[Arr[0]+1,Arr[1]] - 2*s[Arr[0]+1,Arr[1]+1] 
  //   + s[Arr[0]+1,Arr[1]+2];

  // fprintf(stderr, "%g %g %g %g\n", phixx0, phixx1, phixx2, phixx3);
  // dy = phixval[minlist(phixval)]*(dy)*(1.-dy)/2.;
  // fprintf(stderr, "%g\n", phixval[minlist(phixval)]);
  // fprintf(stderr, "####\n" );
  // exit(1);

  return  s[Arr[0],Arr[1]]*coeff[0] + 
          s[Arr[0]+1,Arr[1]]*coeff[1] + 
          s[Arr[0],Arr[1]+1]*coeff[2] + 
          s[Arr[0]+1,Arr[1]+1]*coeff[3];
          // -dx-dy;
}
#endif

/**
#Biquadratic interpolation on a Cartesian grid.

We want to calculate a matrix $A$ such that:
$$\begin{aligned}
f(x,y) &= \sum_{i=0}^2\sum_{j=0}^2a_{ij}x^iy^j \\
       &= \left[\begin{array}{cccc} 
x^2 & x & 1 \end{array}\right]
\left[\begin{array}{cccc} 
a_{2,2} & a_{2,1} & a_{2,0}\\
a_{1,2} & a_{1,1} & a_{1,0}\\
a_{0,2} & a_{0,1} & a_{0,0}\\
\end{array}\right]
\left[\begin{array}{c} 
y^2\\
y\\
1
\end{array}\right]\\
&= X A Y^T
\end{aligned}
$$
Given a set of $3\times 3$ known data points, we have:
$$
\begin{aligned}
F &= BAB^T\\
\left[\begin{array}{ccc}
f(-1,-1) & f(-1,0) & f(-1,1)\\
f(0,-1) & f(0,0) & f(0,1)   \\
f(1,-1) & f(1,0) & f(1,1)   \\
  \end{array}\right]_{F}
  &= 
  \left[\begin{array}{ccc}
  (-1)^2 & -1 & 1\\
  0^2 & 0 & 1\\
  1^2 & 1 & 1\\
  \end{array}\right]_{B}
  \left[\begin{array}{ccc}
  a_{2,2} & a_{2,1} & a_{2,0}\\
  a_{1,2} & a_{1,1} & a_{1,0}\\
  a_{0,2} & a_{0,1} & a_{0,0}\\
  \end{array}\right]_A
  \left[\begin{array}{ccc}
  (-1)^2 & 0 & (1)^2\\
  -1     & 0 & 1    \\
  1      & 1 & 1    \\
  \end{array}\right]_{B^T}
  \end{aligned}
$$
Therefore we have:
$$
\begin{aligned}
F = BAB^T \Rightarrow A &= B^{-1} F \left(B^{T}\right)^{-1}\\
  &= B^{-1} F \left(B^{-1}\right)^{T}\end{aligned}
$$
after some calculations, one gets:
$$
B = \left[\begin{array}{ccc}
 1 & -1 & 1\\
 0 & 0 & 1 \\
 1 & 1 & 1 \\
\end{array}\right] \Rightarrow
B^{-1} = \left[\begin{array}{ccc}
0.5  & -1 & 0.5\\
-0.5 &  0 & 0.5\\
   0 &  1 &   0\\
\end{array}\right]
$$
*/
static inline double mybiquadratic(Point point, scalar s, coord p, int offset){
  double dx = p.x;
  double dy = p.y;
  double X[3] = {dx*dx,dx,1.};
  double Y[3] = {dy*dy,dy,1.};

double Bm1[9] = {0.5, -1, 0.5,
                -0.5,  0, 0.5,
                  0,  1,   0};
double Bm1T[9]= {0.5,-0.5, 0,
                  -1,  0,  1,
                 0.5,  0.5, 0};
double XX[3],YY[3];
  for (int kk=0;kk<=2;kk++){
    XX[kk] = Bm1[kk]*X[0] + Bm1[3+kk]*X[1] + 
      Bm1[6+kk]*X[2];
    YY[kk] = Bm1T[kk*3]*Y[0] + Bm1T[kk*3+1]*Y[1] +
      Bm1T[kk*3+2]*Y[2];
  }


  for (int kk = 0;kk<=2;kk++){
    X[kk] = s[-1,-1+kk,offset]*XX[0] + 
            s[ 0,-1+kk,offset]*XX[1] +
            s[+1,-1+kk,offset]*XX[2];
  }
  return X[0]*YY[0]+X[1]*YY[1]+X[2]*YY[2];
}

// taken from Ghigo

double mytriquadratic(Point point, scalar s, coord p){
/**
We do a triple biquadratic interpolation and then quadratically interpolate the
result.
*/
  double dz = p.z;
  double arr[3];
  arr[0] = mybiquadratic(point, s, p, -1);
  arr[1] = mybiquadratic(point, s, p, 0);
  arr[2] = mybiquadratic(point, s, p, 1);
  return quadratic(dz,arr[0],arr[1],arr[2]);
}

/**
##Bicubic interpolation on a Cartesian grid

Taken from [this link](https://www.ece.mcmaster.ca/~xwu/interp_1.pdf) (slide
16-17).

We want to calculate a matrix $A$ such that:
$$\begin{aligned}
f(x,y) &= \sum_{i=0}^3\sum_{j=0}^3a_{ij}x^iy^j \\
       &= \left[\begin{array}{cccc} 
x^3 & x^2 & x & 1 \end{array}\right]
\left[\begin{array}{cccc} 
a_{3,3} & a_{3,2} & a_{3,1} & a_{3,0}\\
a_{2,3} & a_{2,2} & a_{2,1} & a_{2,0}\\
a_{1,3} & a_{1,2} & a_{1,1} & a_{1,0}\\
a_{0,3} & a_{0,2} & a_{0,1} & a_{0,0}\\
\end{array}\right]
\left[\begin{array}{c} 
y^3\\
y^2\\
y\\
1
\end{array}\right]\\
&= X A Y^T
\end{aligned}
$$
Given a set of $4\times 4$ known data points, we have:
$$
\begin{aligned}
F &= BAB^T\\
\left[\begin{array}{cccc}
f(-1,-1) & f(-1,0) & f(-1,1) & f(-1,2)\\
f(0,-1) & f(0,0) & f(0,1) & f(0,2)\\
f(1,-1) & f(1,0) & f(1,1) & f(1,2)\\
f(2,-1) & f(2,0) & f(2,1) & f(2,2)\\
  \end{array}\right]_{F}
  &= 
  \left[\begin{array}{cccc}
  (-1)^3 & (-1)^2 & -1 & 1\\
  0^3 & 0^2 & 0 & 1\\
  1^3 & 1^2 & 1 & 1\\
  2^3 & 2^2 & 2 & 1\\
  \end{array}\right]_{B}
  \left[\begin{array}{cccc}
  a_{3,3} & a_{3,2} & a_{3,1} & a_{3,0}\\
  a_{2,3} & a_{2,2} & a_{2,1} & a_{2,0}\\
  a_{1,3} & a_{1,2} & a_{1,1} & a_{1,0}\\
  a_{0,3} & a_{0,2} & a_{0,1} & a_{0,0}\\
  \end{array}\right]_A
  \left[\begin{array}{cccc}
  (-1)^3 & 0 & (1)^3 & (2)^3\\
  (-1)^2 & 0 & (1)^2 & (2)^2\\
  -1 & 0 & 1 & 2\\
  1 & 1 & 1 & 1\\
  \end{array}\right]_{B^T}
  \end{aligned}
$$
Therefore we have:
$$
\begin{aligned}
F = BAB^T \Rightarrow A &= B^{-1} F \left(B^{T}\right)^{-1}\\
  &= B^{-1} F \left(B^{-1}\right)^{T}\end{aligned}
$$
after some calculations, one gets:
$$
B = \left[\begin{array}{cccc}
-1 & 1 & -1 & 1\\
0 & 0 & 0 & 1 \\
1 & 1 & 1 & 1 \\
8 & 4 & 2 & 1 \\
\end{array}\right] \Rightarrow
B^{-1} = \left[\begin{array}{cccc}
-1/6 & 1/2 & -1/2 & 1/6\\
1/2 & -1 & 1/2 & 0 \\
-1/3 & -1/2 & 1 & -1/6 \\
0 & 1 & 0 & 0 \\
\end{array}\right]
$$
Now we can do the bicubic interpolation:
$$
f(x,y) = X B^{-1} F \left(B^{-1}\right)^T Y^T = \mathcal{X} F \mathcal{Y}
$$

A simple test case can be seen [here](http://basilisk.fr/sandbox/alimare/1_test_cases/test_bicubic.c).
*/

double bicubic(Point point , scalar s, int Arr[], coord p){
  double dx = p.x - Arr[0];
  double dy = p.y - Arr[1];
  double X[4] = {dx*dx*dx,dx*dx,dx,1.};
  double Y[4] = {dy*dy*dy,dy*dy,dy,1.};


/**
$B^{-1}$ and $(B^{-1})^T$ are hardcoded.
*/
  double Bm1[16] = {-1./6, 1/2., -1/2., 1/6.,
                    1/2., -1   ,  1/2., 0   ,
                    -1./3, -1/2. ,   1  , -1./6,
                    0   ,  1   ,   0  ,  0  };

  double Bm1T[16] = {-1/6., 1/2., -1./3,  0.,
                      1/2.,   -1, -1/2.,  1.,
                     -1/2., 1/2.,    1.,  0.,
                      1/6.,    0, -1./6,  0};
/**
We calculate $\mathcal{X}$ and $\mathcal{Y}$.
*/
  double XX[4],YY[4];
  for (int kk=0;kk<=3;kk++){
    XX[kk] = Bm1[kk]*X[0] + Bm1[4+kk]*X[1] + 
      Bm1[8+kk]*X[2] + Bm1[12+kk]*X[3];
    YY[kk] = Bm1T[kk*4]*Y[0] + Bm1T[kk*4+1]*Y[1] +
      Bm1T[kk*4+2]*Y[2] + Bm1T[kk*4+3]*Y[3];
  }


  for (int kk = 0;kk<=3;kk++){
    X[kk] = s[Arr[0]-1,Arr[1]-1+kk]*XX[0] + 
            s[Arr[0]  ,Arr[1]-1+kk]*XX[1] +
            s[Arr[0]+1,Arr[1]-1+kk]*XX[2] +
            s[Arr[0]+2,Arr[1]-1+kk]*XX[3];
  }

  return X[0]*YY[0]+X[1]*YY[1]+X[2]*YY[2]+X[3]*YY[3];
}

/**
## Cell to node interpolation

Bilinear with quadratic correction (ENO) or bicubic interpolation.
*/
void cell2node(scalar cell_scalar, vertex scalar node_scalar){
  foreach_vertex(){
    coord p_interp = {-0.5, -0.5, -0.5};

    if((point.i-1) > (1 << grid-> depth)){
      p_interp.x = 0.5;
    }

    if((point.j-1) > (1 << grid-> depth)){
      p_interp.y = 0.5;
    }
    
#ifdef BICUBIC    
    assert(dimension == 2); // tricubic not coded so far.
    int Stencil[2] = {-1,-1};
    node_scalar[] = bicubic(point , cell_scalar, Stencil, p_interp);
#elif QUADRATIC
#if dimension == 2
    node_scalar[] = mybiquadratic(point , cell_scalar, p_interp,0);
#else
    if((point.k-1) > (1 << grid-> depth)){
      p_interp.z = 0.5;
    }
    node_scalar[] = mytriquadratic(point , cell_scalar, p_interp);
#endif
#else
    assert(dimension ==2);
    double coeff[4];
    int Stencil[2] = {-1,-1};
    node_scalar[] = mybilin(point , cell_scalar, Stencil, p_interp, coeff);
#endif
  }
  boundary ({node_scalar});
  restriction({node_scalar});

}
