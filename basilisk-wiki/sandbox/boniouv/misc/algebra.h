/** # Definitions and algebra for vectors and tensors*/

/** ## Definition of some usfull macros */
#if dimension == 2
#define POS {x,y} 
#else 
#define POS {x,y,z} 
#endif

const coord coord_null = {0};
const mat3 tens_null = {{0}};
#define Identity\
          {{1,0,0},\
          {0,1,0},\
          {0,0,1}}
const mat3 I = Identity;

/** ## Elementary operations on coord */
coord and_coord(coord a, coord b){
  coord c = {(int)a.x && (int)b.x,(int)a.y &&(int) b.y,(int)a.z &&(int) b.z};
  return c;
}

coord add_coord(coord a, coord b){
  einstein_sum(i){
     a.i += b.i;
  }
  return a;
}

coord diff_coord(coord a, coord b){
  einstein_sum(i){
     a.i -= b.i;
  }
  return a;
}

coord mult_coord(coord a, double d){
  einstein_sum(i){
     a.i *= d;
  }
  return a;
} 

coord div_coord(coord a, double d){
  einstein_sum(i){
     a.i /= d;
  }
  return a;
}

double normL2_coord(coord A){
  double c = 0;
  einstein_sum(i){
     c = A.i*A.i;
  }
  return sqrt(c);
}

//only in 3D for the cross product 
#if dimension == 2
double cross_coord(coord a,coord b){
  double c = a.x*b.y-a.y*b.x;
#elif dimension == 3
coord cross_coord(coord a,coord b){
  coord c = {a.y*b.z-a.z*b.y, 
             a.z*b.x-a.x*b.z, 
             a.x*b.y-a.y*b.x};
#endif
  return c;
}

/** ## Elementary operations on mat3 */
mat3 add_tens(mat3 A,mat3 B){
  einstein_sum(i,j){
     A.i.j += B.i.j;
  }
  return A;
}

mat3 diff_tens(mat3 A,mat3 B){
  einstein_sum(i,j){
     A.i.j -= B.i.j;
  }
  return A;
}

mat3 mult_tens(mat3 A,double d){
  einstein_sum(i,j){
     A.i.j *= d;
  }
  return A;
}

mat3 div_tens(mat3 A,double d){
  einstein_sum(i,j){
     A.i.j /= d;
  }
  return A;
}

double normL2_ten(mat3 A){
  double c = 0;
  einstein_sum(k,l){
     c = A.k.l * A.k.l;
  }
  return sqrt(c);
}

mat3 transpsoe(mat3 A){
  mat3 B;
  einstein_sum(i,j){
    B.i.j = A.j.i;
  }
  return B;
}

mat3 sym(mat3 A){
  mat3 B;
  einstein_sum(i,j){
    B.i.j = 0.5*(A.i.j + A.j.i);
  }
  return B;
}

mat3 coord_prod(coord A,coord B){
  mat3 C;
  einstein_sum(i,j){
    C.i.j = A.i*B.j;
  }
  return C;
}

/** ## Eigen analysis */
coord EigenValue(mat3 A){
  double Delta = pow(A.x.x,2.) - 2.*A.x.x*A.y.y+pow(A.y.y,2.) +4.*A.x.y*A.y.x;
  coord Eig = {0.,0.};
  if(Delta > 0){
    Eig.x = 1./2.*(  A.x.x + A.y.y - sqrt( Delta ));
    Eig.y = 1./2.*(  A.x.x + A.y.y + sqrt( Delta ));
  }
  return Eig;
}

mat3 EigenVectors(mat3 A){
  mat3 V;
  if(A.y.x != 0){
    V.x.x = -1./(2.*A.y.x)*( - A.x.x + A.y.y + sqrt( pow(A.x.x,2.) - 2.*A.x.x*A.y.y+pow(A.y.y,2.) +4.*A.x.y*A.y.x ) ); 
    V.x.y = 1.;
    V.y.x = -1./(2.*A.y.x)*( - A.x.x + A.y.y - sqrt( pow(A.x.x,2.) - 2*A.x.x*A.y.y+pow(A.y.y,2.) +4.*A.x.y*A.y.x ) );
    V.y.y = 1.;
  }else if((A.x.x-A.y.y)!=0){
    V.x.x = 1.;
    V.x.y = 0.;
    V.y.x = -A.x.y/(A.x.x-A.y.y);
    V.y.y = 1.;
  }else{
    V.x.x = 1.;
    V.x.y = 0.;
    V.y.x = 0.;
    V.y.y = 1.;
  }
  V.x = div_coord(V.x,normL2_coord(V.x));
  V.y = div_coord(V.y,normL2_coord(V.y));
  return V;
}

/** ## Distances and positions in periodic coordinates */
/** Set the coordinate of a coord at the positive side of the domain */
coord POS_perio(coord pos,coord per){
  foreach_dimension(){
    if(per.x) pos.x = pos.x>0?pos.x:pos.x+L0;
  }
  return pos;
}

/** compute the periodic distance between two coord */
coord dist_perio_coord(coord a,coord b){
  foreach_dimension(){
    if(a.x>0) b.x = (b.x<a.x-L0/2)*L0+b.x;
    else b.x = (b.x>a.x+L0/2)*(-L0)+b.x;
  }
  return diff_coord(a,b);
}

double dist_perio(coord a,coord b){
  foreach_dimension(){
    if(a.x>0) b.x = (b.x<a.x-L0/2)*L0+b.x;
    else b.x = (b.x>a.x+L0/2)*(-L0)+b.x;
  }
  return normL2_coord(diff_coord(a,b));
}

/** ## Derivative operators */

#if dimension == 2
#define scalar_grad(p,G)\
  {\
    G.x = (p[1,0] - p[-1,0] )/(2.*Delta);\
    G.y = (p[0,1] - p[0,-1] )/(2.*Delta);\
  }

#define scalar_hessian(p,H)\
  {\
    H.x.x = (p[1] + p[-1] - 2*p[])/(Delta*Delta); \
    H.x.y = (p[1,1] - p[-1,1] - p[1,-1] + p[-1,-1])/(2.*Delta*Delta); \
    H.y.x = H.x.y; \
    H.y.y = (p[0,1] + p[0,-1] - 2*p[])/(Delta*Delta); \
  }

#define vec_div(u,D)\
  {\
    D = (u.x[1,0] - u.x[-1,0] )/(2.*Delta) \
      + (u.y[0,1] - u.y[0,-1] )/(2.*Delta);\
  }

#define vec_grad(u,G)\
  {\
    G.x.x =  ( u.x[1,0] - u.x[-1,0] )/(2.*Delta);\
    G.x.y =  ( u.x[0,1] - u.x[0,-1] )/(2.*Delta);\
    G.y.x =  ( u.y[1,0] - u.y[-1,0] )/(2.*Delta);\
    G.y.y =  ( u.y[0,1] - u.y[0,-1] )/(2.*Delta);\
  }

#define vec_lap(u,L)\
  {\
    L.x = (u.x[1,0] + u.x[-1,0] - 2*u.x[]                  \
        + u.x[0,1] + u.x[0,-1] - 2*u.x[] )/(Delta*Delta);\
    L.y = (u.y[1,0] + u.y[-1,0] - 2*u.y[]                  \
        + u.y[0,1] + u.y[0,-1] - 2*u.y[] )/(Delta*Delta);\
  }

#elif dimension == 3

#define scalar_grad(p,G)\
  {\
    G.x = (p[1,0,0] - p[-1,0,0] )/(2.*Delta);\
    G.y = (p[0,1,0] - p[0,-1,0] )/(2.*Delta);\
    G.z = (p[0,0,1] - p[0,0,-1] )/(2.*Delta);\
  }

#define scalar_hessian(p,H)\
  {\
    H.x.x = (p[1] + p[-1] - 2*p[])/(Delta*Delta); \
    H.x.y = (p[1,1] - p[-1,1] - p[1,-1] + p[-1,-1])/(4.*Delta*Delta); \
    H.x.z = (p[1,0,1] - p[-1,0,1] - p[1,0,-1] + p[-1,0,-1])/(4.*Delta*Delta); \
    H.y.x = H.x.y; \
    H.y.y = (p[0,1] + p[0,-1] - 2*p[])/(Delta*Delta); \
    H.y.z = (p[0,1,1] - p[0,-1,1] - p[0,1,-1] + p[0,-1,-1])/(4.*Delta*Delta); \
    H.z.x = H.x.z; \
    H.z.y = H.y.z; \
    H.z.z = (p[0,0,1] + p[0,0,-1] - 2*p[])/(Delta*Delta); \
  }

#define vec_div(u,D)\
  {\
    D = (u.x[1,0,0] - u.x[-1,0,0] )/(2.*Delta) \
      +  (u.y[0,1,0] - u.y[0,-1,0] )/(2.*Delta) \
      +  (u.z[0,0,1] - u.z[0,0,-1] )/(2.*Delta);\
  }

#define vec_grad(u,G)\
  {\
    G.x.x = (u.x[1]     - u.x[-1]    )/(2.*Delta);\
    G.x.y = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);\
    G.x.z = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);\
    G.y.x = (u.y[1]     - u.y[-1]    )/(2.*Delta);\
    G.y.y = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);\
    G.y.z = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);\
    G.z.x = (u.z[1]     - u.z[-1]    )/(2.*Delta);\
    G.z.y = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);\
    G.z.z = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);\
  }

#define vec_lap(u,L)\
  {\
    L.x = (u.x[1] + u.x[-1] - 2*u.x[] \
        + u.x[0,1] + u.x[0,-1] - 2*u.x[] \
        + u.x[0,0,1] + u.x[0,0,-1] - 2*u.x[] )/(Delta*Delta);\
    L.y = (u.y[1] + u.y[-1] - 2*u.y[] \
        + u.y[0,1] + u.y[0,-1] - 2*u.y[] \
        + u.y[0,0,1] + u.y[0,0,-1] - 2*u.y[] )/(Delta*Delta);\
    L.z = (u.z[1] + u.z[-1] - 2*u.z[] \
        + u.z[0,1] + u.z[0,-1] - 2*u.z[] \
        + u.z[0,0,1] + u.z[0,0,-1] - 2*u.z[] )/(Delta*Delta);\
  }
#endif