/** # This file include mathematical operator for coord struct */
/** declaration of the null coord*/

const coord Coord_nul = {0};
const mat3 tens_nul = {{0}};
#define Identity\
          {{1,0,0},\
          {0,1,0},\
          {0,0,1}}
const mat3 I = Identity;

/** and operator on coord */
coord and_coord(coord a,coord b){
  coord c = {(int)a.x && (int)b.x,(int)a.y &&(int) b.y,(int)a.z &&(int) b.z};
  return c;
}
coord add_coord(coord a,coord b){
  #if dimension == 2
  coord c = {a.x+b.x, a.y+b.y};
  #elif dimension == 3
  coord c = {a.x+b.x, a.y+b.y, a.z+b.z};
  #endif
  return c;
}
coord diff_coord(coord a,coord b){
  #if dimension == 2
  coord c = {a.x-b.x,a.y-b.y};
  #elif dimension == 3
  coord c = {a.x-b.x, a.y-b.y, a.z-b.z};
  #endif
  return c;
}
coord mult_coord(coord a,double b){
  #if dimension == 2
  coord c = {a.x*b,  a.y*b};
  #elif dimension == 3
  coord c = {a.x*b,  a.y*b,  a.z*b};
  #endif
  return c;
} 
coord div_coord(coord a,double b){
  #if dimension == 2
  coord c = {a.x/b,  a.y/b};
  #elif dimension == 3
  coord c = {a.x/b,  a.y/b,  a.z/b};
  #endif
  return c;
}
mat3 div_tens(mat3 a,double b){
  einstein_sum(i,j){
     a.i.j /= b;
  }
  return a;
}
mat3 mult_tens(mat3 a,double b){
  einstein_sum(i,j){
     a.i.j *= b;
  }
  return a;
}
mat3 add_tens(mat3 A,mat3 B){
  einstein_sum(i,j){
     A.i.j += B.i.j;
  }
  return A;
}
double tens_norm(mat3 A){
  double norm = 0;
  einstein_sum(k,l){
     norm = A.k.l * A.k.l;
  }
  norm = sqrt(norm);
  return norm;
}
double normL2_coord(coord a){
  #if dimension == 2
  double c = sqrt(pow(a.x,2)+pow(a.y,2));
  #elif dimension == 3
  double c = sqrt(pow(a.x,2)+pow(a.y,2)+pow(a.z,2));
  #endif
  return c;
}
double normL2_coord_sq(coord a){
  #if dimension == 2
  double c = sq(a.x)+sq(a.y);
  #elif dimension == 3
  double c = sq(a.x)+sq(a.y)+sq(a.z);
  #endif
  return c;
}
double normL2_ten(mat3 a){
  #if dimension == 2
  double c = sqrt(pow(a.x.x,2)+pow(a.y.y,2)+pow(a.x.y,2)+pow(a.y.x,2));
  #elif dimension == 3
  double c = sqrt(pow(a.x.x,2)+pow(a.y.y,2)+pow(a.z.z,2)+pow(a.x.y,2)+pow(a.y.x,2)+pow(a.x.z,2)+pow(a.z.x,2)+pow(a.z.y,2)+pow(a.y.z,2));
  #endif
  return c;
}
//only in 3D for the cross product 
#if dimension == 3
coord cross_coord(coord a,coord b){
  coord c = {a.y*b.z-a.z*b.y, 
            a.z*b.x-a.x*b.z, 
            a.x*b.y-a.y*b.x};
#elif dimension == 2
double cross_coord(coord a,coord b){
  double c = a.x*b.y-a.y*b.x;
#endif
  return c;
}

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

/** Set the coordinate of a coord at the positive side of the domain */
coord POS_perio(coord pos,coord per){
  foreach_dimension(){
    if(per.x) pos.x = pos.x > 0 ? pos.x : pos.x+Ls;
  }
  return pos;
}
/** compute the periodic distance between two coord */
coord dist_perio_coord(coord a,coord b){
  foreach_dimension(){
    if(a.x>0) b.x = (b.x<a.x-Ls/2)*Ls+b.x;
    else b.x = (b.x>a.x+Ls/2)*(-Ls)+b.x;
  }
  return diff_coord(a,b);
}


double dist_perio(coord a,coord b){
  foreach_dimension(){
    if(a.x>0) b.x = (b.x<a.x-Ls/2)*Ls+b.x;
    else b.x = (b.x>a.x+Ls/2)*(-Ls)+b.x;
  }
  return normL2_coord(diff_coord(a,b));
}
/** compute the tranpsoe of a tensor */ 
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
  B. i.j = 0.5*(A.i.j + A.j.i);
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
/** compute the symetric part of a tensor */