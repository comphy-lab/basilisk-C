#include "common.h"

#define sym_tensor_equal(a,b) a.x.x[] = b.x.x[]; a.x.y[] = b.x.y[]; a.x.z[] = b.x.z[]; \
                          a.y.y[] = b.y.y[]; a.y.z[] = b.y.z[]; \
                          a.z.z[] = b.z.z[]; 

#define sym_tensor_const_mult(a,b,c) a.x.x[] = b * c.x.x[]; a.x.y[] = b * c.x.y[]; a.x.z[] = b * c.x.z[]; \
                                 a.y.y[] = b * c.y.y[]; a.y.z[] = b * c.y.z[]; \
                                 a.z.z[] = b * c.z.z[]; 

#define sym_tensor_sum_product(a,b) (a.x.x[]*b.x.x[] + a.y.y[]*b.y.y[] + a.z.z[]*b.z.z[] + 2*(a.x.y[]*b.x.y[] + a.x.z[]*b.x.z[] + \
				a.y.z[]*b.y.z[]))

// a = b * c
#define vector_scalar_mult(a,b,c) a.x[] = b[] * c.x[]; a.y[] = b[] * c.y[]; a.z[] = b[] * c.z[];


// a = b * c
#define sym_tensor_scalar_mult(a,b,c) a.x.x[] = b[] * c.x.x[]; a.x.y[] = b[] * c.x.y[]; a.x.z[] = b[] * c.x.z[]; \
                                 a.y.y[] = b[] * c.y.y[]; a.y.z[] = b[] * c.y.z[]; \
                                 a.z.z[] = b[] * c.z.z[]; 

// a = b + c
#define vector_add(a,b,c) a.x[] = b.x[] + c.x[]; a.y[] = b.y[] + c.y[]; a.z[] = b.z[] + c.z[];

// a = b - c
#define vector_subtract(a,b,c) a.x[] = b.x[] - c.x[]; a.y[] = b.y[] - c.y[]; a.z[] = b.z[] - c.z[];

// a = b + c
#define sym_tensor_add(a,b,c) a.x.x[] = b.x.x[] + c.x.x[]; a.x.y[] = b.x.y[] + c.x.y[]; a.x.z[] = b.x.z[] + c.x.z[]; \
                          a.y.y[] = b.y.y[] + c.y.y[]; a.y.z[] = b.y.z[] + c.y.z[]; \
                          a.z.z[] = b.z.z[] + c.z.z[];

// a = b - c
#define sym_tensor_subtract(a,b,c) a.x.x[] = b.x.x[] - c.x.x[]; a.x.y[] = b.x.y[] - c.x.y[]; a.x.z[] = b.x.z[] - c.x.z[]; \
                               a.y.y[] = b.y.y[] - c.y.y[]; a.y.z[] = b.y.z[] - c.y.z[]; \
                               a.z.z[] = b.z.z[] - c.z.z[];



void filter_sym_tensor(tensor a, tensor b) {
  boundary((scalar *) {b});
  foreach(){
      double fxx=0., fxy=0., fxz=0., fyy=0., fyz=0., fzz = 0.;
      foreach_neighbor(1) {
        fxx += b.x.x[];
        fxy += b.x.y[];
        fxz += b.x.z[];
        fyy += b.y.y[];
        fyz += b.y.z[];
        fzz += b.z.z[];
      }
      a.x.x[] = fxx/27.;
      a.x.y[] = fxy/27.;
      a.x.z[] = fxz/27.;
      a.y.y[] = fyy/27.;
      a.y.z[] = fyz/27.;
      a.z.z[] = fzz/27.;
  }
  boundary((scalar *){a});
}


void filter_tensor(tensor a, tensor b) {
  boundary((scalar *) {b});
  foreach(){
      double fxx=0., fxy=0., fxz=0., fyx=0., fyy=0., fyz=0., fzx=0., fzy=0., fzz = 0.;
      foreach_neighbor(1) {
        fxx += b.x.x[];
        fxy += b.x.y[];
        fxz += b.x.z[];
        fyx += b.y.x[];
        fyy += b.y.y[];
        fyz += b.y.z[];
        fzx += b.z.x[];
        fzy += b.z.y[];
        fzz += b.z.z[];
      }
      a.x.x[] = fxx/27.;
      a.x.y[] = fxy/27.;
      a.x.z[] = fxz/27.;
      a.y.x[] = fyx/27.;
      a.y.y[] = fyy/27.;
      a.y.z[] = fyz/27.;
      a.z.x[] = fzx/27.;
      a.z.y[] = fzy/27.;
      a.z.z[] = fzz/27.;
  }
  boundary((scalar *){a});
}



void filter_vector(vector a, vector b) {
  boundary((scalar *) {b});
  foreach(){
      double fx=0., fy=0., fz=0.;
      foreach_neighbor(1) {
        fx += b.x[];
        fy += b.y[];
        fz += b.z[];
      }
      a.x[] = fx/27.;
      a.y[] = fy/27.;
      a.z[] = fz/27.;
  }
  boundary((scalar *){a});
}

void filter_scalar(scalar a, scalar b) {
  boundary((scalar *) {b});
  foreach(){
      double f=0.;
      foreach_neighbor(1) {
        f += b[];
      }
      a[] = f/27.;
  }
  boundary({a});
}



                             
void eddyviscosity(vector u, scalar Evis){	

  symmetric tensor M[];
  symmetric tensor S[];
  symmetric tensor SS[]; // abs(S)*S
  symmetric tensor SSf[]; // abs(S)*S filtered
  tensor ug[]; // grad(u)

  vector u_f[]; // filtered

  scalar Snorm[];
  scalar Cs_sq[];

  double filter_ratio_sq = 6.; // explicit filter to grid filter ratio squared

  foreach(){
    // compute dui_dxj
    ug.x.x[] = (u.x[1, 0, 0]-u.x[-1, 0, 0])/2/Delta;  
    ug.x.y[] = (u.x[0, 1, 0]-u.x[0, -1, 0])/2/Delta;   
    ug.x.z[] = (u.x[0, 0, 1]-u.x[0, 0, -1])/2/Delta; 
    ug.y.x[] = (u.y[1, 0, 0]-u.y[-1, 0, 0])/2/Delta;  
    ug.y.y[] = (u.y[0, 1, 0]-u.y[0, -1, 0])/2/Delta;  
    ug.y.z[] = (u.y[0, 0, 1]-u.y[0, 0, -1])/2/Delta; 
    ug.z.x[] = (u.z[1, 0, 0]-u.z[-1, 0, 0])/2/Delta; 
    ug.z.y[] = (u.z[0, 1, 0]-u.z[0, -1, 0])/2/Delta;  
    ug.z.z[] = (u.z[0, 0, 1]-u.z[0, 0, -1])/2/Delta;

    double S_tr = (ug.x.x[] + ug.y.y[] + ug.z.z[])/3.; // remove trace; divergence free
    S.x.x[] = ug.x.x[] - S_tr;
    S.y.y[] = ug.y.y[] - S_tr;
    S.z.z[] = ug.z.z[] - S_tr;
    S.x.y[] = 0.5*(ug.x.y[] + ug.y.x[]);
    S.x.z[] = 0.5*(ug.x.z[] + ug.z.x[]);
    S.y.z[] = 0.5*(ug.y.z[] + ug.z.y[]);

    // make copy of |S|S in M
    Snorm[] = sqrt(2*(S.x.x[]*S.x.x[] + S.y.y[]*S.y.y[] + \
                      S.z.z[]*S.z.z[] + 2*(S.x.y[]*S.x.y[] + S.x.z[]*S.x.z[] + S.y.z[]*S.y.z[])));
    sym_tensor_scalar_mult(M, Snorm, S);

    

    // copy to |S|S
    sym_tensor_scalar_mult(SS, Snorm, S);
  }

  // Filter |S|S and u

  filter_sym_tensor(SSf, SS);
  filter_vector(u_f, u);

  // compute M
  foreach() {
    sym_tensor_subtract(M, SSf, filter_ratio_sq*SS);
    sym_tensor_const_mult(M, 2*Delta*Delta, M);
  }

  // compute filtered u_i*u_j.  We are done with S, so use for storage.
  foreach(){
    S.x.x[] = u.x[]*u.x[];
    S.x.y[] = u.x[]*u.y[];
    S.x.z[] = u.x[]*u.z[];
    S.y.y[] = u.y[]*u.y[];
    S.y.z[] = u.y[]*u.z[];
    S.z.z[] = u.z[]*u.z[];
  }

  // filter u_i*u_j
  filter_tensor(SSf, S);

  // and store L_ij in S
  foreach(){
    S.x.x[] = SSf.x.x[] - u_f.x[]*u_f.x[];
    S.x.y[] = SSf.x.y[] - u_f.x[]*u_f.y[];
    S.x.z[] = SSf.x.z[] - u_f.x[]*u_f.z[];
    S.y.y[] = SSf.y.y[] - u_f.y[]*u_f.y[];
    S.y.z[] = SSf.y.z[] - u_f.y[]*u_f.z[];
    S.z.z[] = SSf.z.z[] - u_f.z[]*u_f.z[];
  }


  // repeatedly smooth Lij and M - common procedure - makes for more stable solution

  for (int i = 0; i < 3; i++) {
    filter_tensor(SS, M);
    foreach() {
      sym_tensor_equal(M, SS);
    }
  }

  for (int i = 0; i < 3; i++) {
    filter_tensor(SS, S);
    foreach() {
      sym_tensor_equal(S, SS);
    }
  }

  foreach(){       
    Cs_sq[] = sym_tensor_sum_product(S, M)/(sym_tensor_sum_product(M, M) + 1E-16);
    if (Cs_sq[]<0) {
      Cs_sq[] = 0;  // restrict to positive values for stability
      }
    Evis[] = Cs_sq[]*Delta*Delta*Snorm[];
  }
  boundary({Evis});

}




