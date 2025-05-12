/**
# Mesh adaptation and metric-optimisation

This code is a toolbox which provides elements to use mesh adaptation based on 
metric-optimisation. This method consists in the optimisation of the metric 
(*i.e.* the mesh) by solving a constrained optimisation problem : find the mesh 
which minimise the $\mathcal{L}^p$ error for a fixed complexity (*i.e.* a fixed 
number of nodes).

The $L^p$-norm interpolation error is estimated busing the metric-based theory and expresses as $|u - \pi_\mathcal{M} u(x)| = c_n \ 
trace(\mathcal{M}^{-\frac{1}{2}}(x) |H(x)| \mathcal{M}^{-\frac{1}{2}}(x))$, and 
the gobal error is 	
$||u - \pi_\mathcal{M} u||_{\mathcal{L}^p}\ = c_n\ (\int_{\Omega}\ 
(trace(\mathcal{M}^{-\frac{1}{2}}(x)\ |H(x)|\ \mathcal{M}^{-\frac{1}{2}}(x)))^p\ 
dx\ )^{\frac{1}{p}}$, where :

- $\mathcal{M}(x)$ is the metric at the point $x$ ;
- $|H(x)|$ is the "absolute hessian", *i.e.* the hessian reconstructed with the 
absolute values of the eigen values of the hessian ;
- $c_n$ is a constant coefficient which only depend on the dimension of the 
problem for triangular elements, but not for square/cubic elements ;
- $n$ is the dimension ;
- $p$ is the norm.

After solving an optimization problem, one can find that the optimal interpolation error on a quad/octree grid containing square/cubic elements (error on the optimal mesh which minimizes the interpolation error) is $$|u-\pi_\mathcal{M} u|(x) = C_n\, N^{-\frac{2}{n}} \left( \int_\Omega  \left(tr(|H|(x)) \right)^{\frac{np}{2p+n}} dx    \right)^{\frac{2}{n}}  \left(tr(|H|(x)) \right)^{\frac{n}{2p+n}}$$	
, where $N$ is the complexity (number of elements) fixed  for the problem and $C_n=C_2=C_3=1/12$.

This error estimation only depends on the hessian of 
the field $u$ on which is based the mesh adaptation.

For the mesh adaptation, we use $||u - \pi_\mathcal{M} u||_{\mathcal{L}^p, K}$ as a refinement criterion in each element $K$.
In our method, we show the interest to add a maximum level constraint to help the total error remaining local and modelizable by the interpolation error.
To estimate the maximum level constraint, we compute the optimal prefactors predicted by the theory.

Note that, in practice, we allow to compute an error on multiple fields $u_i$ as
$$ E = \sum_i \omega_i ||u_i - \pi_\mathcal{M} u_i||_{\mathcal{L}^p} $$
where $\omega_i$ are user-defined weights.

*/


/**
"[Linear algebra tools](./lin_alg.h)" : define some useful functions and types.
*/

#include "./linalg.h"


/**
# Hessian reconstruction methods

We first need to reconstruct the hessian.
   
## Centered difference

**global variable** : user_hessian

This variable allows to choose the method to reconstruct the hessian 
- 1 for centered difference ;
- 2 for double projection L2 ;
- 3 for Green formula based method.    *Only 2D for now*
- 4 Green simplification.	*Only 2D for now*

*/

int user_hessian=1;

/** Computation of the centered-gradient-based hessian.

*/


#define center_gradient_large(a) center_gradient(a);


void compute_bndry_der (Point point, scalar u, pseudo_v * g){

   
   if (x>=X0+L0-Delta/2.)
      g->x = (u[] - u[-1]) / Delta;
   else if (x<=X0+Delta/2.)
      g->x = (u[1] - u[]) / Delta;
   else if (y>=Y0+L0-Delta/2.)
      g->y = (u[] - u[0,-1]) / Delta;
   else if (y<=Y0+Delta/2.)
      g->y = (u[0,1] - u[]) / Delta;
      
}


   
pseudo_v local_gradient_centered (Point point, scalar u) {

   pseudo_v g;

   foreach_dimension()
      g.x = center_gradient_large(u);

   compute_bndry_der(point,u,&g);
   
   return g;
}



void compute_bndry_hessian (Point point, vector g, pseudo_t * h){

   
   if (x>=X0+L0-Delta/2.){
      h->x.x = (g.x[] - g.x[-1]) / Delta;
      h->y.x = (g.y[] - g.y[-1]) / Delta;
   }
   else if (x<=X0+Delta/2.){
      h->x.x = (g.x[1] - g.x[]) / Delta;
      h->y.x = (g.y[1] - g.y[]) / Delta;
   }
   else if (y>=Y0+L0-Delta/2.){
      h->x.y = (g.x[] - g.x[0,-1]) / Delta;
      h->y.y = (g.y[] - g.y[0,-1]) / Delta;
   }
   else if (y<=Y0+Delta/2.){
      h->x.y = (g.x[0,1] - g.x[]) / Delta;
      h->y.y = (g.y[0,1] - g.y[]) / Delta;
   }

}

pseudo_t compute_hessian_centered_local (Point point, vector g) {

   pseudo_t h; 

   foreach_dimension ()
      h.x = local_gradient_centered (point, g.x);

   //h = local_symmetrize_pseudo_t(h);

   compute_bndry_hessian(point, g, &h);
   
   return h;
}




/**
## Double projection **$L^2$**

[source](https://www.theses.fr/2008PA066622)

### Elements of theory

This method is mathematicaly based on the Clément interpolation. It follows two 
steps : first, it computes the gradient of the considered field $u$, and then 
the hessian as the symetrical gradient of the gradient of $u$.

Let $P_k$ be a node of the mesh $\mathcal{T}$. The ball around $P_k$ constitued 
by the elements which contain $P_k$ is called $S_k$. The elements of $S_k$ are 
called $K_l$. The Clément interpolation allows to compute the gradient on each 
element $K_l$ of $S_k$. Then, the nodal gradient in $P_k$ is reconstructed by 
$\nabla_R u(P_k)\ =\ \frac{\sum_{K_l\ \in\ S_k}\ |K_l|\ \nabla u(K_l)}{|S_k|}$, 
where $|.|$ denotes the volume. Physicaly, the gradient reconstructed with this 
method corresponds to the volume-averaged value of the gradients of the elements 
surronding the current point. 

This procedure is then repreated with each field of the gradient to build the 
hessian : 
$H_R(u) = \frac{H^*(u) + H^*(u)^t}{2}$, where $H^*(u) = \nabla_R (\nabla_R 
u(P_k)^t)$.

### Application to basilisk

The structure of the mesh in basilisk is such that it can be consider that, 
localy, a cell is surrounding by 8 cells (in 2D, and 26 in 3D) with the same 
dimensions. Moreover, the grid is cell-centered, and the cells are squares.
So we want to compute the $\mathcal{L}^2$ gradient (and then the $\mathcal{L}^2$ 
hessian) at the center of each cell. The $S_k$ ball is centered on the current 
cell, and is delimited by the center of its neighbors. That leads naturaly to 
consider "fictive" elements (see image below), on which the elements gradient are calculated, and 
which allow to reconstruct the gradient at the current point. These two steps are simplified analyticaly and below we directly compute the gradient at the current point.

Let $u_{i,j,k}$ be the scalar field considered, and $\Delta$ the length of the cell.

*/



pseudo_v local_gradient_L2 (Point point, scalar u) {

   pseudo_v g;

#if dimension==1
   /** In 1D, this gives the centered gradient 
       $\frac{\partial u}{\partial x} = \frac{u_{i-1} - u_{i+1}}{2\Delta}$.
   */
   g.x = center_gradient_large(u);
#endif
#if dimension==2
   /** In 2D, this gives 
       $\frac{\partial u}{\partial x} = \frac{2u_{i-1,j} + u_{i-1,j-1} + u_{i-1,j+1} 
       - 2u_{i+1,j} - u_{i+1,j-1} - u_{i+1,j+1}}{8\Delta}$.

   */
   double facteur; //normalisation factor used in gradient calculation

   facteur = 1./( 8.*Delta );

   foreach_dimension() {
#if EMBED
      if (fs.x[] && fs.x[1] && fs.x[-1] && fs.x[2]) {
#endif

	 g.x = 2.*u[-1,0] + u[-1,-1] + u[-1,1] - 2.*u[1,0] - u[1,-1] - u[1,1];
	 g.x *= facteur;

#if EMBED
      }
      else // si on est au bord du domaine
         g.x = -center_gradient_large(u);
#endif

   }
#endif
#if dimension==3
   /** In 3D, that gives 
       $\frac{\partial u}{\partial x} = 
      \frac{4u_{i-1,j,k} + 2u_{i-1,j-1,k} + 2u_{i-1,j+1,k} + 2u_{i-1,j,k-1} + 2u_{i-1,j,k+1} 
      + u_{i-1,j-1,k-1} + u_{i-1,j+1,k-1} + u_{i-1,j-1,k+1} + u_{i-1,j+1,k+1}
      - 4u_{i+1,j,k} - 2u_{i+1,j-1,k} - 2u_{i+1,j+1,k} - 2u_{i+1,j,k-1} - 2u_{i+1,j,k+1} 
      - u_{i+1,j-1,k-1} - u_{i+1,j+1,k-1} - u_{i+1,j-1,k+1} - u_{i+1,j+1,k+1}}
      {32\Delta}$.

   */
   double facteur; //normalisation factor used in gradient calculation

   facteur = 1./( 32.*Delta );

   foreach_dimension() {
#if EMBED
      if (fs.x[] && fs.x[1]) {
#endif
	 g.x = 4.*u[-1,0,0] - 4.*u[1,0,0]                                     \
	    + 2.*u[-1,-1,0] + 2.*u[-1,1,0] + 2.*u[-1,0,-1] + 2.*u[-1,0,1]  \
	    - 2.*u[1,-1,0]  - 2.*u[1,1,0]  - 2.*u[1,0,-1]  - 2.*u[1,0,1]   \
	    + u[-1,-1,-1] + u[-1,-1,1] + u[-1,1,-1] + u[-1,1,1]            \
	    - u[1,-1,-1]  - u[1,-1,1]  - u[1,1,-1]  - u[1,1,1];
	 g.x *= facteur;
#if EMBED
      }
      else // si on est au bord du domaine
         g.x = center_gradient_large(u);
#endif

   }
#endif

   compute_bndry_der(point,u,&g);
   
   return g;
}


pseudo_t compute_hessian_L2_local (Point point, vector g) {
   /* Compute the L2 hessian of the vector g_L2 containing the L2 gradient.
    */

   // computation of $H^*(u) = \nabla_R (\nabla_R u(P_k)^t)$

   pseudo_t h; 

   foreach_dimension ()
      h.x = local_gradient_L2 (point, g.x);

   //h = local_symmetrize_pseudo_t(h);


   compute_bndry_hessian(point, g, &h);
   
   return h;
}




/**
## Weak formulation based on Green formula


[source](https://www.theses.fr/2008PA066622)

*/



pseudo_v local_second_derivative_1D (Point point, scalar u) {

   pseudo_v g;
   
   foreach_dimension() {
#if EMBED
      if (cs[]) 
#endif
	 g.x = (u[-1] - 2*u[] + u[1]) / sq(Delta);
#if EMBED
      else
         g.x = 0.;
#endif
   }

   return g;
}

pseudo_t compute_hessian_Green_local (Point point, scalar u, vector g) {

   pseudo_t h;

#if dimension == 2
   foreach_dimension () {
#if EMBED
      if (cs[]) {
#endif
	 h.x.x = 1./4. * (g.x[0,1] + 2*g.x[] + g.x[0,-1]);
	 h.x.y = (- u[-1,1] + u[-1,-1] + u[1,1] - u[1,-1]) * 1./ (4.*sq(Delta));
#if EMBED
      }
      else {
         h.x.x = 0.;
         h.x.y = 0.;
      }
#endif
   }
#endif


   
   compute_bndry_hessian(point, g, &h);
   
	  
   return h;        
}

/**

Simplification cartesian stencil

*/

pseudo_t compute_hessian_Green_local_simple (Point point, scalar u, vector g) {

   pseudo_t h;

#if dimension == 2
   foreach_dimension () {
#if EMBED
      if (cs[]) {
#endif
	 h.x.x = (u[-1] - 2*u[] + u[1]) / sq(Delta);
	 h.x.y = (- u[-1,1] + u[-1,-1] + u[1,1] - u[1,-1]) * 1./ (4.*sq(Delta));
#if EMBED
      }
      else {
         h.x.x = 0.;
         h.x.y = 0.;
      }
#endif
   }
#endif

   compute_bndry_hessian(point, g, &h);

   
	  
   return h;        
}



/**
Choose the hessian reconstruction method

*/


void compute_gradient (int hessian_type, scalar u, vector g) {



   foreach() {
      pseudo_v gra;
      if (hessian_type==1 || hessian_type==4)
	 gra = local_gradient_centered (point, u);
      else if (hessian_type==2)
	 gra = local_gradient_L2 (point, u);
      // gradient computation not needed for Green hessian reconstruction
      // but we initilise another field instead
      else if (hessian_type==3)
	 gra = local_second_derivative_1D (point, u);
 
      foreach_dimension()
	 g.x[] = gra.x;
   }

/** FIXME: remove??? */   
   g.n[left] = neumann (0.); 
   g.n[right] = neumann (0.);  
#if dimension > 1
   g.n[top] = neumann (0.);  
   g.n[bottom] = neumann (0.); 
#endif
#if dimension > 2
   g.n[front] = neumann (0.);  
   g.n[back] = neumann (0.); 
#endif

   boundary((scalar *){g});
}


#if EMBED
void compute_embed_bdry_gradient (scalar u, vector g) {
   foreach() {
      foreach_dimension()
         g.x[] = center_gradient_large(u);
   }
   boundary((scalar *){g});
}
#endif


pseudo_t compute_hessian_local (Point point, int hessian_type, scalar u, vector g) {


   if (hessian_type==1) 
      return compute_hessian_centered_local (point, g); 
   else if (hessian_type==2) 
      return compute_hessian_L2_local (point, g);
   else if (hessian_type==3) 
      return compute_hessian_Green_local (point, u, g);
   else if (hessian_type==4) 
      return compute_hessian_Green_local_simple (point, u, g);
   else {
      fprintf(stderr,"WARNING : hessian_type value does not exists");
      abort();
   }
}



/**
# Error computation

When the hessian is known, we are able to compute the interpolation error.

We define a structure to store some intermediate variables 
useful for error computation.

*/

typedef struct {
   pseudo_t abs_hess;     // absolute hessian
   pseudo_t abs_eig_val;  // diagonal pseudo_t of eigen values of absolute hessian
   pseudo_t eig_vect;     // rotation matrix of hessian
   double det_abs_hess;   // determinant of absolute hessian
   pseudo_t metric_powed; // M^-1/2
} metric_error_var;

pseudo_t troncate_eigen(pseudo_t t);

metric_error_var compute_element_for_metric (pseudo_t hess) {
   /* 
      Compute some element necessary to get the metric, and which we can obtain 
      within the same foreach loop : 

      - Compute the absolute hessian of the hessian, i.e. the hessian 
      reconstructed with the absolutes values of the hessian. 

      - The determinant of the absolute hessian is already compute.

      - idem for the D_Lp factor

   */

   metric_error_var mev;
   pseudo_t temp;
   pseudo_v eig_val;

   // *** Compute |H| *** //
   // compute eigenvalues and vectors
   tensor_eigen(hess, &eig_val, &temp);
   mev.eig_vect = temp;
   mev.abs_eig_val = diagonal_pseudo_t_from_pseudo_v(eig_val);

   // take absolute values of the eigenvalues : eigen values of |H|
   foreach_dimension()
      mev.abs_eig_val.x.x = fabs(mev.abs_eig_val.x.x);

   //troncate eigen values
   mev.abs_eig_val = troncate_eigen(mev.abs_eig_val);

   // compute |H| = R*\Lambda*R^-1
   //--- Notice that R (the eigen vector matrix) is orthonormale 
   //=> R^-1 = R^T. It simplifies a lot.
   temp = multiply_local_tensor (mev.eig_vect, mev.abs_eig_val);
   mev.abs_hess = multiply_local_tensor (temp, transpose_pseudo_t(mev.eig_vect));

   // *** Compute det(|H|) *** //
   // compute det(|H|) = product of eigen values of |H|
   mev.det_abs_hess = 1.;
   foreach_dimension() 
      mev.det_abs_hess *= mev.abs_eig_val.x.x;

   return mev;
}


/** 
Troncate eigen values of |H| : $\lambda_i=max(|\lambda_i|,\varepsilon)$, with 
$\varepsilon=10^{-10}$. It avoids $det(|H|)= 0$, and is more correct for 
numerical computations. This is a standard in codes.
*/


pseudo_t troncate_eigen(pseudo_t t) {

   double epsilon = 1e-10;
   pseudo_t tmp;

   foreach_dimension() 
      tmp.x.x = max(t.x.x,epsilon);

   return tmp;
}



/**
## Assemble the whole local error computation in one function.


"Skeleton" function, which compute the local error from field of interest.
In other words : 

> *One void to compute them all,*

> *And in the code, bind them.*

*/


/**

We compute the local error $$|u-\pi_\mathcal{M} u|(x) = C_n tr(|H|) h^2 $$ and store it in the input scalar local_error


We allow to compute an error on multiple fields $u_i$ as
$$ E = \sum_i \omega_i ||u_i - \pi_\mathcal{M} u_i||_{\mathcal{L}^p} $$
where $\omega_i$ are user-defined weights.

*/

struct Adapt {
  scalar * slist; // list of scalars
  double * max;   // tolerance for each scalar
  int maxlevel;   // maximum level of refinement
  int minlevel;   // minimum level of refinement (default 1)
  scalar * list;  // list of fields to update (default all)
};


void compute_metric_error_isotropic (int Lp, struct Adapt p, scalar local_error) {

  double cn = 1./12.; 
   foreach () 
     local_error[] = 0.;

   // hessian, metric and error computation
   int i = 0;

   for (scalar psi in p.slist) {

     vector grad[];    // gradient
     compute_gradient(user_hessian, psi, grad);

     double wi;
     if (p.max != NULL)
       wi = p.max[i++];
     else
       wi = 1.;   // by default, weight=1

     foreach() {  

       pseudo_t hessl = compute_hessian_local (point, user_hessian, psi, grad); // local hessian
       metric_error_var mev = compute_element_for_metric (hessl);	

       double Htrace = 0.;  
       foreach_dimension() 
         Htrace += mev.abs_eig_val.x.x;

//       double errortmp = pow(Htrace,dimension/(2.*Lp+dimension));
//       local_error[] += wi*pow(pow(errortmp,Lp)*dv(),1./Lp); // optimal local error (without prefactor)
       local_error[] += wi*cn*Htrace*sq(Delta)*pow(dv(),1./Lp); // local error 
     }
     
   }

   boundary({local_error});
  
   return;
}


/**

We estimate the prefactors of the global interpolation error on optimal mesh and uniform mesh.

*/

struct PreFactorData {
  double copt, cuniform; 
};

struct etaData {
  int norm;
  scalar * slist;
  double * wil;
};


struct PreFactorData compute_prefactors ( struct etaData p) {

   int Lp = p.norm;
   struct PreFactorData cd;
   cd.copt = 0.;
   cd.cuniform = 0.;
   double cn = 1./12.; // fixme: verify if it is 1/12 also in the axi case

   int i = 0;
   for (scalar psi in p.slist) {

     scalar error[];
     vector grad[];    // gradient
     compute_gradient(user_hessian, psi, grad);

     double wi;
     if (p.wil != NULL)
       wi = p.wil[i++];
     else
       wi = 1.;

     double derr_uniform=0.;  // error on uniform mesh estimated via current metric 
     double derr=0.;  // optimal global error

     foreach( reduction(+:derr) reduction(+:derr_uniform) ) {  

       pseudo_t hessl = compute_hessian_local (point, user_hessian, psi, grad); // local hessian
       metric_error_var mev = compute_element_for_metric (hessl);	

       double Htrace = 0.;  // for optimal metric
       foreach_dimension() 
         Htrace += mev.abs_eig_val.x.x;

       derr += pow(Htrace,Lp*dimension/(2.*Lp+dimension)) * dv();

       double Htrace_cur = 0.;  
       foreach_dimension() 
         Htrace_cur += fabs(hessl.x.x);

       derr_uniform+= pow(Htrace,Lp)*dv();
     }


     cd.copt += wi*cn*pow(derr,(2.*Lp+dimension)/(Lp*dimension));  // prefactor of optimal error
   
     cd.cuniform += wi*cn*pow(derr_uniform,1./Lp)*sq(L0);  // prefactor of uniform error
     
   }

    return cd; 
}

/** 
Some cases see a total error explosion.
We show in [ARTICLE and thesis](https://theses.hal.science/tel-03966961) that it can be avoided by imposing a minimum size criterion at the cost of a small increase in the interpolation error.
 
For now, this function must be used by the user if necessary*/

double estimate_eta_opt ( struct etaData p ) {
   
   struct PreFactorData cd = compute_prefactors (p);

   return pow(cd.copt/cd.cuniform, dimension/2.)/pow(L0,dimension); //fixme: 3D? axisymmetric?
}



double compute_eta_min( struct etaData p) {

  int Lp = p.norm;
  double num = 0.;
  double den = 0.;
  double voltot = 0.;

   int i = 0;
   for (scalar psi in p.slist) {

     scalar error[];
     vector grad[];    // gradient
     compute_gradient(user_hessian, psi, grad);

     double wi;
     if (p.wil != NULL)
       wi = p.wil[i++];
     else
       wi = 1.;

     foreach( reduction(+:num) reduction(+:voltot) reduction(max:den) ) {  

       pseudo_t hessl = compute_hessian_local (point, user_hessian, psi, grad); // local hessian
       metric_error_var mev = compute_element_for_metric (hessl);	

       double Htrace = 0.;  // for optimal metric
       foreach_dimension() 
         Htrace += mev.abs_eig_val.x.x;

       voltot += wi*dv();
       num    += wi*pow(Htrace,Lp/(2.*Lp+dimension)) * dv();
       den     = wi*max(pow(Htrace,Lp/(2.*Lp+dimension)), den);

     }
   }

    return num/voltot/den; 
}



/**

# Note on the axi-symmetric case

If we write $\gamma_i$ the eigen values of the hessian, the error model is
$$ |u-\pi_\mathcal{M} u| = C_n h^2 \sum_{i=1}^n \gamma_i $$
with $C_n = C_3$ and $\gamma_3 = 0$.

The change of variable which interest us is, in non-axi cases
$$ d = h^{-n} $$
But in axi-symmetric case, the density is 
$$ d = h^{-2} $$
even if an axi-symmetric case is a 3D case.
The minimization problem writes then
$$ \text{find} \text{ min}_\mathcal{M} \left( \int_\Omega \left(|u-\pi_\mathcal{M} u|\right)^p \, dV \right) $$
$$ \text{under the constraint } \int_\Omega d \, dV = N $$

After the resolution, all the $n$ which appears as exponents in the solution come from the original change of variable $d = h^{-n}$.
In other words, in axi-symmetric, this $n$ is the dimension of the grid and not the dimension of the problem.

As a consequence, the functions we wrote are correct, as the variable "dimension" in Basilisk contains the dimension of the grid.
Moreover, the volume integration element $dv()$ is correct in axi.
Thus, we have nothing to change for the main components of the optimal error estimation.

The only remaining quenstion is:
is the constant $C_n$ also equal to $1/12$ in axi, as in 2D and 3D ? 
On instinct I would say yes, but that should be checked one day.

*/




