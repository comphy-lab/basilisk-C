/**
# Compute the curvature of the mesh in 2D and 3D

*/
#if dimension < 3

/**
## Two-dimensional curvature computation

In 2D, we compute the curvature by fitting a 4th-degree polynomial to the mesh using two neighbours on each side of the node of interest
*/

void compute_2d_curvature(lagMesh*) {
  bool up; // decide if we switch the x and y axes
  lagNode* cn; // current node
  for(int i=0; i<mesh->nln; i++) {
    cn = &(mesh->nodes[i]);
    up = (fabs(cn->normal.y) > fabs(cn->normal.x)) ? true : false;
    coord p[5]; // store the coordinates of the current node and of its
                // neighbors'
    for(int j=0; j<5; j++) {
      int index = (mesh->nln + i - 2 + j)%mesh->nln;
      foreach_dimension() p[j].x = up ? mesh->nodes[index].pos.x :
        mesh->nodes[index].pos.y;
    }
    /** If one of the neighboring nodes is across a periodic boundary, we
correct its position */
    for(int j=0; j<5; j++) {
      if (j!=2) {
        foreach_dimension() {
          if (ACROSS_PERIODIC(p[j].x,p[2].x)) {
            p[j].x += (ACROSS_PERIODIC(p[j].x + L0, p[2].x)) ? -L0 : L0;
          }
        }
      }
    }

/** Since the bending force will involve taking the laplacian of the curvature,
we seek a cuvrature to fourth order accuracy, so we need to interpolate the
membrane with a fourth-degree polynomial, which we need to differentiate twice
to get the curvature:
$$P_4(x) = \sum_{j=i-2}^{i+2} y_j \prod_{k \neq j} \frac{x - x_j}{x_k - x_j} $$

$$P'_4(x) = \sum_{j=i-2}^{i+2} y_j \left( \prod_{k \neq j} \frac{1}{x_k - x_j}
\right) \sum_{l \neq j}\prod_{m \neq j, m \neq l} x - x_m $$

$$P''_4(x) = \sum_{j=i-2}^{i+2} y_j \left( \prod_{k \neq j} \frac{1}{x_k - x_j}
\right) \sum_{l \neq j}\sum_{m \neq j, m \neq l}\sum_{n \neq j, n \neq l, n
\neq m}x - x_n $$
*/
    double dy = 0.;
    double ddy = 0.;
    for(int j=0; j<5; j++) {
      double b1 = 0.; double b2 = 0.;
      for(int l=0; l<5; l++) {
        if (l!=j) {
          double c1 = 1.; double c2 = 0.;
          for(int m=0; m<5; m++) {
            if (m!=j && m!=l) {
              double d2 = 1.;
              for(int n=0; n<5; n++) {
                if (n!=j && n!=l && n!= m) {
                  d2 *= p[2].x - p[n].x;
                }
              }
              c1 *= p[2].x - p[m].x;
              c2 += d2;
            }
          }
          b1 += c1;
          b2 += c2;
        }
      }
      for(int k=0; k<5; k++) {
        if (k!=j) {
          b1 /= (p[k].x - p[j].x);
          b2 /= (p[k].x - p[j].x);
        }
      }
      dy += b1*p[j].y;
      ddy += b2*p[j].y;
    }

    /** The formula for the signed curvature of a function y(x) is
$$ \kappa = \frac{y''}{(1 + y'^2)^{\frac{3}{2}}}. $$
The sign is dertemined from a parametrization of the curve: walking
anticlockwise along the curve, if we turn left the curvature is positive. This
statement can be easily written as a dot product between edge i's normal
vector and edge (i+1)'s direction vector.*/
    coord a, b;
    foreach_dimension() {
      a.x = mesh->edges[mesh->nodes[i].edge_ids[0]].normal.x;
      b.x = mesh->edges[mesh->nodes[i].edge_ids[1]].normal.x;
    }
    int s = (a.x*b.x + a.y*b.y > 0) ? 1 : -1;
    cn->curv = s*fabs(ddy)/cube(sqrt(1 + sq(dy)));
  }
}

/**
## Three-dimensional curvature computation

In 3D things are a little more complicated. Fitting a 4-th degree polynomial to
the two-ring neighbors is too cumbersome and instead we fit a second-degree
paraboloid to the first-ring neighbors, following the approach of Yazdani and
Bagchi.

In most cases the node of interest has six neighbors and the system is
overconstrained: we use the least-square method which we
perform several times using a new normal vector computed from the paraboloid
equation.

Since this method involves inverting small $5 \times 5$ matrices, we include the
header file [matrix_toolbox.h](matrix_toolbox.h) containing helper functions to
perform this operation using the LU-decomposition.
*/

#else // dimension == 3
#include "matrix-toolbox.h"

/**
The function below fits a second-degree paraboloid to the one-ring neighbors
of a node of the membrane, using the [ordinary least-squares method](https://en.wikipedia.org/wiki/Ordinary_least_squares#Matrix/vector_formulation).
In the rare (exactly 12) cases where the node has only 5 neighbours, the least
squares method reduces to a simple 5x5 matrix inversion.

Input:

* $XX$ is the matrix of regressors, of size $5 \times 6$
* $\beta$ is the $5 \times 1$ vector of unknowns parametes
* $yy$ is the $6 \times 1$ vector of response variables
* perform_least_squares is a boolean indicating if the system is over-determined

Output:

The vector $\beta$ is populated with the parameters minimizing the sum of the
square of the error
*/
void fit_paraboloid(double** XX, double* beta, double* yy, bool perform_least_squares) {
  if (!perform_least_squares) {
    int* P = malloc(6*sizeof(int));
    LUPDecompose(XX, 5, 1.e-10, P);
    LUPSolve(XX, P, yy, 5, beta);
    free(P);
  }
  else {
    int* P = malloc(7*sizeof(int));
    /** We compute $\bm{X^T} \bm{X}$ and its inverse */
    double** AA = malloc(5*sizeof(double*));
    for(int i=0; i<5; i++) AA[i] = malloc(5*sizeof(double));
    double** IAA = malloc(5*sizeof(double*));
    for(int i=0; i<5; i++) IAA[i] = malloc(5*sizeof(double));
    for(int i=0; i<5; i++) {
      for(int j=0; j<5; j++) {
        AA[i][j] = 0.;
        for(int k=0; k<6; k++) AA[i][j] += XX[k][j]*XX[k][i];
      }
    }
    LUPDecompose(AA, 5, 1.e-10, P);
    LUPInvert(AA, P, 5, IAA);
    /** And then we apply the ordinary linear least squares differences:
    $$ \bm{\beta} = \left( \bm{X^T} \bm{X} \right)^{-1} \bm{X^T}\bm{y} $$
    */
    for(int i=0; i<5; i++) {
      beta[i] = 0.;
      for(int j=0; j<5; j++) {
        double Xty = 0.;
        for(int k=0; k<6; k++) Xty += XX[k][j]*yy[k];
        beta[i] += IAA[i][j]*Xty;
      }
    }
    free(P);
    for(int i=0; i<5; i++) {
      free(AA[i]);
      free(IAA[i]);
    }
    free(AA);
    free(IAA);
  }
}

/**
The function below computes the surface Laplacian (also knows a Laplace-Beltrami
operator) at node $i$ of either the membrane (diff_curv = false) or the
curvature of the membrane (diff_curv = true).

If diff_curv is set to false, this
function populates the curv and gcurv attribute of node i with the mean
curvature $\kappa = (c_1 + c_2)/2$ and the Gaussian curvature $\kappa_g = c_1
c_2$, with $c_1$ and $c_2$ the principal curvatures at node $i$. In this case
the return value is not relevant and set to 0.

If diff_curv is set to true, this function instead returns the surface Laplacian
of the (previously computed) curvature.
*/
double laplace_beltrami(lagMesh* mesh, int i, bool diff_curv) {
  /** If we wish to compute $\Delta_{LB} \kappa$, the result is stored in the
  variable $lbcurv$ and returned by the function. Otherwise the $lbcurv$ is
  not used. */
  double lbcurv = 0.; // lbcurv for "Laplace-Beltrami of the curvature"
  lagNode* cn = &(mesh->nodes[i]); // cn for "current node"
  int* ngb = cn->neighbor_ids;

  /** The following three variables are the data structures that will be
  fed to the least squares method */
  double** XX = malloc(6*sizeof(double*));
  for(int j=0; j<6; j++) XX[j] = malloc(5*sizeof(double));
  double* yy = malloc(6*sizeof(double));
  double* beta = malloc(5*sizeof(double));

  /** We successively fit a paraboloid to the one-ring neighbors until the
  vector normal to the paraboloid at the center node is converged. */
  double normal_change = 1.;
  int nb_fit_iterations = 0;
  int max_iterations = diff_curv ? 1 : 10;
  while (normal_change > 1.e-10 && nb_fit_iterations < max_iterations) {
    nb_fit_iterations++;
    /** Create local frame of reference and the linear system to solve */
    coord ex, ey, ez;
    foreach_dimension() ez.x = cn->normal.x;
    foreach_dimension() ey.x = GENERAL_1DIST(mesh->nodes[ngb[0]].pos.x,
      cn->pos.x);
    double eydez, ney, nex;
    eydez = 0.; ney = 0.; nex = 0.;
    foreach_dimension() eydez += ey.x*ez.x; // eydez for "ey dot ez"
    foreach_dimension() ey.x -= eydez*ez.x;
    ney = cnorm(ey);
    foreach_dimension() ey.x /= ney;
    foreach_dimension() ex.x = ey.y*ez.z - ey.z*ez.y;
    nex = cnorm(ex);
    foreach_dimension() ex.x /= nex;
    double M[3][3];
    M[0][0] = ex.x; M[1][0] = ex.y; M[2][0] = ex.z;
    M[0][1] = ey.x; M[1][1] = ey.y; M[2][1] = ey.z;
    M[0][2] = ez.x; M[1][2] = ez.y; M[2][2] = ez.z;

    double ipngb[6][3]; // ipngb for "initial position of neighbors"
    double rpngb[6][3]; // rpngb for "rotated position of neighbors"
    for(int j=0; j<cn->nb_neighbors; j++) {
      ipngb[j][0] = GENERAL_1DIST(mesh->nodes[ngb[j]].pos.x, cn->pos.x);
      ipngb[j][1] = GENERAL_1DIST(mesh->nodes[ngb[j]].pos.y, cn->pos.y);
      ipngb[j][2] = GENERAL_1DIST(mesh->nodes[ngb[j]].pos.z, cn->pos.z);
      for(int k=0; k<3; k++) {
        rpngb[j][k] = 0;
        for(int l=0; l<3; l++) {
          rpngb[j][k] += M[l][k]*ipngb[j][l];
        }
      }
    }

    /** Store the coordinates of the neighboring nodes in the appropriate
    data structure before using the ordinary least squares method.*/
    for(int j=0; j<cn->nb_neighbors; j++) {
      yy[j] = (diff_curv) ? (mesh->nodes[ngb[j]].curv -
        mesh->nodes[ngb[j]].ref_curv - cn->curv + cn->ref_curv) : rpngb[j][2];
      XX[j][0] = sq(rpngb[j][0]);
      XX[j][1] = rpngb[j][0]*rpngb[j][1];
      XX[j][2] = sq(rpngb[j][1]);
      XX[j][3] = rpngb[j][0];
      XX[j][4] = rpngb[j][1];
    }

    /** Solve the linear system directly (in case of 5 neighbors) or with
    the least-squares method (in case of 6 neighbors) */
    bool perform_least_squares = (cn->nb_neighbors > 5);
    fit_paraboloid(XX, beta, yy, perform_least_squares);

    if (!diff_curv) {
      /** Compute the mean and Gaussian curvatures, as well as a refined
      normal vector, from the fitted paraboloid. */
      cn->curv = -(beta[0] + beta[2] + beta[0]*sq(beta[4])
        + beta[2]*sq(beta[3]) - beta[1]*beta[3]*beta[4])
        /sqrt(cube(1 + sq(beta[3]) + sq(beta[4])));
      cn->gcurv = (4*beta[0]*beta[2] - sq(beta[1]))
        /sq(1 + sq(beta[3]) + sq(beta[4]));

      coord buff, prev_n;
      foreach_dimension() prev_n.x = cn->normal.x;
      buff.x = -beta[3]; buff.y = -beta[4]; buff.z = 1;
      double nn = cnorm(buff);
      foreach_dimension() buff.x /= nn;
      cn->normal.x = M[0][0]*buff.x + M[0][1]*buff.y + M[0][2]*buff.z;
      cn->normal.y = M[1][0]*buff.x + M[1][1]*buff.y + M[1][2]*buff.z;
      cn->normal.z = M[2][0]*buff.x + M[2][1]*buff.y + M[2][2]*buff.z;
      normal_change = 0.;
      foreach_dimension() {
        buff.x = cn->normal.x - prev_n.x;
        if (fabs(buff.x) > normal_change) normal_change = fabs(buff.x);
      }
      cn->nb_fit_iterations = nb_fit_iterations;
    }
    else {
      lbcurv = 2*(beta[0] + beta[2]);
    }
  }
  for(int j=0; j<6; j++) free(XX[j]);
  free(XX);
  free(yy);
  free(beta);
  return lbcurv;
}
#endif

/**
## General computation of the nodal curvatures

The function below computes the signed curvature of the Lagrangian mesh at each
node.
*/
void comp_curvature(lagMesh* mesh) {
  if (!mesh->updated_curvatures) {
    comp_normals(mesh);
    #if dimension < 3
      compute_2d_curvature(mesh);
    #else
      for(int i=0; i<mesh->nln; i++) {
        laplace_beltrami(mesh, i, false);
      }
    #endif
    mesh->updated_curvatures = true;
  }
}

void initialize_refcurv_onecaps(lagMesh* mesh) {
  #if REF_CURV
    comp_curvature(mesh);
  #endif
  for(int j=0; j<mesh->nln; j++) {
    #if (REF_CURV)
      #if GLOBAL_REF_CURV
        mesh->nodes[j].ref_curv = C0;
      #else
        mesh->nodes[j].ref_curv = mesh->nodes[j].curv;
      #endif
    #else
      mesh->nodes[j].ref_curv = 0.;
    #endif
  }
}

void initialize_refcurv() {
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive)
      initialize_refcurv_onecaps(&CAPS(i));
  }
}

event init (i = 0) {
  #if (RESTART_CASE == 0)
    initialize_refcurv();
  #endif
}
