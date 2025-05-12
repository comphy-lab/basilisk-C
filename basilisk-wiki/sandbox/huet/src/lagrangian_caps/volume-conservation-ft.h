/**
# Enforcing volume conservation

In this file, we follow the method of Sigüenza et al. to enforce exact volume 
conservation of the capsule. The details can be found in appendix A of 
[\[1\]](#siguenza2016validation).
*/

/**
We will need to find the smallest real root of a cubic equation, which we
implement in a separate file.
*/
#include "smallest_root_cubic.h"

/**
The function below computes the quantity $\bm{\alpha_m}$ as specified in Eq.\,(41) of [1](#siguenza2016validation).
*/
trace
coord compute_alpha_m(lagMesh* mesh, int m) {
  coord alpha = {0., 0., 0.};

  for(int i=0; i<mesh->nodes[m].nb_triangles; i++) {
    int tid = mesh->nodes[m].triangle_ids[i];
    int local_id = 0;
    while (mesh->triangles[tid].node_ids[local_id] != m && local_id < 3) 
      local_id++;
    int nids[2];  // `nids` for "neighbor ids"
    nids[0] = mesh->triangles[tid].node_ids[(local_id + 1)%3];
    nids[1] = mesh->triangles[tid].node_ids[(local_id + 2)%3];

    foreach_dimension()
      alpha.x += -periodic_friendly_cross_product_x(
        mesh->nodes[nids[0]].pos, mesh->nodes[nids[1]].pos, 
        mesh->centroid)/12;
  }

  return alpha;
}

/**
The function below solves an optimization problem to enforce the
conservation of volume of a capsule. The cost function to minimize is the following:
\begin{equation}
  J_\Lambda(\bm{\delta X}) = \sum_{i = 1}^{n} \|\bm{\delta X}_i\|^2 - \Lambda \left( V(\bm{X} + \bm{\delta X}) - V_0 \right),
\end{equation}
where $\bm{X} = \left[\bm{X}_1 \cdots \bm{X}_n \right]$ is a tensor 
storing the coordinates $\bm{X}_i$ of the $n$ capsule nodes, 
$\bm{\delta X} = \left[\bm{\delta X}_1 \cdots \bm{\delta X}_n \right]$ 
is a tensor storing their displacements, $V(\bm{X})$ is the capsule 
volume in configuration $\bm{X}$, $V_0$ is the initial volume and 
$\Lambda$ is a Lagrange multiplier. Sigüenza et al. show that $\Lambda$
is a (real) root of a third order polynomial which coefficients are
given by Eqs.\,(43-46) in [1](#siguenza2016validation). An assumption 
for their derivation is $\| \bm{\delta X} \| \ll \| \bm{X} \|$, so the
solution to the optimization will not quite be exact, but this 
approximation is asymptotically true as $\Delta t$ tends to zero.
*/
trace
void enforce_optimal_volume_conservation(lagMesh* mesh) {
  /** First, the $\bm{\alpha_m}$ are computed and stored in an array 
  of size the number of Lagrangian nodes.*/
  coord* alpha = malloc(mesh->nln*sizeof(coord));
  for(int m=0; m<mesh->nln; m++)
    alpha[m] = compute_alpha_m(mesh, m);

  /** Once all the $\bm{\alpha_m}$ are known, the coefficients of the
  polynomial in $\Lambda$ can be computed:
  \begin{equation}
    A \Lambda^3 + B \Lambda^2 + C \Lambda + D = 0
  \end{equation}
  with
  \begin{equation}
    A = \frac{1}{18} \sum_{i=1}^{N_t} \[ 
      \bm{\alpha_{i,1}} \cdot (\bm{\alpha_{i,2}} \cross \bm{\alpha_{i,3}}) + 
      \bm{\alpha_{i,2}} \cdot (\bm{\alpha_{i,3}} \cross \bm{\alpha_{i,1}}) +
      \bm{\alpha_{i,3}} \cdot (\bm{\alpha_{i,1}} \cross \bm{\alpha_{i,2}}) \],
  \end{equation}
  \begin{equation}
    B = \frac{1}{6} \sum_{i=1}^{N_t} \[ 
      \bm{x_{i,1}} \cdot (\bm{\alpha_{i,2}} \cross \bm{\alpha_{i,3}}) + 
      \bm{x_{i,2}} \cdot (\bm{\alpha_{i,3}} \cross \bm{\alpha_{i,1}}) +
      \bm{x_{i,3}} \cdot (\bm{\alpha_{i,1}} \cross \bm{\alpha_{i,2}}) \],
  \end{equation}
  \begin{equation}
    C = \frac{1}{6} \sum_{i=1}^{N_t} \[ 
      \bm{\alpha_{i,1}} \cdot (\bm{x_{i,2}} \cross \bm{x_{i,3}}) +
      \bm{\alpha_{i,2}} \cdot (\bm{x_{i,3}} \cross \bm{x_{i,1}}) +
      \bm{\alpha_{i,3}} \cdot (\bm{x_{i,1}} \cross \bm{x_{i,2}})
      \],
  \end{equation}
  and 
  \begin{equation}
    D = V - V_0
  \end{equation}
  where $\star_{i,1}$ is the quantity $\star$ at vertex 1 of triangle $i$, $V$ is the capsule volume and $V_0$ is the initial (conserved)
  volume.
  In practice, we normalize the polynomial so that its largest coefficient is 1.
  */
  double coeff_polynomial[4];
  for(int k=0; k<4; k++) coeff_polynomial[k] = 0;
  for(int j=0; j<mesh->nlt; j++) {
    int tn[3]; // `tn` for "triangle nodes"
    for(int k=0; k<3; k++) 
      tn[k] = mesh->triangles[j].node_ids[k];
    for(int k=0; k<3; k++) {
      coord cp0; coord cp2; // `cp` for "cross-product"
      foreach_dimension() {
        cp0.x = alpha[tn[(k+1)%3]].y*alpha[tn[(k+2)%3]].z
          - alpha[tn[(k+1)%3]].z*alpha[tn[(k+2)%3]].y;
        cp2.x = periodic_friendly_cross_product_x(
          mesh->nodes[tn[(k+1)%3]].pos, mesh->nodes[tn[(k+2)%3]].pos,
          mesh->centroid);
      }
      coeff_polynomial[3] += cdot(alpha[tn[k]], cp0);
      coeff_polynomial[2] += periodic_friendly_dot_product(
      mesh->nodes[tn[k]].pos, cp0, mesh->centroid);
      coeff_polynomial[1] += cdot(alpha[tn[k]], cp2);
    }
  }
  coeff_polynomial[3] /= 18;
  coeff_polynomial[2] /= 6;
  coeff_polynomial[1] /= 6;
  coeff_polynomial[0] = mesh->volume - mesh->initial_volume;
  double normalize_factor = 
    max(max(fabs(coeff_polynomial[3]), coeff_polynomial[2]), 
    max(fabs(coeff_polynomial[1]), coeff_polynomial[0]));
  coeff_polynomial[3] /= normalize_factor;
  coeff_polynomial[2] /= normalize_factor;
  coeff_polynomial[1] /= normalize_factor;
  coeff_polynomial[0] /= normalize_factor;
  /**
  Once the coefficients are found, we find its smallest real root (in
  absolute value) thanks to a separate routine.
  */
  double lambda = find_smallest_real_root(coeff_polynomial);

  /**
  Finally, the correction applied to each node $i$ is
  \begin{equation}
    \bm{x_i} \leftarrow \bm{x_i} + \Lambda \bm{\alpha_m} 
  \end{equation}
  */
  for(int i=0; i<mesh->nln; i++) {
    foreach_dimension() 
      mesh->nodes[i].pos.x += lambda*alpha[i].x;
  }
  free(alpha);
  correct_lag_pos(mesh);
}

/**
## References

~~~bib
@Article{siguenza2016validation,
  author  = {Sig{\"u}enza, Julien and Mendez, Simon and Ambard, Dominique and Dubois, Fr{\'e}d{\'e}ric and Jourdan, Franck and Mozul, R{\'e}my and Nicoud, Franck},
  journal   = {Journal of Computational Physics},
  title   = {Validation of an immersed thick boundary method for simulating fluid--structure interactions of deformable membranes},
  year    = {2016},
  pages   = {723--746},
  volume  = {322},
  file    = {:1-s2.0-S0021999116302662-main.pdf:PDF},
  publisher = {Elsevier},
}

~~~
*/