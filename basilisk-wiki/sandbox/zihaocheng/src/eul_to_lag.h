/**
## Introduction
Two IBM-related operations are defined: the velocity interpolations from the ambient Eulerian grid cells to the Lagrangian markers, and updates to the translation and angular position of each particle after one time step. 

## Variable list
* `vel_x,y,z` & `vel_x,y,z_sum`: velocity of each Lagrangian marker for MPI implementations
* `dist_x,y,z`: distance between the cell center to the Lagrangian marker
* `weight_x,y,z`: kernels
*/

#include "ibm.h"

#if _MPI
double vel_x[p_n][tn_node];
double vel_y[p_n][tn_node];
double vel_x_sum[p_n][tn_node];
double vel_y_sum[p_n][tn_node];
#if dimension == 3
double vel_z[p_n][tn_node];
double vel_z_sum[p_n][tn_node];
#endif
#endif

/**
## Velocity interpolation
We first reset the velocity of all Lagrangian markers to zero:
*/
event interpolation(i++)
{
    if (i > 0)
    {
        for (int k = 0; k < p_n; k++)
        {
            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                foreach_dimension()
                    particle[k].node[n].vel.x = 0;
                #if _MPI
                vel_x[k][n] = 0.;
                vel_y[k][n] = 0.;
                vel_x_sum[k][n] = 0.;
                vel_y_sum[k][n] = 0.;
                #if dimension == 3
                vel_z[k][n] = 0.;
                vel_z_sum[k][n] = 0.;
                #endif
                #endif
            }
        }

/**
We then calculate the distance between the center of grid cells within the chosen stencil and the target Lagrangian marker, evaluate the corresponding kernels, and eventually interpolate the velocity of each marker:
$$\dot{\mathbf{X}}_j(t)=\sum_{\mathbf{x}} \Delta x^3\mathbf{u}(\mathbf{x},t)\it{\Delta}\left(\mathbf{X}_j(t),\mathbf{x}\right)$$
*/
        for (int k = 0; k < p_n; k++)
        {
            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                foreach_cache(particle[k].node[n].my_stencil)
                {
                    const double dist_x = particle[k].node[n].pos.x - x;
                    const double dist_y = particle[k].node[n].pos.y - y;
                    #if dimension == 3
                    const double dist_z = particle[k].node[n].pos.z - z;
                    #endif

                    #if dimension == 2
                    #if (IBM_stencil == 1 || IBM_stencil == 11)
                    if (fabs(dist_x) <= 1.5 && fabs(dist_y) <= 1.5)
                    #else // (IBM_stencil == 4 || IBM_stencil == 14)
                    if (fabs(dist_x) <= 2.5 && fabs(dist_y) <= 2.5)
                    #endif
                    #elif dimension == 3
                    #if (IBM_stencil == 1 || IBM_stencil == 11)
                    if (fabs(dist_x) <= 1.5 && fabs(dist_y) <= 1.5 && fabs(dist_z) <= 1.5)
                    #else // (IBM_stencil == 4 || IBM_stencil == 14)
                    if (fabs(dist_x) <= 2.5 && fabs(dist_y) <= 2.5 && fabs(dist_z) <= 2.5)
                    #endif
                    #endif
                    {
                        const double weight_x = stencil(dist_x);
                        const double weight_y = stencil(dist_y);
                        #if dimension == 3
                        const double weight_z = stencil(dist_z);
                        #endif

                        #if _MPI
                        #if dimension == 2
                        vel_x[k][n] += (u.x[] + 0.5 * force.x[] / rho[]) * weight_x * weight_y;
                        vel_y[k][n] += (u.y[] + 0.5 * force.y[] / rho[]) * weight_x * weight_y;
                        #elif dimension == 3
                        vel_x[k][n] += (u.x[] + 0.5 * force.x[] / rho[]) * weight_x * weight_y * weight_z;
                        vel_y[k][n] += (u.y[] + 0.5 * force.y[] / rho[]) * weight_x * weight_y * weight_z;
                        vel_z[k][n] += (u.z[] + 0.5 * force.z[] / rho[]) * weight_x * weight_y * weight_z;
                        #endif
                        #else
                        #if dimension == 2
                        particle[k].node[n].vel.x += (u.x[] + 0.5 * force.x[] / rho[]) * weight_x * weight_y;
                        particle[k].node[n].vel.y += (u.y[] + 0.5 * force.y[] / rho[]) * weight_x * weight_y;
                        #elif dimension == 3
                        particle[k].node[n].vel.x += (u.x[] + 0.5 * force.x[] / rho[]) * weight_x * weight_y * weight_z;
                        particle[k].node[n].vel.y += (u.y[] + 0.5 * force.y[] / rho[]) * weight_x * weight_y * weight_z;
                        particle[k].node[n].vel.z += (u.z[] + 0.5 * force.z[] / rho[]) * weight_x * weight_y * weight_z;
                        #endif
                        #endif
                    }
                }
            }
        }
    }
}

/**
## Advecting the Lagrangian markers
We obtain the new position of each marker by the explicit forward Euler method:
$$\mathbf X_j(t+\Delta t)=\mathbf X_j(t)+\dot{\mathbf X}_j(t)\Delta t$$ 
Note that in the present feedback IBM, the markers are allowed to be "pulled" away from their reference positions with different displacements, i.e. particles are not perfectly rigid with this scheme.
*/
event update_particle_position(i++)
{
    if (i > 0)
    {
        for (int k = 0; k < p_n; k++)
        {
            foreach_dimension()
                particle[k].center.pos.x = 0.;

            #if _MPI
            MPI_Reduce(vel_x[k], vel_x_sum[k], particle[k].num_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Bcast(vel_x_sum[k], particle[k].num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Reduce(vel_y[k], vel_y_sum[k], particle[k].num_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Bcast(vel_y_sum[k], particle[k].num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            #if dimension == 3
            MPI_Reduce(vel_z[k], vel_z_sum[k], particle[k].num_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Bcast(vel_z_sum[k], particle[k].num_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            #endif
            #endif

            // update node and center positions
            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                #if _MPI
                particle[k].node[n].pos.x += vel_x_sum[k][n];
                particle[k].node[n].pos.y += vel_y_sum[k][n];
                #if dimension == 3
                particle[k].node[n].pos.z += vel_z_sum[k][n];
                #endif
                #else
                foreach_dimension()
                    particle[k].node[n].pos.x += particle[k].node[n].vel.x;
                #endif

                foreach_dimension()
                    particle[k].center.pos.x += particle[k].node[n].pos.x / particle[k].num_nodes;
            }

/**
For the free motion of particles, we need to evaluate the reference positions of the markers according to the momentum equations. First, we sum up the force and torque experienced by each particle:
*/
            #if motion_type == 1
            foreach_dimension()
            {
                particle[k].F_tot.x = 0.;
                particle[k].T_tot.z = 0.;
            }

            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                foreach_dimension()
                {
                    particle[k].F_tot.x += particle[k].node[n].lag_F.x;
                    particle[k].T_tot.z += particle[k].node[n].lag_T.z;
                }
            }

/**
Then, we evaluate the hydrodynamic force on each particle by subtracting the buoyancy and "internal fluid" inertia from the total force, and similarly for torque. We eventually have [\[1\]](#Feng2009)[\[2\]](#Cheng2022):
$$\mathbf U_s^{n+1} = \left(1+\frac{\rho_f}{\rho_s}\right)\mathbf U_s^n - \frac{\rho_f}{\rho_s}\mathbf U_s^{n-1} + \frac{\left(-\mathbf F_{tot}+\left(\rho_s-\rho_f\right)V_s\mathbf g\right)\Delta t}{\rho_sV_s}$$
$$\mathbf \Omega_s^{n+1} = \left(1+\frac{\rho_f}{\rho_s}\right)\mathbf \Omega_s^n - \frac{\rho_f}{\rho_s}\mathbf \Omega_s^{n-1} - \frac{\mathbf T_{tot}\Delta t}{I_s}$$
*/
            foreach_dimension()
            {
                particle[k].center.vel.x = particle[k].center.vel_pre.x + particle[k].center.acl.x / particle[k].rho_ratio 
                                           - particle[k].F_tot.x / (particle_volume * particle[k].rho_ratio)  + (1. - 1. / particle[k].rho_ratio) * pa[k].x;
                particle[k].center.acl.x = particle[k].center.vel.x - particle[k].center.vel_pre.x;
                particle[k].center.pos_ref.x += particle[k].center.vel_pre.x + particle[k].center.acl.x / 2.;
                particle[k].center.vel_pre.x = particle[k].center.vel.x;

                particle[k].center.agl_vel.z = particle[k].center.agl_vel_pre.z + particle[k].center.agl_acl.z / particle[k].rho_ratio 
                                               - particle[k].T_tot.z / (moment_inertia * particle[k].rho_ratio);
                particle[k].center.agl_acl.z = particle[k].center.agl_vel.z - particle[k].center.agl_vel_pre.z;
                particle[k].angle.z += particle[k].center.agl_vel_pre.z + particle[k].center.agl_acl.z / 2.;
                particle[k].center.agl_vel_pre.z = particle[k].center.agl_vel.z;
            }

/**
Finally, we perform rotation on each particle:
*/
            #if dimension == 2
            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                double xref = particle[k].radius * cos(2. * M_PI * (double)n / particle[k].num_nodes);
                double yref = particle[k].radius * sin(2. * M_PI * (double)n / particle[k].num_nodes);
                particle[k].node[n].pos_ref.x = particle[k].center.pos_ref.x + xref * cos(particle[k].angle.z) - yref * sin(particle[k].angle.z);
                particle[k].node[n].pos_ref.y = particle[k].center.pos_ref.y + xref * sin(particle[k].angle.z) + yref * cos(particle[k].angle.z);
            }
            #elif dimension == 3
            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                double indice = (double)n + 0.5;
                double phi = acos(1. - 2. * indice / particle[k].num_nodes);
                double theta = M_PI * (1. + sqrt(5.)) * indice;
                double xref = particle[k].radius * cos(theta) * sin(phi);
                double yref = particle[k].radius * sin(theta) * sin(phi);
                double zref = particle[k].radius * cos(phi);
                particle[k].node[n].pos_ref.x = particle[k].center.pos_ref.x 
                                            + xref * cos(particle[k].angle.z) * cos(particle[k].angle.y) 
                                            + yref * (cos(particle[k].angle.z) * sin(particle[k].angle.y) * sin(particle[k].angle.x) - sin(particle[k].angle.z) * cos(particle[k].angle.x)) 
                                            + zref * (cos(particle[k].angle.z) * sin(particle[k].angle.y) * cos(particle[k].angle.x) + sin(particle[k].angle.z) * sin(particle[k].angle.x));
                particle[k].node[n].pos_ref.y = particle[k].center.pos_ref.y 
                                            + xref * sin(particle[k].angle.z) * cos(particle[k].angle.y) 
                                            + yref * (sin(particle[k].angle.z) * sin(particle[k].angle.y) * sin(particle[k].angle.x) + cos(particle[k].angle.z) * cos(particle[k].angle.x)) 
                                            + zref * (sin(particle[k].angle.z) * sin(particle[k].angle.y) * cos(particle[k].angle.x) - cos(particle[k].angle.z) * sin(particle[k].angle.x));
                particle[k].node[n].pos_ref.z = particle[k].center.pos_ref.z 
                                            - xref * sin(particle[k].angle.y) 
                                            + yref * cos(particle[k].angle.y) * sin(particle[k].angle.x) 
                                            + zref * cos(particle[k].angle.y) * cos(particle[k].angle.x);
            }
            #endif // dimension
            #endif // motion_type
        }
    }
}

/**
## Reference
~~~bib
@Article{Feng2009,
  author   = {Z. G. Feng and E. E. Michaelides},
  journal  = {Computers \& Fluids},
  title    = {Robust treatment of no-slip boundary condition and velocity updating for the lattice-{B}oltzmann simulation of particulate flows},
  year     = {2009},
  issn     = {0045-7930},
  number   = {2},
  pages    = {370-381},
  volume   = {38},
  doi      = {https://doi.org/10.1016/j.compfluid.2008.04.013},
  file     = {:LBM/1-s2.0-S0045793008001060-main.pdf:PDF},
  url      = {https://www.sciencedirect.com/science/article/pii/S0045793008001060},
}

@article{Cheng2022,
title = {An immersed boundary/multi-relaxation time lattice {B}oltzmann method on adaptive octree grids for the particle-resolved simulation of particle-laden flows},
journal = {Journal of Computational Physics},
volume = {471},
pages = {111669},
year = {2022},
issn = {0021-9991},
doi = {https://doi.org/10.1016/j.jcp.2022.111669},
url = {https://www.sciencedirect.com/science/article/pii/S002199912200732X},
author = {Zihao Cheng and Anthony Wachs},
keywords = {Lattice Boltzmann method, Immersed boundary method, Adaptive mesh refinement, Parallel computing, Particle-laden flow, Fluid-solid interaction},
}
~~~
*/