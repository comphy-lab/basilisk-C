/**
## Introduction
Two IBM-related operations are defined: the evaluations of the "penalty" force and torque on each particle, and the force distributions from the Lagrangian markers to their ambient Eulerian grid cells.

## Variable list
* `area`: average distance between two adjacent markers in 2D and area occupied by each marker in 3D
* `indice, phi, theta`: parameters for evenly distributing markers on a spherical surface by the golden spiral method (see also [this link](http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/))
* `dist_x,y,z`: distance between the cell center to the Lagrangian marker
* `weight_x,y,z`: kernels
*/

#include "ibm.h"

/**
## Evaluation of force and torque
We first calculate the average spacing (area) of each marker:
*/
event compute_particle_forces(i++)
{
    if (i > 0)
    {
        for (int k = 0; k < p_n; k++)
        {
            #if dimension == 2
            const double area = 2. * M_PI * particle[k].radius / particle[k].num_nodes;
            #else // dimension == 3
            const double area = 4. * M_PI * particle[k].radius * particle[k].radius / particle[k].num_nodes;
            #endif
          
/**
For particles with the precribed motions, we give the reference positions of Lagrangian markers at each time step:
*/
            for (int n = 0; n < particle[k].num_nodes; n++)
            {
                #if motion_type == 2
                foreach_dimension()
                    particle[k].center.pos_ref.x = pc[k].x;
                #if dimension == 2
                particle[k].node[n].pos_ref.x = particle[k].center.pos_ref.x + particle[k].radius * cos(2. * M_PI * (double)n / particle[k].num_nodes);
                particle[k].node[n].pos_ref.y = particle[k].center.pos_ref.y + particle[k].radius * sin(2. * M_PI * (double)n / particle[k].num_nodes);
                #elif dimension == 3
                double indice = (double)n + 0.5;
                double phi = acos(1. - 2. * indice / particle[k].num_nodes);
                double theta = M_PI * (1. + sqrt(5.)) * indice;
                particle[k].node[n].pos_ref.x = particle[k].center.pos_ref.x + particle[k].radius * cos(theta) * sin(phi);
                #endif // dimension
                #endif // motion_type

/**
We evaluate the "penalty" force and torque on each marker based on the deviation to the reference position:
$$\mathbf f_j(t)=-\kappa\frac{d}{\Delta x}\left[\mathbf X_j(t)-\mathbf X_j^{(0)}(t)\right]$$
*/
                foreach_dimension()
                    particle[k].node[n].lag_F.x = -particle[k].stiffness * (particle[k].node[n].pos.x - particle[k].node[n].pos_ref.x) * area;
            
                #if dimension == 2
                particle[k].node[n].lag_T.z = (particle[k].node[n].pos.x - particle[k].center.pos.x)*particle[k].node[n].lag_F.y
                                              - (particle[k].node[n].pos.y - particle[k].center.pos.y)*particle[k].node[n].lag_F.x;
                #elif dimension == 3
                foreach_dimension()
                    particle[k].node[n].lag_T.z = (particle[k].node[n].pos.x - particle[k].center.pos.x)*particle[k].node[n].lag_F.y 
                                              - (particle[k].node[n].pos.y - particle[k].center.pos.y)*particle[k].node[n].lag_F.x;
                #endif // dimension
            }
        }
    }
}

/**
## Force distribution
We first reset the particle-fluid force to zero. Note that here we re-arrange the global indices of interpolation stencils due to AMR.
*/
event spread(i++)
{
    if (i > 0)
    {
        Rearrange_indices(particle);
        reset({force}, 0.);

/**
We then calculate the distance between the center of grid cells within the chosen stencil and the target Lagrangian marker, evaluate the corresponding kernels, and eventually distribute the force on markers to their surrounding grid cells:
$$\mathbf F(\mathbf x,t)=\sum_{j}\mathbf f_j(t)\it{\Delta}\left(\mathbf X_j(t),\mathbf x\right)$$
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

                        foreach_dimension()
                            #if dimension == 2
                            force.x[] += (particle[k].node[n].lag_F.x * weight_x * weight_y);
                            #else // dimension == 3
                            force.x[] += (particle[k].node[n].lag_F.x * weight_x * weight_y * weight_z);
                            #endif
                    }
                }
            }
        }
    }
}