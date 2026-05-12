// This is the headfile of computing contact line velocity in 3D 

// 003 with complex interpolation

/**
We will need basic functions for volume fraction computations. */
//#include "grid/octree.h"   // Required for grid management and foreach loops
#include "vof.h" 
//#include "fractions.h"     // Required if f is a VOF (Volume of Fluid) field


//#ifndef F1_ERR
#define F1_ERR 1e-3 // Tolerance on the interfacial cells
//#endif

// 2nd-order interpolation function
double quatral_interpolation(double fx, double vala, double valb, double valc, double fa, double fb, double fc) {
    return vala*(fx-fb)*(fx-fc)/(fa-fb)/(fa-fc) + valb*(fx-fa)*(fx-fc)/(fb-fa)/(fb-fc) + valc*(fx-fb)*(fx-fa)/(fc-fb)/(fc-fa);
}

// 1st-order linear interpolation function
double linear_interpolation(double fx, double vala, double valb, double fa, double fb) {
    return vala + (vala-valb)/(fa-fb)*(fx-fa);
}

// Gradient macros
#define center_gradient_x(a) ((a[1,0,0] - a[-1,0,0])/(2.*Delta));
#define center_gradient_z(a) ((a[0,0,1] - a[0,0,-1])/(2.*Delta));


//void contact_line_velocity_in_xz_plane(scalar f, vector u, double width, int level, double *output_values) 

void contact_line_velocity_and_normal_vector(scalar f, vector u, scalar ucl2)
{
    scalar ff[]; // Face vector for fluxes
    vector n[];  // Interface normal vector


    // Reconstruct the interface
    reconstruction(f, n, ff);

    foreach(){
        if (f[] > F1_ERR && f[] < 1. - F1_ERR && y <= 1. * Delta ) {

            // Gradient calculation
            double gradient_fx = center_gradient_x(f);
            double gradient_fz = center_gradient_z(f);
            double interpolated_u_x =0.0, interpolated_u_z =0.0;

            // Initialize interpolation points
            double fp_a =0., fp_b =0., up_ax=0., up_az=0., up_bx=0., up_bz=0.; // fp_a: f[P-1]; fp_b: f[P+1];
            double dx1 = fabs(n.x[]/n.z[]), dz1=fabs(n.z[]/n.x[]);

            if (fabs(n.x[])<fabs(n.z[])) { 
                // if |n.x|<|n.z|, we interpolate P-1 and P+1 along the x direction
                if (n.x[]*n.z[]>0.){ 
                    // Positive slope
                    fp_a = f[0,0,-1] - dx1*(f[-1,0,-1]-f[0,0,-1]);
                    up_ax = u.x[0,0,-1] - dx1*(u.x[-1,0,-1]-u.x[0,0,-1]);
                    up_az = u.z[0,0,-1] - dx1*(u.z[-1,0,-1]-u.z[0,0,-1]);
                    
                    fp_b = f[0,0,1] + dx1*(f[1,0,1]-f[0,0,1]);
                    up_bx = u.x[0,0,1] + dx1*(u.x[1,0,1]-u.x[0,0,1]);
                    up_bz = u.z[0,0,1] + dx1*(u.z[1,0,1]-u.z[0,0,1]);
                }
                else{
                    // Negative slope
                    fp_a = f[0,0,-1] - dx1*(f[1,0,-1]-f[0,0,-1]);
                    up_ax = u.x[0,0,-1] - dx1*(u.x[1,0,-1]-u.x[0,0,-1]);
                    up_az = u.z[0,0,-1] - dx1*(u.z[1,0,-1]-u.z[0,0,-1]);
                    
                    fp_b = f[0,0,1] + dx1*(f[-1,0,1]-f[0,0,1]);
                    up_bx = u.x[0,0,1] + dx1*(u.x[-1,0,1]-u.x[0,0,1]);
                    up_bz = u.z[0,0,1] + dx1*(u.z[-1,0,1]-u.z[0,0,1]);

                }
            }
            else { 
                // if |n.x|>|n.z|, we interpolate P-1 and P+1 along the z direction
                if (n.x[]*n.z[]>0.){
                    fp_a = f[-1,0,0] - dz1*(f[-1,0,-1]-f[-1,0,0]);
                    up_ax = u.x[-1,0,0] - dz1*(u.x[-1,0,-1]-u.x[-1,0,0]);
                    up_az = u.z[-1,0,0] - dz1*(u.z[-1,0,-1]-u.z[-1,0,0]);
                    
                    fp_b = f[1,0,0] + dz1*(f[1,0,1]-f[1,0,0]);
                    up_bx = u.x[1,0,0] + dz1*(u.x[1,0,1]-u.x[1,0,0]);
                    up_bz = u.z[1,0,0] + dz1*(u.z[1,0,1]-u.z[1,0,0]);
                }
                else{
                    fp_a = f[-1,0,0] - dz1*(f[-1,0,1]-f[-1,0,0]);
                    up_ax = u.x[-1,0,0] - dz1*(u.x[-1,0,1]-u.x[-1,0,0]);
                    up_az = u.z[-1,0,0] - dz1*(u.z[-1,0,1]-u.z[-1,0,0]);
                    
                    fp_b = f[1,0,0] + dz1*(f[1,0,-1]-f[1,0,0]);
                    up_bx = u.x[1,0,0] + dz1*(u.x[1,0,-1]-u.x[1,0,0]);
                    up_bz = u.z[1,0,0] + dz1*(u.z[1,0,-1]-u.z[1,0,0]);

                }

            }

            // Determine interpolation direction
                        if ((f[] >= 0.5 && ((gradient_fx < 0. && gradient_fz <= 0.) || (gradient_fx < 0. && gradient_fz > 0.))) ||
    (f[] < 0.5 && ((gradient_fx >= 0. && gradient_fz >= 0.) || (gradient_fx >= 0. && gradient_fz <= 0.)))) {
                // forwards
                interpolated_u_x = linear_interpolation(0.5, u.x[], up_bx, f[], fp_b);
                interpolated_u_z = linear_interpolation(0.5, u.z[],  up_bz, f[], fp_b);

            } 
            else {
                //backwards
                interpolated_u_x = linear_interpolation(0.5, u.x[], up_ax, f[], fp_a);
                interpolated_u_z = linear_interpolation(0.5, u.z[], up_az, f[], fp_a);
            }



            double dot_product_xz = 0.0;
            

            // Compute the velocity in the xz plane (u_x and u_z contributions)
            double cos_alpha = n.x[] /sqrt(n.x[] * n.x[] + n.z[] * n.z[]);
            double sin_alpha = n.z[] /sqrt(n.x[] * n.x[] + n.z[] * n.z[]);
            double uline_xz = cos_alpha * interpolated_u_x + sin_alpha * interpolated_u_z;
            
            
             // Dot product
            dot_product_xz = (n.x[] * uline_xz * cos_alpha + n.z[] * uline_xz * sin_alpha);
            ucl2[] = dot_product_xz >= 0.0 ? fabs(uline_xz) : -fabs(uline_xz);
        }
        else{
            ucl2[]=0.0;
        }


    }

   
}

