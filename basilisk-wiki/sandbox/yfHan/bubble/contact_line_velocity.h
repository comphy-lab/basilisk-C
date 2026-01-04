/**
This headfile computes the contact line velocity in 2D domain
*/


/**
We will need basic functions for volume fraction computations. */
#include "vof.h" 

#ifndef F1_ERR
#define F1_ERR 1e-6 //tolerance on the interfacial cells
#endif
/**
F_ERR is the accepted error over f to avoid division by zero; it is also used to identify interfacial cells.
*/



void contact_line_velocity_and_normal_vector(scalar f, vector u, double width, int level , double *output_values)
{

 scalar ff[]; // Face vector for fluxes
 vector n[];       // Interface normal vector

 /** Initialize VOF and reconstruct the interface 
 */
reconstruction(f, n, ff );


double uline = 0.0;
double count = 0.0;
double dot_product = 0.0;
  
  foreach(reduction(+:uline) reduction(+:count) reduction(+:dot_product)) {
    if (f[] > F1_ERR && f[] < 1. - F1_ERR && x <= 1. * width / pow(2., level)) {
      uline += u.y[];
      count++;
  /** Compute the dot product of normal vector and contact line velocity */
      dot_product = dot_product + (n.x[]*u.x[] + n.y[]*u.y[]);
      //dot_product =  n.y[]*u.y[];
      
    }
  }
  
  /** 
  Compute the mean value of contact line velocity within the interfacial cells
  */
  if (count > 0)
  {  uline /= count;
    output_values[0]=uline;
    output_values[1]=dot_product;
    
  }
  else
    printf("The interface does not contact on the wall!");



}