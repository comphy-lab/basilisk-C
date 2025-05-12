/**
## Introduction
In this file, you can find the discrete delta functions used for the force spreading and velocity interpolation in the feedback Immersed Boundary Method (IBM). Here I provide four different types of delta function:

* Regular 2-point (IBM_stencil=1) and 4-point (IBM_stencil=4) delta functions
* Smoothed 2-point (IBM_stencil=11) and 4-point delta functions  (IBM_stencil=14) proposed by Yang et, al [\[1\]](#Yang2009).

## Variable list
* `stencil`: delta function value
* `IBM_stencil`: delta function type
* `dist`:    distance between the cell center and the Lagrangian node on boundary.
*/

#ifndef IBM_stencil
#define IBM_stencil (1)
#endif

double stencil(double dist)
{

    double stencil = 0.;

    if (IBM_stencil == 1)
    {
        if (fabs(dist) <= 1.)
        {
            stencil = 1. - fabs(dist);
        }
        else
        {
            stencil = 0.;
        }   
    }
    else if (IBM_stencil == 4)
    {
        if (fabs(dist) < 1.)
        {
            stencil = 0.125 * (3. - 2 * fabs(dist) + sqrt(1. + 4 * fabs(dist) - 4 * sq(dist)));
        }
        else if (fabs(dist) <= 2. && fabs(dist) > 1.)
        {
            stencil = 0.125 * (5. - 2 * fabs(dist) - sqrt(-7. + 12 * fabs(dist) - 4 * sq(dist)));
        }
    }
    else if (IBM_stencil == 11)
    {
        if (fabs(dist) <= 0.5)
        {
            stencil = 3./4. - sq(dist);
        }
        else if (fabs(dist) > 0.5 && fabs(dist) <= 1.5)
        {
            stencil = 9./8. - 3.*fabs(dist)/2.+ sq(dist)/2.;
        }
        else
        {
            stencil = 0.;
        }
    }
    else if (IBM_stencil == 14)
    {
        if (fabs(dist) <= 0.5)
        {
            stencil = 3./8. + M_PI/32. - sq(dist)/4.;
        }
        else if (fabs(dist) > 0.5 && fabs(dist) <= 1.5)
        {
            stencil = 1./4. + (1-fabs(dist))*sqrt(-2.+8.*fabs(dist)-4*sq(dist))/8. - asin(sqrt(2.)*(fabs(dist)-1.))/8.;
        }
        else if (fabs(dist) > 1.5 && fabs(dist) <= 2.5)
        {
            stencil = 17./16.-M_PI/64.-3.*fabs(dist)/4.+sq(dist)/8.+(fabs(dist)-2.)*sqrt(-14.+16.*fabs(dist)-4.*sq(dist))/16.+asin(sqrt(2.)*(fabs(dist)-2.))/16.;
        }
        else
        {
            stencil = 0.;
        }
    }

    return stencil;
}

/**
## References 

~~~bib
@Article{Yang2009,
  author   = {X. Yang and X. Zhang and Z. Li and G. He},
  journal  = {Journal of Computational Physics},
  title    = {A smoothing technique for discrete delta functions with application to immersed boundary method in moving boundary simulations},
  year     = {2009},
  number   = {20},
  pages    = {7821-7836},
  volume   = {228},
  doi      = {https://doi.org/10.1016/j.jcp.2009.07.023},
  file     = {:Oscillating_Cylinder/1-s2.0-S0021999109004136-main.pdf:PDF},
  url      = {https://www.sciencedirect.com/science/article/pii/S0021999109004136},
}
~~~
*/
