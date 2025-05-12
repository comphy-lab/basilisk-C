/**
---
header-includes:
  - \usepackage{algorithm}
  - \usepackage{algpseudocode}
---

# Motivations

*Non-closed, moving interfaces* are two-phase interfaces 
that extend to infinity, therefore that are prolonging 
beyond the simulation box.
For non-closed, moving interfaces, **curvature computation** at 
intersecting points with borders is hence a challenge, 
as it is critical for the stability of the interface. 

Since curvature computation is based on *height functions*, 
and height functions depend on volume fractions belonging to 
the column in which they are computed, we do need 
to assess precisely volume fractions values in *ghost cells*. 

The present file details an algorithm able to compute 
efficiently the volume fractions in ghost cells of an 
intersected boundary, when the normal at the interface 
crossing this boundary is already known.



# Algorithm

An algorithm for computing volume fraction of ghost cells on 
a boundary crossed by a non-closed, moving interface, is 
designed, based on subroutines $\mathcal{V}$ and $\mathcal{V}^{-1}$
developed by [Scardovelli \& Zaleski, (2000)](#scardovelli2000).

It follows the notations of the sketch below (example for the top border):

<img src="img_vel_tan/sketch_algo_vel_tan.png" width="400" class="center"/>

________________
**Update of volume fraction value in ghost cell $\mathcal{C}^*$**

________________

**foreach** in-domain cell $\mathcal{C}$ on a boundary, **do**

> **if** $0 < C[\,] < 1$ **then** *(current volume fraction)*
    
>> Define normal vector $\mathbf{n}$ of cell $\mathcal{C}$;

>> $\alpha \gets \mathcal{V}^{-1}(\mathbf{n},C[\,])$; 

>> **return** $C^*[\,] \gets \mathcal{V}(\mathbf{n}, \alpha,
    \text{coords of bottom-left/top-right corners of 
    $\mathcal{C}^*$}$ rel. to $\mathcal{C}$).

> **end if**

> **if** $0 < C[1] < 1$ **then** *(downstream volume fraction)*

>> Define normal vector $\mathbf{n}$ of cell $\mathcal{C}^+$;

>> $\alpha \gets \mathcal{V}^{-1}(\mathbf{n},C[1])$; 

>> \textbf{return} $C^*[\,] \gets \mathcal{V}(\mathbf{n}, \alpha,
    \text{coords of bottom-left/top-right corners of  
    $\mathcal{C}^*$}$ rel. to $\mathcal{C}^+$).

> **end if**

> **if** $0 < C[-1] < 1$ **then** *(upstream volume fraction)*

>> Define normal vector $\mathbf{n}$ of cell $\mathcal{C}^-$;

>> $\alpha \gets \mathcal{V}^{-1}(\mathbf{n},C[-1])$; 

>> \textbf{return} $C^*[\,] \gets \mathcal{V}(\mathbf{n}, \alpha,
    \text{coords of bottom-left/top-right corners of
    $\mathcal{C}^*$}$ rel. to $\mathcal{C}^-$).

> **end if**

> \textbf{return} $C^*[\,] \gets C[\,]$.

**end for**


## References

~~~bib
@article{scardovelli2000,
title = {Analytical Relations Connecting Linear Interfaces and Volume Fractions in Rectangular Grids},
journal = {Journal of Computational Physics},
volume = {164},
number = {1},
pages = {228-237},
year = {2000},
issn = {0021-9991},
doi = {https://doi.org/10.1006/jcph.2000.6567},
url = {https://www.sciencedirect.com/science/article/pii/S0021999100965677},
author = {Ruben Scardovelli and Stephane Zaleski},
abstract = {The computational-geometric problems arising when a linear interface cuts a cube are considered. They are of interest in particular for the calculation of volume fractions or interface positions in three-dimensional interface calculations in the Volume of Fluid (VOF) methods. Typically, the normal vector is known. One then wants to compute the volume fraction knowing the interface position, or conversely the interface position knowing the volume fraction. Explicit expressions of general use are given, and the algorithms used to search for solutions are described in detail. Explicit formulas for cubic roots are found to be less than two thirds as time consuming as Newton–Raphson iterations.}
}
~~~


# Code for a straight interface crossing left and bottom borders

To sum up, the main idea for computing the correct volume fraction 
`f[border]` in the `border` ghost cell is quite simple; if we know:

  * the normal at the interface crossing the `border` boundary;
  * a non-zero volume fraction `ffX` in the neighbourhood of the cell covered;

then an intercept can be evaluated with the help of 
[appropriated functions](http://basilisk.fr/src/geometry.h).

Thanks to this intercept, the normal, and the use of the coordinates 
of `f[border]`'s cell regarding the position of `ffX`, we are finally able to 
compute the volume fraction in the ghost cell.

So, as a first assumption, we need to know analytically the normal 
to the interface $\mathbf{n}$. 
For a straight liquid interface crossing left and bottom borders 
(liquid phase occupying the bottom-left corner), 
the normal (directed from the liquid phase towards the gas one) 
reads in 2D:

$$
n_x = \dfrac{\sin \beta_0}{|\cos \beta_0| + |\sin \beta_0|}
\quad ; \quad
n_y = \dfrac{\cos \beta_0}{|\cos \beta_0| + |\sin \beta_0|}
$$
where $\beta_0 = \arctan(y/x)$ and using the $L_1-$norm convention.

Then, the rest of the code is straightforward thanks to the 
details given in the [*Algorithm*](#algorithm) section.


*/

double f_BC_left(double ff0, double ff1, double ff2){
  // ff0, ff1, ff2 = f[], f[0,1], f[0,-1] respectively
  if ((ff0 > 0.0) && (ff0 < 1.0)){
    coord ntemp = {
      sin(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0))), 
      cos(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)))
    };
    double alpha_temp = plane_alpha (ff0, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-1.5,-0.5}, (coord){-0.5,0.5});
  }

  if ((ff1 > 0.0) && (ff1 < 1.0)){
    coord ntemp = {
      sin(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0))), 
      cos(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)))
    };
    double alpha_temp = plane_alpha (ff1, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-1.5,-1.5}, (coord){-0.5,-0.5});
  }

  if ((ff2 > 0.0) && (ff2 < 1.0)){
    coord ntemp = {
      sin(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0))), 
      cos(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)))
    };
    double alpha_temp = plane_alpha (ff2, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-1.5,0.5}, (coord){-0.5,1.5});
  }

  return ff0;
}

double f_BC_bottom(double ff0, double ff1, double ff2){
  // ff0, ff1, ff2 = f[], f[1,0], f[-1,0] respectively
  if ((ff0 > 0.0) && (ff0 < 1.0)){
    coord ntemp = {
      sin(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0))), 
      cos(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)))
    };
    double alpha_temp = plane_alpha (ff0, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-0.5,-1.5}, (coord){0.5,-0.5});
  }

  if ((ff1 > 0.0) && (ff1 < 1.0)){
    coord ntemp = {
      sin(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0))), 
      cos(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)))
    };
    double alpha_temp = plane_alpha (ff1, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-1.5,-0.5}, (coord){-0.5,-0.5});
  }

  if ((ff2 > 0.0) && (ff2 < 1.0)){
    coord ntemp = {
      sin(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0))), 
      cos(BETA_0)/(fabs(sin(BETA_0)) + fabs(cos(BETA_0)))
    };
    double alpha_temp = plane_alpha (ff2, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){0.5,-1.5}, (coord){1.5,-0.5});
  }

  return ff0;
}
