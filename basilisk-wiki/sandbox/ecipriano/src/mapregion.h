/**
# Map Region

Map a region of the domain with a Heaviside function.
This function sets to 1 the cells containing the interface
and the gas phase (the liquid phase if inverse).
The number of layers defines how many layers of cells close to
the interface, and belonging to the liquid phase, must
be included (0 by default).

* *H*: Heaviside function to be filled
* *f*: Volume fraction field
* *nl*: Number of layers to be included
*/

#ifndef MAPREGION_H
#define MAPREGION_H

void mapregion (
  scalar H,             // heaviside function
  scalar f,             // vof field (f = 1 if liquid)
  int nl = 0,           // number of additional layers (default 0, optional 1 or 2)
  int inverse = 0,      // the vof field if = 1 if gas (default false)
  int narrow = 0,       // map just a narrow band around the interface (default false)
  int nointerface = 0,  // if false heaviside set to zero at the interface
  double tol = 1e-10    // tolerance which defines the interfacial cells
)
{
  scalar fc[];
  foreach()
    fc[] = (!inverse) ? f[] : 1. - f[];

  foreach() {
    if (narrow)
      H[] = (fc[] > tol && fc[] < 1.-tol) ? 1. : 0.;
    else
      H[] = (fc[] < 1.-tol) ? 1. : 0.;
    if (fc[] > 1.-tol && nl > 0) {
      bool lightup = false;
      foreach_neighbor(nl) {
        if (fc[] > tol && fc[] < 1.-tol) {
          lightup = true;
          break;
        }
      }
      H[] = lightup ? 1. : 0.;
    }
    if (nointerface)
      H[] = (fc[] > tol && fc[] < 1.-tol) ? 0. : H[];
  }
}

#endif
