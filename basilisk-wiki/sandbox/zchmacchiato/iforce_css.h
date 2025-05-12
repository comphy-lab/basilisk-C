/**
# Interfacial forces

We assume that the interfacial acceleration can be expressed as
$$
\phi\mathbf{n}\delta_s/\rho
$$
with $\mathbf{n}$ the interface normal, $\delta_s$ the interface Dirac
function, $\rho$ the density and $\phi$ a generic scalar field. Using
a CSF/Peskin-like approximation, this can be expressed as
$$
\phi\nabla f/\rho
$$
with $f$ the volume fraction field describing the interface.

The interfacial force potential $\phi$ is associated to each VOF
tracer. This is done easily by adding the following [field
attributes](/Basilisk C#field-attributes). */

attribute {
  scalar phi;
}

/**
Interfacial forces are a source term in the right-hand-side of the
evolution equation for the velocity of the [centered Navier--Stokes
solver](navier-stokes/centered.h) i.e. it is an acceleration. If
necessary, we allocate a new vector field to store it. */

event defaults (i = 0) {  
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
  }
}

/**
The calculation of the acceleration is done by this event, overloaded
from [its definition](navier-stokes/centered.h#acceleration-term) in
the centered Navier--Stokes solver. */

event acceleration (i++)
{
  
  /**
  We check for all VOF interfaces for which $\phi$ is allocated. The
  corresponding volume fraction fields will be stored in *list*. */

  scalar * list = NULL;
  for (scalar f in interfaces)
    if (f.phi.i) {
      list = list_add (list, f);

      /**
      To avoid undeterminations due to round-off errors, we remove
      values of the volume fraction larger than one or smaller than
      zero. */

      foreach()
	f[] = clamp (f[], 0., 1.);
    }

  /**
  On trees we need to make sure that the volume fraction gradient
  is computed exactly like the pressure gradient. This is necessary to
  ensure well-balancing of the pressure gradient and interfacial force
  term. To do so, we apply the same prolongation to the volume
  fraction field as applied to the pressure field. */
  
#if TREE
  for (scalar f in list) {
    f.prolongation = p.prolongation;
    f.dirty = true; // boundary conditions need to be updated
  }
#endif

  /**
  Finally, for each interface for which $\phi$ is allocated, we
  compute the interfacial force acceleration
  $$
  \phi\mathbf{n}\delta_s/\rho \approx \alpha\phi\nabla f
  $$ 
  */

  face vector ia = a;
  scalar CSS11[],CSS12[],CSS22[];
  face vector CSS[];
  
    foreach(){
    CSS11[]=0.0;
    CSS12[]=0.0;
    CSS22[]=0.0;
    
  }
  foreach_face(){
    CSS.x[]=0;
  }


  CSS11[left] = neumann(0);  
  CSS12[left] = neumann(0);  
  CSS22[left] = neumann(0);  

  CSS11[top] = neumann(0);  
  CSS12[top] = neumann(0);  
  CSS22[top] = neumann(0);  

  CSS.n[left] = neumann(0);
  CSS.n[top] = neumann(0);
  /**
  2-D case: Declare tensor component. xx-> CSS11, xy,yx->CSS12, yy->CSS22.
  face vector CSS. Initilize, and give the boundary condition.
  */
  
for (scalar f in list){

  foreach(){
    double mole=sqrt(sq((f[0,0] - f[0,-1])/Delta)+sq((f[0,0] - f[-1,0])/Delta)+1e-6);
    CSS11[]=  sq((f[0,0] - f[0,-1])/Delta)/mole;
    CSS12[]= - (f[0,0] - f[0,-1])*(f[0,0] - f[-1,0])/sq(Delta)/mole;
    CSS22[]=  sq((f[0,0] - f[-1,0])/Delta)/mole;
  }
  foreach_face(){

  CSS.x[]= ((CSS11[1,0]-CSS11[-1,0])+
           (CSS12[0,1]-CSS12[0,-1]))/2.0/Delta;

  CSS.y[]= ((CSS12[1,0]-CSS12[-1,0])+
           (CSS22[0,1]-CSS22[0,-1]))/2.0/Delta;
  }
  boundary({CSS,CSS11,CSS12,CSS22});
  }
/**
Compute the tensor components and the divergence.

$\nabla f$ is computed at face. So, the divergence $\nabla \cdot \Pi_s$ is then computed by central difference, and it is still located at face.
*/
  
  foreach_face()
    for (scalar f in list)
      if (f[] != f[-1] && fm.x[] > 0.) {

	/**
	We need to compute the potential *phif* on the face, using its
	values at the center of the cell. If both potentials are
	defined, we take the average, otherwise we take a single
	value. If all fails we set the potential to zero: this should
	happen only because of very pathological cases e.g. weird
	boundary conditions for the volume fraction. */
	
	scalar phi = f.phi;
	double phif =
	  (phi[] < nodata && phi[-1] < nodata) ?
	  (phi[] + phi[-1])/2. :
	  phi[] < nodata ? phi[] :
	  phi[-1] < nodata ? phi[-1] :
	  0.;
/** we need to change the CSF model to CSS model;
The equation of the CSF is  
$$\phi\mathbf{n}\delta_s/\rho \approx \alpha\phi\nabla f $$
The equation of the CSS is 
$$\sigma/\rho\nabla\cdot \left( (\mathbf{I}-\mathbf{n}\otimes\mathbf{n})\delta_s\right) \approx \sigma/\rho\nabla\cdot
\left[\left(|\nabla f|^2 \mathbf{I}-\nabla f\otimes \nabla f \right)/|\nabla f|\right]$$

The inner stress tensor is exactly

$$
\Pi_s=\left( \begin{array}{cc}
n_y \partial_y f & - n_x \partial_y f \\
-n_y \partial_x f & n_x \partial_x f  
\end{array}\right)
$$

Then the divergence of the stress tensor is

$$
\nabla \cdot \Pi_s=\left(\begin{array}{c}
\partial_x(n_y \partial_y f)-\partial_y(n_x \partial_y f)\\
\partial_y(n_x \partial_x f)-\partial_x(n_y \partial_x f)
\end{array}\right)
$$

Finally, the equation can be expressed as
$$
a=\sigma/\rho\left(\begin{array}{c}
\partial_x(n_y \partial_y f)-\partial_y(n_x \partial_y f)\\
\partial_y(n_x \partial_x f)-\partial_x(n_y \partial_x f)
\end{array}\right)
$$
*/      ia.x[] += alpha.x[]/fm.x[]*phif*CSS.x[];
/**
Finally, the acceleration is updated.
*/
	//ia.x[] += alpha.x[]/fm.x[]*phif*(f[] - f[-1])/Delta;
      }

  /**
  On trees, we need to restore the prolongation values for the
  volume fraction field. */
  
#if TREE
  for (scalar f in list) {
    f.prolongation = fraction_refine;
    f.dirty = true; // boundary conditions need to be updated
  }
#endif
  
  /**
  Finally we free the potential fields and the list of volume
  fractions. */

  for (scalar f in list) {
    scalar phi = f.phi;
    delete ({phi});
    f.phi.i = 0;
  }
  free (list);
}

/**
## References

See Section 3, pages 8-9 of:

~~~bib
@hal{popinet2018, hal-01528255}
~~~
*/
