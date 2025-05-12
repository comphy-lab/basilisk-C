/**
# Law of the Wall
Here we implement a closure for turbulent flow over a rough surface according to the so-called ["log-law of the wall"](https://en.wikipedia.org/wiki/Law_of_the_wall). The idea is that fully resolving a turbulent boundary layer near a rough surface within a larger domain may be unfeasable for very high Reynodlsnumber flows. This is expected to be the case when the grid resolution is much larger than the typical size the roughness elements. Rather ad hoc, the effect of a wall can be summarized as that is exerts a drag force (momentum flux) on the flow. 

We use the well established result for a profile of a horizontal velocity profile ($u(y)$) in a wall-bounded turbulent flow in the so-called log layer:

$$u(y) = \frac{u_{\tau}}{k}\text{ln}\left( \frac{y}{y_0}\right),$$

where $\text{ln}\left( x \right)$ is the natural logarithm of a dummy variable $x$, $u_\tau$ is the friction velocity, $k = 0.4$ the Von Karmann constant, $y$ is the height above the roughness elements that are characterezed by a so-called 'roughness length' $y_0$, i.e the height at which the velocity vanishes. 

Using this, we write for the discretized horizontal velocity in the cell ($U_n = \sqrt{u_x^2 + u_z^2}$) adjacent to the rough and no-slip boundary,

$$U_n = \frac{\int_{y_0}^{y_0+\Delta} u \text{d}y}{\Delta} = \frac{u_\tau}{k\Delta}\left[ y\left( \text{ln}\left( \frac{y}{y_0}\right) -1\right) \right]_{y_0}^{y_0+\Delta}. $$

When we assume $\Delta \gg y_0$, we can rewrite this result as,

$$u_\tau\approx \frac{kU_n}{\text{ln}\left( \frac{\Delta }{y_0} \right)-1 }. $$

Using Pythagoras, 

$$u_{\tau,x}^2 + u_{\tau,z}^2 = u_\tau^2,$$

we may write,

$$u_{\tau,i}^2 = \left( \frac{u_i}{U_n} \right) ^2 u_\tau^2$$

which closes the problem, as we have expressed the stress as a function of the known values and constants $u_x, u_z, y_0, \Delta$ and $k$

## implementation

We assume there is no heterogeneity for $y_0$ and in order to prevent an excessive number of calculations of the logaritmic function, we make use of a look-up table (`lut`) containing values for $\left( \frac{k}{\text{ln} \left( \frac{\Delta}{y_0} \right) - 1} \right) ^ 2$ that can be used for adaptive and non-adaptive simulations using upto 19 levels of refinement. Noting that it is the user's responsibility not to apply the closure when the assumptions do not make sense. I.e. e.g. in [the Roughness sublayer](http://glossary.ametsoc.org/wiki/Roughness_sublayer).
*/
#ifndef Y_0
#define Y_0  (0.01)
#endif
double lut[20];

event init (t = 0){
  double k = 0.4;
  lut[0] = 0.; //level > 0
  for (int m = 1; m <= 19; m++)
    lut[m] = sq(k/(log(L0/(((double)(1<<m))*Y_0)) - 1.));
}
/**
Now we use the [midpoint method](https://en.wikipedia.org/wiki/Midpoint_method) to advance in time as we take into account the surface drag at the `bottom` surface.  
*/
#define dvdt(u) (sign(u)*lut[level]*sq(u)/Delta)
event law_of_the_wall(i++){
  double ui; //scratch for the mid-point-value estimate 
  foreach_boundary(bottom){
    ui = u.x[] - (dt*dvdt(u.x[])/2.); 
    u.x[] -= dt*dvdt(ui);
    ui = u.z[] - (dt*dvdt(u.z[])/2.);
    u.z[] -= dt*dvdt(ui);
  }
}


/**
## Tests

* [The concept applied in a SCM for the Neutral Trubulent Ekman layer (SCM)](turbulentekman.c)
* There should be more and convincing tests
*/ 