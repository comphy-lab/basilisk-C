/**
# Filter for grid-scale oscillations

**Note that this is obsolete and is kept only for historical interest.**

Because it is colocated, the [layered solver](hydro.h) can be prone to
grid-scale oscillations. This noise is usually small but can be
significant when gradients of bathymetry/wave fields are
under-resolved.

The filter below can be used to reduce this noise. It is inspired from
the work of [Engquist et al, 1989](#engquist1989).

The strength of the filter is controlled by the filtering timescale
`filter`: the smaller the value, the stronger the filter. This can be
interpreted as the timescale during which a disturbance must not move
to be seen by the filter. */

double filter = 0.;

event viscous_term (i++)
{
  if (filter > 0.) {
    foreach()
      foreach_dimension() {
        double Hm = 0., H = 0., Hp = 0.;
	foreach_layer()
	  Hm += h[-1], H += h[], Hp += h[1];
        if (Hm > dry && H > dry && Hp > dry) {
	  double Dp = eta[1] - eta[], Dm = eta[] - eta[-1];

	  /**
	  The filter is only applied to extrema in the free-surface
	  height $\eta$ (first condition) contained within two
	  inflection points (second condition): this effectively
	  selects only grid-scale oscillations and avoid filtering
	  smooth extrema. */
	  
	  if (Dp*Dm < 0. && ((eta[2] + eta[] - 2.*eta[1])*
			     (eta[-1] + eta[1] - 2.*eta[]) < 0. ||
			     (eta[-1] + eta[1] - 2.*eta[])*
			     (eta[-2] + eta[] - 2.*eta[-1]) < 0.)) {

	    /**
	    We then compute the shift in the value of $\eta$ necessary
	    to smooth the extremum, see Algorithm 2.1 in [Engquist et
	    al, 1989](#engquist1989). */ 
	    
	    double dp, dm;
	    if (fabs(Dp) > fabs(Dm)) {
	      dp = fabs(Dp);
	      dm = fabs(Dm);
	    }
	    else {
	      dp = fabs(Dm);
	      dm = fabs(Dp);
	    }
	    double d = min(dm, dp/2.);
	    double a = Dp > 0. ? 1. : -1.;

	    /**
	    We apply only part of the correction, weighted by the
	    timescale. */
	    
	    eta[] += min(dt/filter, 1.)*a*d;
	    double Hnew = eta[] - zb[];
	    if (Hnew > dry) {
	      foreach_layer()
		h[] *= Hnew/H;
	    }
	    else {
	      for (scalar s in tracers)
		foreach_layer()
		  s[] = 0.;
	    }
	  }
	}
      }
  }
}

/**
## References

~~~bib
@article{engquist1989,
  title={Nonlinear filters for efficient shock computation},
  author={Engquist, Bj{\"o}rn and L{\"o}tstedt, Per and 
          Sj{\"o}green, Bj{\"o}rn},
  journal={Mathematics of Computation},
  volume={52},
  number={186},
  pages={509--537},
  year={1989},
  url={https://www.ams.org/journals/mcom/1989-52-186/S0025-5718-1989-0955750-9/S0025-5718-1989-0955750-9.pdf}
}
~~~
*/
