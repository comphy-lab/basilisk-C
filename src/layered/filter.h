/**
# Fourth-order filter

This is not ready for general consumption. Just kept for reference. */

double filter = 1./16.; // maximum filtering for fourth-order filter, see Klemp

event viscous_term (i++)
{
  // fourth-order filter
#if dimension == 1
  foreach()
    foreach_layer()
      for (scalar s in {w,u})
	s[] -= (s[2] + s[-2] - 4.*(s[1] + s[-1]) + 6.*s[])*filter;
#else // dimension == 2
  // this is the 4th-order diagonal-term operator of [Klemp, 2017](#klemp2017)
  foreach()
    if (x > 114 && x < 169 && y > 10 && y < 65)
    foreach_layer()
      for (scalar s in {w,u})
        s[] -= ((s[2,2] + s[2,-2] + s[-2,-2] + s[-2,2])/16. +
		(s[-1,2] + s[1,2] + s[-1,-2] + s[1,-2] +
		 s[2,-1] + s[2,1] + s[-2,-1] + s[-2,1])/4. +
		3.*(s[2,0] + s[-2,0] + s[0,2] + s[0,-2])/8. -
		(s[1,1] + s[1,-1] + s[-1,-1] + s[-1,1]) -
		5.*(s[1,0] + s[-1,0] + s[0,1] + s[0,-1])/2. +
		41.*s[]/4.)*filter;
#endif // dimension == 2
}

/**
## References

~~~bib
@article{klemp2017,
  title={Damping characteristics of horizontal Laplacian diffusion filters},
  author={Klemp, Joseph B},
  journal={Monthly Weather Review},
  volume={145},
  number={11},
  pages={4365--4379},
  year={2017}
}
~~~
*/
