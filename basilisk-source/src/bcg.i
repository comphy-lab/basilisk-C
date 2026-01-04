%{
  extern void tracer_fluxes (const scalar f, 
			     const face vector u,
			     face vector flux,
			     double dt,
			     scalar src);
  
  struct Advection {
    scalar * tracers;
    face vector u;
    double dt;
    scalar * src; // optional
  };
  
  extern void advection (struct Advection p);
%}
