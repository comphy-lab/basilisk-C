%include "run.i"
%include "timestep.i"

%{
  extern double (* gradient) (double, double, double);
  extern scalar * tracers;
%}

%include "tracer.i"
