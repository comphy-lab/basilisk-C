%include "predictor-corrector.i"

%{
  extern double G;
  extern double dry;
  extern void conserve_elevation (void);
%}

extern void conserve_elevation (void);
