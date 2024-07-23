%{
  struct _interpolate {
    scalar v;
    double x, y, z;
  };
  extern double interpolate (struct _interpolate p);
%}

%inline %{
  typedef struct {
    double avg, rms, max, area;
  } norm;
  extern norm normf (scalar f);

  typedef struct {
    double min, max, sum, stddev, area;
  } stats;
  extern stats statsf (scalar f);

  extern void vorticity (const vector u, scalar omega);
%}

%apply (double * IN_ARRAY1, int DIM1) {(double * x, int len1)};
%apply (double * ARGOUT_ARRAY1, int DIM1) {(double * val, int len)};
%inline %{
  void _interpolate1D (scalar v, double * x, int len1, double * val, int len) {
    int i;
    for (i = 0; i < len; i++) {
      struct _interpolate p = {v, 0.9999999999*x[i]};
      val[i] = interpolate (p);
    }
  }
%}

%apply (double * IN_ARRAY2, int DIM1, int DIM2) {
  (double * x, int len3, int len4),
  (double * y, int len5, int len6)
}
%apply (double * INPLACE_ARRAY2, int DIM1, int DIM2) {
  (double * val, int len1, int len2)
}
%inline %{
  void _interpolate2D (scalar v, 
                       double * x, int len3, int len4,
                       double * y, int len5, int len6,
                       double * val, int len1, int len2) {
    int i;
    for (i = 0; i < len1*len2; i++) {
      struct _interpolate p = {v, 0.9999999999*x[i], 0.9999999999*y[i]};
      val[i] = interpolate (p);
    }
  }
%}
