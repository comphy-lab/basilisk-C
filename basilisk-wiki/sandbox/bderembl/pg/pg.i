%include "predictor-corrector.i"
%include "poisson.i"
%include "timestep.i"

%apply (double * IN_ARRAY1, int DIM1) {
  (double * val1, int len1)
}
%inline %{
  void pyset_field (int ifield, double * val1, int len1);
%}
%apply (double * INPLACE_ARRAY1, int DIM1) {
  (double * val2, int len2)
}
%inline %{
  void pyget_field (int ifield, double * val2, int len2);
  void pystep ( double * val1, int len1,
                double * val2, int len2);
%}
%apply (double * IN_ARRAY1, int DIM1) {
  (double * val3, int len3)
}
%inline %{
void pystep_lt ( double * val1, int len1,
		 double * val2, int len2,
		 double * val3, int len3);
%}

%inline %{
  void pyset_vars ();
  void pytrash_vars ();
  void pyset_contpar(int pycontpar);
  void pyadjust_contpar(double contpar_val);
  void pyinit_vertgrid(int pynl);
%}