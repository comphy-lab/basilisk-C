#define TWO_PHASE_VARPROP 1

scalar rho1v[], rho2v[], mu1v[], mu2v[];

event defaults (i = 0) {
  rho1v.nodump = rho2v.nodump = true;
  mu1v.nodump = mu2v.nodump = true;
}

extern scalar f, rhov;
extern face vector alphav;
#if FILTERED
extern scalar sf;
#else
# define sf f
#endif

#define aavg(f,val1,val2) (clamp(f,0,1)*val1 + (1.-clamp(f,0,1))*val2)

event properties (i++)
{
  foreach_face() {
    double ff = face_value (sf, 0);
    double rho1f = face_value (rho1v, 0) + 1e-10;
    double rho2f = face_value (rho2v, 0) + 1e-10;
    alphav.x[] = fm.x[]/aavg (ff, rho1f, rho2f);

    double mu1f = face_value (mu1v, 0);
    double mu2f = face_value (mu2v, 0);
    face vector muv = mu;
    muv.x[] = fm.x[]*aavg (ff, mu1f, mu2f);
  }

  foreach()
    rhov[] = cm[]*aavg (sf[], rho1v[], rho2v[]);
}

