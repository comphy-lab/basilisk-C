//ACC(acc declare copyin(theta))
/* Compile minmod2 for device */
ACC(acc routine seq)
double mm2 (double s0, double s1, double s2, double theta)
{
  double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
  double result = 0.;
  if (s0 < s1 && s1 < s2) {
    if (d2 < d1) d1 = d2;
    //    return min(d1, d3);
    result = min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    if (d2 > d1) d1 = d2;
    //    return max(d1, d3);
    result = max(d1, d3);
  }
  //  return 0.;
  return result;
}


void gradientsfast (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
  scalar s; vector v;
  for (s,v in f,g) {
    if (s.gradient) {
      /* foreach(@) { */
      /*   foreach_dimension() */
      /*     v.x[] = minmod2 (s[-1], s[], s[1])/Delta; */
      /* } */

    size_t ngrid = ((Cartesian *)grid)->n;
    size_t nall = list_len(all)+1;

#pragma acc parallel loop gang num_workers(4) vector_length(32) independent present(acc_dataarr)
    for (size_t iacc = 1; iacc <= ngrid; iacc++) {
#pragma acc loop worker vector independent
      for (size_t jacc = 1; jacc <= ngrid; jacc++) {
	double theta = 1.3;
	Point pt;
	pt.n = ngrid;
	pt.i = iacc;
	pt.j = jacc;
	double Delta = L0*(1./pt.n);
	Delta *= Radius*3.14159265358979/180.;
	acc_dataarr[(pt.n + 2)*(pt.n + 2)*v.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = mm2 (acc_dataarr[(pt.n + 2)*(pt.n + 2)*s.i + (pt.i + -1)*(pt.n + 2) + (pt.j + 0)], acc_dataarr[(pt.n + 2)*(pt.n + 2)*s.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)], acc_dataarr[(pt.n + 2)*(pt.n + 2)*s.i + (pt.i + 1)*(pt.n + 2) + (pt.j + 0)], theta)/Delta;
	acc_dataarr[(pt.n + 2)*(pt.n + 2)*v.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = mm2 (acc_dataarr[(pt.n + 2)*(pt.n + 2)*s.i + (pt.i + 0)*(pt.n + 2) + (pt.j + -1)], acc_dataarr[(pt.n + 2)*(pt.n + 2)*s.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)], acc_dataarr[(pt.n + 2)*(pt.n + 2)*s.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 1)], theta)/Delta;
      }
    }
  }
    else // centered
      foreach(@) {
        foreach_dimension()
          v.x[] = (s[1] - s[-1])/(2.*Delta);
      }
  }
  boundary ((scalar *) g);
}


/**
# A solver for the Saint-Venant equations

*/

scalar zb[], h[], eta[];
vector u[];

double G = 1.;
double dry = 1e-10;

#include "predictor-corrector.h"

scalar * evolving = {h, u};

trace
static void advance_saint_venant (scalar * output, scalar * input, 
				  scalar * updates, double dt)
{
  // recover scalar and vector fields from lists
  scalar hi = input[0], ho = output[0], dh = updates[0];
  vector ui = vector(input[1]), uo = vector(output[1]), dhu = vector(updates[1]);

  // new fields in ho[], uo[]
  /* foreach(@) { */
  /*   double hold = hi[]; */
  /*   ho[] = hold + dt*dh[]; */
  /*   eta[] = zb[] + ho[]; */
  /*   if (ho[] > dry) */
  /*     foreach_dimension() */
  /* 	uo.x[] = (hold*ui.x[] + dt*dhu.x[])/ho[]; */
  /*   else */
  /*     foreach_dimension() */
  /* 	uo.x[] = 0.; */
  /* } */

  /* Inserted autgenerated code here */
  size_t ngrid = ((Cartesian *)grid)->n;
  size_t ngridall = (ngrid + 2)*(ngrid + 2);

  size_t uix_baseidx = ngridall*ui.x.i;
  size_t uiy_baseidx = ngridall*ui.y.i;

  size_t dhux_baseidx = ngridall*dhu.x.i;
  size_t dhuy_baseidx = ngridall*dhu.y.i;

  size_t hi_baseidx = ngridall*hi.i;
  size_t ho_baseidx = ngridall*ho.i;

  size_t dh_baseidx = ngridall*dh.i;
  size_t eta_baseidx = ngridall*eta.i;

  size_t zb_baseidx = ngridall*zb.i;

  size_t uox_baseidx = ngridall*uo.x.i;
  size_t uoy_baseidx = ngridall*uo.y.i;

#pragma acc parallel loop gang num_workers(16) vector_length(32) independent present(acc_dataarr)
  for (size_t iacc = 1; iacc <= ngrid; iacc++) {
#pragma acc loop worker vector independent
    for (size_t jacc = 1; jacc <= ngrid; jacc++) {

      // Compute general array index in variable block
      size_t idx = iacc*(ngrid + 2) + jacc;

      // Compute field indices
      size_t uix_idx = uix_baseidx + idx;
      size_t uiy_idx = uiy_baseidx + idx;

      size_t dhux_idx = dhux_baseidx + idx;
      size_t dhuy_idx = dhuy_baseidx + idx;

      size_t hi_idx = hi_baseidx + idx;
      size_t ho_idx = ho_baseidx + idx;

      size_t dh_idx = dh_baseidx + idx;
      size_t eta_idx = eta_baseidx + idx;

      size_t zb_idx = zb_baseidx + idx;

      size_t uox_idx = uox_baseidx + idx;
      size_t uoy_idx = uoy_baseidx + idx;

      // Get data
      double uix = acc_dataarr[uix_idx];
      double uiy = acc_dataarr[uiy_idx];

      double dhux = acc_dataarr[dhux_idx];
      double dhuy = acc_dataarr[dhuy_idx];

      // Original loop body
      double hold = acc_dataarr[hi_idx];

      double ho = hold + dt*acc_dataarr[dh_idx];
      acc_dataarr[eta_idx] = acc_dataarr[zb_idx] + ho;
      if (ho > dry)
	{
	  acc_dataarr[uox_idx] = (hold*uix + dt*dhux)/ho;
	  acc_dataarr[uoy_idx] = (hold*uiy + dt*dhuy)/ho;
	}
      else
	{
	  acc_dataarr[uox_idx] = 0.;
	  acc_dataarr[uoy_idx] = 0.;
	}
      acc_dataarr[ho_idx] = ho;
    }
  }

  // fixme: on trees eta is defined as eta = zb + h and not zb +
  // ho in the refine_eta() and coarsen_eta() functions below
  boundary ({ho, eta, uo});
}

#include "riemann.h"

trace
double update_saint_venant (scalar * evolving, scalar * updates, double dtmax)
{

  scalar h = evolving[0];
  vector u = vector(evolving[1]);

  face vector Fh[], S[];
  tensor Fq[];

  vector gh[], geta[];
  tensor gu[];
  for (scalar s in {gh, geta, gu}) {
    s.gradient = zero;
  }
  gradientsfast ({h, eta, u}, {gh, geta, gu});

  double _dtreduce = dtmax;
  int ig = 0, jg = 0;
  size_t ngrid = ((Cartesian *)grid)->n;
  size_t nall = list_len(all)+1;

#pragma acc update device(N, X0, Y0, L0, G, CFL)
#pragma acc parallel loop gang num_workers(8) vector_length(8) independent  present(acc_dataarr)  reduction(min:dtreduce)
    for (size_t iacc = 1; iacc <= ngrid + 1; iacc++) {
#pragma acc loop vector worker independent
      for (size_t jacc = 1; jacc <= ngrid + 1; jacc++) {
	Point pt;
	pt.n = ngrid;
	pt.i = iacc;
	pt.j = jacc;

	double Delta = L0*(1./pt.n);
	if ((pt.j <= pt.n)) {
	  Delta *= Radius*3.14159265358979/180.;
	  double hi = acc_dataarr[(pt.n + 2)*(pt.n + 2)*h.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)], hn = acc_dataarr[(pt.n + 2)*(pt.n + 2)*h.i + (pt.i + -1)*(pt.n + 2) + (pt.j + 0)], dt = _dtreduce;
	  if (hi > dry || hn > dry) {
	    double dx = Delta/2.;
	    double zi = acc_dataarr[(pt.n + 2)*(pt.n + 2)*eta.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] - hi;
	    double zl = zi - dx*(acc_dataarr[(pt.n + 2)*(pt.n + 2)*geta.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*gh.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]);
	    double zn = acc_dataarr[(pt.n + 2)*(pt.n + 2)*eta.i + (pt.i + -1)*(pt.n + 2) + (pt.j + 0)] - hn;
	    double zr = zn + dx*(acc_dataarr[(pt.n + 2)*(pt.n + 2)*geta.x.i + (pt.i + -1)*(pt.n + 2) + (pt.j + 0)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*gh.x.i + (pt.i + -1)*(pt.n + 2) + (pt.j + 0)]);
	    double zlr = ((zl) > (zr) ? (zl) : (zr));
	    double hl = hi - dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gh.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)];
	    double up = acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] - dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gu.x.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)];
	    double hp = ((0.) > (hl + zl - zlr) ? (0.) : (hl + zl - zlr));
	    double hr = hn + dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gh.x.i + (pt.i + -1)*(pt.n + 2) + (pt.j + 0)];
	    double um = acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.x.i + (pt.i + -1)*(pt.n + 2) + (pt.j + 0)] + dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gu.x.x.i + (pt.i + -1)*(pt.n + 2) + (pt.j + 0)];
	    double hm = ((0.) > (hr + zr - zlr) ? (0.) : (hr + zr - zlr));
	    double fh, fu, fv;
	    kurganov (hm, hp, um, up, Delta*acc_dataarr[(pt.n + 2)*(pt.n + 2)*cm.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]/acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)], &fh, &fu, &dt);
	    fv = (fh > 0. ? acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.y.i + (pt.i + -1)*(pt.n + 2) + (pt.j + 0)] + dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gu.y.x.i + (pt.i + -1)*(pt.n + 2) + (pt.j + 0)] : acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] - dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gu.y.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)])*fh;
	    double sl = G/2.*(((hp)*(hp)) - ((hl)*(hl)) + (hl + hi)*(zi - zl));
	    double sr = G/2.*(((hm)*(hm)) - ((hr)*(hr)) + (hr + hn)*(zn - zr));
	    acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fh.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*fh;
	    acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.x.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*(fu - sl);
	    acc_dataarr[(pt.n + 2)*(pt.n + 2)*S.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*(fu - sr);
	    acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.y.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*fv;
	  }
	  else
	    acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fh.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.x.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*S.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.y.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = 0.;
	  _dtreduce = ((_dtreduce) < (dt) ? (_dtreduce) : (dt));
	}

	//
      }
    }
    if (_dtreduce < dtmax) dtmax = _dtreduce;
    _dtreduce = dtmax;

#pragma acc parallel loop gang num_workers(8) vector_length(8) independent  present(acc_dataarr)  reduction(min:dtreduce)
    for (size_t iacc = 1; iacc <= ngrid + 1; iacc++) {
#pragma acc loop vector worker independent
      for (size_t jacc = 1; jacc <= ngrid + 1; jacc++) {

	Point pt;
	pt.n = ngrid;
	pt.i = iacc;
	pt.j = jacc;

	//

	double Delta = L0*(1./pt.n);
	if ((pt.i <= pt.n)) {
	  Delta *= Radius*3.14159265358979/180.;
	  double hi = acc_dataarr[(pt.n + 2)*(pt.n + 2)*h.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)], hn = acc_dataarr[(pt.n + 2)*(pt.n + 2)*h.i + (pt.i + 0)*(pt.n + 2) + (pt.j + -1)], dt = _dtreduce;
	  if (hi > dry || hn > dry) {
	    double dx = Delta/2.;
	    double zi = acc_dataarr[(pt.n + 2)*(pt.n + 2)*eta.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] - hi;
	    double zl = zi - dx*(acc_dataarr[(pt.n + 2)*(pt.n + 2)*geta.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*gh.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]);
	    double zn = acc_dataarr[(pt.n + 2)*(pt.n + 2)*eta.i + (pt.i + 0)*(pt.n + 2) + (pt.j + -1)] - hn;
	    double zr = zn + dx*(acc_dataarr[(pt.n + 2)*(pt.n + 2)*geta.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + -1)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*gh.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + -1)]);
	    double zlr = ((zl) > (zr) ? (zl) : (zr));
	    double hl = hi - dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gh.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)];
	    double up = acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] - dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gu.y.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)];
	    double hp = ((0.) > (hl + zl - zlr) ? (0.) : (hl + zl - zlr));
	    double hr = hn + dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gh.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + -1)];
	    double um = acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + -1)] + dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gu.y.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + -1)];
	    double hm = ((0.) > (hr + zr - zlr) ? (0.) : (hr + zr - zlr));
	    double fh, fu, fv;
	    kurganov (hm, hp, um, up, Delta*acc_dataarr[(pt.n + 2)*(pt.n + 2)*cm.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]/acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)], &fh, &fu, &dt);
	    fv = (fh > 0. ? acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + -1)] + dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gu.x.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + -1)] : acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] - dx*acc_dataarr[(pt.n + 2)*(pt.n + 2)*gu.x.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)])*fh;
	    double sl = G/2.*(((hp)*(hp)) - ((hl)*(hl)) + (hl + hi)*(zi - zl));
	    double sr = G/2.*(((hm)*(hm)) - ((hr)*(hr)) + (hr + hn)*(zn - zr));
	    acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fh.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*fh;
	    acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.y.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*(fu - sl);
	    acc_dataarr[(pt.n + 2)*(pt.n + 2)*S.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*(fu - sr);
	    acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.x.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*fv;
	  }
	  else
	    acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fh.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.y.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*S.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.x.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = 0.;
	  _dtreduce = ((_dtreduce) < (dt) ? (_dtreduce) : (dt));
	}
      }
    }
    if (_dtreduce < dtmax) dtmax = _dtreduce;

  boundary_flux ({Fh, S, Fq});
  
  scalar dh = updates[0];
  vector dhu = vector(updates[1]);

  /* foreach(@) { */
  /*   dh[] = (Fh.x[] + Fh.y[] - Fh.x[1,0] - Fh.y[0,1])/(cm[]*Delta); */
  /*   foreach_dimension() */
  /*     dhu.x[] = (Fq.x.x[] + Fq.x.y[] - S.x[1,0] - Fq.x.y[0,1])/(cm[]*Delta); */

  /*   double dmdl = (fm.x[1,0] - fm.x[])/(cm[]*Delta); */
  /*   double dmdt = (fm.y[0,1] - fm.y[])/(cm[]*Delta); */
  /*   double fG = u.y[]*dmdl - u.x[]*dmdt; */
  /*   dhu.x[] += h[]*(G*h[]/2.*dmdl + fG*u.y[]); */
  /*   dhu.y[] += h[]*(G*h[]/2.*dmdt - fG*u.x[]); */
  /* } */

#pragma acc parallel loop gang num_workers(4) vector_length(32) independent present(acc_dataarr)
  for (size_t iacc = 1; iacc <= ngrid; iacc++) {
#pragma acc loop worker vector independent
    for (size_t jacc = 1; jacc <= ngrid; jacc++) {
      Point pt;
      pt.n = ngrid;
      pt.i = iacc;
      pt.j = jacc;
      double Delta = L0*(1./pt.n);
      Delta *= Radius*3.14159265358979/180.;
      acc_dataarr[(pt.n + 2)*(pt.n + 2)*dh.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = (acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fh.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] + acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fh.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fh.x.i + (pt.i + 1)*(pt.n + 2) + (pt.j + 0)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fh.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 1)])/(acc_dataarr[(pt.n + 2)*(pt.n + 2)*cm.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*Delta);
      acc_dataarr[(pt.n + 2)*(pt.n + 2)*dhu.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = (acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.x.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] + acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.x.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*S.x.i + (pt.i + 1)*(pt.n + 2) + (pt.j + 0)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.x.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 1)])/(acc_dataarr[(pt.n + 2)*(pt.n + 2)*cm.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*Delta);
      acc_dataarr[(pt.n + 2)*(pt.n + 2)*dhu.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] = (acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.y.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] + acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.y.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*S.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 1)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*Fq.y.x.i + (pt.i + 1)*(pt.n + 2) + (pt.j + 0)])/(acc_dataarr[(pt.n + 2)*(pt.n + 2)*cm.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*Delta);
      double dmdl = (acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.x.i + (pt.i + 1)*(pt.n + 2) + (pt.j + 0)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)])/(acc_dataarr[(pt.n + 2)*(pt.n + 2)*cm.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*Delta);
      double dmdt = (acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 1)] - acc_dataarr[(pt.n + 2)*(pt.n + 2)*fm.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)])/(acc_dataarr[(pt.n + 2)*(pt.n + 2)*cm.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*Delta);
      double fG = acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*dmdl - acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*dmdt;
      acc_dataarr[(pt.n + 2)*(pt.n + 2)*dhu.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] += acc_dataarr[(pt.n + 2)*(pt.n + 2)*h.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*(G*acc_dataarr[(pt.n + 2)*(pt.n + 2)*h.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]/2.*dmdl + fG*acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]);
      acc_dataarr[(pt.n + 2)*(pt.n + 2)*dhu.y.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)] += acc_dataarr[(pt.n + 2)*(pt.n + 2)*h.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]*(G*acc_dataarr[(pt.n + 2)*(pt.n + 2)*h.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]/2.*dmdt - fG*acc_dataarr[(pt.n + 2)*(pt.n + 2)*u.x.i + (pt.i + 0)*(pt.n + 2) + (pt.j + 0)]);
    }
  }

  return dtmax;
}

event defaults (i = 0)
{

  advance = advance_saint_venant;
  update = update_saint_venant;
}

event init (i = 0)
{
  foreach(@)
    eta[] = zb[] + h[];
  boundary (all);
}

#include "elevation.h"
