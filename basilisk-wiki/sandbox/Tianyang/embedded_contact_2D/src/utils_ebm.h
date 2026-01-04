norm normf_no_cm (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
  foreach(reduction(max:max) reduction(+:avg) 
	  reduction(+:rms) reduction(+:volume)) 
    if (f[] != nodata) {
      double v = fabs(f[]);
      if (v > max) max = v;
      volume += sq(Delta);
      avg    += sq(Delta)*v;
      rms    += sq(Delta)*sq(v);
    }
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}

stats statsf_no_cm (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
	  reduction(max:max) reduction(min:min)) 
    if (dv() > 0. && f[] != nodata) {
      volume += sq(Delta);
      sum    += sq(Delta)*f[];
      sum2   += sq(Delta)*sq(f[]);
      if (f[] > max) max = f[];
      if (f[] < min) min = f[];
    }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}

coord normal_contact (coord ns, coord nf, double angle)
{
  assert (dimension == 2); // fixme: 2D only for the moment
  coord n;
  if (- ns.x*nf.y + ns.y*nf.x > 0) { // 2D cross product
    n.x = - ns.x*cos(angle) + ns.y*sin(angle);
    n.y = - ns.x*sin(angle) - ns.y*cos(angle);
  }
  else {
    n.x = - ns.x*cos(angle) - ns.y*sin(angle);
    n.y =   ns.x*sin(angle) - ns.y*cos(angle);
  }
  return n;
}