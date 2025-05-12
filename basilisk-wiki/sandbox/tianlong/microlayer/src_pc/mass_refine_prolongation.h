void mass_refine (Point point, scalar s)
{
  scalar f = s.c;
  double s_parent = s[];
  foreach_child() {
  if (f[] > F_ERR && f[] < 1. - F_ERR)
    s[] = s_parent;
    else
    s[] = 0.;
  }
}

void refinement_avg (Point point, scalar s) {
  scalar f = s.c;
  double sp = s[];
  foreach_child() {
    if (f[] > F_ERR && f[] < 1. - F_ERR)
      s[] = sp;
      else
      s[] = 0.;
  }
}