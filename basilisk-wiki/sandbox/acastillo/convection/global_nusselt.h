#if dimension == 2
#define dA Delta
#define dV sq(Delta)
#endif

#if dimension == 3
#define dA sq(Delta)
#define dV pow(Delta,3)
#endif

double nusselt_top (scalar T)
{
  double nutop = 0.;
  foreach_boundary (top,reduction(+:nutop))
    nutop += dA*(T[] - T[ghost])/Delta;
  return nutop;
}

double nusselt_bot (scalar T)
{
  double nubot=0.;
  foreach_boundary (bottom,reduction(+:nubot))
    nubot += dA*(T[ghost] - T[])/Delta;
    return nubot;
}

double nusselt_right (scalar T)
{
  double nuright=0.;
  foreach_boundary (right,reduction(+:nuright))
    nuright += dA*(T[] - T[ghost])/Delta;
  return nuright;
}

double nusselt_left (scalar T)
{
  double nuleft=0.;
  foreach_boundary (left,reduction(+:nuleft))
    nuleft += dA*(T[ghost] - T[])/Delta;
  return nuleft;
}

double nusselt_vol (scalar T, vector u)
{
  double nuvol=0.;
  foreach(reduction(+:nuvol))
    nuvol += dV*(sqrt(Ra)*u.y[]*T[] - (T[0,1] - T[0,-1])/(2*Delta));
  return nuvol;
}
