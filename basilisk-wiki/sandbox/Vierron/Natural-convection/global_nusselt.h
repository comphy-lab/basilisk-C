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
  //double Ttop = 0.;
  //double Tbot = 0.;

  foreach_boundary (top,reduction(+:nutop))
    nutop += dA*(T[] - T[ghost])/Delta;
#if 0
  foreach_boundary (top,reduction(+:Ttop))
    Ttop += (T[] + T[ghost])/2.;

  foreach_boundary (bottom,reduction(+:Tbot))
    Tbot += (T[] + T[ghost])/2.;

  //return nutop/L0*(Tbot-Ttop);
#endif
  return nutop/L0;
}


double nusselt_bot (scalar T)
{
  double nubot=0.;
  //double Ttop = 0.;
  //double Tbot = 0.;

  foreach_boundary (bottom,reduction(+:nubot))
    nubot += dA*(T[ghost] - T[])/Delta;
#if 0
  foreach_boundary (top,reduction(+:Ttop))
    Ttop += (T[] + T[ghost])/2.;

  foreach_boundary (bottom,reduction(+:Tbot))
    Tbot += (T[] + T[ghost])/2.;

  //return nubot/L0*(Tbot-Ttop);
#endif
  return nubot/L0;
}

double nusselt_vol (scalar T, vector u)
{
  double nuvol=0.;
  foreach(reduction(+:nuvol))
    nuvol += dV*(sqrt(Ra)*u.y[]*T[] - (T[0,1] - T[0,-1])/(2*Delta));
  return nuvol/L0;
}
