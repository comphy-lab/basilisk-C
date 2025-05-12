double nusselt_top (scalar T)
{
  double nutop=0.;
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


double nusselt_top2 (scalar T)
{
  double nutop2 =0.;
  foreach_boundary (top,reduction(+:nutop2))
    nutop2 += dA*(T[] - Ttop)/(Delta/2.);
  return nutop2;
}


double nusselt_bot2 (scalar T)
{
  double nubot2=0.;
  foreach_boundary (bottom,reduction(+:nubot2))
    nubot2 += dA*(Tbot - T[])/(Delta/2.);
  return nubot2;
}

double nusselt_top3 (scalar T)
{
  double nutop3=0.;
  foreach_boundary (top,reduction(+:nutop3))
    nutop3 += dA*(4*T[0,-1]-3*T[0,0]-T[0,-2])/(2*Delta);
  return nutop3;
}

double nusselt_bot3 (scalar T)
{
  double nubot3=0.;
  foreach_boundary (bottom,reduction(+:nubot3))
    nubot3 += -dA*(4*T[0,1]-3*T[0,0]-T[0,2])/(2*Delta);
  return nubot3;
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