#ifndef _MY_BC_H
#define _MY_BC_H

#define myneumann(val) myneumann_(point, neighbor, _s, val)

foreach_dimension()
static double myneumann_x (Point point, scalar s, double val, int shift)
{
  int ilayer = abs(shift);
  if(ilayer > 2 || ilayer < 1)
  {
    printf("erro in boundary conditions!\n");
    exit(1);
  }
  int sig = sign(shift);
  int id_val_r = 1 * sig - shift;
  double dn = ilayer == 1 ? Delta : 3.0 * Delta;
  double val_r = s[id_val_r];
  return val_r + sig * val * dn;
}

double myneumann_(Point point, Point neighbor, scalar s, double val)
{
  if (neighbor.i != point.i)
    return myneumann_x (point, s, val, neighbor.i - point.i);
  if (neighbor.j != point.j)
    return myneumann_y (point, s, val, neighbor.j - point.j);
#if(dimension == 3)
  if (neighbor.k != point.k)
    return myneumann_z (point, s, val, neighbor.k - point.k);
#endif
  assert (false); // not reached
  return 0.;
}

#define mydirichlet(val) mydirichlet_(point, neighbor, _s, val)

foreach_dimension()
static double mydirichlet_x (Point point, scalar s, double val, int shift)
{
  int ilayer = abs(shift);
  if(ilayer > 2 || ilayer < 1)
  {
    printf("erro in boundary conditions!\n");
    exit(1);
  }
  int sig = sign(shift);
  int id_val_r = 1 * sig - shift;
  double dn = ilayer == 1 ? Delta : 3.0 * Delta;
  double val_r = s[id_val_r];
  return 2.0 * val  - val_r;
}

double mydirichlet_(Point point, Point neighbor, scalar s, double val)
{
  if (neighbor.i != point.i)
    return mydirichlet_x (point, s, val, neighbor.i - point.i);
  if (neighbor.j != point.j)
    return mydirichlet_y (point, s, val, neighbor.j - point.j);
#if(dimension == 3)
  if (neighbor.k != point.k)
    return mydirichlet_z (point, s, val, neighbor.k - point.k);
#endif
  assert (false); // not reached
  return 0.;
}

#define myOutletP() myOutletP_(point, neighbor, _s)

foreach_dimension()
static double myOutletP_x (Point point, scalar s, int shift)
{
  int ilayer = abs(shift);
  if(ilayer > 2 || ilayer < 1)
  {
    printf("erro in boundary conditions!\n");
    exit(1);
  }
  int sig = sign(shift);
  int i1 = 0;
  int i2 = 0;
  if(sig > 0)
    i2 = -1;
  else
    i2 = 1;
  
  double ghost1 = 2.0 * s[i1] - s[i2];
  double ghost2 = 2.0 * ghost1 - s[i1];

  double final = ilayer == 1 ? ghost1 : ghost2;
  
  return final;
}

double myOutletP_(Point point, Point neighbor, scalar s)
{
  if (neighbor.i != point.i)
    return myOutletP_x (point, s, neighbor.i - point.i);
  if (neighbor.j != point.j)
    return myOutletP_y (point, s, neighbor.j - point.j);
#if(dimension == 3)
  if (neighbor.k != point.k)
    return myOutletP_z (point, s, neighbor.k - point.k);
#endif
  assert (false); // not reached
  return 0.;
}


#endif