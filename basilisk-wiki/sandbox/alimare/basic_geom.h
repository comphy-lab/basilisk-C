double rectangle(double x, double y, coord center, coord size) {

  double P1_Plus = x - size.x/2.  - center.x;
  double P1_Minus = x + size.x/2. - center.x;
  double P1 = max (P1_Plus, -P1_Minus);

  double P2_Plus = y - size.y/2.  - center.y;
  double P2_Minus = y + size.y/2. - center.y;
  double P2 = max (P2_Plus, -P2_Minus);

  double c = max ( P1,P2 );
  return c;
}

double circle(double x, double y,  coord center, double radius) {
  double R2  =  sq(x - center.x) + sq (y - center.y);
  return ( sqrt(R2) - radius);
}

double ellipse(double x, double y,  coord center, double radius,
               double a, double b) {
  double R2  =  sq((x - center.x)/a) + sq ((y - center.y)/b);
  return ( sqrt(R2) - radius);
}

double sphere(double x, double y,  double z, coord center, double radius) {
  double R2  =  sq(x - center.x) + sq (y - center.y) + sq (z - center.z);
  return ( sqrt(R2) - radius);
}


// TODO: cylinder