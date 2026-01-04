/**
Tests function declarations and pointers */

double func (coord q)
{
  double c = q.z*5, d = c*q.y;
  return c + d;

  return c - d;
}

double func (coord p);

int main()
{{ interpreter_verbosity (4);
    
  coord p, p2 = {0,0,5};
  p.z = 3;
  
  double a = 3, b = 4, c = a*p.z, d = a*p2.z;

  void (* pfunc) (coord) = func;
  
  double q = pfunc (p2);

  a = 3*q;

  double r[3];  

  p = p2;
  a = p.z*3;

  printf ("%g\n", a);

  {
    double * u = r + 1;
    
    u += 2;
    u[-2] = (a++)*b;
    
    a = a + u[-2];
    
    c = 2*a;
  }
  
  {
    double * u = r + 2;
    u[-2] = a*b;
    a = 2*u[-2];
    
    double * p = &a;
    
    b = 1*p[0];
  }
}}
