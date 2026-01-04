/**
Checks minmod2() */

#include "utils.h"

int main()
{
  {
    double a = 1.[1], b = 2.[1], c = 3.[1];
    unset_double (&a);
    unset_double (&b);
    unset_double (&c);
    display_value (minmod2 (a, b, c));
    
    a = 1.[2], b = 2.[2], c = 3.[2];
    unset_double (&a);
    unset_double (&b);
    unset_double (&c);
    display_value (minmod2 (a, b, c));
  }
}
