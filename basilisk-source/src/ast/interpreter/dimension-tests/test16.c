/**
# Dimensions in constant expressions */

int main()
{

  /**
  Constant expressions can be dimensioned. */
  
  double a = (1. + 2.)[1];
  display_value (a);
  double b = 3. + a;

  /**
  Constants in constant expressions can also be dimensioned (but this
  is silly and confusing). */
  
  double c = 4.[1] + 5.;
  double d = 6. + c;
}
