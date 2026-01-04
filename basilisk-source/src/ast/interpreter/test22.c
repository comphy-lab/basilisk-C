/**
# Checks that logical operations are correctly unset

This is a bit tricky because of "partial evaluation". */

int main()
{
  int b = 1, a;
  display_value (b && a);
  display_value (b || a);
  display_value (a && b);
  display_value (a || b);
  b = 0;
  display_value (b && a);
  display_value (b || a);
  display_value (a && b);
  display_value (a || b);
}
