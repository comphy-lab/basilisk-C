/**
# Checks the initialization of (postfix) (multidimensional) arrays */

void myfunc1 (double d[2])
{
  for (int i = 0; i < 2; i++)
    display_value (d[i]);
}

void myfunc (double a[2][2])
{
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      display_value (a[i][j]);
}

int main()
{{  interpreter_verbosity (2);
    double x1 = 1, y1 = 2, x2 = 3, y2 = 4;
    myfunc1 ((double[2]){x1, y1}); // this is a postfix unidimensional array
    myfunc ((double[2][2]) {{x1, y1},{x2, y2}}); // this is a postfix multidimensional array
    double b[2][2] = {{x1, y1},{x2, y2}}; // this is an init_declarator multidimensional array
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
	display_value (b[i][j]);
}}
