/**
# No warning for Point point

There should be a warning or even an error according to [/src/NEWS](). 

*This is actually OK. The error will occur only if one tries to access fields, coordinates etc. after this call.* 

*/

int main()
{
  init_grid (1);
  Point point = locate (0);
  double a = x;
}
