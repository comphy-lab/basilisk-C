/**
Checks that the interpreter correctly converts doubles/floats to
bool. */

int main()
{
  double a = 0.1;
  if (a)
    display_value (1); // must display 1
  if (0.1 || 0.2)
    display_value (1); // must display 1
}
