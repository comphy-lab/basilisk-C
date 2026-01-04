/**
Checks that indiviual structure members are properly initialized. */

int main()
{
  {
    interpreter_verbosity (4);
    coord a = {0}, c = {1};
    double b = 0;
    if (b)
      a = c;
    if (a.y); // must be defined and zero
  }
}
