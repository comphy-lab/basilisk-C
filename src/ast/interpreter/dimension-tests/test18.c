/**
# Check the consistency of the dimensions of return values */

double func1()
{
  double ret;
  if (a)
    ret = 1.;
  if (ret == 1.)
    ret = 2.;
  return ret;
}

double func2()
{
  double ret;
  if (a)
    ret = 1.;
  if (ret == 1.)
    return 2.;
  return ret;
}

coord func3()
{
  coord ret;
  foreach_dimension() {
    if (a)
      ret.x = 1.;
    if (ret.x == 1.)
      ret.x = 2.;
  }
  return ret;
}

int main()
{
  {
    interpreter_verbosity (3);
    func1() == 3. [1];
    func2() == 3. [2];
    foreach_dimension()
      func3().x == 3. [3];
  }
}
