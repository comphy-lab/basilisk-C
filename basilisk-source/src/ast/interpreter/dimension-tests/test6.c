/**
Checks dimension of returned constants. */

double func (double x)
{
  return 1.;
}

int main()
{
  {
    double a = 1 [1] + func (1);
  }
}
