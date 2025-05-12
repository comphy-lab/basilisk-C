/**
# Local constants can have different dimensions */

double func (double b) {
  double a = 1.; // this local constant has the dimension of b
  return a + b;
}

int main()
{
  {
    interpreter_verbosity (4);
    func (1. [1]);
    func (1. [2]);
  }
}
