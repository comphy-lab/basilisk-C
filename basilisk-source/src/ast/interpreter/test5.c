/**
(Probably) checks for potential scope troubles */

int b = 1;

void func (int a)
{
  a = b;
}

int main()
{
  {
    interpreter_verbosity (4);
    struct {
      double x;
    } b;
    func (1);
  }
}
