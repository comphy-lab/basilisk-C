/**
# Loop tests for the interpreter */

void func (int i) {}

int main()
{
{
  interpreter_verbosity (4);
  
  int i = 0;
  while (i < 3) {
    func (i);
    i++;
  }

  do {
    i--;
    func (i);
  } while (i > 0);
  
  for (int j = 0; j < 3; j++)
    func (j);

  /**
  And some structure tests which have nothing to do with loops. */
  
  func (sizeof (scalar));

  scalar s;
  s = (scalar){1};
  func (s.i);
 }
}
