/**
# Variable scopes */

double c;

double func (double x)
{
  display_value (x);
  display_value (c);
}

int main()
{ // 2
  func (c);
  { // 3
    double a;
    display_value (a);
    func (a);
    { // 4
      double b;
      display_value (b);
      display_value (c);
      func (c);
    }
    display_value (a);
  }
}
