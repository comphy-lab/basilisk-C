/**
Check the dimensions of implicit zeros in structures. */ 

struct Func {
  double a; // explicit (1)
  double b; // implicit (0)
};

int main()
{
  { interpreter_verbosity (2);
    struct Func p = {1};
    display_value (p.a);
    display_value (p.b);
  }
}
