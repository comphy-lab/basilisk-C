/**
# Careful with dimensions depending on undefined branches */

int main()
{{ interpreter_verbosity (3);

    double a = 1, b = 3 [1];
    if (x)     // This is an undefined branch
      a *= b;  // [a] = [1] + [3]
}}
    
/**
The dimensions must be the same whether the branch has been taken
or not, this implies that:
[1] (the branch has not been taken)
must be equal to:
[1] + [3] (the branch has been taken)
and thus:
[3] == [0]
which is not possible. 

The solution is simple:

~~~c
double a, b = 3 [1];
if (x)
  a = b;
else
  a = 1;
~~~
*/
