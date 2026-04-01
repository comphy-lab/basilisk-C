#include <stdio.h>

/**
# Define an enum with name `Foo`
We directly define the enum type `Foo` using typedef similarly to the way structs are usually defined in C.
This is valid C99 code and is accepted without warning by gcc, clang, and the intel C compiler.
Enums are defined this way for example in the HDF5 library.
Using the HDF5 library will therefore cause the warnings to appear.
*/
typedef enum { FOO, BAR, BAZ } Foo;

/**
# Use the enum
This will produce a warning from qcc.
The code for this minimal example still does the correct thing.
*/
static void show(Foo f) {
  switch (f) {
    case FOO: puts("FOO"); break;
    case BAR: puts("BAR"); break;
    case BAZ: puts("BAZ"); break;
  }
}

int main(void) {
  show(FOO);
  show(BAR);
  show(BAZ);
}

/**
Note: Defining the enum like this will not produce a warning:
```C
enum Foo { FOO, BAR, BAZ };
typedef enum Foo Foo;
```
*/
