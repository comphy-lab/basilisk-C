/**
# Test writing to arrays with an "unset" index 

This does not work yet. See the UNSET_ARRAY flag in [interpreter.c](). */

int main()
{
  { interpreter_verbosity (4);
    int i;

    double b[] = {4,5};
    b[i]; // fixme: must be "unset"
    b[0]; // must not be "unset"
    b[i] = 6; // must not be (unset)
    b[0]; // must be "unset"

    struct { int x, y; } a[] = {{1,2}};

    a[i] = (struct { int x, y; }){3,4};
    a[0].x; // fixme: must be (unset)
    a[0].y; // fixme: must be (unset)

    a[0] = (struct { int x, y; }){5,6};
    a[0].x; // must not be (unset)
    a[0].y; // must not be (unset)

    a[i].x; // fixme: should be (unset)
    
    a[i].x = 7; // must not be (unset)
    a[0].x; // fixme: should be (unset)
    a[0].y; // must not be (unset)
  }
}
