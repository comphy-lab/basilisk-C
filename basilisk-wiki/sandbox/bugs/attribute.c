/**
qcc fails for the code below, because the type system is too
rudimentary i.e. it does not recognize that f.rho is indeed a scalar. */

attribute {
  scalar rho;
}

scalar f[];

event init (i = 0) {
  f.rho = new scalar;
}

/**
Changing the code to

~~~literatec
event init (i = 0) {
  scalar a = new scalar;
  f.rho  = a;
}
~~~

works. */
