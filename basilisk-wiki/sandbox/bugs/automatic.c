/**
# Threads may get stuck

The following example illustrates a scenerio where not all MPI treads
try to set (automatic) boundary conditions. It may wait forever or
crashes occur when other threads continue to set boundary conditions
in other contexts.

<div class="message">
Stéphane Popinet replies:

As soon as you "branch" execution depending on process number, i.e. here

~~~c
if (pid() == 0)
  do_something();
~~~

you need to be sure that "do_something();" is completely independent from what is happening in other processes (otherwise this thread will get stuck forever waiting for other processes).

This is not specific to Basilisk.</div>
*/

int main() {
  init_grid(N);
  scalar s[];
  foreach()
    s[] = x;
  if (pid() == 0)
    interpolate (s,  0.1, 0.1);
  printf ("%d\n", pid());
}