/**
Run this using something like

~~~bash
make restore.tst
cd restore
./restore
~~~

should give

~~~
0 0
1 1
2 2
3 3
4 4
5 5
6 6
~~~

then

~~~bash
cp dump restart
./restore
~~~

should give

~~~
2 2
3 3
4 4
5 5
6 6
~~~

but this does this instead

~~~
3 6 1
~~~

*/

#include "layered/hydro.h"

int main()
{
  DT = 1;
  run();
}

event init (i = 0) {
  restore ("restart-restore.dump");
}

event dumping (t = 2)
  dump();

event ending (t = 6);

event logfile (i++)
  fprintf (stderr, "%d %g %g\n", i, t, dt);
