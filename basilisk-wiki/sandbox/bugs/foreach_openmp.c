/**
MWE showing unintended behaviour of the qcc compiler to identify race conditions
when OpenMP is enabled.

Code should **fail** compilation.

Compilation using

```qcc test.c -fopenmp -lm``` 

succeeds.
*/

#include "grid/quadtree.h"
#include "run.h"

scalar f1[], f2[];
scalar *l1 = {f1};
scalar *l2 = {f2};
/**
Dummy case setup. Not necessary since we are not *running* anything.
*/
int main ()
{
	//domain size
	init_grid( 1<<7 );

	run();
}

/**
```scalar``` ```curVar, oldVar``` are defined outside the foreach loop. Code should fail to compile
with -fopenmp because there is a race condition.
Same behaviour is seen when ```scalar``` is replaced by ```double```/```int```.
See [Google Groups Discussion](https://groups.google.com/g/basilisk-fr/c/dZ4wIC0DwsU) for more details.
*/
event update(i< 2; i++)
{
	scalar curVar, oldVar;
	foreach()
	{
		for (curVar, oldVar in l1, l2)
		{
			oldVar[] = curVar[];
		}
	}
	//fields_stats();
}

event updateDouble(i< 2; i++)
{
	double var;
	foreach()
	{
		var = sq(f1[]);
		f1[] = var;
	}
	//fields_stats();
}

event updateInt(i< 2; i++)
{
	int var;
	foreach()
	{
		var = sq(f1[]);
		f1[] = var;
	}
	//fields_stats();
}
