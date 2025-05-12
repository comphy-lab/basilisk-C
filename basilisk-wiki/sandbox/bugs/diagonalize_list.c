/**
# qcc crashes when using the diagonalize macro with a list of scalars
 */

/**
The following macro:
*/
macro testdiag (scalar a, double expr)
{
    foreach() {
        double _d;
        diagonalize(a)
            _d = expr;
        a[] *= _d ;
    }
}
/**
can be called with success with, e.g.:
*/
/**
~~~c
testdiag (u, u[1] + u[-1] - 2.*u[]);
~~~
*/
/**
However, writing:
*/
macro testdiag_list (scalar * list, double * expr)
{
    foreach() {
        for (int i = 0; i < list_len (list); i++) {
            scalar a = list[i];
            double _d;
            diagonalize(a) // comment this line to avoid qcc segfault
                _d = expr[i];
            a[] *= _d ;
        }
    }
}
/**
results in `qcc` segfaulting when called with:
*/
/**
~~~c
scalar * list_scalars = {u};
testdiag_list (list_scalars, (double[]){u[1] + u[-1] - 2.*u[]});
~~~
*/

int main()
{
    init_grid (16);
    periodic (right);
    scalar u[];

    foreach()
        u[] = 1.;

    scalar * list_scalars = {u};

    testdiag (u, u[1] + u[-1] - 2.*u[]);
    testdiag_list (list_scalars, (double[]){u[1] + u[-1] - 2.*u[]});
}
