#include "grid/cartesian.h"
#include "run.h"
#include "timestep.h"
#include "bcg.h"
#include "./TVD_split.h"
#include "./TVD_nosplit.h"

scalar f1[], f2[], f3[];
face vector uf[];

/**
Initial condition: Gaussian blob centered at (0.5, 0.5), nosplit TVD fails in this test
*/
double initial(double x, double y) {
    return exp(-200. * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
}

int main() {
    origin(-0.5, -0.5);
    L0 = 1.;
    N = 256;
    DT = 0.5;
    periodic(right);
    periodic(top);

    FILE *output_file = fopen("moving_output", "w");
    if (!output_file) {
        fprintf(stderr, "Error opening output file\n");
        return 1;
    }

    run();
    fclose(output_file);
}

event init(i = 0) {
    /**
    Face velocity defined separately for convenience in modification
    */
    foreach_face(x) uf.x[] = 1.;
    foreach_face(y) uf.y[] = 1.;

    foreach() {
        f1[] = initial(x, y);
        f2[] = initial(x, y);
        f3[] = initial(x, y);
    }
}

/**
Flag `last` arranges event at end of sequence rather than calling order
*/
event timestep(i++, last) {
    dt = dtnext(timestep(uf, DT));
}

event advection(i++, last) {
    scalar *tracer1 = {f1};
    scalar *tracer2 = {f2};
    scalar *tracer3 = {f3};
    TVD_advection(tracer1, uf, dt, i);
    advection(tracer2, uf, dt);
    TVD_nosplit_advection(tracer3, uf, dt);
}

event dtprint(i += 10, last) {
    fprintf(stderr, "i = %d \t dt = %f\n", i, dt);
}

/**
Output error metrics
*/
event output_event(t = {0, 0.5, 1}) {
    double maxerr1 = 0., maxerr2 = 0., maxerr3 = 0.;
    foreach() {
        double err1 = fabs(f1[] - initial(x - t, y));
        if (err1 > maxerr1) maxerr1 = err1;

        double err2 = fabs(f2[] - initial(x - t, y));
        if (err2 > maxerr2) maxerr2 = err2;

        double err3 = fabs(f3[] - initial(x - t, y));
        if (err3 > maxerr3) maxerr3 = err3;
    }

    FILE *output_file = fopen("moving_output", "a");
    if (output_file) {
        fprintf(output_file, "t=%g TVD maxerr1=%g BCG maxerr2=%g TVD_nosplit maxerr3=%g\n", t, maxerr1, maxerr2, maxerr3);
        fclose(output_file);
    } else {
        fprintf(stderr, "Error writing to output file\n");
    }
}

event stop(t = 1);