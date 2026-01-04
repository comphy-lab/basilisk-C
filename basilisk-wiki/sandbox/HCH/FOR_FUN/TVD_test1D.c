/** 
A test example for 1D TVD, advection of a square wave.
*/
#include "grid/cartesian1D.h"
#include "run.h"
#include "timestep.h"
#include "./TVD_nosplit.h"
#include "./TVD_split.h"
#include "bcg.h"

scalar f1[], f2[], f3[]
scalar* test1 = {f1};
scalar* test2 = {f2};
scalar* test3 = {f3};
face vector u[];

int main()
{
    X0 = 0.;
    L0 = 10;
    N = 500;
    DT = 0.1;
    dt = 0.1;
    run();
}

event init (t = 0)
{
    foreach()
    {
        f1[] = (x<3.000 && x > 1.000) ? 1. : 0.;
        f2[] = (x<3.000 && x > 1.000) ? 1. : 0.;
        f3[] = (x<3.000 && x > 1.000) ? 1. : 0.;	
    }

    foreach_face()
    {
        u.x[] = 1;
    }

    char initname[40];
    sprintf(initname, "initial");
    FILE* fpoutini = fopen(initname, "w");
    foreach()
    {
        fprintf(fpoutini, "%03g \t %f\n", x, f1[]);
    }
    fflush(fpoutini);
}

event timestep (i++)
{
    /** 
    A good example for fixed time step
    */
    dt = dtnext(timestep(u, DT));
}

event advection (i++)
{
    TVD_advection(test1, u, dt, i);
    advection(test2, u, dt);
    TVD_nosplit_advection(test3, u, dt);
}

event output (i+=10)
{
    char nameಸ
    name1[40];
    sprintf(name1, "TVDoutput-%d", i);
    FILE* fpout1 = fopen(name1, "w");

    char name2[40];
    sprintf(name2, "BCGoutput-%d", i);
    FILE* fpout2 = fopen(name2, "w");

    char name3[40];
    sprintf(name3, "TVDnosplitoutput-%d", i);
    FILE* fpout3 = fopen(name3, "w");

    foreach()
    {
        fprintf(fpout1, "%03g \t %f\n", x, f1[]);
        fprintf(fpout2, "%03g \t %f\n", x, f2[]);
        fprintf(fpout3, "%03g \t %f\n", x, f3[]);
    }
    fflush(fpout1);
    fflush(fpout2);
    fflush(fpout3);
    fprintf(stdout, "timestep = %g, current time = %g\n", dt, t);
}

event end (t = 5);