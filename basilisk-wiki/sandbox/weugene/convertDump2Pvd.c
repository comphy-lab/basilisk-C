/** 
# Aim
As it is known Basilisk can save data in dump format and in pvtu-vtu or pvd-pvtu-vtu files. Saving in dump files (binary foramt) is more reasonable due to its relatively small size in comparison with Paraview format(mixed binary and ASCII format). Primarily thanks to its pure binary format and saving all the data which can be subsequently postprocessed.

##In short about pvd file.
*.pvd file is an extension of Øystein Lande's module "output_vtu_foreach.h" which can be found in [his sandbox](http://basilisk.fr/sandbox/oystelan/), where a bunch of pvtu (header file, outputted by master processor) and vtu (outputted by each processor) files are outputted.
In my version, I am saving time of outputing in one pvd file (header files refers to all outputted pvtu files and save corresponding time).

    test.pvd
    res/test_0_0000.pvtu
    res/test_0000_n000.vtu
    res/test_0000_n001.vtu
    res/test_0_0001.pvtu
    res/test_0001_n000.vtu
    res/test_0001_n001.vtu


You can find my version of modification in my sandbox [here](http://basilisk.fr/sandbox/weugene/)


The present code helps to 

    1) Postprocess computations vorticity, lambda2, etc.,
    2) Cut not interesting regions, reducing the size of Paraview files,
    3) Output to the pvd file (Detailed information can be found in my sandbox). 

This module can be run in parallel in order to distribute a single dump file into vtu files between processors. It can accelerate calculations if you are using  paraview scripting tools, such as pvbatch: 

    mpirun -np 4 pvbatch ./my_script.py
*/
#include "grid/octree.h"
#include "run.h"
#include "timestep.h"
#include "utils.h"
#include "lambda2.h"
#include "../src_local/output_vtu_foreach.h"

#include <stdio.h>
#include <stdlib.h>
#include <wordexp.h>
#include <ctype.h>
/**
List all interesting fields.
*/
scalar fs[], omega[], l2[], l[];
scalar f[], * interfaces = {f};
vector u[];
/**
Parameters for reducing refinement in not interesting regions using unrefine function.
*/
double length_min = 1e+30, length_max = -1e+30, length = 1;
/**
Convert each i_take-th dump files
*/
int i_take = 1;
double myt=0;
/**
List of all dump files
*/
wordexp_t fp;
char **w, dump_name[30];

//get time from dump name (dump-1.1234 => 1.1234)
double get_double(const char *str)
{
    /* First skip non-digit characters */
    /* Special case to handle negative numbers and the `+` sign */
    while (*str && !(isdigit(*str) || ((*str == '-' || *str == '+') && isdigit(*(str + 1)))))
        str++;
    /* The parse to a double */
    return strtod(str, NULL);
}

int main (int argc, char * argv[]) {
    // set which dump files will be converted: each $(i_take)th
    // by default each dump will be converted
    if (argc > 1)
        i_take = atoi (argv[1]);

    wordexp("dump-*", &fp, 0);
    w = fp.we_wordv;
    for (int i = 0; i < fp.we_wordc; i++) fprintf(ferr, "All dump files: %s\n", w[i]);

    for (int i = 0; i < fp.we_wordc; i += i_take) {
        myt =  fabs(get_double(w[i]));
        strcpy(dump_name, w[i]);
        fprintf(ferr, "reading dump file: %s at time t= %g\n", dump_name, myt);
        run();
    }
    wordfree(&fp);
}

event init (t = 0) {
    bool success = restore (file = dump_name);
    fprintf(ferr, "file has been read: L0=%g\n", L0);

    if (!success) {
        fprintf(ferr, "can't open the file %s. Missing this file, go to the next file\n", dump_name);
        return 0;
    }
}

event calculate_aux_fields(i++){
    foreach() l[] = level;
    vorticity (u, omega);
    lambda2 (u, l2);
}
/**
Here region far away from a single bubble will be unrefined. Here we find out the center of the bubble, then shift by -5 and 4 to find out a high resolution box. And specify a high resolution cylinder parameters.
*/
event coarsen_grid(i++){
    double xcg = 0, dvtmp, volume = 0, volumeg = 0 ;
    length_min = 1e+30, length_max = -1e+30, length = 0;
    foreach( reduction(+:xcg) reduction(+:volume) reduction(+:volumeg)) {
        if (fs[]<1){
            dvtmp = (1.0 - f[])*(1.0 - fs[])*dv(); // gas volume
            volumeg += dvtmp;//gas liquid
            volume += (1.0 - fs[])*dv();//channel volume
            xcg   += x*dvtmp;// Along x
        }
    }
    xcg /= volumeg;
    length_min = xcg - 5;
    length_max = xcg + 4;
    length = length_max - length_min;

    fprintf (ferr, "x= %g length_min= %g length_max= %g length= %g it_fp= %d\n"
                   "volume= %g volumeg= %g\n",
             xcg, length_min, length_max, length, iter_fp,
             volume, volumeg);
    unrefine ( (x < length_min || x > length_max) && level >= 1);
    unrefine ( (sq(y) + sq(z) > sq(0.55)) && level >= 1);
}



event vtk_file (i++)
{
  /** 
  There is a two-line calling of the function for output vtu files where myt is time of outputting, e.g. myt = t + dt
  */
    char subname[80]; sprintf(subname, "dump2pvd_compressed");
    output_vtu_MPI( subname, myt, (scalar *) {fs, f, l, l2, omega}, (vector *) {u});
}

event stop(t = 100){ // t = 100 should  be sufficiently big in order to reach this event
    return 0;
};


