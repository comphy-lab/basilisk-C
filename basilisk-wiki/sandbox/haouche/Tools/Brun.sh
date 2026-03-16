# Here is a way to use MPI with Basilisk:

# Here you just have to change the name of your file "file_name.c"
# and change the value of "N" for example, if you want to 
#use 4 processor remove N and set up 4

qcc -source -grid=quadtree -D_MPI=1 file_name.c 
mpicc -std=c99 -D_XOPEN_SOURCE=700 _file_name.c -lm -O2 -Wall -o parallel mpiexec -np N ./parallel

# And in the terminal write this command to execute your code:
>> sh Brun.sh

# It's compatible with output_vtu.h, just instead to use this:
//int pt = 0;
//event interface(t += time_step; t <= final_time) {
//    scalar *list;
//    bool linear;
//    char name[80];
//    sprintf(name, "snap-%d.vtu",pt);
//    FILE *ptr = fopen(name,"w");
//    output_vtu_bin_foreach((scalar *) {f, p}, (vector *) {u}, ptr, false);
//    fclose(ptr);
//
//    pt = pt + 1;
//}

# Replace in your code by this one:
//int pt = 0;
//event interface(t += time_step; t <= final_time) {
//    char name[80];
//    scalar* list = {f, p, u};
//    sprintf (name, "snap-%d", pt);
//    p.nodump = false;
//    dump(name, list);
//
//    pt = pt + 1;
//}

# And don't forget to remove #include "output_vtu.h" 
# in your file_name.c

# Thanks to Yashika and Anil for sharing this code.