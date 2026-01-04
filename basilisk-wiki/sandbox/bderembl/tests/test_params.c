/**
# Test for parameter reading in extra.h

Here is the typical skeleton file that I use when I need namelists.
   
The 3 key steps in main are (i) param declaration, (ii) adding parameters,
(iii) cleanup.

~~~literatec
params = array_new();
add_param ("N", &N, "int");
array_free (params);
~~~

*/
#include "../libs/extra.h"

/**
   Declare two global variable, nu and dh, which will be modified by the
   namelist.
*/
double nu = 0.;
double dh[2] = {0., 1.};

int main(int argc,char* argv[]) {

/**
   The type of parameters that we can read in a namelist are "int", "double", and "array"
*/
  params = array_new();
  add_param ("N", &N, "int");
  add_param ("nu", &nu, "double");
  add_param ("dh", &dh[0], "array");


  // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param,argv[1]); // default: params.in

  read_params(file_param);

/**
   Create an output directory named outdir_0001 or use a higher number if 0001
   already exists
*/

  create_outdir();

/**
   copy the configuration file in that directory
*/
  backup_config(file_param);

/**
   Free the params variable
*/
  array_free (params);
}
