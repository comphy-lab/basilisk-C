/**
# How to use `output_vtu.h` with MPI (Brun.sh)

Here, I show you how to extract your data when using MPI.  
We start by creating an `extract_file.c` and including the necessary libraries:
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "output_vtu.h" 

/**
In this example, our `main_file.c` uses:

- `#include "navier-stokes/centered.h"`
- `#include "two-phase.h"`

More generally, you should include in `extract_file.c` the same headers you used in your `main_file.c`.

Then, simply copy the code below. Here, we create a folder named `vtufile` to store the output:
*/

int pt = 0;
int main(int argc, char *argv[]){	
	if (argc > 1) pt = atoi(argv[1]);
	system("mkdir -p vtufile");
	run();	
}

/**
We now use `output_vtu.h`, and specifically the function `output_vtu_bin_foreach`, to save our data (e.g., the volume fraction $f$, the velocity field $\mathbf{u}$ and the pressure field $p$).
*/

event init (t = 0.0){
	char filename[80];
	while(1){
		sprintf(filename, "snap-%d",pt);
		printf("extracting %s\n",filename);
		fflush(stdout);
		if(!restore (file = filename)){
			printf("files ended");
			exit(0);
		}
	
		scalar *list;		
		bool linear;
		char name[80];
		sprintf(name, "vtufile/snap-%d.vtu",pt);
		FILE *ptr = fopen(name,"w");
		output_vtu_bin_foreach((scalar *) {f, p}, (vector *) {u}, ptr, false);
		fclose(ptr);
	
		pt = pt + 1;
	}
}

/**
This is how you can extract your simulation data when using `Brun.sh`. 

To summerize:

- in your `main_file.c` remove `#include output_vtu.h`
- in your `extract_file.c` include all the same headers as in your `main file.c`, plus `#include output_vtu.h`
- In the `output_vtu_bin_foreach` call, choose the variables you want to save: volume fraction, velocity field, pressure, etc.

Special thanks to Yashika and Anil for sharing this method!
*/