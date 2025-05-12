/** 
The following code is a simple script illustrating the observed behaviour. 
We start by initializing a small 4 cell domain with four layers
*/
#include "grid/multigrid1D.h"
#include "layered/hydro.h"

int main() {
    N = 4;
    nl = 4;
    L0 = 4;
    run();
}

/**
In the init event, we construct a scalar field s and fill it with zeroes. Then we print the scalar field illustrating that all cell values are zero. We then change the value of a single cell from zero to one. Now, when we print the scalar field, a total of nl cell values have been changed to one! Finally we illustrate the computed sum of all cell values to be nl rather than one. 
*/
event init(i=0){

  /**
  The supposed bug was due to a misunderstanding of the way multilayer
  fields are allocated. See the [documentation](/Basilisk%20C#layers). */
  
    // Initializing empty scalar field
#if 0
    scalar s[]; // this is not a multilayer field allocation...
#else
    scalar s = new scalar[nl]; // this is a multilayer field allocation
#endif
    
    foreach(){
        foreach_layer(){
            s[] = 0;
        }
    }

    FILE * fp = fopen ("out.txt", "w");
  
    // Printing value of scalar field in each cell
    fprintf(fp, "Initial field:\n");
    foreach(){
        foreach_layer(){
            fprintf(fp, "s[] = %g\n", s[]);
        }
    }
    
    // Changing value of s in a single cell. 
    int inject_id = 2;
    int k = 0;
    foreach(){
        foreach_layer(){
            if (k == inject_id){
              s[] = 1;  // Same behaviour observed when using s[]++;
            }
            k++;
        }
    }

    // Printing changed scalar field 
    fprintf(fp, "\nChanged field:\n");
    int expected_sum = 1;
    int computed_sum = 0;
    foreach(){
        foreach_layer(){
            fprintf(fp, "s[] = %g\n", s[]);
            computed_sum += s[];
        }
    }
    fprintf(fp, "expected: %d    computed: %d\n", expected_sum, computed_sum);
    fclose (fp);
}

/**
The output is not what is expected:

![](scalar_field_multilayer_potential_bug/out.txt)
*/
