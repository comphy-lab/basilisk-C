#include "grid/quadtree.h"
#include "input.h"
#include "output.h"

#define LEVEL 3

scalar D11[], D12[], D21[];

int main (int argc, char * argv[])
{
  size (1.);
  N = 1 << LEVEL;
  init_grid (1 << LEVEL);
    
    foreach()
    {
       D11[] = 1.*x;
       D12[] = 1.*y;       
       D21[] = 2.*x;
//       D22[] = 2.*y;
    }
    
    // Last item D22 is a constant scalar
    const scalar D22[] = 0.5; 
    char outfile[200];
    
    scalar * LLIST[] = {{D11,D12},
                        {D21,D22}};

    for(int row=0;row<2;row++)
    {
      int col=0;
      for((const) scalar s in LLIST[row])   // (const) mandatory?
      {
        sprintf(outfile,"D%d%d.txt",row+1,col+1);
        scalar tmp[];
        foreach()
          tmp[] = s[];
        FILE * fp = fopen(outfile,"wt");
        output_field({tmp},fp,n = N,linear = true);
        // output_field({s},fp);   // => fails  
        fclose(fp);
        col++;
      }
    }
    
    return(0);
}

/**
![D22.txt](test_lists/D22.txt)
*/