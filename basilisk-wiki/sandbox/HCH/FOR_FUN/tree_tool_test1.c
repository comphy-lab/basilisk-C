/**
This file is to test mesh tools including is_active, is_local, is refined, is_leaf and their combination occured in tree headerfile.
 */
#include "grid/quadtree.h"
#include "utils.h"

#define WIDTH 16
void meshdraw();
void foreachdraw();

int main()
{
        size(WIDTH);
        origin (-8,-8);
        init_grid (4);

/**
Three layers of mesh is created whose centroid located at (0,0).
 */
        refine(x>=0&&y>=0&&level<3);
        refine(x>=3&&y>=3&&level<4);

        meshdraw();
        foreachdraw();
        
        return 0;
}

void meshdraw()
{
        char name[20];
        sprintf(name, "meshout");
        FILE*meshfp = fopen(name, "w");
        output_cells(meshfp);
        fflush(meshfp);
}

void foreachdraw()
{
        char name[20];
        sprintf(name, "foreach_cell");
        FILE*fp1 = fopen(name, "w");

        sprintf(name, "is_active");
        FILE*fp2 = fopen(name, "w");

        sprintf(name, "is_leaf");
        FILE*fp3 = fopen(name, "w");

        sprintf(name, "is_local");
        FILE*fp4 = fopen(name, "w");

        sprintf(name, "is_refined");
        FILE*fp5 = fopen(name, "w");

        sprintf(name, "leaf_child");
        FILE*fp6 = fopen(name, "w");

        sprintf(name, "leaf_local");
        FILE*fp7 = fopen(name, "w");

        int i;

        foreach_cell()
        {
                i = 1;
                if(true)
                {
                        fprintf(fp1,"%d\t%f\t%f\n",i,x,y);
                        i++;
                }

                i = 1;
                if(is_active(cell))
                {
                        fprintf(fp2,"%d\t%f\t%f\n",i,x,y);
                        i++;
                }

                i = 1;
                if(is_leaf(cell))
                {
                        fprintf(fp3,"%d\t%f\t%f\n",i,x,y);
                        foreach_child()
                        {
                                fprintf(fp6,"%d\t%f\t%f\n",i,x,y);
                        }
                        i++;
                }

                i = 1;
                if(is_local(cell))
                {
                        fprintf(fp4,"%d\t%f\t%f\n",i,x,y);
                        i++;
                }

                i = 1;
                if(is_refined(cell))
                {
                        fprintf(fp5,"%d\t%f\t%f\n",i,x,y);
                        i++;
                }
        }

/**
I especially test the following case which frequently shows up in tree concerned file. 
 */
        foreach_cell()
        {
                i = 1;
                if(is_leaf(cell))
                {
                        continue;
                }
                if(is_local(cell))
                {
                        fprintf(fp7,"%d\t%f\t%f\n",i,x,y);
                        i++;
                }
        }

        fflush(fp1);
        fflush(fp2);
        fflush(fp3);
        fflush(fp4);
        fflush(fp5);
        fflush(fp6);
        fflush(fp7);
}