/** 
Two-phase Navier-Stokes solver is included to do the post-processing
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

int main()
{
	run();
}

event init(i=0)
{
	char filename[100];
	double tbg =0;
	double tend =20.00;
	double t;
	double tstep = 0.01;
/**
this for loop reads the dumpfiles and extract the information out.
*/
	for(t=tbg;t<=tend;t+=tstep)
	{ 
		sprintf(filename,"dumpfile_%.2f",t);
		restore(file = filename);
		vector h[];
		heights(f,h);
		double xMax = -HUGE;
		double y1 = 0;
		foreach()
		{
		    if(h.x[]!=nodata)
		    {
			if(x+Delta*height(h.x[])>xMax)
			{
			    xMax = x + Delta*height(h.x[]);
			    y1 = y;
			}
		    }
		}
		fprintf(stderr,"%f %f %f\n",t,xMax,y1);
    }
}
event end(i=0)
{

}
