/** 
This script is used to do post-processing, it can output the coordinate of the air sheet's tip at evering saving time step*/
#include "navier-stokes/centered.h"
#define FILETERED 1
#include "two-phase.h"
#include "tension.h"

int main()
{
	run();
}

event init(i=0)
{
	char filename[100];
	double tbg = 0.00;
	double tend = 100.00;
	double t;
	double tstep = 0.01;
	for(t=tbg;t<=tend;t+=tstep)
	{ 
		sprintf(filename,"dumpfile_%.2f",t);
		restore(file = filename);
		vector h[];
		heights(f,h);
		double yMin = HUGE;
		double x1 = 0;
		foreach()
		{
			if(h.y[]!=nodata)
			{
				if(y+Delta*height(h.y[])<yMin)
				{
					yMin = y + Delta*height(h.y[]);
					x1 = x;
				}
			}
		}
		double yMax = -HUGE;
		double x2 = 0;
		foreach()
		{
			if(h.y[]!=nodata)
			{
				if(y+Delta*height(h.y[])>yMax)
				{
					yMax = y + Delta*height(h.y[]);
				    x2 = x;
				}
			}
		}
		double xMax = -HUGE;
		double y3 = 0;
		foreach()
		{
		    if(h.x[]!=nodata)
		    {
			if(x+Delta*height(h.x[])>xMax)
			{
			    xMax = x + Delta*height(h.x[]);
			    y3 = y;
			}
		    }
		}
		fprintf(stderr,"%f %f %f %f %f %f %f\n",t,xMax,y3,x1,yMin,x2,yMax);		/** By using the coordinates of these characteristic points, you can obtain the corresponding thickness, Reynolds number, local Weber number, local Ohnesorge number and aspect ratio.*/
    }
}
event end(i=0)
{
	

}
