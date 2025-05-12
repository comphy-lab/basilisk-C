/** 
#WAVEPROBES

   This function uses the built in height function in basilisk to calculate the surface
   elevation at a given point x (or x and y in 3D). The probes lower and upper elevation must be specified in the probe_heightlimits array.
   Calculated heights from heights(f,h) is only available in cells close to the free surface, hence the wave probe is discretized using <ssize> 
   as the number of points for which the distance to the surface is calculated. the specifiec number of points are evenly distributed between 
   the lower and upper values specified in <probe_heightlimits>. Finally the function will return the value from the point closest to the free surface.

   For 3D, the sea surface is assumed to be in the xy-plane. i.e., the probes are assumed to be aligned with the z-axis
   
   Author: Oystein Lande
   Date: 2017-09-30
   Updated: 2018-09-05 corrected some unfortunate bugs and made it more userfriendly.
*/
#include "heights.h"
vector h[];

/** Start off by defining a helpful function*/
// Get index of minimum value
int getIndexMinProbe(double arr[], int size) {

   int i;
   int iii = 0;
   double sum = arr[0];

   for (i = 0; i < size; ++i) {
      if (arr[i] < sum){
        sum = arr[i];
        iii = i;  
        }
   }
   return iii;
}

/* Waveprobe function for 2D simulations*/
#if dimension == 2
double wprobe(double xcoord,double probe_heightlimits[],int ssize){
	double yi_dist[ssize];
    double yi[ssize];
    double ycoord;

	  //int NN = sizeof(yi_dist) / sizeof(double);
	for (int ccc=0;ccc<ssize;ccc++){
			ycoord = probe_heightlimits[0]+((probe_heightlimits[1]-probe_heightlimits[0])/(ssize-1.0))*ccc; 
            Point point = locate(xcoord,ycoord);
		if (h.y[] != nodata) {  
			yi[ccc] = y + height(h.y[])*Delta;

			yi_dist[ccc] = abs(y - yi[ccc]);
		}
		else{
			yi_dist[ccc] = 999;
      yi[ccc] = -999;
		}
		      
	}    
    return yi[getIndexMinProbe(yi_dist,ssize)];    
}
#endif
/** Waveprobe function for 3D simulations */
#if dimension == 3
double wprobe(double xcoord,double ycoord,double probe_heightlimits[],int ssize){
	double zi_dist[ssize];
    double zi[ssize];
    double zcoord;
	  //int NN = sizeof(yi_dist) / sizeof(double);
	for (int ccc=0;ccc<ssize;ccc++){
		zcoord = probe_heightlimits[0]+((probe_heightlimits[1]-probe_heightlimits[0])/(ssize-1.0))*ccc; 
        Point point = locate(xcoord,ycoord,zcoord); 
		if (h.y[] != nodata) {  
			zi[ccc] = z + height(h.z[])*Delta;
			zi_dist[ccc] = abs(z - zi[ccc]);
		}
		else{
			zi_dist[ccc] = 999;
      		zi[ccc] = -999;
		}
		      
	}    
    return zi[getIndexMinProbe(zi_dist,ssize)];    
}

#endif

/**
## Example of usage given below
// find wave elevation at position x = 0, and write it to file, every 0.05sec
event waveprobe (t+=0.05;t<=MAXTIME) {
        heights(
	static FILE * fp0 = fopen("waveprobe0.dat", "w");
	double ycoords[2]  = {-0.05,0.05}; // defines the minimum and maximum y-value 
	double yMax0 = wprobe(0.0,ycoords,5); // The last value in wprobe specifies how many discretization points which should be used in the given ycoords range. 	
	// update file
	fprintf(fp0, "%g %g\n",t,yMax0);
	fflush(fp0);
}
*/