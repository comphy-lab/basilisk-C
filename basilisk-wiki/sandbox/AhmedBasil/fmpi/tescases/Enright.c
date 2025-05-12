/**
LeVeque
https://epubs.siam.org/doi/10.1137/0733033

*/

#define Tend 1.
#define tout (Tend/20.)

#include "grid/octree.h"
#include "run.h"

scalar f[];

int main() 	{
#if dimension == 2
	fprintf(stdout, "\nCode for 3d only. Exiting");
	return 0;
#endif
	size(1.0);
	//origin(-1.0,-1.0,-1.0);
  init_grid (1 << 5);
	
	//dump("box1x1x1.dmp");
	
	run();		
	return 0;
}

#include "../../front-common.h"	
#include "../../front-out.h"	
//#include "../fmpi-rebalance.h"
#include "../fmpi-redistribute.h"

//#include "front-quality.h"

event init (i = 0) {
	Front * fr = get_front();

	add_sphere(fr, 0.5, 0.75, 0.75, 0.15);
#if _FDEBUG	
	fmpi_debug(fr);
#endif
	dump("dumpdirectory/Enright_start.dmp");
	dump_front(file = "dumpdirectory/start.dmpf", fmpi_ = true);
}

/**
Advancing front
*/
event stability(i++)
{
	fprintf(stdout,"\nt %g", t);
	if (dt>0.0001)
		dt = 0.0001;
}
//event advance(i++, t += dt) {
event advance(i++; i<100) {
//event advance(i++; i<1) {
	Front *fr = fr_sphere;
	double xm, ym, zm, um, vm, wm, T = Tend;

	foreach_frontpoint(fr) {
		xm = frontpoint->x[0];
		ym = frontpoint->x[1];
		zm = frontpoint->x[2];
		Point point = locate(xm, ym, zm);	
		assert(level >= 0);

		//Enright
		um = 2*sin(pi*xm)*sin(pi*xm)*sin(2*pi*ym)*sin(2*pi*zm)*cos(pi*t/T); 
		vm = -sin(2*pi*xm)*sin(pi*ym)*sin(pi*ym)*sin(2*pi*zm)*cos(pi*t/T); 
		wm = -sin(2*pi*xm)*sin(2*pi*ym)*sin(pi*zm)*sin(pi*zm)*cos(pi*t/T); 
double dmin = L0/(1<<grid->maxdepth);
assert(um*dt < dmin);
assert(vm*dt < dmin);
assert(zm*dt < dmin);
		frontpoint->x[0] +=  um*dt;	
		frontpoint->x[1] +=  vm*dt;	
		frontpoint->x[2] +=  wm*dt;	
	}

	redistribute(fr);
	char file[100]; sprintf(file, "dumpdirectory/Enright-%d.dmp", i); dump_front(file = file, fmpi_ = true);
}

/*
event regrid(i+=10){
event regrid(i++){
	Front *fr = get_front();
	regrid_front(fr);
}
*/
/*
event output(t+=tout) {
	char f2[40];
	sprintf(f2, "dumpdirectory/Enright_dumpf_t%0.2f", t);
	dump_front(f2);
}

event laststep(t=Tend) {
}

event end(t=end) {
}
*/
