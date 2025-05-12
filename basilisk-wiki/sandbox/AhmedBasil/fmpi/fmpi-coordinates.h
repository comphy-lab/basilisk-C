#if _MPI
#include "../front-common.h"
#include "fmpi.h"

void fmpi_coordinates_fn1 (Array * a, Array * b){
	if (a->len == 0) return;
	Frontpoint * frontpoints = get_front()->frontpoints, * frontpoint;
	size_t size = dimension*sizeof(double); 
	int * ptid = (int *) (a->p); 
	for(int p=0; p<(int) (a->len/sizeof(int)); ++p){
		frontpoint = frontpoints + ptid[p];
		assert(local_point());
		array_append(b, frontpoint->x, size);
	}
}

void fmpi_coordinates_fn2 (Array * a, Array * b){
	array_append (b, a->p, a->len);
}

void fmpi_coordinates(Front *fr, Fmpi * f){

	reset_fmpi_arrays(f);			

	int pmap, len = 1;
	size_t ssize = sizeof(int), rsize = dimension*sizeof(double); 
	foreach_frontpoint_all(fr){
		if(local_point() == 0) {
			pmap = f->mappid[frontpoint->pid]; //fixme : slow
			assert(pmap != -1);
			array_append( f->sr[pmap].q, &(frontpoint->alias), ssize);
			array_append( f->sr[pmap].temp, &ptid, ssize);
			len = 0;
		}
	}

	fmpi(f, fmpi_coordinates_fn1, fmpi_coordinates_fn2);

	if(len) return;
	
	Frontpoint * frontpoints = fr->frontpoints;
	foreach_fmpi_sr(f){
		int * ptid = (int *) sr->temp->p;
		double * x = (double *) sr->r->p;
assert((ssize*sr->r->len) == (rsize*sr->q->len));	
		for(int i=0; i< (int) (sr->temp->len/ssize); ++i, x += dimension)	
			memcpy(frontpoints[ptid[i]].x, x, rsize);
	}

}

#endif