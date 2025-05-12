#if _MPI
#include "../front-common.h"
#include "fmpi.h"

void fmpi_remap_frontpoints_fn1 (Array * a, Array * b){
	if (a->len == 0) return;

	Frontpoint * frontpoints = get_front()->frontpoints, * frontpoint;
	size_t size = sizeof(int);
	int * ptid = (int *) (a->p);
	for(int p=0; p< (int) (a->len/size); ++p, ++ptid){
		frontpoint = frontpoints + (*ptid);
		array_append(b, &(frontpoint->pid), size);
		array_append(b, &(frontpoint->alias), size);
	}
}

void fmpi_remap_frontpoints_fn2 (Array * a, Array * b){
	array_append (b, a->p, a->len);
}

void fmpi_remap_frontpoints(Front *fr, Fmpi * f){

	reset_fmpi_arrays(f);			

	int map;
	size_t ssize = sizeof(int), rsize = sizeof(gID); 
	bool len = 1;
	foreach_frontelement(fr){
		gID * g = frontelement->p;
		for(int i=0; i<dimension; ++i, ++g) {	
			if( (g->pid) != pid() ) {
				map = f->mappid[g->pid]; //fixme: not fast
				assert(map != -1);
assert(f->sr[map].pid == g->pid);
				array_append(f->sr[map].q,    &(g->alias), ssize);
				array_append(f->sr[map].temp, &elid, ssize);
				array_append(f->sr[map].temp, &i, ssize);
				len = 0;
			}
			else{
				g->pid 	 = frontpoints[g->alias].pid;
				g->alias = frontpoints[g->alias].alias;
				frontelement->v[i] =  (g->pid != pid()) ? -1 : g->alias;
			}
		}
	}

	fmpi(f, fmpi_remap_frontpoints_fn1,
	        fmpi_remap_frontpoints_fn2);

	if(len) return;
/*
	foreach_frontelement(fr){
		gID * g = frontelement->p;
		for(int i=0; i<dimension; ++i, ++g) {	
			if( (g->pid) == pid() ) {
				g->pid   = frontpoints[g->alias].pid;
				g->alias = frontpoints[g->alias].alias;
				frontelement->v[i] =  (g->pid != pid()) ? -1 : g->alias;
			}
		}	
	}
*/

	Frontelement * frontelements = fr->frontelements, * frontelement;
	foreach_fmpi_sr(f){
		int * elid  = (int *) sr->temp->p;
		gID * g  = (gID *) sr->r->p, * gr;
assert((ssize*sr->r->len) == (rsize*sr->q->len));	
		for(int i=0; i< (int) (sr->temp->len/(2*ssize)); ++i, elid += 2, ++g)	{
			frontelement = frontelements + elid[0];
			gr = frontelement->p + (elid[1]); 
assert( (gr->pid) == sr->pid );
			memcpy(gr, g, rsize);
			frontelement->v[elid[1]] =  (g->pid != pid()) ? -1 : g->alias;
		}
	}

}

#endif