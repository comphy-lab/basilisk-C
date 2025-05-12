#if _MPI
#include "../front-common.h"
#include "fmpi.h"

void fmpi_remap_frontelements_fn1 (Array * a, Array * b){
	if (a->len == 0) return;

	Front * fr = get_front();
	Frontelement * frontelements = fr->frontelements, * frontelement;
	size_t size = sizeof(int);
	int * elid = (int *) (a->p);
	for(int p=0; p< (int) (a->len/size); ++p, ++elid){
		frontelement = frontelements + (*elid);
		array_append(b, &(frontelement->pid), size);
		array_append(b, &(frontelement->alias), size);
	}
}

void fmpi_remap_frontelements_fn2 (Array * a, Array * b){
	array_append (b, a->p, a->len);
}

void fmpi_remap_frontelements(Front *fr, Fmpi * f){

	reset_fmpi_arrays(f);			

	int map;
	size_t ssize = sizeof(int), rsize = sizeof(gID); 
	bool len = 1;
	foreach_frontelement(fr){
		gID * g = frontelement->e;
		for(int i=0; i<dimension; ++i, ++g) {	
			if( (g->pid) != pid() ) {
				map = f->mappid[g->pid]; //fixme: not fast
				assert(map != -1);
				array_append(f->sr[map].q,    &(g->alias), ssize);
				array_append(f->sr[map].temp, &elid, ssize);
				array_append(f->sr[map].temp, &i, ssize);
				len = 0;
			}
			else{
				g->pid   = frontelements[g->alias].pid;
				g->alias = frontelements[g->alias].alias;
				frontelement->n[i] =  (g->pid != pid()) ? -1 : g->alias;
			}
		}
	}

	fmpi(f, fmpi_remap_frontelements_fn1,
	        fmpi_remap_frontelements_fn2);

	if(len) return;
/*
	foreach_frontelement(fr){
		gID * g = frontelement->e;
		for(int i=0; i<dimension; ++i, ++g) {	
			if( (g->pid) == pid() ) {
				g->pid   = frontelements[g->alias].pid;
				g->alias = frontelements[g->alias].alias;
				frontelement->n[i] =  (g->pid != pid()) ? -1 : g->alias;
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
			gr = frontelement->e + (elid[1]); 
			memcpy(gr, g, rsize);
			frontelement->n[elid[1]] =  (gr->pid != pid()) ? -1 : gr->alias;
		}
	}

}

#endif