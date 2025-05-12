#if _MPI
#include "../front-common.h"
#include "fmpi.h"

void fmpi_connectivity_fn1 (Array * a, Array * b){
	if (a->len == 0) return;

	Front * fr = get_front();
	Frontelement * frontelements = fr->frontelements, * frontelement;
	size_t size = sizeof(int), bsize = dimension*sizeof(gID);
	int * elid = (int *) (a->p);
	for(int p=0; p< (int) (a->len/size); ++p, ++elid){
		frontelement = frontelements + (*elid);
		array_append(b, frontelement->p, bsize);
		array_append(b, frontelement->e, bsize);
	}
}

void fmpi_connectivity_fn2 (Array * a, Array * b){
	array_append (b, a->p, a->len);
}

void fmpi_connectivity(Front *fr, Fmpi * f){

	reset_fmpi_arrays(f);			

	int map;
	size_t ssize = sizeof(int), size = sizeof(gID);
	size_t rsize = 2*dimension*sizeof(gID); 
	bool len = 1;
	foreach_frontelement_all(fr){
		if(local_element()) break;
		map = f->mappid[frontelement->pid];
		assert(map != -1);
		array_append(f->sr[map].q,    &(frontelement->alias), ssize);
		array_append(f->sr[map].temp, &elid, ssize);
		len = 0;
	}

	fmpi(f, fmpi_connectivity_fn1,
	        fmpi_connectivity_fn2);

	if(len) return;

	Frontelement * frontelements = fr->frontelements, * frontelement;
	foreach_fmpi_sr(f){
		int * elid  = (int *) sr->temp->p;
		gID * g  = (gID *) sr->r->p, * g2;
assert((ssize*sr->r->len) == (rsize*sr->q->len));	
		for(int i=0; i< (int) (sr->temp->len/(ssize)); ++i, ++elid)	{
			frontelement = frontelements + elid[0];
			g2 = frontelement->p;
			for(int j=0; j<dimension; ++j, ++g2) {
				memcpy(g2, g, size); 
				frontelement->v[j] =  (g->pid != pid()) ? -1 : g->alias;
				++g;
			}
			g2 = frontelement->e;
			for(int j=0; j<dimension; ++j, ++g2) {
				memcpy(g2, g, size); 
				frontelement->n[j] =  (g->pid != pid()) ? -1 : g->alias;
				++g;
			}
		}
	}

}

#endif