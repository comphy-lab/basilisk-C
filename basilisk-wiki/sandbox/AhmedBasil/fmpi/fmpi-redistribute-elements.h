#if _MPI
#include "../front-common.h"
#include "fmpi.h"

void fmpi_redistribute_frontelements_fn1 (Array * a, Array * b){
	if (a->len == 0) return;

	Front * fr = get_front();
	Frontelement * frontelements = fr->frontelements, * frontelement;
	size_t sizea = sizeof(double)*dimension, sizeb = dimension*sizeof(gID);
	size_t size = sizea + 2*sizeb; 
	char * c = (char *) (a->p);
	double x[dimension];
	int newe;
	for(int p=0; p< (int) (a->len/size); ++p){
    memcpy (x, c, sizea); c += sizea;
		assert(locate_coord(x).level >= 0);
		add_local_element(fr, newe);
		frontelement = frontelements + newe;
    memcpy (frontelement->p, c, sizeb); c += sizeb;
    memcpy (frontelement->e, c, sizeb); c += sizeb;
		for(int i=0; i<dimension; ++i){
			(frontelement->v)[i] =-1; (frontelement->n)[i] =-1;
		}
		array_append(b, &newe, sizeof(int));
	}
}

void fmpi_redistribute_frontelements_fn2 (Array * a, Array * b){
	array_append (b, a->p, a->len);
}

void fmpi_redistribute_frontelements(Front *fr, Fmpi * f){

	reset_fmpi_arrays(f);			

	int rpid;
	size_t ssizea = sizeof(double)*dimension, ssizeb = dimension*sizeof(gID);
	size_t ssize = ssizea + 2*ssizeb, rsize = sizeof(int); 
	Array * a; 
	bool len = 1;
	foreach_frontelement(fr){
		rpid = locate_elementrank();
		if(rpid == pid()) continue;
		{
			//fixme : expensive
			int p = f->mappid[rpid];
			assert(p != -1);
			a = f->sr[p].q;
			array_append(a, frontelement->c, ssizea); //fixme : not required. just to assert(); 
			array_append(a, frontelement->p, ssizeb);
			array_append(a, frontelement->e, ssizeb);
			array_append(f->sr[p].temp, &elid, sizeof(int));
			len = 0;
		}
	}

	fmpi(f, fmpi_redistribute_frontelements_fn1,
	        fmpi_redistribute_frontelements_fn2);

	if(len) return;

	Frontelement * frontelements = fr->frontelements;
	foreach_fmpi_sr(f){
		int * elid  = (int *) sr->temp->p, * alias = (int *) sr->r->p;
		int pid = sr->pid;
assert((ssize*sr->r->len) == (rsize*sr->q->len));	
		for(int i=0; i< (int) (sr->temp->len/sizeof(int)); ++i, alias++, elid++)	{
			//memcpy(&(frontelements[*elid].alias), alias, rsize);
			frontelements[*elid].alias =  *alias;
			frontelements[*elid].pid = pid;
		}
	}

}

#endif