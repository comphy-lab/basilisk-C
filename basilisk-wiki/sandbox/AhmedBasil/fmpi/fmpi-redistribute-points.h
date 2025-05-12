#if _MPI
#include "../front-common.h"
#include "fmpi.h"

void fmpi_redistribute_frontpoints_fn1 (Array * a, Array * b){
	if (a->len == 0) return;

	Front * fr = get_front();
	Frontpoint * frontpoints = fr->frontpoints, * frontpoint;
	size_t size = dimension*sizeof(double);
	double * x = (double *) (a->p);
	int newp;
	for(int p=0; p< (int) (a->len/size); ++p, x += dimension){
		add_local_point(fr, newp);
		frontpoint = frontpoints + newp;
    memcpy (frontpoint->x, x, size);
		assert(locate_coord(x).level >= 0);
		array_append(b, &newp, sizeof(int));
	}
}

void fmpi_redistribute_frontpoints_fn2 (Array * a, Array * b){
	array_append (b, a->p, a->len);
}

void fmpi_redistribute_frontpoints(Front *fr, Fmpi * f){

	reset_fmpi_arrays(f);			

	int rpid;
	size_t ssize = dimension*sizeof(double), rsize = sizeof(int); 
	Array * a; 
	bool len = 1;
	foreach_frontpoint(fr){
		rpid = locate_pointrank();
		if(rpid == pid()) continue;
		{
			//fixme : expensive
			int p = f->mappid[rpid];
			assert(p != -1);
			a = f->sr[p].q;
			array_append(a, frontpoint->x, ssize);
			array_append(f->sr[p].temp, &ptid, sizeof(int));
			len = 0;
		}
	}

	fmpi(f, fmpi_redistribute_frontpoints_fn1,
	        fmpi_redistribute_frontpoints_fn2);

	if(len) return;

	Frontpoint * frontpoints = fr->frontpoints;
	foreach_fmpi_sr(f){
		int * ptid  = (int *) sr->temp->p, * alias = (int *) sr->r->p;
		int pid = sr->pid;
assert((ssize*sr->r->len) == (rsize*sr->q->len));	
		for(int i=0; i< (int) (sr->temp->len/sizeof(int)); ++i, alias++, ptid++)	{
			//memcpy(&(frontpoints[*ptid].alias), alias, rsize);
			frontpoints[*ptid].alias =  *alias;
			frontpoints[*ptid].pid = pid;
		}
	}

}

#endif