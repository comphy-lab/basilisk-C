#include "front-common.h"

#if dimension == 3

#define MAX_EDGES 10000 
typedef struct{
	int E1, E2, P1, P2;
}Frontedge;

typedef struct{
	int NS;
	Frontedge edges[MAX_EDGES];
}Frontedges;

Frontedges * default_edge_list = NULL;

Frontedges * get_frontedge(Front * fr){
	if(default_edge_list == NULL)
		default_edge_list = (Frontedges *)malloc(sizeof(Frontedges));

	Frontedges * fs = default_edge_list;
	Frontedge * edges = fs->edges;
	fs->NS = 0;
	int nbr;
	foreach_frontelement(fr){
		for(int n=0; n<3; ++n){
			nbr = frontelement->n[n];
			if(nbr > elid) {
				assert(fs->NS < MAX_EDGES);
				edges[fs->NS].E1 =  elid;
				edges[fs->NS].E2 =  nbr;
				edges[fs->NS].P1 =  frontelement->v[n];
				edges[fs->NS].P2 =  frontelement->v[(n+1)%3];
				fs->NS++;
			}			
		}
	}
	return (fs);
}

@def foreach_frontedge(fr, fs) 
	{
		Frontedge * edges = fs->edges; 
		Frontedge * edge;
		Frontelement * frontelements = fr->frontelements;
		Frontpoint * frontpoints =  fr->frontpoints;
		NOT_UNUSED(edges);
		NOT_UNUSED(edge);
		NOT_UNUSED(frontelements);
		NOT_UNUSED(frontpoints);
		for (int ie=0; ie < fs->NS; ++ie){
			edge = edges+ie;
@
@def end_foreach_frontedge()
		}
	}
@

event end(t= end) {
	if(default_edge_list)
		free(default_edge_list);
	default_edge_list = NULL;
}

#endif

