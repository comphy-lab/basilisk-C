#if _MPI
#include "../front-common.h"
#include "fmpi.h"

#include "fmpi-coordinates.h"
#include "fmpi-redistribute-elements.h"
#include "fmpi-remap-elements.h"
#include "fmpi-remap-points.h"
#include "fmpi-connectivity.h"
#include "fmpi-redistribute-points.h"

void redistribute(Front *fr){

	if(npe() == 1)  return;

	Fmpi * f = fmpi_new();

	fmpi_coordinates(fr, f); 

	foreach_frontelement(fr) {
		frontelement->c[0] = centroidx();
		frontelement->c[1] = centroidy();
#if dimension == 3
		frontelement->c[2] = centroidz();
#endif
	}
	//B) delete NL frontelems
	foreach_frontelement_all(fr) {
		if (local_element() == 0) 
			delete_element(fr, elid);
	}
	//B) delete NL frontpoints
	foreach_frontpoint_all(fr) {
		if(local_point() == 0)
			delete_point(fr, ptid);
	}

	//A) derefernce all NL corners
	foreach_frontelement(fr) {
		for(int i=0; i<dimension; ++i) 
			if((frontelement->p)[i].pid != pid()) 
				(frontelement->v)[i] = -1; //derefernce //fixme:commented
			
		for(int i=0; i<dimension; ++i) 
			if((frontelement->e)[i].pid != pid()) 
				frontelement->n[i] = -1; //derefernce //fixme:commented
			
	}

	//C) redistribute frontpoints
	//redistribute_frontpoints(fr, redistributep);
	fmpi_redistribute_frontpoints(fr, f);

	//D) redistribute frontelements
	//redistribute_frontelements(fr, redistributee);
	fmpi_redistribute_frontelements(fr, f);

	//E) remap/update gID of corners
		//E.1) query NL corners(points)
		//E.2) remap corners (using local_remap and remap_query_return)
	//remap_frontpoints(fr, remap);
	fmpi_remap_frontpoints(fr, f);

	//F) remap/update gID of nbrs
		//F.1) query NL nbrs(elems)
		//F.2) remap nbrs (using local_remap and remap_query_return)
	//remap_frontelements(fr, remap);
	fmpi_remap_frontelements(fr, f);

	//G) delete NL frontpoints(MO)
	foreach_frontpoint_all(fr) {
		if(local_point() == 0) 
			delete_point(fr, ptid);
	}

	//H) delete NL frontelements(MO)
	foreach_frontelement_all(fr) {
		if(local_element() == 0)
			delete_element(fr, elid);
	}

	//I) create local instances of NL nbrs (rereference)
	int newe;
	foreach_frontelement(fr) {
		for(int i=0; i<dimension; ++i) {
			if((frontelement->e)[i].pid != pid()){
				add_element(fr, newe);
				frontelements[newe].pid = (frontelement->e)[i].pid;
				frontelements[newe].alias = (frontelement->e)[i].alias;
				frontelement->n[i] = newe;//rereferncing
			}
		}
	}

	//J) updates connectivity of NL nbrs
	fmpi_connectivity(fr, f);

	//K) create local instances of NL corners (rereference)
	int newp;
	foreach_frontelement_all(fr) {
		for(int i=0; i<dimension; ++i) {
			if((frontelement->p)[i].pid != pid()){
				add_point(fr, newp);
				frontpoints[newp].pid = (frontelement->p)[i].pid;
				frontpoints[newp].alias = (frontelement->p)[i].alias;
				(frontelement->v)[i] = newp;//rereferncing
			}
		}
	}

	//L) update coorndinates of NL corners	
	fmpi_coordinates(fr, f); 

	fmpi_free(f);
}
#endif