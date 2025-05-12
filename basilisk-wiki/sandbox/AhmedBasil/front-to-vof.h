#include "front-common.h" 

scalar volfrac[];

struct polyvertex{
	struct polyvertex * n, * p;
	double x[3];
	bool inside;
};

typedef struct polyvertex PolyVertex;
	
typedef struct {
	int nv;
	PolyVertex * p;
}Polygon;

PolyVertex * insert_vertex(Polygon * P, PolyVertex * p) {

	P->nv++;
	PolyVertex * newp = (PolyVertex *)malloc(sizeof(PolyVertex));
	PolyVertex * pn = p->n;

	pn->p = newp;
	p->n = newp;
	newp->n = pn;
	newp->p = p;

	return newp;
}

void delete_vertex(Polygon * P, PolyVertex * p) {
	P->nv--;
	PolyVertex * pp = p->p;
	PolyVertex * pn = p->n;

	P->p = NULL;
	free(p);
	if (P->nv == 0)
		return;

	P->p = pn;
	pp->n = pn;
	pn->p = pp;
}

Polygon * polygon_triangle(Front * fr, int elid) {

	Polygon * P = (Polygon *)malloc(sizeof(Polygon));
	P->nv = 3;
	PolyVertex * a = (PolyVertex *)malloc(sizeof(PolyVertex));
	PolyVertex * b = (PolyVertex *)malloc(sizeof(PolyVertex));
	PolyVertex * c = (PolyVertex *)malloc(sizeof(PolyVertex));

	a->n = b; b->n = c; c->n = a;
	a->p = c; c->p = b; b->p = a;

	P->p = a;
	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelement =&(fr->frontelements[elid]);
	for(int i=0;i<3; ++i) {
		a->x[i] = frontpoints[frontelement->v[0]].x[i];
		b->x[i] = frontpoints[frontelement->v[1]].x[i];
		c->x[i] = frontpoints[frontelement->v[2]].x[i];
	}
	return P;
}

@def foreach_polyvertices(P) 
	{
		PolyVertex * p = P->p, * n;
		int nv = P->nv;
		while (nv--)
		{
			n = p->n;
@
@def end_foreach_polyvertices() 
			p = n;
		}
	}
@

void free_polygon(Polygon * P) {
	foreach_polyvertices(P) {
		free(p); p = NULL;
	}
	free(P);
}


bool cut_the_polygon(Polygon * P, Point point) {
	foreach_polyvertices(P) {
		if(p->x[0] < x-0.5*Delta || p->x[0] > x+0.5*Delta) 
			p->inside = false;
		else
			p->inside = true;
	}
	foreach_polyvertices(P) {
		if((p->x[0] - x + 0.5*Delta) * (n->x[0] - x + 0.5*Delta) < 0.) {
			PolyVertex * newp = insert_vertex	(P, p);
			newp->inside = true;
			newp->x[0] = x - 0.5*Delta;
			double r = (newp->x[0] - p->x[0]) / (n->x[0] - p->x[0]); 
			newp->x[1] = p->x[1] + r*(n->x[1] - p->x[1]);
			newp->x[2] = p->x[2] + r*(n->x[2] - p->x[2]);
		}
		if((p->x[0] - x - 0.5*Delta) * (n->x[0] - x - 0.5*Delta) < 0.) {
			PolyVertex * newp = insert_vertex	(P, p);
			newp->inside = true;
			newp->x[0] = x + 0.5*Delta;
			double r = (newp->x[0] - p->x[0]) / (n->x[0] - p->x[0]); 
			newp->x[1] = p->x[1] + r*(n->x[1] - p->x[1]);
			newp->x[2] = p->x[2] + r*(n->x[2] - p->x[2]);
		}
	}
	
	foreach_polyvertices(P) {
		if(!p->inside) 
			delete_vertex(P, p);
	}
	foreach_polyvertices(P) {
		if(p->x[1] < y-0.5*Delta || p->x[1] > y+0.5*Delta) 
			p->inside = false;
		else
			p->inside = true;
	}
	foreach_polyvertices(P) {
		if((p->x[1] - y + 0.5*Delta) * (n->x[1] - y + 0.5*Delta) < 0.) {
			PolyVertex * newp = insert_vertex	(P, p);
			newp->inside = true;
			newp->x[1] = y - 0.5*Delta;
			double r = (newp->x[1] - p->x[1]) / (n->x[1] - p->x[1]); 
			newp->x[0] = p->x[0] + r*(n->x[0] - p->x[0]);
			newp->x[2] = p->x[2] + r*(n->x[2] - p->x[2]);
		}
		if((p->x[1] - y - 0.5*Delta) * (n->x[1] - y - 0.5*Delta) < 0.) {
			PolyVertex * newp = insert_vertex	(P, p);
			newp->inside = true;
			newp->x[1] = y + 0.5*Delta;
			double r = (newp->x[1] - p->x[1]) / (n->x[1] - p->x[1]); 
			newp->x[0] = p->x[0] + r*(n->x[0] - p->x[0]);
			newp->x[2] = p->x[2] + r*(n->x[2] - p->x[2]);
		}
	}
	
	foreach_polyvertices(P) {
		if(!p->inside) 
			delete_vertex(P, p);
	}

	foreach_polyvertices(P) {
		if(p->x[2] < z-0.5*Delta || p->x[2] > z+0.5*Delta) 
			p->inside = false;
		else
			p->inside = true;
	}
	foreach_polyvertices(P) {
		if((p->x[2] - z + 0.5*Delta) * (n->x[2] - z + 0.5*Delta) < 0.) {
			PolyVertex * newp = insert_vertex	(P, p);
			newp->inside = true;
			newp->x[2] = z - 0.5*Delta;
			double r = (newp->x[2] - p->x[2]) / (n->x[2] - p->x[2]); 
			newp->x[0] = p->x[0] + r*(n->x[0] - p->x[0]);
			newp->x[1] = p->x[1] + r*(n->x[1] - p->x[1]);
		}
		if((p->x[2] - z - 0.5*Delta) * (n->x[2] - z - 0.5*Delta) < 0.) {
			PolyVertex * newp = insert_vertex	(P, p);
			newp->inside = true;
			newp->x[2] = z + 0.5*Delta;
			double r = (newp->x[2] - p->x[2]) / (n->x[2] - p->x[2]); 
			newp->x[0] = p->x[0] + r*(n->x[0] - p->x[0]);
			newp->x[1] = p->x[1] + r*(n->x[1] - p->x[1]);
		}
	}
	
	foreach_polyvertices(P) {
		if(!p->inside) 
			delete_vertex(P, p);
	}
	return 1;
}