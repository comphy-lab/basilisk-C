#include "front-utils.h"
#include "vofi.h"

#if _FTERMINAL
#if _fperiodicx || dimension == 3
void THIS_IS_A_NON_COMPILABLE_LINE (3);
#endif
#endif

#if _fperiodicx
#define clip_fperiodic(x)           ((x<0) ? (x+L0) : (x<L0) ? x : (x-L0))
#define apply_fperiodicx(x)         (X0 + clip_fperiodic(x-X0))
#define apply_fperiodic(x,x0)       (X0 + clip_fperiodic(x-x0))
#define apply_fperiodic_cts(x,x0)   (x0 + clip_fperiodic(x-x0))
#define apply_fperiodic_diff(x1,x2,ref)   (clip_fperiodic(x1-ref)-clip_fperiodic(x2-ref))
#define integer_fperiodic(x)        ((x<GHOSTS)  ?  (x + (1<<point.level)) : \
                                    (x<(GHOSTS + (1<<point.level)))  ?  x : \
                                    (x - (1<<point.level)))
#endif

#define FDEBUG  0

#if FDEBUG == 1
FILE * fdbgout = NULL;

void reset_fdbgout() {
	bool opened = false;
	if (fdbgout != NULL) 
		fclose(fdbgout);
	//fdbgout = NULL;
	fdbgout = lfopen("dumpdirectory/fdebug", "w");
}

void reopen_fdbgout() {
	bool opened = false;
	if (fdbgout != NULL) {
		opened = true;
		fclose(fdbgout);
		fdbgout = NULL;
	}
	//MPI_Barrier make sure interrupt from other processors(if any) happens only after the file closes
	MPI_Barrier(MPI_COMM_WORLD);	
	if (opened)
		fdbgout = lfopen("dumpdirectory/fdebug", "a");
	else
		fdbgout = lfopen("dumpdirectory/fdebug", "w");
}

void reopen_fdbgout_local() {
	bool opened = false;
	if (fdbgout) {
		opened = true;
		fclose(fdbgout);
		fdbgout = NULL;
	}
	if (opened)
		fdbgout = lfopen("dumpdirectory/fdebug", "a");
	else
		fdbgout = lfopen("dumpdirectory/fdebug", "w");
}

event init(t = 0) {
	reopen_fdbgout();
}

event end(t=end) {
	if (fdbgout)
		fclose(fdbgout);
}
#endif

#define MAX_FRONTS 10
#define MAX_ELEMS 100000
#define MAX_POINTS 100000

typedef struct{
	double amax;
	double amin;
#if dimension == 3
	double aspmax;
#endif
}Front_Regrid;

Front_Regrid * frregrid = NULL;

void set_regrid_props(int level) {
	/*level: element length is of the order of L0/(1<<level) . ~0.5 L0 /(1<<level)
	 *     Default value, level = grid->maxdepth
	 *     if(level <=0 ||level > grid->maxdepth), Then level = maxdepth; */

	if(frregrid) return;
	if(!grid) {
		fprintf(stdout, "\n Grid not yet initialised.\
		                 \n Regrid properties are relative to grid delta_min	\
		                 \n Initialise grid first. Exiting");
		assert(false); 
	}

	int depth = grid->maxdepth;
#if _MPI
	if (npe() == 1) 
		depth = grid->depth;
#endif
	int h = (level <=0 || level >= depth) ? depth : level;

	double dmin = L0/(1<<h);
	frregrid = (Front_Regrid *) malloc(sizeof(Front_Regrid));
	frregrid->amin = 0.32*dmin; 
	//values taken from paris-code
	frregrid->amax = 0.65*dmin;
	//frregrid->amax = 0.96*dmin;
#if dimension == 3
	//frregrid->aspmax = 1.4;	
	frregrid->aspmax = 1.54;	
#endif
}

scalar c_refer[];

typedef struct{
	//Coordinates
	double x[dimension];
	//Tension
	double T[dimension];
#if _MPI
	//parent processor id whose domain cntains this point
	int pid;
	//point Id in the parent processor
	int alias;
#endif
#if _FTERMINAL
	//ghost or not
	bool ghost;
#endif
}Frontpoint;

typedef struct{
	int pid;
	int alias;
}gID;	//global ID (unique identification of pt/element)

typedef struct{
	//integer references to  the vertices
	int v[dimension];
	//integer ref to neighbour elements
	int n[dimension];
#if _MPI
	double c[dimension]; //fixme : delete this
	gID p[dimension];
	gID e[dimension]; 
	//parent processor id whose domain cntains the centroid of this element
	int pid;
	//element Id in the parent processor (in case frontelement.pid != pid() )
	int alias;
#endif
#if _FTERMINAL
	//ghost or not
	bool ghost;
#endif
}Frontelement;

typedef struct {
	//An array of objects of Frontpoint
	Frontpoint *frontpoints;
	int *PConnectNext;
	int *PConnectPrev;
	bool *PAllocated;
	int *NP, *NPEmpty;
	int *FirstPoint;
	int *FirstEmptyPoint;
#if _FTERMINAL
	int P[2], E[2]; 
	void (*boundary[2]) (bool);
	void (*undo_boundary[2]) (bool);
#endif

	int *point_id; //fixme : rename later

	//An array of objects of Frontelement
	Frontelement *frontelements;
	int *EConnectNext;
	bool *EAllocated;
	int *EConnectPrev;
	int *NE,*NEEmpty;
	int *FirstElement;
	int *FirstEmptyElement;
}Front;


typedef struct {
	//An array of objects of Frontpoint
	Frontpoint *frontpoints;
	int *PConnectNext;
	int *PConnectPrev;
	bool *PAllocated;
	int *NP, *NPEmpty;
	int *FirstPoint;
	int *FirstEmptyPoint;
#if _FTERMINAL
	int P[2], E[2]; 
	void (*boundary[2]) (bool);
	void (*undo_boundary[2]) (bool);
#endif

	int *point_id; //fixme : rename later

	//An array of objects of Frontelement
	Frontelement *frontelements;
	int *EConnectNext;
	bool *EAllocated;
	int *EConnectPrev;
	int *NE,*NEEmpty;
	int *FirstElement;
	int *FirstEmptyElement;
}Membrane;

typedef double * fpscalar;
typedef double * fescalar;

fescalar new_fescalar() {
	fescalar val = (fescalar)malloc(MAX_ELEMS*sizeof(double));
	return val;
}

void free_fescalar(fescalar val){
	free(val);
}

fpscalar new_fpscalar() {
	fpscalar val = (fpscalar)malloc(MAX_POINTS*sizeof(double));
	return val;
}

#define free_fpscalar free_fescalar

typedef struct {
	fpscalar x;
	fpscalar y;
#if dimension == 3
	fpscalar z;
#endif
}fpvector;

typedef struct {
	fescalar x;
	fescalar y;
#if dimension == 3
	fescalar z;
#endif
}fevector;

fpvector * new_fpvector() {
	fpvector * v = (fpvector *)malloc(sizeof(fpvector));
	v->x = new_fpscalar();
	v->y = new_fpscalar();
#if dimension == 3
	v->z = new_fpscalar();
#endif
	return v;
}

void allocate_fpvector(fpvector * v) {
	v->x = new_fpscalar();
	v->y = new_fpscalar();
#if dimension == 3
	v->z = new_fpscalar();
#endif
}

/**
	fixme : There is no provision to make sure..
  ..memory allocation doesn't happen to an pointer which already..
  ..has mem. Also freeing doesn;t happen for a non-pointing pointer. ..
	..Mem leak is possible. 
	Another problem: Allocating+Deallocating huge mem from heap slow down the process*/

void free_fpvector(fpvector * v){
	free_fpscalar(v->x);
	free_fpscalar(v->y);
#if dimension == 3
	free_fpscalar(v->z);
#endif
}

fevector * new_fevector() {
	fevector * v = (fevector *)malloc(sizeof(fevector));
	v->x = new_fescalar();
	v->y = new_fescalar();
#if dimension == 3
	v->z = new_fescalar();
#endif
	return v;
}

void allocate_fevector(fevector * v) {
	v->x = new_fescalar();
	v->y = new_fescalar();
#if dimension == 3
	v->z = new_fescalar();
#endif
}

void free_fevector(fevector * v){
	free_fescalar(v->x);
	free_fescalar(v->y);
#if dimension == 3
	free_fescalar(v->z);
#endif
}

#if dimension==2
Front *fr_circle = NULL;
#elif dimension==3
Front *fr_sphere = NULL;
#endif

/**
foreach_frontpoint(fr) is macro that iterates through each marker points of the front*/
@def foreach_frontpoint_all(fr) 
	{
		int * PrevPoint   = fr->PConnectPrev;
		int * point_id    = fr->point_id;	
		int NP            = *(fr->NP);
		int FirstPoint    = *(fr->FirstPoint);	
		int nextpt;
		NOT_UNUSED(point_id);

		Frontpoint * frontpoints = fr->frontpoints;
		Frontpoint * frontpoint;
		NOT_UNUSED(frontpoint);
		NOT_UNUSED(frontpoints);
		int ptid = FirstPoint;
		for (int ipt=0; ipt<NP; ++ipt, ptid = nextpt) {
			frontpoint = frontpoints + ptid;
			nextpt = PrevPoint[ptid];
		
@
@def end_foreach_frontpoint_all()
		}
	}
@
/**
To see if the point is local */
#if _MPI
#define local_point() (pid()==frontpoint->pid)
#endif

#if _FTERMINAL
#define ghost_point() (frontpoint->ghost)
#endif
/**
Iterate through local elements only*/

@def foreach_frontpoint(fr) 
	foreach_frontpoint_all(fr)
#if _MPI
		if(local_point()) {
#endif
#if _FTERMINAL
			if(ghost_point() == false) {
#endif
@
@def end_foreach_frontpoint()
#if _FTERMINAL
			}
#endif
#if _MPI
		}
#endif
	end_foreach_frontpoint_all()
@

/**
foreach_frontelement(fr) is macro that iterates through each marker elements of the front*/
@def foreach_frontelement_all(fr) 
	{
		int * PrevElement = fr->EConnectPrev;
		int * point_id    = fr->point_id;	
		int NE            = *(fr->NE);
		int FirstElement  = *(fr->FirstElement);	
		int nextel;
		NOT_UNUSED(point_id);

		Frontelement * frontelements = fr->frontelements;
		Frontpoint * frontpoints = fr->frontpoints;
		Frontelement * frontelement;
		NOT_UNUSED(frontelement);
		NOT_UNUSED(frontelements);
		NOT_UNUSED(frontpoints);
		int elid = FirstElement;
		for (int ielem=0; ielem<NE; ++ielem, elid = nextel) {
			frontelement = frontelements + elid;
			nextel = PrevElement[elid];
@
@def end_foreach_frontelement_all()
		}
	}
@
/** 
To check if the element's centroid lies in this proc */
#if _MPI
#define local_element() (pid()==frontelement->pid)
#endif
#if _FTERMINAL
#define ghost_element() (frontelement->ghost)
#endif
/**
Iterate through lcal elements only*/

@def foreach_frontelement(fr) 
	foreach_frontelement_all(fr)
#if _MPI
		if(local_element()) {
#endif
#if _FTERMINAL
			if(ghost_element() == false) {
#endif
@
@def end_foreach_frontelement()
#if _FTERMINAL
			}
#endif
#if _MPI
		}
#endif
	end_foreach_frontelement_all()
@
/**
Macro to add or delete obj*/
@def add_element(fr, newe) 
	add_obj(fr->EConnectNext,  fr->EConnectPrev, &newe, fr->FirstElement,
          fr->FirstEmptyElement, fr->NE, fr->NEEmpty, fr->EAllocated);
#if _FTERMINAL
	fr->frontelements[newe].ghost =  false;
#endif
@

@def add_point(fr, newp) 
	add_obj(fr->PConnectNext, fr->PConnectPrev, &newp, fr->FirstPoint,
	           fr->FirstEmptyPoint, fr->NP, fr->NPEmpty, fr->PAllocated );
#if _FTERMINAL
	fr->frontpoints[newp].ghost =  false;
#endif
@


#if dimension == 2
#if !_MPI
@def foreach_frontelement_cyclic(fr) 
	{
		int * point_id    = fr->point_id;	
		int NE            = *(fr->NE);
		int FirstElement  = *(fr->FirstElement);	
		int nextel;
		NOT_UNUSED(point_id);

		Frontelement * frontelements = fr->frontelements;
		Frontpoint * frontpoints = fr->frontpoints;
		Frontelement * frontelement;
		NOT_UNUSED(frontelement);
		NOT_UNUSED(frontelements);
		NOT_UNUSED(frontpoints);
		int elid = FirstElement;
		for (int ielem=0; ielem<NE; ++ielem, elid = nextel) {
			frontelement = frontelements + elid;
			nextel = frontelement->n[0];
@
@def end_foreach_frontelement_cyclic()
		}
		assert(nextel == FirstElement);
	}
@
#endif
#endif

#if _MPI
#define add_local_element(fr, newe) \
	add_element(fr, newe); \
	frontelements[newe].pid = pid(); \
	frontelements[newe].alias = newe;

#define add_local_point(fr, newp) \
	add_point(fr, newp); \
	frontpoints[newp].pid = pid(); \
	frontpoints[newp].alias = newp;
#endif

#define delete_element(fr, olde) \
	delete_obj(fr->EConnectNext, fr->EConnectPrev, &olde, fr->FirstElement,\
	           fr->FirstEmptyElement, fr->NE, fr->NEEmpty, fr->EAllocated );

#define delete_point(fr, oldp) \
	delete_obj(fr->PConnectNext, fr->PConnectPrev, &oldp, fr->FirstPoint,\
	           fr->FirstEmptyPoint, fr->NP, fr->NPEmpty, fr->PAllocated );

#if dimension == 2
#define centroidx()  (0.5*(frontpoints[frontelement->v[0]].x[0] + frontpoints[frontelement->v[1]].x[0]))
#define centroidy()   (0.5*(frontpoints[frontelement->v[0]].x[1] + frontpoints[frontelement->v[1]].x[1]))
#elif dimension == 3
#define centroidx()   ((frontpoints[frontelement->v[0]].x[0] + frontpoints[frontelement->v[1]].x[0]\
                         + frontpoints[frontelement->v[2]].x[0])/3.0)
#define centroidy()   ((frontpoints[frontelement->v[0]].x[1] + frontpoints[frontelement->v[1]].x[1]\
                         + frontpoints[frontelement->v[2]].x[1])/3.0)
#define centroidz()   ((frontpoints[frontelement->v[0]].x[2] + frontpoints[frontelement->v[1]].x[2]\
                         + frontpoints[frontelement->v[2]].x[2])/3.0)
#endif

//Neighbours
#define Nbri(i,e) (frontelements + e->n[i])
#define Nbr0(e)   (frontelements + e->n[0])
#define Nbr1(e)   (frontelements + e->n[1])
#if dimension == 3
#define Nbr2(e)   (frontelements + e->n[2])
#endif
//Neighbour indices
#define NbrIi(i,e) e->n[i]
#define NbrI0(e)   e->n[0]
#define NbrI1(e)   e->n[1]
#if dimension == 3
#define NbrI2(e)   e->n[2]
#endif
//Vertices
#define Vertexi(i,e) (&frontpoints[e->v[i]])
#define Vertex0(e)   (&frontpoints[e->v[0]])
#define Vertex1(e)   (&frontpoints[e->v[1]])
#if dimension == 3
#define Vertex2(e)   (&frontpoints[e->v[2]])
#endif
//vertex indices
#define VertexIi(i,e) e->v[i]
#define VertexI0(e)   e->v[0]
#define VertexI1(e)   e->v[1]
#if dimension == 3
#define VertexI2(e)   e->v[2]
#endif

#define Vertexij(i,j,e) Vertexi(i,Nbri(j))
#if dimension == 2
#define Vertex00(e)     Vertex0(Nbr0(e))
#define Vertex11(e)     Vertex1(Nbr1(e))
#define VertexI00(e)    VertexI0(Nbr0(e))
#define VertexI11(e)    VertexI1(Nbr1(e))
#endif

#if dimension == 2
		//fixfperiodic
#define distance(x, y) sqrt((x[0]-y[0])*(x[0]-y[0]) \
                       + (x[1]-y[1])*(x[1]-y[1]))
bool front_connectivity_valid(Front * fr){
	int n;
	foreach_frontelement_all(fr) {
		n = frontelement->n[0];
		if(n != -1)
			if(elid != frontelements[n].n[1]) {
				return false;
			}

		n = frontelement->n[1];
		if(n!=-1)
			if(elid != frontelements[n].n[0]) {
				return false;
			}
	}
	return true;
}
#elif dimension == 3
#define distance(x, y) sqrt((x[0]-y[0])*(x[0]-y[0]) \
                       + (x[1]-y[1])*(x[1]-y[1])\
                       + (x[2]-y[2])*(x[2]-y[2]))

bool front_connectivity_valid(Front * fr){
	int pa, pb, pc;
	int nbr;
	bool found;
	foreach_frontelement(fr) {
		if(frontelement->v[0] == frontelement->v[1]) 
			return false;
		if(frontelement->v[2] == frontelement->v[1]) 
			return false;
		if(frontelement->v[0] == frontelement->v[2]) 
			return false;
		if(frontelement->n[0] == frontelement->n[1]) 
			return false;
		if(frontelement->n[2] == frontelement->n[1]) 
			return false;
		if(frontelement->n[0] == frontelement->n[2]) 
			return false;
		for(int i=0; i<3; ++i) {
			pa = frontelement->v[i];
			pb = frontelement->v[(i+1)%3];
			pc = frontelement->v[(i+2)%3];
			nbr = frontelement->n[i];
			found = 0;
			for(int j=0; j<3; ++j) 
				if(frontelements[nbr].n[j] == elid){
					if(frontelements[nbr].v[j] != pb) 
						return false;
					if (frontelements[nbr].v[(j+1)%3] != pa)
						return false;
					if (frontelements[nbr].v[(j+2)%3] == pc)
						return false;
					found = 1;
				}
			if(!found)
				return false;
		}
	}
	return true;
}
#endif

#if dimension == 2
void element_normal(double *n, double *p0, double *p1){
	double norm = distance(p0,p1);
	//fixme : check it again	
	n[0] = -(p1[1] - p0[1])/norm;
	n[1] =   p1[0] - p0[0]/norm;
} 
#elif dimension == 3
void element_normal(double *n, double *p0, double *p1, double *p2){
	double u[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
	double v[3] = {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};
	//fixme : check it again	
	n[0] = u[1]*v[2]-u[2]*v[1];
	n[1] = u[2]*v[0]-u[0]*v[2];
	n[2] = u[0]*v[1]-u[1]*v[0];
	double norm = n[0]*n[0] + n[1]*n[1] +n[2]*n[2];
	n[0] /= norm;
	n[1] /= norm;
	n[2] /= norm;
}
#endif

/**
Create an empty front;
*/
Front *create_empty_front(){
		
	Front *fr = (Front *)malloc(sizeof(Front));

	fr->PConnectNext = (int *)malloc((MAX_POINTS + 1) * sizeof(int));
	fr->PConnectPrev = (int *)malloc((MAX_POINTS + 1) * sizeof(int));
	fr->PConnectNext = fr->PConnectNext + 1;
	fr->frontpoints  = (Frontpoint *)malloc(MAX_POINTS*sizeof(Frontpoint));
	fr->PAllocated  =  (bool *)malloc(MAX_POINTS*sizeof(bool));
	for(int i=0; i<MAX_POINTS; ++i)
		fr->PAllocated[i] = false;

	fr->point_id = (int *)malloc((MAX_POINTS + 1) * sizeof(int));

	fr->EConnectNext = (int *)malloc((MAX_ELEMS + 1) * sizeof(int));
	fr->EConnectPrev = (int *)malloc((MAX_ELEMS + 1) * sizeof(int));
	fr->EConnectNext = fr->EConnectNext + 1;
	fr->frontelements = (Frontelement *)malloc(MAX_ELEMS*sizeof(Frontelement));
	fr->EAllocated  =  (bool *)malloc(MAX_ELEMS*sizeof(bool));
	for(int i=0; i<MAX_ELEMS; ++i)
		fr->EAllocated[i] = false;

	fr->NP = (int *)malloc(sizeof(int));
	fr->NPEmpty = (int *)malloc(sizeof(int));
	fr->FirstPoint = (int *)malloc(sizeof(int));
	fr->FirstEmptyPoint = (int *)malloc(sizeof(int));
	fr->NE = (int *)malloc(sizeof(int));
	fr->NEEmpty = (int *)malloc(sizeof(int));
	fr->FirstElement = (int *)malloc(sizeof(int));
	fr->FirstEmptyElement = (int *)malloc(sizeof(int));

	//Set  default values 	
	*(fr->NP) = 0;
	*(fr->NPEmpty) = MAX_POINTS;
	*(fr->NE) = 0;
	*(fr->NEEmpty) = MAX_ELEMS;

	//Setting default connectivity(Cyclic connection)
	int *connectf,*connectb;

	connectf = fr->PConnectNext;
	connectb = fr->PConnectPrev;
	for (int l=0; l<=MAX_POINTS ; ++l){
		connectf[l-1] = l;
		connectb[l] = l-1;
	}
	*(fr->FirstPoint) = -1;
	*(fr->FirstEmptyPoint) = 0;

	connectf = fr->EConnectNext;
	connectb = fr->EConnectPrev;	
	for (int l=0; l<=MAX_ELEMS ; ++l){
		connectf[l-1] = l;
		connectb[l] = l-1;
	}
	*(fr->FirstElement) = -1;
	*(fr->FirstEmptyElement) = 0;

	return(fr);
}

Front * get_front(){
#if dimension==2
	if (fr_circle == NULL)
		fr_circle = create_empty_front();
	return fr_circle;
#elif dimension==3		
	if (fr_sphere == NULL)
		fr_sphere = create_empty_front();
	return fr_sphere;
#endif
}

/** Reset Connectivity*/
void reset_front(Front *f){
	reset_obj_list(f->PConnectNext, f->PConnectPrev, MAX_POINTS,
	       f->FirstPoint, f->FirstEmptyPoint, f->NP, f->NPEmpty, f->PAllocated);

	reset_obj_list(f->EConnectNext, f->EConnectPrev, MAX_ELEMS,
	       f->FirstElement, f->FirstEmptyElement, f->NE, f->NEEmpty, f->EAllocated);
}

/**
Free memory*/
void free_front(Front *fr){
	if(fr == NULL) return;

	fr->PConnectNext = fr->PConnectNext - 1;
	free(fr->PConnectNext);	
	free(fr->PConnectPrev);	
	free(fr->point_id);
	free(fr->frontpoints);
	fr->EConnectNext = fr->EConnectNext - 1;
	free(fr->EConnectNext);	
	free(fr->EConnectPrev);	
	free(fr->frontelements);
	free(fr->PAllocated);
	free(fr->EAllocated);

	free(fr->NP);
	free(fr->NPEmpty);
	free(fr->FirstPoint);
	free(fr->FirstEmptyPoint);
	free(fr->NE);
	free(fr->NEEmpty);
	free(fr->FirstElement);
	free(fr->FirstEmptyElement);
	
	free(fr);

	fr = NULL;
	
	return;
}

#if _FTERMINAL
#define apply_fsymmetric_bottom(x) (2*Y0 - x)
//candidates for function pointer to applied at the two terminal ends of 
void fsymmetric_bottom(bool terminal){
	Front * fr = get_front();
	Frontpoint * frontpoints =  fr->frontpoints;
	Frontelement * frontelements = fr->frontelements;
	Frontelement * frontelement = frontelements + fr->E[terminal];
	assert(frontelement->n[terminal] == -1); 
	assert(frontelement->v[terminal] == fr->P[terminal]);

	int newp, newe;
	add_point(fr, newp); 
	frontpoints[newp].ghost = true;
	add_element(fr, newe); 
	frontelements[newe].ghost = true;
	
	//connectivity
	frontelements[newe].v[terminal] = newp;
	frontelements[newe].v[!terminal] = fr->P[terminal];
	frontelements[newe].n[!terminal] = fr->E[terminal];
	frontelement->n[terminal] = newe;

	//coordinates
	frontpoints[newp].x[0] = frontpoints[frontelement->v[!terminal]].x[0];
	frontpoints[newp].x[1] = 2*Y0 - frontpoints[frontelement->v[!terminal]].x[1];
}

void undo_fsymmetric_bottom(bool terminal){
	Front * fr = get_front();
	Frontpoint * frontpoints =  fr->frontpoints;
	Frontelement * frontelements = fr->frontelements;
	Frontelement * frontelement = frontelements + fr->E[terminal];
	int	nbr = frontelement->n[terminal];	
	assert(frontelements[nbr].ghost == true);
	int nbrc = frontelements[nbr].v[terminal];  
	assert(frontpoints[nbrc].ghost == true);

	delete_point(fr, nbrc); 
	delete_element(fr, nbr); 
	
	frontelement->n[terminal] = -1;
}

void fsymmetric_st_bottom(bool terminal){
	Front * fr = get_front();
	Frontpoint * frontpoint =  fr->frontpoints + fr->P[terminal];	
	frontpoint->T[0] *= 2.;
	frontpoint->T[1] = 0.;
}	
	
void front_ghost_boundary(Front *fr) {
	assert(fr->boundary[0]);
	fr->boundary[0](0);
	assert(fr->boundary[1]);
	fr->boundary[1](1);
}

void undo_front_ghost_boundary(Front *fr) {
	assert(fr->undo_boundary[0]);
	fr->undo_boundary[0](0);
	assert(fr->undo_boundary[1]);
	fr->undo_boundary[1](1);
}
#endif
 
/**
List containing all fronts.
*/
typedef struct{
	Front **FArray;
	int NF;
}FrontList;

FrontList *initialize_front_list(){
	FrontList *fl = (FrontList *)malloc(sizeof(FrontList *));
	fl->NF = 0;
	fl->FArray	= (Front **)malloc(MAX_FRONTS*sizeof(Front *)); 
	
	return(fl); 
}

/**
Push a front to front_list.
*/
int push_to_front_list(FrontList *front_list, Front *front){
	if (front_list->NF > MAX_FRONTS){
		//ERROR: cannot add one more front
		return(0); 
	}
	create_empty_front(front);
	front_list->FArray[front_list->NF]	=	front;
	front_list->NF += 1;		
	return(1);
}

#if dimension == 2
double elem_length(Front *fr, int elem){
	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelement = &(frontelements[elem]);
	int p0 = VertexI0(frontelement);
	int p1 = VertexI1(frontelement);
#if _fperiodicx
	double xmin = frontpoints[p0].x[0] - 0.25*L0;
	double d = sq(apply_fperiodic_cts(frontpoints[p0].x[0], xmin) - 
	              apply_fperiodic_cts(frontpoints[p1].x[0], xmin)) + 
	           sq(frontpoints[p0].x[1] - frontpoints[p1].x[1]);
#else
	double d = sq(frontpoints[p0].x[0]-frontpoints[p1].x[0]) + sq(frontpoints[p0].x[1]-frontpoints[p1].x[1]);
#endif
	return (sqrt(d));
}
#elif dimension == 3 
double elem_length(Front *fr, int elem, int index){
	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelement = &(frontelements[elem]);
	
	int p0 = VertexIi(index, frontelement);
	int p1 = VertexIi((index+1)%3, frontelement);
	return (distance(frontpoints[p0].x, frontpoints[p1].x));
}

double min_edge_length(Front *fr, int elem){
	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelement = &(frontelements[elem]);
	
	int p0 = VertexI0(frontelement);
	int p1 = VertexI1(frontelement);
	int p2 = VertexI2(frontelement);
	double d0, d1, d2;
	d0 = distance(frontpoints[p0].x, frontpoints[p1].x);
	d1 = distance(frontpoints[p1].x, frontpoints[p2].x);
	d2 = distance(frontpoints[p2].x, frontpoints[p0].x);
	return (min(min(d0, d1), d2));
}

double deletethis(Front *fr, int elem){
	Frontelement * frontelements = fr->frontelements;
	//Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelement = &(frontelements[elem]);
	
	int p0 = VertexI0(frontelement);
	int p1 = VertexI1(frontelement);
	int p2 = VertexI2(frontelement);
		fprintf(stdout, "Nbr %d %d %d", p0, p1, p2);
	return 0.;
}


double max_edge_length(Front *fr, int elem){
	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelement = &(frontelements[elem]);
	
	int p0 = VertexI0(frontelement);
	int p1 = VertexI1(frontelement);
	int p2 = VertexI2(frontelement);
	double d0, d1, d2;
	d0 = distance(frontpoints[p0].x, frontpoints[p1].x);
	d1 = distance(frontpoints[p1].x, frontpoints[p2].x);
	d2 = distance(frontpoints[p2].x, frontpoints[p0].x);
	return (max(max(d0, d1), d2));
}
#endif


/**
circle front
*/

#if !_MPI
#if dimension == 2
int add_circle(Front *fr, double xc, double yc,
                double radin){
	set_regrid_props(-1);	
	double elength = 0.5*(frregrid->amax + frregrid->amin);
	int np = (int) floor(2*pi*radin/elength);

	if(np > MAX_POINTS || np < 10) {
		printf("np should be 10 and MAX_POINTS");
		return(-1);
	}

	double dtheta = 2.0*pi/((double) np);

	reset_front(fr);

	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelements = fr->frontelements;

	*(fr->NP) = np;
	*(fr->NPEmpty) = MAX_POINTS - np;
	*(fr->FirstPoint) = np-1;
	*(fr->FirstEmptyPoint) = np;

	*(fr->NE) = np;
	*(fr->NEEmpty) = MAX_ELEMS - np;
	*(fr->FirstElement) = np-1;
	*(fr->FirstEmptyElement) = np;

	for(int p=0; p<np; ++p){
		frontpoints[p].x[0] = xc + radin*cos(p*dtheta);
		frontpoints[p].x[1] = yc + radin*sin(p*dtheta);

		frontelements[p].v[0] = p;
		frontelements[p].v[1] = p+1;

		frontelements[p].n[0] = p-1;
		frontelements[p].n[1] = p+1;
	}

	frontelements[np-1].v[1] = 0;
	frontelements[np-1].n[1] = 0;
	frontelements[0].n[0] = np-1;

	foreach_frontpoint(fr)
		fr->PAllocated[ptid] = true;

	foreach_frontelement(fr)
		fr->EAllocated[elid] = true;

	return(0);

}

#if _FTERMINAL
int add_semicircle(Front *fr, double xc, double yc_obsolete,
                double radin){
	double yc = Y0;
	set_regrid_props(-1);	
	double elength = 0.5*(frregrid->amax + frregrid->amin);
	int np = (int) floor(pi*radin/elength);

	if(np > MAX_POINTS || np < 10) {
		printf("np should be 10 and MAX_POINTS");
		return(-1);
	}

	double dtheta = pi/((double) (np-1));

	reset_front(fr);

	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelements = fr->frontelements;

	*(fr->NP) = np;
	*(fr->NPEmpty) = MAX_POINTS - np;
	*(fr->FirstPoint) = np-1;
	*(fr->FirstEmptyPoint) = np;

	int ne = np-1;
	*(fr->NE) = ne;
	*(fr->NEEmpty) = MAX_ELEMS - ne;
	*(fr->FirstElement) = ne-1;
	*(fr->FirstEmptyElement) = ne;

	for(int p=0; p<ne; ++p){
		frontpoints[p].x[0] = xc + radin*cos(p*dtheta);
		frontpoints[p].x[1] = yc + radin*sin(p*dtheta);

		frontelements[p].v[0] = p;
		frontelements[p].v[1] = p+1;

		frontelements[p].n[0] = p-1;
		frontelements[p].n[1] = p+1;
	}
	frontpoints[np-1].x[0] = xc - radin;
	frontpoints[np-1].x[1] = yc;
	//overwriting
	frontpoints[0].x[1] = yc;
	frontelements[ne-1].n[1] = -1;

	foreach_frontpoint_all(fr) {
		fr->PAllocated[ptid] = true;
		frontpoint->ghost =  false;
	}

	foreach_frontelement_all(fr) {
		fr->EAllocated[elid] = true; 
		frontelement->ghost =  false;
	}

	fr->E[0] = 0;
	fr->E[1] = ne-1;
	fr->P[0] = 0;
	fr->P[1] = np-1;

//fixme: setting the function pointers
	fr->boundary[0] = fsymmetric_bottom; 
	fr->boundary[1] = fsymmetric_bottom; 
	fr->undo_boundary[0] = undo_fsymmetric_bottom; 
	fr->undo_boundary[1] = undo_fsymmetric_bottom; 
	return(0);

}
#endif

int add_xcosinewave(Front *fr, double yc, double ampl){
	
	if((yc+ampl > Y0+L0) || (yc - ampl < Y0))
		return(-2);

	set_regrid_props(-1);	
	double elength = 0.5*(frregrid->amax + frregrid->amin);
	int np = (int) floor(L0/elength);
	np = np - np%2;

	if(np > MAX_POINTS || np < 10) {
		printf("np should be 10 and MAX_POINTS");
		return(-1);
	}

	double dtheta = L0/np;

	reset_front(fr);

	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelements = fr->frontelements;

	*(fr->NP) = np;
	*(fr->NPEmpty) = MAX_POINTS - np;
	*(fr->FirstPoint) = np-1;
	*(fr->FirstEmptyPoint) = np;

	*(fr->NE) = np;
	*(fr->NEEmpty) = MAX_ELEMS - np;
	*(fr->FirstElement) = np-1;
	*(fr->FirstEmptyElement) = np;

	for(int p=0; p<np; ++p){
		frontpoints[p].x[0] = X0 + p*dtheta;
		frontpoints[p].x[1] = yc + ampl*cos(2*pi*p*dtheta/L0);

		frontelements[p].v[0] = p;
		frontelements[p].v[1] = p+1;

		frontelements[p].n[0] = p-1;
		frontelements[p].n[1] = p+1;
	}

	frontelements[np-1].v[1] = 0;
	frontelements[np-1].n[1] = 0;
	frontelements[0].n[0] = np-1;

	foreach_frontpoint(fr) 
		fr->PAllocated[ptid] = true;

	foreach_frontelement(fr) 
		fr->EAllocated[elid] = true;

	return(np);

}

#elif dimension == 3 
int add_sphere(Front *fr, double xc, double yc, double zc,
                double radin){

	int iip, ist, iie,
	    ia, ib, ic, id, **icp, **ine,
	    iqq, ne, np;
	double _PI, dph, phi, theta, **pt;

	set_regrid_props(-1);
	double elength = 0.5*(frregrid->amax + frregrid->amin);
	int nps = floor(radin/elength * sqrt(2*pi/sqrt(3.)));
	
	if(nps < 5) {
		printf("surface grid resolution will be poor. Increase no of grids per diameter of sphere");
		exit(-1);
	}
	np = 4*nps*nps+2;
	ne = 4*2*nps*nps;
	if(np > MAX_POINTS || ne >MAX_ELEMS) {
		printf("\nRequires larger memory pool to store surface grid");
		exit(-1);
	}

	_PI = 4.0 * atan(1.0);
	dph = 0.5*_PI/((double) nps);

	icp  = (int **)malloc((4*2*nps*nps)*sizeof(int *));
	ine  = (int **)malloc((4*2*nps*nps)*sizeof(int *));
  for(int l = 0; l< 8*nps*nps; ++l){
		icp[l] = (int *)malloc(3*sizeof(int));
		ine[l] = (int *)malloc(3*sizeof(int));
	}
	pt  = (double **)malloc((4*nps*nps+2)*sizeof(double *));
  for(int l = 0; l< 4*nps*nps+2; ++l){
		pt[l] = (double *)malloc(3*sizeof(double));
	}

	pt[0][0] = xc; pt[0][1] = yc; pt[0][2] = zc+radin;
	pt[1][0] = xc; pt[1][1] = yc; pt[1][2] = zc-radin;
	for(int iq = 1; iq<=4; ++iq){
		for(int i2 = 1; i2<=nps; ++i2){
			for(int i1 = 1; i1<=nps; ++i1){

				iip = (iq-1)*nps*nps + (i2-1)*nps + i1 + 2;
				phi = dph * ((double) i1-i2);
				ist = i2-1;
				if(i1-i2 < 0) ist = i1-1;
				theta = 0.50*_PI*( ((double) iq-1) + ((double) ist)/(nps - fabs(i1-i2)) );
				pt[iip - 1][0] = xc + radin*cos(phi)*cos(theta);
				pt[iip - 1][1] = yc + radin*cos(phi)*sin(theta);
				pt[iip - 1][2] = zc + radin*sin(phi);

				iie = 2*i1+2*nps*(i2-1)+2*nps*nps*(iq-1);
				ia = iip;
				ib = iip+nps;
				ic = ib+1;
				id = ia+1;
				if(i1 == nps) {
					iqq=iq;
					if(iqq == 4)	iqq=0;
					ic = 2+iqq*nps*nps+nps+1-i2;
					id = ic+1;
				}
				if(i2 == nps) {
					iqq=iq;
					if(iqq == 4) iqq = 0;
					ib = 2+iqq*nps*nps+(nps+1-i1)*nps+1;
					ic = ib-nps;
				}
				if((i1 == nps) && (i2 == 1)) id=1;
				if((i2 == nps) && (i1 == 1)) ib=2;
				icp[iie-2][0] = ia-1;
				icp[iie-2][1] = ib-1;
				icp[iie-2][2] = ic-1;
				icp[iie-1][0] = ia-1;
				icp[iie-1][1] = ic-1;
				icp[iie-1][2] = id-1;
				ine[iie-2][0] = iie-3;
				ine[iie-2][1] = iie+2*nps-1;
				ine[iie-2][2] = iie-1;
				ine[iie-1][0] = iie-2;
				ine[iie-1][1] = iie;
				ine[iie-1][2] = iie-2*nps-2;
				if(i1 == 1  ) {
					iqq=iq-1;
					if(iqq == 0) iqq=4;
					ine[iie-2][0] = iqq*2*nps*nps-2*i2;
				}
				if(i1 == nps) {
					iqq=iq;
					if(iqq == 4) iqq=0;
					ine[iie-1][1] = iqq*2*nps*nps+2*(nps+1-i2)-1;
				}
				if(i2 == 1  ) {
					iqq = iq-1;
					if(iqq == 0) iqq=4;
					ine[iie-1][2] = iqq*2*nps*nps-2*nps*(i1-1)-1;
				} 
				if(i2 == nps) {	
					iqq = iq;
					if(iqq == 4) iqq=0;
					ine[iie-2][1] = iqq*2*nps*nps+2*nps*(nps-i1);
				}

			}
		}
	}

	reset_front(fr);

	//equivalent to adding np points and ne elements
	*(fr->NP) = np;
	*(fr->NPEmpty) = MAX_POINTS - np;
	*(fr->FirstPoint) = np - 1;
	*(fr->FirstEmptyPoint) = np;
	*(fr->NE) = ne;
	*(fr->NEEmpty) = MAX_ELEMS - ne;
	*(fr->FirstElement) = ne - 1;
	*(fr->FirstEmptyElement) = ne;

	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;

	for(int pid=0; pid<np; ++pid){
		frontpoints[pid].x[0] = pt[pid][0];
		frontpoints[pid].x[1] = pt[pid][1];
		frontpoints[pid].x[2] = pt[pid][2];
	} 

	for(int eid=0; eid<ne; ++eid){
		frontelements[eid].v[0] = icp[eid][0];
		frontelements[eid].v[1] = icp[eid][1];
		frontelements[eid].v[2] = icp[eid][2];
	}

	for(int eid=0; eid<ne; ++eid){
		frontelements[eid].n[0] = ine[eid][0];
		frontelements[eid].n[1] = ine[eid][1];
		frontelements[eid].n[2] = ine[eid][2];
	}

	//free mem
	for(int i=0; i<8*nps*nps; ++i) {
		free(icp[i]);
		free(ine[i]);
	}
	free(icp);
	free(ine);

	for(int i=0; i<4*nps*nps+2; i++){
		free(pt[i]);
	}
	free(pt);

	foreach_frontpoint(fr)
		fr->PAllocated[ptid] = true;

	foreach_frontelement(fr)
		fr->EAllocated[elid] = true;

	//assert(front_connectivity_valid(fr));

	return(0);

}
#endif 	//if dimension
#endif	//if !_MPI


/**
Routines for regrid*/



/**
Finding the quality of a grid.
Following function returns 0 if all the three 
quality standards are satisfied. Returns 1 otherwise
*/
#if dimension == 2 
bool front_quality(Front *fr){
	
	double elength;
	bool q = true, qg;
	foreach_frontelement(fr) {
		//fixme : fperiodic
		elength = elem_length(fr, elid);	
		if(elength < frregrid->amin || elength > frregrid->amax)  {
			q = false;
			break;
		}
	}
#if _MPI	
	MPI_Allreduce(&q, &qg, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD); 
#else
	qg = q;
#endif

	return qg;
}
#elif dimension == 3 
double elem_area(Front * fr, int elid){
	Frontpoint * fp = fr->frontpoints;
	Frontelement * fe = fr->frontelements;
	double * p0 = fp[fe[elid].v[0]].x;
	double * p1 = fp[fe[elid].v[1]].x;
	double * p2 = fp[fe[elid].v[2]].x;
	double area, x0, y0, z0, x1, y1, z1;
	x0 = p1[0] - p0[0]; 		
	y0 = p1[1] - p0[1]; 		
	z0 = p1[2] - p0[2]; 		
	x1 = p2[0] - p0[0]; 		
	y1 = p2[1] - p0[1]; 		
	z1 = p2[2] - p0[2]; 		

	//Area  = 0.5 | u x v |	
	area = 0.50 * sqrt((y0*z1 - y1*z0)*(y0*z1 - y1*z0) + 
	                   (x1*z0 - x0*z1)*(x1*z0 - x0*z1) +
	                   (x0*y1 - x1*y0)*(x0*y1 - x1*y0) );
	return (area);
}

double aspect_ratio(Front * fr, int elid){
	double s = elem_length(fr, elid, 0) +
	           elem_length(fr, elid, 1) +
	           elem_length(fr, elid, 2);
	s /= 3.;
	double area = elem_area(fr, elid);
	return (0.25*sqrt(3.)*s*s/area);
}

int front_quality(Front *fr){
	double amin = sq(frregrid->amin)*1.7/4.;
	double amax = sq(frregrid->amax)*1.7/4.;
	double earea;
	foreach_frontelement(fr) {
		//fixfperiodic
		earea = elem_area(fr, elid);	
		if(earea < amin) return(1);
		if(earea > amax) return(1);
		//fixme: include aspect ratio constraint also
	}
	return 0;
	//fixme : not yet implemented
}
#endif
//-------------------------------------------------------------------

#define nbr_int(elemid, indxid)  fr->frontelements[elemid].n[indxid]
#define vertex_int(elemid, indxid)  fr->frontelements[elemid].v[indxid]

//-----------------------------------------------------------------------

#if  !_MPI

#if dimension == 2
int delete_front_element(Front *fr, int elid){
#if _FTERMINAL
	assert(elid != fr->E[0] && elid != fr->E[1]);
#endif
	fprintf(stdout, "-");
	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelement = &(frontelements[elid]);
	int p0, p1, p00, p11;
	int n0, n1;

	n0 = NbrI0(frontelement);
	n1 = NbrI1(frontelement);

	p0  = VertexI0(frontelement);
	p1  = VertexI1(frontelement);
	p00 = VertexI00(frontelement);
	p11 = VertexI11(frontelement);
	
	//Conditions to abandon deleting.
	//I don't know if any exists.
	
#if _fperiodicx
	double xp = frontpoints[p00].x[0] - 0.25*L0;
	frontpoints[p0].x[0] = ( - apply_fperiodic_cts( frontpoints[p00].x[0], xp) +
	                     9.0 * apply_fperiodic_cts( frontpoints[p0].x[0],  xp) + 
	                     9.0 * apply_fperiodic_cts( frontpoints[p1].x[0],  xp) - 
	                           apply_fperiodic_cts( frontpoints[p11].x[0], xp))/16.;
	frontpoints[p0].x[0] = apply_fperiodicx(frontpoints[p0].x[0]);
#else
	frontpoints[p0].x[0] = ( - frontpoints[p00].x[0] +  9.0 *(frontpoints[p0].x[0] + frontpoints[p1].x[0]) - frontpoints[p11].x[0] )/16.;
#endif
	frontpoints[p0].x[1] = ( - frontpoints[p00].x[1] +  9.0 *(frontpoints[p0].x[1] + frontpoints[p1].x[1]) - frontpoints[p11].x[1] )/16.;


	//Delete Element elid and p1
	delete_obj(fr->EConnectNext, fr->EConnectPrev, &elid, fr->FirstElement,
	           fr->FirstEmptyElement, fr->NE, fr->NEEmpty, fr->EAllocated);
	delete_obj(fr->PConnectNext, fr->PConnectPrev, &p1, fr->FirstPoint,
	           fr->FirstEmptyPoint, fr->NP, fr->NPEmpty, fr->PAllocated );

	//Resetting neighbour elements around deleted elements
	frontelements[n0].n[1] = n1;
	frontelements[n1].n[0] = n0;
	
	//Reset vertices of all elements having p1 as a vertex from p1 to p0
	frontelements[n1].v[0] = p0;
		
	return(1);	
}

/**
Add an element.
*/
int add_front_element(Front *fr, int elid){
#if _FTERMINAL
	assert(elid != fr->E[0] && elid != fr->E[1]);
#endif
	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;
	int p0, p1, p00, p11;
	Frontelement * frontelement = frontelements+elid;
	int n1;
	int newp, newe;

	//n0 = NbrI0(frontelement);
	n1 = NbrI1(frontelement);

	p0  = VertexI0(frontelement);
	p1  = VertexI1(frontelement);
	p00 = VertexI00(frontelement);
	p11 = VertexI11(frontelement);

	NOT_UNUSED(newp);
	NOT_UNUSED(newe);

	//Create Element newelem, and Point newpt
	add_element(fr, newe);
	add_point(fr, newp);

#if _fperiodicx
	double xp = frontpoints[p00].x[0] - 0.25*L0;
	frontpoints[newp].x[0] = ( - apply_fperiodic_cts( frontpoints[p00].x[0], xp) +
	                     9.0 * apply_fperiodic_cts( frontpoints[p0].x[0],  xp) + 
	                     9.0 * apply_fperiodic_cts( frontpoints[p1].x[0],  xp) - 
	                           apply_fperiodic_cts( frontpoints[p11].x[0], xp))/16.;
	frontpoints[newp].x[0] = apply_fperiodicx(frontpoints[newp].x[0]);
#else
	frontpoints[newp].x[0] = ( - frontpoints[p00].x[0] +  9.0 *(frontpoints[p0].x[0] + frontpoints[p1].x[0]) - frontpoints[p11].x[0] )/16.;
#endif
	frontpoints[newp].x[1] = ( - frontpoints[p00].x[1] +  9.0 *(frontpoints[p0].x[1] + frontpoints[p1].x[1]) - frontpoints[p11].x[1] )/16.;

	//Resetting Corners and Nbrs
	frontelements[elid].n[1] = newe;
	frontelements[newe].n[0] = elid;
	frontelements[newe].n[1] = n1;
	frontelements[n1].n[0]   = newe;

	frontelements[elid].v[1] = newp;
	frontelements[newe].v[0] = newp;
	frontelements[newe].v[1] = p1;
	return(1);		
}

#if _FTERMINAL
void regrid_front_terminals(Front * fr) {
	front_ghost_boundary(fr);

	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;
	int p0, p1, pm, p00, p11;
	int n11;
	double a, b, r;
	int elem;
	Frontelement * frontelement;

	for(int i=0; i<2; ++i) {
		int j = !i;
		elem  = fr->E[i];
		frontelement =  frontelements+elem;
		assert(frontelement->n[i]  != -1);
		p00 = frontelements[frontelement->n[i]].v[i];
		p0  = frontelement->v[i];
		pm  = frontelement->v[j];
		p1  = frontelements[frontelement->n[j]].v[j];
		n11 = frontelements[frontelement->n[j]].n[j];
		p11 = frontelements[n11].v[j];
		a =  elem_length(fr, elem);
		b =  elem_length(fr, frontelement->n[j]);
		r =  b/a;
		if (r < 1.) {
			if (a < frregrid->amax) continue;
			frontpoints[pm].x[0] = polynomial_p4_eq(frontpoints[p00].x[0], frontpoints[p0].x[0],
		                           frontpoints[pm].x[0], frontpoints[p1].x[0], 1.+0.5*(1.+r));
			frontpoints[pm].x[1] = polynomial_p4_eq(frontpoints[p00].x[1], frontpoints[p0].x[1],
		                           frontpoints[pm].x[1], frontpoints[p1].x[1], 1.+0.5*(1.+r));
		}
		else {
			if (a > frregrid->amin) continue;
			frontpoints[pm].x[0] = polynomial_p4_eq(frontpoints[p0].x[0], frontpoints[pm].x[0],
		                           frontpoints[p1].x[0], frontpoints[p11].x[0], 0.5*(1.+r));
			frontpoints[pm].x[1] = polynomial_p4_eq(frontpoints[p0].x[1], frontpoints[pm].x[1],
		                           frontpoints[p1].x[1], frontpoints[p11].x[1], 0.5*(1.+r));
		}
	}

	undo_front_ghost_boundary(fr);
}	
#endif

void (* regrid_front_scheme)(Front *) = NULL;
void regrid_front(Front *fr){
	//fprintf(stdout, "\n[");
	int *ConnectPrev,
			*NE, *FirstElement,
	    elem, ielem;
	double  elength;

	ConnectPrev = fr->EConnectPrev;
	NE = fr->NE;
	FirstElement = fr->FirstElement;	

	for(int ipass = 0; ipass<6; ++ipass){
#if _FTERMINAL
		regrid_front_terminals(fr);
#endif

		//Adding elements if needed
		ielem = 0;
		elem = *FirstElement;
		while(ielem < *NE) {
//fixme:terminal elements cannot split/merge
#if _FTERMINAL
			if(elem == fr->E[0] || elem == fr->E[1]) {
				ielem++;
				elem = ConnectPrev[elem];
				continue;
			}
#endif
			elength = elem_length(fr, elem);
		  if( elength > frregrid->amax ) { 
				add_front_element(fr, elem);
				//if an element added, start the loop from the beginning
				ielem = 0;
				elem = *FirstElement;
				continue;					
			}
			ielem++;
			elem = ConnectPrev[elem];
		}

		//Deleting elements if needed
		ielem = 0;
		elem = *FirstElement;
		while(ielem < *NE) {
#if _FTERMINAL
			if(elem == fr->E[0] || elem == fr->E[1]) {
				ielem++;
				elem = ConnectPrev[elem];
				continue;
			}
#endif
			elength = elem_length(fr, elem);
		  if( elength < frregrid->amin ) { 
				delete_front_element(fr, elem);
						//After deleting, start the loop from the beginning
				ielem = 0;
				elem = *FirstElement;
				continue;									
			}
			ielem++;
			elem = ConnectPrev[elem];
		}
		
		if(front_quality(fr) == 0) {
			break;
		}

	}
	//fprintf(stdout, "]");
}


#elif dimension == 3 
//Fixme : Algorithm here is computationally very costly
int nflocal(Front *fr, int *m, int *n0, int *n1, int *n2,
            int *mc, int *nc0, int *nc1, int *nc2, int *p){
	Frontelement * fe = fr->frontelements;
	
	*n1 = (*n0 + 1)%3;
	*n2 = (*n0 + 2)%3;
	p[0] = fe[*m].v[*n0];  
	p[1] = fe[*m].v[*n1];  
	p[2] = fe[*m].v[*n2];  

	*mc = fe[*m].n[*n0];
	for(int i=0 ; i<3 ; ++i){
		if(fe[*mc].v[i] == p[0]) 
			*nc0= i;
	}
	*nc1 = (*nc0 + 1)%3;
	*nc2 = (*nc0 + 2)%3;
	p[3] = fe[*mc].v[*nc1];

	int vertex;
	for (int i=0 ; i<3; ++i) {
		vertex = fe[fe[*m].n[*n1]].v[i];
		if( vertex != p[1] && vertex != p[2])		p[4] = vertex; 
		vertex = fe[fe[*m].n[*n2]].v[i];
		if( vertex != p[2] && vertex != p[0])		p[5] = vertex; 	
	}

	return(1);
}

/**
Delete an element (And also one of its neighbour).
*/
int delete_front_element(Front *fr, int elem_id, int vertex_id){
	Frontelement * fe = fr->frontelements;
	Frontpoint * fp = fr->frontpoints;
	int m, n0, n1, n2, 
	    mc, nc0, nc1, nc2,
	    m1, m2, mc0, mc1, p[8],
	    mtemp, ntemp,
	    vertex;

	m = elem_id;
	n0 = vertex_id;

	nflocal(fr, &m, &n0, &n1, &n2,
	         &mc, &nc0, &nc1, &nc2, p);
	
	m1  = fe[m].n[n1];
	m2  = fe[m].n[n2];
	mc0 = fe[mc].n[nc0];
	mc1 = fe[mc].n[nc1];

	for(int i=0; i<3; ++i){
		vertex = fe[mc0].v[i];
		if(vertex != p[0] && vertex != p[3]) p[6] = vertex;
		vertex = fe[mc1].v[i];
		if(vertex != p[3] && vertex != p[1]) p[7] = vertex;
	}
	
	//Conditions to abandon deleting.
	for(int i=0; i<3; ++i){
		for (int j=0; j<3; ++j){
			if( fe[m1].n[i] == fe[m2].n[j] && fe[m1].n[i] != m) {
				//i.e if m1 and m2 have a common neighbour other than m
				return(0); 
			}
			if( fe[mc0].n[i] == fe[mc1].n[j] && fe[mc0].n[i] != mc) {
				//i.e if mc0 and mc1 have a common neighbour other than mc
				return(0);
			}
		}
	}


	for(int i=0; i<3; ++i){
		fp[p[0]].x[i] = 0.50 *( fp[p[0]].x[i] + fp[p[1]].x[i]) -
		                 0.125 *( fp[p[2]].x[i] + fp[p[3]].x[i]) +
		                 0.0625 *( fp[p[4]].x[i] + fp[p[5]].x[i] + fp[p[6]].x[i] + fp[p[7]].x[i]);
		
	}
	//Delete Elements m,mc and Point p[1]
	delete_element(fr, m);
	delete_element(fr, mc);
	delete_point(fr, p[1]);

	//Reset corner of all elements having P[1] as a corner from P[1] to P[0]
	mtemp = m1;
	ntemp = 0;
	bool found;
	int stop = 0;
	bool print = false;
	while(++stop){
		for(int i=0; i<3; ++i){
			if(fe[mtemp].v[i] == p[1]) {
				ntemp = i;
			}
		}
		if(fe[mtemp].n[ntemp] == mc) break;
		mtemp = fe[mtemp].n[ntemp]; 
		if(stop == 100) {
			print = true;
			break;
		}
	}

if(print) {

fprintf(stdout,"\nm");
fprintf(stdout,"%d	%d-%d %d-%d %d-%d", m,  fe[m].v[0], fe[m].n[0], fe[m].v[1], fe[m].n[1], fe[m].v[2], fe[m].n[2]); 
fprintf(stdout,"\nm1 ");
fprintf(stdout,"%d	%d-%d %d-%d %d-%d", m1, fe[m1].v[0], fe[m1].n[0], fe[m1].v[1], fe[m1].n[1], fe[m1].v[2], fe[m1].n[2]); 
fprintf(stdout,"\nm2 ");
fprintf(stdout,"%d 	%d-%d %d-%d %d-%d", m2, fe[m2].v[0], fe[m2].n[0], fe[m2].v[1], fe[m2].n[1], fe[m2].v[2], fe[m2].n[2]); 
fprintf(stdout,"\nmc ");
fprintf(stdout,"%d	%d-%d %d-%d %d-%d", mc, fe[mc].v[0], fe[mc].n[0], fe[mc].v[1], fe[mc].n[1], fe[mc].v[2], fe[mc].n[2]); 
fprintf(stdout,"\nmc0 ");
fprintf(stdout,"%d	%d-%d %d-%d %d-%d", mc0, fe[mc0].v[0], fe[mc0].n[0], fe[mc0].v[1], fe[mc0].n[1], fe[mc0].v[2], fe[mc0].n[2]); 
fprintf(stdout,"\nmc1 ");
fprintf(stdout,"%d	%d-%d %d-%d %d-%d", mc1, fe[mc1].v[0], fe[mc1].n[0], fe[mc1].v[1], fe[mc1].n[1], fe[mc1].v[2], fe[mc1].n[2]); 
fprintf(stdout,"\nLOOP mc%d p[1]%d", mc, p[1]);

if(print) fprintf(stdout,"\nLOOP mc%d p[1]%d", mc, p[1]);
}
stop = 0;
mtemp = m1;
	while(++stop){
		found = 0;

if(print)	fprintf(stdout,"\n\t%d(%d-%d %d-%d %d-%d)", mtemp,  fe[mtemp].v[0], fe[mtemp].n[0], fe[mtemp].v[1], fe[mtemp].n[1], fe[mtemp].v[2], fe[mtemp].n[2]); 
		for(int i=0; i<3; ++i){
			if(fe[mtemp].v[i] == p[1]) {
				fe[mtemp].v[i] = p[0];
				ntemp = i;
				found = true;
			}
		}
		assert(found);// To avoid infinite looping
		if(fe[mtemp].n[ntemp] == mc) break;
		mtemp = fe[mtemp].n[ntemp]; 
		assert(stop<100);
	}

	//Resetting neighbour elements around deleed elements
	for(int i=0; i<3; ++i){
		if(fe[m1].n[i] == m) fe[m1].n[i] = m2;
		if(fe[m2].n[i] == m) fe[m2].n[i] = m1;
		if(fe[mc0].n[i] == mc) fe[mc0].n[i] = mc1;
		if(fe[mc1].n[i] == mc) fe[mc1].n[i] = mc0;
	}
	
	if(front_connectivity_valid(fr)==0){
fprintf(stdout,"\nm %d mc %d m1 %d m2 %d, mc0 %d, mc1 %d",m,mc,m1,m2,mc0,mc1);
fprintf(stdout,"\np");
for(int i=0;i<8;++i)
fprintf(stdout," %d", p[i]);

fprintf(stdout,"\nm");
fprintf(stdout,"%d	%d-%d %d-%d %d-%d", m,  fe[m].v[0], fe[m].n[0], fe[m].v[1], fe[m].n[1], fe[m].v[2], fe[m].n[2]); 
fprintf(stdout,"\nm1 ");
fprintf(stdout,"%d	%d-%d %d-%d %d-%d", m1, fe[m1].v[0], fe[m1].n[0], fe[m1].v[1], fe[m1].n[1], fe[m1].v[2], fe[m1].n[2]); 
fprintf(stdout,"\nm2 ");
fprintf(stdout,"%d 	%d-%d %d-%d %d-%d", m2, fe[m2].v[0], fe[m2].n[0], fe[m2].v[1], fe[m2].n[1], fe[m2].v[2], fe[m2].n[2]); 
fprintf(stdout,"\nmc ");
fprintf(stdout,"%d	%d-%d %d-%d %d-%d", mc, fe[mc].v[0], fe[mc].n[0], fe[mc].v[1], fe[mc].n[1], fe[mc].v[2], fe[mc].n[2]); 
fprintf(stdout,"\nmc0 ");
fprintf(stdout,"%d	%d-%d %d-%d %d-%d", mc0, fe[mc0].v[0], fe[mc0].n[0], fe[mc0].v[1], fe[mc0].n[1], fe[mc0].v[2], fe[mc0].n[2]); 
fprintf(stdout,"\nmc1 ");
fprintf(stdout,"%d	%d-%d %d-%d %d-%d", mc1, fe[mc1].v[0], fe[mc1].n[0], fe[mc1].v[1], fe[mc1].n[1], fe[mc1].v[2], fe[mc1].n[2]); 
	}
	return(1);	//Delete successful		

}
/**
Add an element.
*/
int add_front_element(Front *fr, int elem_id, int vertex_id){
	Frontelement * fe = fr->frontelements;
	Frontpoint * fp = fr->frontpoints;
	int m, n0, n1, n2, 
	    mc, nc0, nc1, nc2,
	    m2, mc0, mc1, p[8],
	    newpt, newm, newmc,
	    vertex;

	m = elem_id;
	n0 = vertex_id;

	nflocal(fr, &m, &n0, &n1, &n2,
	         &mc, &nc0, &nc1, &nc2, p);
	
	//m1 = Nbr[m][n1];
	m2  = fe[m].n[n2];
	mc0 = fe[mc].n[nc0];
	mc1 = fe[mc].n[nc1];
	for(int i=0; i<3; ++i){
		vertex = fe[mc0].v[i];
		if(vertex != p[0] && vertex != p[3]) p[6] = vertex;
		vertex = fe[mc1].v[i];
		if(vertex != p[3] && vertex != p[1]) p[7] = vertex;
	}

	//Create Elements newm, newmc and Point newpt
	add_element(fr, newm);
	add_element(fr, newmc);
	add_point(fr, newpt);
	//Coordinates of new point
	for(int i=0; i<3; ++i){
		fp[newpt].x[i] = 0.50 *( fp[p[0]].x[i] + fp[p[1]].x[i]) -
		                 0.125 *( fp[p[2]].x[i] + fp[p[3]].x[i]) +
		                 0.0625 *( fp[p[4]].x[i] + fp[p[5]].x[i] + fp[p[6]].x[i] + fp[p[7]].x[i]);
	}

	//Resetting Corners and Nbrs
	fe[m].v[n0] = newpt;
	fe[mc].v[nc0] = newpt;
	
	fe[newm].v[0] = newpt;
	fe[newm].v[1] = p[2];
	fe[newm].v[2] = p[0];

	fe[newmc].v[0] = newpt;
	fe[newmc].v[1] = p[0];
	fe[newmc].v[2] = p[3];
	
	fe[m].n[n2] = newm;
	fe[mc].n[nc0] = newmc;

	fe[newm].n[0] = m;
	fe[newm].n[1] = m2;
	fe[newm].n[2] = newmc;
		
	fe[newmc].n[0] = newm;
	fe[newmc].n[1] = mc0;
	fe[newmc].n[2] = mc;

	for (int i=0; i<3; ++i){
		if(fe[m2].n[i] == m)    fe[m2].n[i] = newm;
		if(fe[mc0].n[i] == mc)    fe[mc0].n[i] = newmc;
	}

	//assert(front_connectivity_valid(fr));

	return(1);		
}
void regrid_front(Front *fr){
	int *ConnectPrev,	*NE, *FirstElement,
	    elem, ielem, p0, p1, p2, n0, n1, n2, ntemp;
	double  s0, s1, s2, stemp;

	ConnectPrev = fr->EConnectPrev;
	Frontpoint * fp = fr->frontpoints;
	Frontelement * fe = fr->frontelements;
	NE = fr->NE;
	FirstElement = fr->FirstElement;	

	for(int ipass = 0; ipass<6; ++ipass){

		//Adding elements if needed
		ielem = 0;
		elem = *FirstElement;
		while(ielem < *NE) {

			p0 = fe[elem].v[0];
			p1 = fe[elem].v[1];
			p2 = fe[elem].v[2];
			s0 = distance(fp[p0].x, fp[p1].x);
			s1 = distance(fp[p1].x, fp[p2].x);
			s2 = distance(fp[p2].x, fp[p0].x);
			n0 = 0; n1 = 1; n2 = 2;		

			//Sorting sidelengths in increasing order s0 <= s1 <= s2
			if(s0 > s1){
				stemp = s0;	s0 = s1; s1 = stemp;
				ntemp = n0;	n0 = n1; n1 = ntemp;
			}
			if(s1 > s2){
				stemp = s1;	s1 = s2; s2 = stemp;
				ntemp = n1;	n1 = n2; n2 = ntemp;
			}
			if(s0 > s1){
				stemp = s0;	s0 = s1; s1 = stemp;
				ntemp = n0;	n0 = n1; n1 = ntemp;
			}
				
			//fixme  : some conditions switched off. see next line which is commented
		  if( (s2 > frregrid->amax) || ((s0 > frregrid->amin ) && (aspect_ratio(fr, elem) > frregrid->aspmax) ) ) {
		  //if( s2 > frregrid->amax) {
	fprintf(stdout, "+");

					add_front_element(fr, elem, n2);

					//if an element added, start the loop from the beginning
					ielem = 0;
					elem = *FirstElement;
					continue;					
			}

			ielem++;
			elem = ConnectPrev[elem];
		}

		ielem = 0;
		elem = *FirstElement;
		while(ielem < *NE) {

			p0 = fe[elem].v[0];
			p1 = fe[elem].v[1];
			p2 = fe[elem].v[2];
			s0 = distance(fp[p0].x, fp[p1].x);
			s1 = distance(fp[p1].x, fp[p2].x);
			s2 = distance(fp[p2].x, fp[p0].x);
			n0 = 0; n1 = 1; n2 = 2;		

			//Sorting sidelengths in increasing order s0 <= s1 <= s2
			if(s0 > s1){
				stemp = s0;	s0 = s1; s1 = stemp;
				ntemp = n0;	n0 = n1; n1 = ntemp;
			}
			if(s1 > s2){
				stemp = s1;	s1 = s2; s2 = stemp;
				ntemp = n1;	n1 = n2; n2 = ntemp;
			}
			if(s0 > s1){
				stemp = s0;	s0 = s1; s1 = stemp;
				ntemp = n0;	n0 = n1; n1 = ntemp;
			}
				
			if( s0 < frregrid->amin) {
					if (delete_front_element(fr, elem, n0)) {
	fprintf(stdout, "-");
						ielem = 0;
						elem = *FirstElement;
						continue;
					}
			}

			ielem++;
			elem = ConnectPrev[elem];
		}

		//if(front_quality(fr) == 0) 
		//	break;

	}
}
#endif
#endif // !_MPI


//---------------------------locate_modified
//returns rank of the cell containing point{x,y,z}
int locate_rank (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = { .level = l };
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + GHOSTS;
#if dimension >= 2
    point.j = (p.y - Y0)/L0*n + GHOSTS;
#endif
#if dimension >= 3
    point.k = (p.z - Z0)/L0*n + GHOSTS;
#endif
    if (point.i >= 0 && point.i < n + 2*GHOSTS
#if dimension >= 2
	&& point.j >= 0 && point.j < n + 2*GHOSTS
#endif
#if dimension >= 3
	&& point.k >= 0 && point.k < n + 2*GHOSTS
#endif
	) {
      if (allocated(0) && is_leaf(cell))
	return cell.pid;
    }
    else
      break;
  }
	assert(false);
}

/**
Functions required by 
1)event vof()
	1.1)Advance marker points
      name vof() is used because ..

*/

/**
Advancing front*/
void (* advance_markers_scheme) (vector U, double d_t) = NULL;

#if dimension == 2
void advance_markers(vector U, double d_t) {
	Front * fr = fr_circle;
	int     is, js;

	double xm, ym, um, vm, dx, dy;
	
	foreach_frontpoint(fr) {
		
		xm = frontpoint->x[0];
		ym = frontpoint->x[1];
		//Locate the cell that contain the marker point	
		Point point = locate(xm,ym);

		dx = (xm - x)/Delta; 	
		dy = (ym - y)/Delta; 	
		is = (dx > 0.) ;
		js = (dy > 0.) ;

		um = 0.0; vm = 0.0;
		for (int ii = 0; ii<4; ++ii) {
			for (int jj = 0; jj<4; ++jj) {
				um += w_ij(dx-ii-is+2, dy-jj-js+2) *U.x[ii+is-2,jj+js-2];
				vm += w_ij(dx-ii-is+2, dy-jj-js+2) *U.y[ii+is-2,jj+js-2];
			}	
		}
	
		//Update location of marker point	
		frontpoint->x[0] +=  um*d_t;
#if _FTERMINAL	
		if (ptid == fr->P[0] || ptid == fr->P[1]) continue;
#endif
		frontpoint->x[1] +=  vm*d_t;	

#if _fperiodicx
		frontpoint->x[0] = apply_fperiodicx(frontpoint->x[0]);
#endif
	}	
}
#elif dimension == 3 
void advance_markers(vector U, double d_t) {
	Front *fr = fr_sphere;
	int     is, js, ks;

	double xm, ym, zm, um, vm, wm, dx, dy, dz;
	
	//Only local points
	foreach_frontpoint(fr) {
		
		xm = (frontpoint->x)[0];
		ym = (frontpoint->x)[1];
		zm = (frontpoint->x)[2];
		//Locate the cell that contain the marker point	
		//fixme : computationally expensive
		Point point = locate(xm, ym, zm);

		dx = (xm - x)/Delta; 	
		dy = (ym - y)/Delta; 	
		dz = (zm - z)/Delta; 	
		is = (dx > 0.) ;
		js = (dy > 0.) ;
		ks = (dz > 0.) ;

		um = 0.0; vm = 0.0; wm = 0.0;
		//fixme : how bad is using local initialaizing in for loop?
		for (int ii = 0; ii<4; ++ii) {
			for (int jj = 0; jj<4; ++jj) {
				for (int kk = 0; kk<4; ++kk) {
					um += w_ijk(dx-ii-is+2, dy-jj-js+2, dz-kk-ks+2) *U.x[ii+is-2,jj+js-2, kk+ks-2];
					vm += w_ijk(dx-ii-is+2, dy-jj-js+2, dz-kk-ks+2) *U.y[ii+is-2,jj+js-2, kk+ks-2];
					wm += w_ijk(dx-ii-is+2, dy-jj-js+2, dz-kk-ks+2) *U.z[ii+is-2,jj+js-2, kk+ks-2];
				}
			}	
		}
	
		//Update location of marker point	
		(frontpoint->x)[0] +=  um*d_t;	
		(frontpoint->x)[1] +=  vm*d_t;	
		(frontpoint->x)[2] +=  wm*d_t;	

	}
	
	//fprintf(stdout, "\n	EVENT(%f) :Advance Front ---", t);
}
#endif

#if dimension == 2
fpscalar front_curvature(){
#if _MPI
	assert(false);//fixme: NOt implmented yet
#endif
#if _fperiodicx
	assert(false);//fixme: NOt implmented yet
#endif
#if _FTERMINAL
	assert(false);//fixme: NOt implmented yet
#endif
	fpscalar kappa = new_fpscalar();
	fpscalar length = new_fpscalar();
	Front * fr = get_front();
	int p0, p1, p00, p11;
	double v[2], v2;

	foreach_frontpoint_all(fr){
		frontpoint->T[0] = 0.;
		frontpoint->T[1] = 0.;
		kappa[ptid] = 0.;
		length[ptid] = 0.;
	}

	foreach_frontelement_all(fr){
		p0  = VertexI0(frontelement);
		p1  = VertexI1(frontelement);
		p00 = VertexI00(frontelement);
		p11 = VertexI11(frontelement);
		v[0] = 27. *(frontpoints[p1].x[0] - frontpoints[p0].x[0]) - 
		            (frontpoints[p11].x[0] - frontpoints[p00].x[0]);
		v[1] = 27. *(frontpoints[p1].x[1] - frontpoints[p0].x[1]) - 
		            (frontpoints[p11].x[1] - frontpoints[p00].x[1]);

		v2 = sqrt(v[0]*v[0] + v[1]*v[1]);
		v[0] /= v2;
		v[1] /= v2;

		frontpoints[p0].T[0] += v[0];
		frontpoints[p0].T[1] += v[1];
		frontpoints[p1].T[0] -= v[0];
		frontpoints[p1].T[1] -= v[1];

		length[p0] += distance(frontpoints[p0].x,frontpoints[p1].x)/2.;
		length[p1] += distance(frontpoints[p0].x,frontpoints[p1].x)/2.;
	}

	foreach_frontpoint(fr) {
/**
\kappa \Delta s \textbf{n}
*/
		double k_n[2] = {frontpoint->T[0], frontpoint->T[1]}; 
		kappa[ptid] = sqrt(k_n[0] *k_n[0] + k_n[1] *k_n[1] )/length[ptid];
		//fprintf(stdout, " (%g %g %g %g)",k_n[0], k_n[1], length[ptid], kappa[ptid]);
	}
	free_fpscalar(length);
	return kappa;
} 


void front_to_vof_alternative(Front * fr, scalar c, bool inverse){

	vector norm[], vol[];	
	scalar facets[];

	foreach() {
		facets[] = 0.;
		foreach_dimension() {
			norm.x[] = 0.;
			vol.x[] = 0;
		}
	}

	int p0, p1;
	double t[2]; 
	coord  a, b, o;

	foreach_frontelement(fr) {
		p0  = VertexI0(frontelement);
		p1  = VertexI1(frontelement);

		a.x = frontpoints[p0].x[0];	a.y = frontpoints[p0].x[1];
		b.x = frontpoints[p1].x[0];	b.y = frontpoints[p1].x[1];
	
		Point point = locate(a.x, a.y);
	
		for(int I_ = -1; I_ <= 1; ++I_) {
			o.x = x + I_*Delta;
			for(int J_ = -1; J_ <= 1; ++J_) {
				o.y = y + J_*Delta;
				if( clip_line (o, Delta, a, b, t)) {
					//assert( 0.<=t[0]<=1. && 0.<=t[1]<=1.); 
					facets[I_, J_] += 1.;
					norm.x[I_, J_] += normalx(a,b)/tangentn(a,b);
					norm.y[I_, J_] += normaly(a,b)/tangentn(a,b);
/**
	$Area = \frac{1}{2} b h$
	$h = 0.5 ( X(t_0) + X(t_1)) - x_{i-1/2, j}$
	$b = |Y(t_0) - Y(t_1)|	
*/
					vol.x[I_, J_] += (a.x + 0.5*(t[0] + t[1])*(b.x - a.x) - (x + (I_ - 0.5)*Delta) )*fabs((t[1] - t[0])*(b.y - a.y))/Delta/Delta;
					vol.y[I_, J_] += (a.y + 0.5*(t[0] + t[1])*(b.y - a.y) - (y + (J_ - 0.5)*Delta) )*fabs((t[1] - t[0])*(b.x - a.x))/Delta/Delta;
				}
			}
		}
		//fixme : An alternative algo can be impl.
		//instead of going through 3x3
	}

	double v, n;

	foreach() {
		if (facets[] >= 1.) {
			v = (fabs(norm.x[]) > fabs(norm.y[])) ? vol.x[]: vol.y[];
			n = (fabs(norm.x[]) > fabs(norm.y[])) ? norm.x[] : norm.y[]; assert(fabs(n) > 0.);
			
			c[] = ((n > 0) != inverse) ? (v) : (1. - v);
		}
		else {
			c[] = c[] > 0.5; 
		}
		assert(c[] >= 0. && c[] <= 1.);
	}
/**
	fixme : 
	1)Bug arises when a cell of l < Max_level have front passing through.
	Maybe avoided by interpolation(restriction) [fraction wise].
	But not sure(?)
	2)periodic and mpi can't be handled
	
	//dump(file = "dumpdirectory/volfrac.dmp", list = (scalar *){norm, vol, c, f, facets});
*/
	
}

#endif