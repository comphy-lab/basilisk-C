#if _MPI
#include "front-common.h"

#define FMPI_BUFF_SIZE 2000
#define FMPI_TAG_DEFAULT 1

//fixme : deletethis . This pgm intended for 4 procs
#define locate_point(x, y) (2*(x > 0.5) + (y>0.5))

//Some macros
//iterate through NL points
@def foreach_frontpoint_nl(fr) 
	foreach_frontpoint_all(fr)
		if(!local_point()) {
@
@def end_foreach_frontpoint_nl()
		}
	end_foreach_frontpoint_all()
@

//iterate through NL elems
@def foreach_frontelement_nl(fr) 
	foreach_frontelement_all(fr)
		if(!local_element()) {
@
@def end_foreach_frontelement_nl()
		}
	end_foreach_frontelement_all()
@

//iterate through NL points with break
@def foreach_frontpoint_nlb(fr) 
	foreach_frontpoint_all(fr)
		if(!local_point()) {
@
@def end_foreach_frontpoint_nlb()
		}
		else
			break;
	end_foreach_frontpoint_all()
@

//iterate through NL elems with break
@def foreach_frontelement_nlb(fr) 
	foreach_frontelement_all(fr)
		if(!local_element()) {
@
@def end_foreach_frontelement_nlb()
		}
		else
			break;
	end_foreach_frontelement_all()
@


#if FDEBUG == 1
#include "front-debug.h"
#endif

#define FDEBUGB 0
#if FDEBUGB == 1
#include "fmpi-debug.h"
event init(t=0) {
	init_fdebugb();
}
event end(t=end) {
	finalise_fdebugb();
}
#endif

//types of packet data that are send through MPI tunnel
typedef struct {
	double x[dimension];
}Fmpi_ppacket;

typedef struct {
	double c[dimension];
	gID    p[dimension];
	gID    e[dimension];
}Fmpi_epacket;

#if dimension == 2
typedef struct {
	double x[2];
	gID elid;
	gID n, v;
	int dir;   //fixme: only 2 possible values 0/1.'int' can be replaced by 'bool'
}Fmpi_ghostpacket;

void pack_ghostelement(Front * fr, Frontelement * frontelement, Fmpi_ghostpacket * packet) {
	assert(frontelement->pid != pid());

	Frontelement * frontelements = fr->frontelements; 
	Frontpoint * frontpoints = fr->frontpoints;

	packet->elid.alias = frontelement->alias;
	packet->elid.pid = frontelement->pid;
	for(int dir=0; dir <=1; ++dir){
		int e = frontelement->n[dir];
		if(e >= 0) {
			packet->dir = dir;
			assert(frontelement->e[dir].pid == pid());
			packet->n.pid = pid();	
			packet->n.alias = e;
			
			int p = frontelements[e].v[dir];
			assert(p > -1); //fixme: Not reqd
			packet->v.pid = frontpoints[p].pid;	
			packet->v.alias = frontpoints[p].alias;	
			packet->x[0] = frontpoints[p].x[0];	
			packet->x[1] = frontpoints[p].x[1];	

			//assert(frontelement->n[!dir] == -1);
			return;
//fixme: there are cases where both sides of NL elems are L
		}
	}
}


#if FDEBUG == 1
void print_ghostelement(Fmpi_ghostpacket * packet) {
	fprintf(fdbgout, "\n\t RECV");
	if(packet->dir == 0)
	fprintf(fdbgout, "	(%d %d)->",
					packet->n.pid,
					packet->n.alias);
	fprintf(fdbgout, " (%d %d)",
					packet->elid.pid,
					packet->elid.alias);
	if(packet->dir == 1)
	fprintf(fdbgout, " ->(%d %d)",
					packet->n.pid,
					packet->n.alias);
}
#endif

void unpack_ghostelement(Front * fr, Fmpi_ghostpacket * packet) {
	Frontelement * frontelements = fr->frontelements; 
	int elid = packet->elid.alias;
	bool dir = (bool) packet->dir;
	Frontelement * frontelement = frontelements + elid;
	if( frontelement->e[dir].pid == pid() ){
		elid = frontelement->n[dir];
		frontelement = &frontelements[elid];
	}
	assert(frontelement->pid == pid()); 
	assert(packet->elid.pid == pid()); 
	//assert(frontelement->e[dir].pid == packet->n.pid); //fixme: extra precaution
	frontelement->e[dir].alias = packet->n.alias;
	Frontpoint * frontpoints = fr->frontpoints;

	int e = frontelement->n[dir];
	frontelements[e].alias = packet->n.alias; 
	frontelements[e].p[dir].alias = packet->v.alias; 
	frontelements[e].p[dir].pid = packet->v.pid; 
	
	int p = frontelements[e].v[dir];
	assert(p > -1); //fixme: Not reqd
	frontpoints[p].pid = packet->v.pid;	
	frontpoints[p].alias =  packet->v.alias;	
	frontpoints[p].x[0] = packet->x[0];	
	frontpoints[p].x[1] = packet->x[1];	

}

#elif dimension == 3
typedef struct {
	double x[dimension];
	gID elid;
	gID n, v;
	int dir;  //0, 1 or 2 
}Fmpi_ghostpacket;

void pack_ghostelement(Front * fr, Frontelement * frontelement, Fmpi_ghostpacket * packet) {
	assert(frontelement->pid != pid());

	Frontelement * frontelements = fr->frontelements; 
	Frontpoint * frontpoints = fr->frontpoints;

	packet->elid.alias = frontelement->alias;
	packet->elid.pid = frontelement->pid;

	int index = -1;

	for(int dir=0; dir <=2; ++dir){
		int e = frontelement->n[dir];
		if(e >= 0) {
			packet->dir = dir;
			for(int i=0; i<3; ++i) {
				if(frontelements[e].n[i] == frontelement->alias)  {
					index = i;
					assert(frontelement->p[dir].alias == frontelements[e].v[(i+1)%3]);
					break;
				}
			}
			assert(index >= 0);
	
			assert(frontelement->e[dir].pid == pid());
			packet->n.pid = pid();	
			packet->n.alias = e;
			
			int p = frontelements[e].v[dir];
			assert(p > -1); //fixme: Not reqd
			packet->v.pid = frontpoints[p].pid;	
			packet->v.alias = frontpoints[p].alias;	
			packet->x[0] = frontpoints[p].x[0];	
			packet->x[1] = frontpoints[p].x[1];	

			//assert(frontelement->n[!dir] == -1);
			return;
//fixme: there are cases where both sides of NL elems are L
		}
	}
}

void unpack_ghostelement(Front * fr, Fmpi_ghostpacket * packet) {
	Frontelement * frontelements = fr->frontelements; 
	int elid = packet->elid.alias;
	bool dir = (bool) packet->dir;
	Frontelement * frontelement = frontelements + elid;
	if( frontelement->e[dir].pid == pid() ){
		elid = frontelement->n[dir];
		frontelement = &frontelements[elid];
	}
	assert(frontelement->pid == pid()); 
	assert(packet->elid.pid == pid()); 
	//assert(frontelement->e[dir].pid == packet->n.pid); //fixme: extra precaution
	frontelement->e[dir].alias = packet->n.alias;
	Frontpoint * frontpoints = fr->frontpoints;

	int e = frontelement->n[dir];
	frontelements[e].alias = packet->n.alias; 
	frontelements[e].p[dir].alias = packet->v.alias; 
	frontelements[e].p[dir].pid = packet->v.pid; 
	
	int p = frontelements[e].v[dir];
	assert(p > -1); //fixme: Not reqd
	frontpoints[p].pid = packet->v.pid;	
	frontpoints[p].alias =  packet->v.alias;	
	frontpoints[p].x[0] = packet->x[0];	
	frontpoints[p].x[1] = packet->x[1];	

}

#endif


void pack_frontpoint(Frontpoint * frontpoint, Fmpi_ppacket * packet){
	packet->x[0] = (frontpoint->x)[0]; 
	packet->x[1] = (frontpoint->x)[1];
#if dimension == 3 
	packet->x[2] = (frontpoint->x)[2]; 
#endif
}
	
void unpack_frontpoint(Frontpoint * frontpoint, Fmpi_ppacket * packet){
	(frontpoint->x)[0] = packet->x[0]; 
	(frontpoint->x)[1] = packet->x[1]; 
#if dimension == 3 
	(frontpoint->x)[2] = packet->x[2]; 
#endif
}

void pack_frontelement(double xc, double yc, 
#if dimension ==3
	                     double zc,
#endif
	                     Frontelement * frontelement, Fmpi_epacket * packet){
	packet->c[0] = xc;
	packet->c[1] = yc;
#if dimension == 3 
	packet->c[2] = zc;
#endif

	packet->p[0].pid = (frontelement->p)[0].pid;
	packet->p[0].alias = (frontelement->p)[0].alias;
	packet->p[1].pid = (frontelement->p)[1].pid;
	packet->p[1].alias = (frontelement->p)[1].alias; 
#if dimension == 3 
	packet->p[2].pid = (frontelement->p)[2].pid;
	packet->p[2].alias = (frontelement->p)[2].alias;
#endif

	packet->e[0].pid = (frontelement->e)[0].pid;
	packet->e[0].alias = (frontelement->e)[0].alias;
	packet->e[1].pid = (frontelement->e)[1].pid;
	packet->e[1].alias = (frontelement->e)[1].alias; 
#if dimension == 3 
	packet->e[2].pid = (frontelement->e)[2].pid;
	packet->e[2].alias = (frontelement->e)[2].alias;
#endif
}
	
void unpack_frontelement(Frontelement * frontelement, Fmpi_epacket * packet){
	(frontelement->p)[0].pid = packet->p[0].pid;
	(frontelement->p)[0].alias = packet->p[0].alias;
	(frontelement->v)[0] =-1;
	(frontelement->p)[1].pid = packet->p[1].pid;
	(frontelement->p)[1].alias = packet->p[1].alias;
	(frontelement->v)[1] =-1;
#if dimension == 3
	(frontelement->p)[2].pid = packet->p[2].pid;
	(frontelement->p)[2].alias = packet->p[2].alias; 
	(frontelement->v)[0] =-1;
#endif

	(frontelement->e)[0].pid = packet->e[0].pid;
	(frontelement->e)[0].alias = packet->e[0].alias;
	frontelement->n[0] =-1;
	(frontelement->e)[1].pid = packet->e[1].pid;
	(frontelement->e)[1].alias = packet->e[1].alias;
	frontelement->n[1] =-1;
#if dimension == 3
	(frontelement->e)[2].pid = packet->e[2].pid;
	(frontelement->e)[2].alias = packet->e[2].alias;
	frontelement->n[0] =-1;
#endif
}
//tunnel to communicate any type of front data bw procs
typedef struct {
	//proc to be communicatted
	MPI_Request request;
	int pid;
	int len;
	int id[FMPI_BUFF_SIZE];
	void * qbuff;
	void * rbuff;
}Fmpi_comm_sr;

typedef struct{
	int npid;
	int * mappid;
	Fmpi_comm_sr * s, * r;
}Fmpi_comm;

Fmpi_comm * create_Fmpi_comm(size_t qsize, size_t rsize) {
	Fmpi_comm * c = (Fmpi_comm *)malloc(sizeof(Fmpi_comm));

	assert(qsize%sizeof(int) == 0);
	assert(rsize%sizeof(int) == 0);
		
	MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
	RcvPid * rcvpid = (mpi->mpi_level).snd;
	int npid = rcvpid->npid;

	c->npid = npid;
	c->s = (Fmpi_comm_sr *)malloc(npid*sizeof(Fmpi_comm_sr));
	c->r = (Fmpi_comm_sr *)malloc(npid*sizeof(Fmpi_comm_sr));
	Fmpi_comm_sr *s, *r;
	for(int p=0; p<npid; ++p){
		s = &(c->s[p]);
		s->pid = (rcvpid->rcv[p]).pid;
		s->rbuff = (void *)malloc(FMPI_BUFF_SIZE*rsize);
		s->qbuff = (void *)malloc(FMPI_BUFF_SIZE*qsize);

		r = &(c->r[p]);
		r->pid = (rcvpid->rcv[p]).pid;
		r->rbuff = (void *)malloc(FMPI_BUFF_SIZE*rsize);
		r->qbuff = (void *)malloc(FMPI_BUFF_SIZE*qsize);
	}

	c->mappid = (int *)malloc(npe()*sizeof(int));
	for(int p=0; p<npe(); p++) 
		c->mappid[p] = -1;
	for(int p=0; p<npid; ++p)
		c->mappid[(rcvpid->rcv[p]).pid] = p;	

	return c;
}

void  empty_Fmpi_comm(Fmpi_comm * c){
	free(c->s->rbuff);
	free(c->r->rbuff);
	free(c->s->qbuff);
	free(c->r->qbuff);
	free(c->s);
	free(c->r);
	free(c->mappid);
	free(c);
}

//{
	Fmpi_comm * redistributepc;
	Fmpi_comm * redistributeec;
	Fmpi_comm * remapc;
	Fmpi_comm * connectivityc;
	Fmpi_comm * coordinatec;
//}
event init(t=0) {
	redistributepc = create_Fmpi_comm(sizeof(Fmpi_ppacket), sizeof(int));
	redistributeec = create_Fmpi_comm(sizeof(Fmpi_epacket), sizeof(int));
	remapc = create_Fmpi_comm(sizeof(int), sizeof(gID));
	connectivityc = create_Fmpi_comm(sizeof(int), sizeof(gID));
	coordinatec = create_Fmpi_comm(sizeof(int), sizeof(Fmpi_ppacket));
}
event end(t=end) {
	empty_Fmpi_comm(redistributepc);
	empty_Fmpi_comm(redistributeec);
	empty_Fmpi_comm(remapc);
	empty_Fmpi_comm(connectivityc);
	empty_Fmpi_comm(coordinatec);
}

void sendrecv_ncomm(Fmpi_comm * c){
	Fmpi_comm_sr *s, *r;
	int tag = FMPI_TAG_DEFAULT;

	for(int p = 0; p<c->npid; ++p){
		s = &(c->s[p]);
		MPI_Isend(&(s->len), 1, MPI_INT, s->pid, tag, MPI_COMM_WORLD, &(s->request)); 	
	}

	for(int p = 0; p<c->npid; ++p){
		r = &(c->r[p]);
		MPI_Irecv(&(r->len), 1, MPI_INT, r->pid, tag, MPI_COMM_WORLD, &(r->request)); 	
	}

	for(int p = 0; p<c->npid; ++p){
		r = &(c->r[p]);
		MPI_Wait(&(r->request), MPI_STATUS_IGNORE); 	
	}

	for(int p = 0; p<c->npid; ++p){
		s = &(c->s[p]);
		MPI_Wait(&(s->request), MPI_STATUS_IGNORE); 	
	}
}

void sendrecv_commlist(Fmpi_comm * c, size_t  qsize){

	Fmpi_comm_sr *s, *r;
	int tag = FMPI_TAG_DEFAULT;
	int size = qsize/sizeof(int);

	for(int p = 0; p<c->npid; ++p){
		s = &(c->s[p]);
		if(s->len) 
			MPI_Isend(s->qbuff, (s->len)*size, MPI_INT, s->pid, tag, MPI_COMM_WORLD, &(s->request));
	}

	for(int p = 0; p<c->npid; ++p){
		r = &(c->r[p]);
		if(r->len)
			MPI_Irecv(r->qbuff, (r->len)*size, MPI_INT, r->pid, tag, MPI_COMM_WORLD, &(r->request));
	}

	for(int p = 0; p<c->npid; ++p){
		r = &(c->r[p]);
		if(r->len)
			MPI_Wait(&(r->request), MPI_STATUS_IGNORE); 	
	}

	for(int p = 0; p<c->npid; ++p){
		s = &(c->s[p]);
		if(s->len)
			MPI_Wait(&(s->request), MPI_STATUS_IGNORE); 	
	}
}

void sendrecv_commlist_return(Fmpi_comm * c, size_t rsize){

	Fmpi_comm_sr *s, *r;
	int tag = FMPI_TAG_DEFAULT;
	int size =  rsize/sizeof(int);

	for(int p = 0; p<c->npid; ++p){
		s = &(c->r[p]);
		if(s->len) 
			MPI_Isend(s->rbuff, (s->len)*size, MPI_INT, s->pid, tag, MPI_COMM_WORLD, &(s->request)); 	
	}

	for(int p = 0; p<c->npid; ++p){
		r = &(c->s[p]);
		if(r->len) 
			MPI_Irecv(r->rbuff, (r->len)*size, MPI_INT, r->pid, tag, MPI_COMM_WORLD, &(r->request)); 	
	}

	for(int p = 0; p<c->npid; ++p){
		r = &(c->s[p]);
		if(r->len) 
			MPI_Wait(&(r->request), MPI_STATUS_IGNORE); 	
	}
	
	for(int p = 0; p<c->npid; ++p){
		s = &(c->r[p]);
		if(s->len)
			MPI_Wait(&(s->request), MPI_STATUS_IGNORE); 	
	}
}

#if dimension == 2 //fixme:reconnection NOT yet done in 3D
void reconnect_front(Front *fr, Fmpi_comm * comm){

	Fmpi_comm_sr * r;
	Fmpi_ghostpacket * packet;
	
	for(int p=0; p<comm->npid; ++p)
		comm->s[p].len = 0;
	
	foreach_frontelement_all(fr){
		if(local_element()) continue;
		int p = comm->mappid[frontelement->pid];
		r = &(comm->s[p]);
		assert(r->len < FMPI_BUFF_SIZE);
		packet = &(((Fmpi_ghostpacket *) r->qbuff)[r->len]);
		pack_ghostelement(fr, frontelement, packet);
		r->len++;
	}

	sendrecv_ncomm(comm);
	//for(int p=0; p<comm->npid; ++p)
		//comm->r[p].len = comm->s[p].len;

	sendrecv_commlist(comm, sizeof(Fmpi_ghostpacket));

	for(int p = 0; p<comm->npid; ++p){
		r = &(comm->r[p]);
		packet = (Fmpi_ghostpacket *) r->qbuff;
		for(int i=0; i<r->len; ++i) {
			unpack_ghostelement(fr, packet+i);		
		}
	}
	
}
#endif

int add_circle(Front *fr, double xc, double yc,
                double radin){

#if dimension == 3
	assert(false);
#endif

	set_regrid_props(-1);	
	double elength = 0.5*(frregrid->amax + frregrid->amin);
	int np = (int) floor(2*pi*radin/elength);

	if(np > MAX_POINTS || np < 10) {
		printf("np should be b/w 10 and MAX_POINTS");
		return(-1);
	}
#if FDEBUG == 1
	reopen_fdbgout();
	fprintf(fdbgout, "\nInit Circle in MPI");
	fprintf(fdbgout, "\nnp %d", np);	
	reopen_fdbgout();
#endif

	double dtheta = 2.0*pi/((double) np);

	reset_front(fr);
	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelements = fr->frontelements;
	int newp, newe;
	//fixme : This routine is very costly
	//But happens only once
	double x0, y0;
	bool * local_pt = (bool *)malloc(np*sizeof(bool));
	bool * local_el = (bool *)malloc(np*sizeof(bool));
	int * local_pid = (int *)malloc(np*sizeof(int));
	int * local_eid = (int *)malloc(np*sizeof(int));

	//add all local points
	for (int p=0; p<np; ++p) {

		x0 = xc + radin*cos(p*dtheta);
		y0 = yc + radin*sin(p*dtheta);
		
		if(locate(x0, y0).level>=0) {
			add_point(fr, newp);

			frontpoints[newp].x[0] = x0;
			frontpoints[newp].x[1] = y0;

			frontpoints[newp].pid = pid();
			frontpoints[newp].alias = newp;

			local_pt[p] = 1;
			local_pid[p] = newp;
		}
		else {
			local_pt[p] = 0;
			local_pid[p] = -1;
		}
	}
#if FDEBUG == 1
	reopen_fdbgout();
	fprintf(fdbgout, "\nLocal points");
	for (int p=0; p<np; ++p) 
		if(local_pt[p])
			fprintf(fdbgout, "\n\t%d)%d", p, local_pid[p]);
	reopen_fdbgout();
#endif

	//add all local elements
	for (int e=0; e<np; ++e) {
		x0 = xc + 0.5*radin*(cos(e*dtheta) + cos((e+1)*dtheta));
		y0 = yc + 0.5*radin*(sin(e*dtheta) + sin((e+1)*dtheta));

		if(locate(x0, y0).level>=0) {
			add_element(fr, newe);

			//corners & nbrs not yet fixed
			frontelements[newe].pid = pid();
			frontelements[newe].alias = newe;

			local_el[e] = 1;
			local_eid[e] = newe;
		}
		else {
			local_el[e] = 0;
			local_eid[e] = -1;
		}
	}
#if FDEBUG == 1
	reopen_fdbgout();
	fprintf(fdbgout, "\nLocal elems");
	for (int p=0; p<np; ++p) 
		if(local_el[p])
			fprintf(fdbgout, "\n\t%d)%d", p, local_eid[p]);
	reopen_fdbgout();
#endif

	for (int e=0; e<np; ++e) {
	//search for non-local nbrs elements
		if(local_el[e] == 0 && local_eid[e] == -1) {
			if(local_el[(e+1)%np] || local_el[(e-1+np)%np]) {
				add_element(fr, newe);

				frontelements[newe].pid = -1;

				local_eid[e] = newe;
			}
		}
	//search for non-local corner points
		if(local_pt[e] == 0 && local_pid[e] == -1) {
			if(local_el[e] || local_el[(e-1+np)%np]) {
				add_point(fr, newp);

				frontpoints[newp].pid = -1;

				local_pid[e] = newp;
			}
		}
	}

#if FDEBUG == 1
	reopen_fdbgout();
	fprintf(fdbgout, "\nNLE");
	for (int e=0; e<np; ++e) {
		if(!local_el[e] && local_eid[e] >=0 )
			fprintf(fdbgout, "\n\t%d) %d", e, local_eid[e]);
	}
	fprintf(fdbgout, "\nNLP");
	for (int p=0; p<np; ++p) {
		if(!local_pt[p] && local_pid[p] >=0 )
			fprintf(fdbgout, "\n\t%d) %d", p, local_pid[p]);
	}
	reopen_fdbgout();
#endif

	Fmpi_comm * query = create_Fmpi_comm(sizeof(int), sizeof(int));
	if(npe()>1) {
		query_frontelements(np, local_eid, fr, query);
		query_frontpoints(np, local_pid, fr, query);
	}
	empty_Fmpi_comm(query);

	//set neighbours and corners
	int P, E;
	for (int e=0; e<np; ++e) {
		if(local_el[e]) {
			assert(local_eid[e] >= 0);

			P = e;
			assert(local_pid[P] >= 0);
			frontelements[local_eid[e]].v[0]       = local_pid[P];
			frontelements[local_eid[e]].p[0].pid   = frontpoints[local_pid[P]].pid;
			frontelements[local_eid[e]].p[0].alias = frontpoints[local_pid[P]].alias;
		
			P = (e+1)%np;
			assert(local_pid[P] >= 0);
			frontelements[local_eid[e]].v[1]       = local_pid[P];
			frontelements[local_eid[e]].p[1].pid   = frontpoints[local_pid[P]].pid;
			frontelements[local_eid[e]].p[1].alias = frontpoints[local_pid[P]].alias;

			E = (e-1+np)%np;
			assert(local_eid[E] >= 0);
			if(local_el[E])
				frontelements[local_eid[e]].n[0] = local_eid[E];
			else
				frontelements[local_eid[e]].n[0] = -1;
			frontelements[local_eid[e]].e[0].pid   = frontelements[local_eid[E]].pid;
			frontelements[local_eid[e]].e[0].alias = frontelements[local_eid[E]].alias;

			E = (e+1+np)%np;
			assert(local_eid[E] >= 0);
			if(local_el[E])
				frontelements[local_eid[e]].n[1] = local_eid[E];
			else
				frontelements[local_eid[e]].n[1] = -1;
			frontelements[local_eid[e]].e[1].pid   = frontelements[local_eid[E]].pid;
			frontelements[local_eid[e]].e[1].alias = frontelements[local_eid[E]].alias;
		}
	}

	//delete_all NonLocal elements
	foreach_frontelement_all(fr) {
		if(local_element() == 0)
			delete_element(fr, elid);
	}

	free(local_pt);
	free(local_el);
	free(local_pid);
	free(local_eid);


	Fmpi_comm * coordinate = create_Fmpi_comm(sizeof(int), sizeof(Fmpi_ppacket));
	coordinate_frontpoints(fr, coordinate);
	empty_Fmpi_comm(coordinate);
	 
	return(0);

}

#if dimension == 3 
int add_sphere(Front *fr, double xc, double yc, double zc,
                double radin){

	int iip, ist, iie,
	    ia, ib, ic, id, **icp, **ine,
	    iqq, ne, np;
	bool * local_pt, * local_el; 
	int * local_pid, * local_eid;
	double _PI, dph, phi, theta, **pt;

	set_regrid_props(-1);
	double elength = 0.5*(frregrid->amax + frregrid->amin);
	int nps = floor(radin/elength * sqrt(2*pi/sqrt(3.)));
	
	printf("\nelength%g dmin%g nps%d", elength, L0/(1<<grid->maxdepth), nps);
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
	local_el  = (bool *)malloc((4*2*nps*nps)*sizeof(bool)); //
	local_eid  = (int *)malloc((4*2*nps*nps)*sizeof(int));//
  for(int l = 0; l< 8*nps*nps; ++l){
		icp[l] = (int *)malloc(3*sizeof(int));
		ine[l] = (int *)malloc(3*sizeof(int));
	}
	pt  = (double **)malloc((4*nps*nps+2)*sizeof(double *));
	local_pt  = (bool *)malloc((4*nps*nps+2)*sizeof(bool)); //
	local_pid  = (int *)malloc((4*nps*nps+2)*sizeof(int));//
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

	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;
	int nbr, corner, newp, newe; //
	double xC, yC, zC; //centroid of elem

	//add local points
	for(int p=0; p<np; ++p){
		if(locate(pt[p][0], pt[p][1], pt[p][2]).level >=0){
			add_point(fr, newp);

			frontpoints[newp].x[0] = pt[p][0];
			frontpoints[newp].x[1] = pt[p][1];
			frontpoints[newp].x[2] = pt[p][2];

			frontpoints[newp].pid = pid();
			frontpoints[newp].alias = newp;

			local_pt[p] = 1;
			local_pid[p] = newp;
		}
		else {
			local_pt[p] = 0;
			local_pid[p] = -1;
		}
	}
	
	int npt;
	MPI_Reduce(fr->NP, &npt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
	if(pid() == 0)
		assert(npt == np);

	for(int e=0; e<ne; ++e){
		xC = (pt[icp[e][0]][0] + pt[icp[e][1]][0] + pt[icp[e][2]][0])/3.;
		yC = (pt[icp[e][0]][1] + pt[icp[e][1]][1] + pt[icp[e][2]][1])/3.;
		zC = (pt[icp[e][0]][2] + pt[icp[e][1]][2] + pt[icp[e][2]][2])/3.;
		if(locate(xC, yC, zC).level >= 0) {
			add_element(fr, newe);

			//corners & nbrs not yet fixed
			frontelements[newe].pid = pid();
			frontelements[newe].alias = newe;

			local_el[e] = 1;
			local_eid[e] = newe;
		}
		else {
			local_el[e] = 0;
			local_eid[e] = -1;
		}
	}

	MPI_Reduce(fr->NE, &npt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
	if(pid() == 0)
		assert(npt == ne);

	for (int e=0; e<ne; ++e) {
//fprintf(stdout, "\n	e%d[%d|%d]", e, local_el[e], local_eid[e]);
		if(!local_el[e]) continue;
		for(int n=0; n<3; ++n){
			//search for non-local nbrs elements
			nbr = ine[e][n];
			if(local_eid[nbr] == -1) {
		 		add_element(fr, newe);
				frontelements[newe].pid = -1;
				local_eid[nbr] = newe;
			}
//fprintf(stdout, "	n%d[%d|%d]", nbr, local_el[nbr], local_eid[nbr]);
			//search for non-local corner points
			corner = icp[e][n];
			if(local_pid[corner] == -1) {
				add_point(fr, newp);
				frontpoints[newp].pid = -1;
				local_pid[corner] = newp;

				frontpoints[newp].x[0] = pt[corner][0];
				frontpoints[newp].x[1] = pt[corner][1];
				frontpoints[newp].x[2] = pt[corner][2];
			}
		}
	}

	if(npe()>1) {
		Fmpi_comm * query = create_Fmpi_comm(sizeof(int), sizeof(int));
		query_frontelements(ne, local_eid, fr, query);
		query_frontpoints(np, local_pid, fr, query);
		empty_Fmpi_comm(query);
	}

	for (int e=0; e<ne; ++e) {
		if(local_el[e]) {
//fprintf(stdout,"\ne%d", local_eid[e]);
			assert(local_eid[e] >= 0);
			
			for(int n=0; n<3; ++n) {
				//set neighbours
				nbr = ine[e][n];
				assert(local_eid[nbr] >= 0);
				if(local_el[nbr])
					frontelements[local_eid[e]].n[n] = local_eid[nbr];
				else
					frontelements[local_eid[e]].n[n] = -1;
				frontelements[local_eid[e]].e[n].pid = frontelements[local_eid[nbr]].pid;
				frontelements[local_eid[e]].e[n].alias = frontelements[local_eid[nbr]].alias;
//fprintf(stdout,"	n/%d/%d/%d/%d", nbr, frontelements[local_eid[e]].n[n], frontelements[local_eid[e]].e[n].pid, frontelements[local_eid[e]].e[n].alias);
	
				//set corners
				corner = icp[e][n];
				assert(local_pid[corner] >= 0);
				frontelements[local_eid[e]].v[n]  = local_pid[corner];
				frontelements[local_eid[e]].p[n].pid = frontpoints[local_pid[corner]].pid;
				frontelements[local_eid[e]].p[n].alias = frontpoints[local_pid[corner]].alias;
			}


		}
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


	//assert(front_connectivity_valid(fr));

	return(0);

}
#endif 	//if dimension

//modified functions of tree-mpi.h

face vector mypid[];
	//fixme : This algorithm works if listf has length == 1

static void apply_bc_return (Rcv * rcv, scalar * list, scalar * listv,
		      vector * listf, int l, MPI_Status s)
{
	//fixme : This algorithm works if listf has length == 1
  double * b = rcv->buf;
  foreach_cache_level(rcv->halo[l], l) {
    for (scalar s in list)
      s[] += *b++;
    for (vector v in listf)
      foreach_dimension() {

				if(((int) mypid.x[]) != rcv->pid)  {
					v.x[] += *b;
					mypid.x[] = (double) rcv->pid;
				}

				b++;

				if (*b != nodata && allocated(1)) {
					if(((int) mypid.x[1]) != rcv->pid)  {
	  				v.x[1] += *b;//fixme : source of bug
						mypid.x[1] = (double) rcv->pid;
					}
				}

				b++;
      }
    for (scalar s in listv) {
      for (int i = 0; i <= 1; i++)
	for (int j = 0; j <= 1; j++)
#if dimension == 3
	  for (int k = 0; k <= 1; k++)
#endif
	    {
	      if (*b != nodata && allocated(i,j,k))
		s[i,j,k] = *b;
	      b++;
	    }
    }
  }
  size_t size = b - (double *) rcv->buf;
  free (rcv->buf);
  rcv->buf = NULL;

  int rlen;
  MPI_Get_count (&s, MPI_DOUBLE, &rlen);
  if (rlen != size) {
    fprintf (stderr,
	     "rlen (%d) != size (%ld), %d receiving from %d at level %d\n"
	     "Calling debug_mpi(NULL)...\n"
	     "Aborting...\n",
	     rlen, size, pid(), rcv->pid, l);
    fflush (stderr);
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -2);
  }
}


static void rcv_pid_receive_return (RcvPid * m, scalar * list, scalar * listv,
			     vector * listf, int l)
{
  if (m->npid == 0)
    return;

	//fixme : This algorithm works if listf has length == 1
	foreach_face()
		mypid.x[] = -1.0;
  
  int len = list_len (list) + 2*dimension*vectors_len (listf) +
    (1 << dimension)*list_len (listv);

  MPI_Request r[m->npid];
  Rcv * rrcv[m->npid]; // fixme: using NULL requests should be OK
  int nr = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = malloc (sizeof (double)*rcv->halo[l].n*len);
      MPI_Irecv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
			           BOUNDARY_TAG(l), MPI_COMM_WORLD, &r[nr]);
      rrcv[nr++] = rcv;
    }
  }

  /* non-blocking receives (does nothing when using blocking receives) */
  if (nr > 0) {
    int i;
    MPI_Status s;
    mpi_waitany (nr, r, &i, &s);
    while (i != MPI_UNDEFINED) {
      Rcv * rcv = rrcv[i];
      assert (l <= rcv->depth && rcv->halo[l].n > 0);
      assert (rcv->buf);
      apply_bc_return (rcv, list, listv, listf, l, s);
      mpi_waitany (nr, r, &i, &s);
    }
  }
}

static void rcv_pid_sync_reverse (SndRcv * m, scalar * list, int l)
{
  scalar * listr = NULL, * listv = NULL;
  vector * listf = NULL;
  for (scalar s in list)
    if (!is_constant(s)) {
      if (s.face)
	listf = vectors_add (listf, s.v);
      else if (s.restriction == restriction_vertex)
	listv = list_add (listv, s);
      else
	listr = list_add (listr, s);
    }
  rcv_pid_send (m->rcv, listr, listv, listf, l);
  rcv_pid_receive_return (m->snd, listr, listv, listf, l);
  rcv_pid_wait (m->rcv);
  free (listr);
  free (listf);
  free (listv);
}

//static void mpi_boundary_level_reverse (const Boundary * b, scalar * list, int l)
void mpi_boundary_level_reverse (const Boundary * b, scalar * list, int l)
{
/*#if FDEBUG == 1
	fprintf(fdbgout, "\nmpi_boundary_level_reverse() start");
	reopen_fdbgout();
#endif
  */
	MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync_reverse (&m->mpi_level, list, l);
/*#if FDEBUG == 1
	fprintf(fdbgout, "\nmpi_boundary_level_reverse() ends");
	reopen_fdbgout();
#endif
*/
}
/**
Regrid
*/

#if dimension == 2
int delete_front_element(Front *fr, int elid){

	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;
	Frontelement * frontelement = &(frontelements[elid]);

	if(frontelement->e[0].pid != pid())
		return 0;
	if(frontelement->e[1].pid != pid())
		return 0;

	fprintf(stdout, "-");
	int p0, p1, p00, p11;
	int n0, n1;

	n0 = NbrI0(frontelement);
	n1 = NbrI1(frontelement);

	p0  = VertexI0(frontelement);
	p1  = VertexI1(frontelement);
	p00 = VertexI00(frontelement);
	p11 = VertexI11(frontelement);
	
	frontpoints[p0].x[0] = ( - frontpoints[p00].x[0] +  9.0 *(frontpoints[p0].x[0] + frontpoints[p1].x[0]) - frontpoints[p11].x[0] )/16.;
	frontpoints[p0].x[1] = ( - frontpoints[p00].x[1] +  9.0 *(frontpoints[p0].x[1] + frontpoints[p1].x[1]) - frontpoints[p11].x[1] )/16.;


	//Delete Element elid and p1
	delete_obj(fr->EConnectNext, fr->EConnectPrev, &elid, fr->FirstElement,
	           fr->FirstEmptyElement, fr->NE, fr->NEEmpty, fr->EAllocated);
	delete_obj(fr->PConnectNext, fr->PConnectPrev, &p1, fr->FirstPoint,
	           fr->FirstEmptyPoint, fr->NP, fr->NPEmpty, fr->PAllocated );

	//Resetting neighbour elements around deleted elements
	frontelements[n0].n[1] = n1;
	frontelements[n0].e[1].alias = n1;
	frontelements[n1].n[0] = n0;
	frontelements[n1].e[0].alias = n0;
	
	//Reset vertices of all elements having p1 as a vertex from p1 to p0
	frontelements[n1].v[0] = p0;
	frontelements[n1].p[0].alias = frontpoints[p0].alias;
	frontelements[n1].p[0].pid = frontpoints[p0].pid;
		
	return(1);	
}

/**
Add an element.
*/
int add_front_element(Front *fr, int elid){


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
	add_local_element(fr, newe);
	add_local_point(fr, newp);

	frontpoints[newp].x[0] = ( - frontpoints[p00].x[0] +  9.0 *(frontpoints[p0].x[0] + frontpoints[p1].x[0]) - frontpoints[p11].x[0] )/16.;
	frontpoints[newp].x[1] = ( - frontpoints[p00].x[1] +  9.0 *(frontpoints[p0].x[1] + frontpoints[p1].x[1]) - frontpoints[p11].x[1] )/16.;

	//Resetting Corners and Nbrs
	frontelements[elid].n[1] = newe;
	frontelements[elid].e[1].alias = newe;
	frontelements[elid].e[1].pid = pid();

	frontelements[newe].n[0] = elid;
	frontelements[newe].e[0].alias = elid;
	frontelements[newe].e[0].pid = pid();

	frontelements[newe].n[1] = n1;
	frontelements[newe].e[1].alias = frontelements[n1].alias;
	frontelements[newe].e[1].pid = frontelements[n1].pid;

	frontelements[n1].n[0]   = newe;
	frontelements[n1].e[0].alias  = newe;
	//frontelements[n1].e[0].pid   = pid();

	frontelements[elid].v[1] = newp;
	frontelements[elid].p[1].alias = newp;
	frontelements[elid].p[1].pid = pid();

	frontelements[newe].v[0] = newp;
	frontelements[newe].p[0].alias = newp;
	frontelements[newe].p[0].pid = pid();

	frontelements[newe].v[1] = p1;
	frontelements[newe].p[1].alias = frontpoints[p1].alias;
	frontelements[newe].p[1].pid = frontpoints[p1].pid;

	return(1);		
}


void regrid_front(Front *fr){
	Fmpi_comm * coordinate = create_Fmpi_comm(sizeof(int), sizeof(Fmpi_ppacket));	
	coordinate_frontpoints(fr, coordinate); 
	empty_Fmpi_comm(coordinate);

	Fmpi_comm * reconn = create_Fmpi_comm(sizeof(Fmpi_ghostpacket), sizeof(int));
	//fprintf(stdout, "\n[");
	int *ConnectPrev,
			*NE, *FirstElement,
	    elem, ielem;
	double  elength;
	Frontelement * frontelements = fr->frontelements;

	ConnectPrev = fr->EConnectPrev;
	NE = fr->NE;
	FirstElement = fr->FirstElement;	


	for(int ipass = 0; ipass<6; ++ipass){

		//Adding elements if needed
		ielem = 0;
		elem = *FirstElement;
		while(ielem < *NE) {
//fixme:terminal elements cannot split/merge
			elength = elem_length(fr, elem);
		  if( elength > frregrid->amax && frontelements[elem].pid == pid()) { 
				add_front_element(fr, elem);
				//if an element added, start the loop from the beginning
				ielem = 0;
				elem = *FirstElement;
				continue;					
			}
			ielem++;
			elem = ConnectPrev[elem];
		}

/*
#if FDEBUG == 1
fprintf(fdbgout, " - ");
reopen_fdbgout();
#endif
		//Deleting elements if needed
		ielem = 0;
		elem = *FirstElement;
		while(ielem < *NE) {
			elength = elem_length(fr, elem);
		  if( elength < frregrid->amin && frontelements[elem].pid == pid()) { 
				if( delete_front_element(fr, elem) ){
deleteme2++;
						//After deleting, start the loop from the beginning
					ielem = 0;
					elem = *FirstElement;
					continue;
				}									
			}
			ielem++;
			elem = ConnectPrev[elem];
		}
if(deleteme2){
		if(!front_connectivity_valid(fr)) {
				fprintf(fdbgout, "\nFOUND_A_PROBLEM");
		}
				fprintf(fdbgout, "\ndelete_fr");
	foreach_frontelement_all(fr){
		fprintf(fdbgout, "\n\t<%d %d %d>L%d <%f %f>", frontelement->n[0], elid, frontelement->n[1],local_element(), frontpoints[frontelement->v[0]].x[0], frontpoints[frontelement->v[1]].x[0]);
	}
}
valid = front_connectivity_valid(fr);
MPI_Allreduce(&valid, &validg, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
reopen_fdbgout();
assert(validg);
*/

		reconnect_front(fr, reconn);
	
		if(front_quality(fr)) 
			break;

	}
	empty_Fmpi_comm(reconn);
}

#elif dimension == 3
@def local_valence(eid, indx) 
	{
		if(frontelements[eid].e[indx].pid != pid())
			return false; 
		int p = frontelements[eid].v[indx];
		int e = frontelements[eid].n[indx];
		int j = -1;
		while(1){
			for(int i=0; i<3; ++i){
				if(frontelements[e].v[i] == p) {
					j = i;
				}
			}
			if(j==-1)
				return false;
			if(frontelements[e].e[j].pid != pid())
				return false;
			if(frontelements[e].n[j] == eid) break;
			e = frontelements[e].n[j];
			j = -1; 
		}
	}
@
//Fixme : Algorithm here is computationally very costly
bool nflocal(Front *fr, int *m, int *n0, int *n1, int *n2,
            int *mc, int *nc0, int *nc1, int *nc2, int *p){
	Frontelement * frontelements = fr->frontelements;
	
	*n1 = (*n0 + 1)%3;
	*n2 = (*n0 + 2)%3;

	p[0] = frontelements[*m].v[*n0];  
	p[1] = frontelements[*m].v[*n1];  
	p[2] = frontelements[*m].v[*n2];  

	*mc = frontelements[*m].n[*n0];
	for(int i=0 ; i<3 ; ++i){
		if(frontelements[*mc].v[i] == p[0]) 
			*nc0= i;
	}
	*nc1 = (*nc0 + 1)%3;
	*nc2 = (*nc0 + 2)%3;
	p[3] = frontelements[*mc].v[*nc1];

	int vertex;
	for (int i=0 ; i<3; ++i) {
		vertex = frontelements[frontelements[*m].n[*n1]].v[i];
		if( vertex != p[1] && vertex != p[2])		p[4] = vertex; 
		vertex = frontelements[frontelements[*m].n[*n2]].v[i];
		if( vertex != p[2] && vertex != p[0])		p[5] = vertex; 	
	}

	return true;
}

/**
Delete an element (And also one of its neighbour).
*/
bool delete_front_element(Front *fr, int elem_id, int vertex_id){

	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;
	int m, n0, n1, n2, 
	    mc, nc0, nc1, nc2,
	    m1, m2, mc0, mc1, p[8],
	    mtemp, ntemp,
	    vertex;

	m = elem_id;
	n0 = vertex_id;

	local_valence(elem_id, n0);
	local_valence(elem_id, (n0+1)%3);

	nflocal(fr, &m, &n0, &n1, &n2,
	         &mc, &nc0, &nc1, &nc2, p);
	
	m1  = frontelements[m].n[n1];
	m2  = frontelements[m].n[n2];
	mc0 = frontelements[mc].n[nc0];
	mc1 = frontelements[mc].n[nc1];

	for(int i=0; i<3; ++i){
		vertex = frontelements[mc0].v[i];
		if(vertex != p[0] && vertex != p[3]) p[6] = vertex;
		vertex = frontelements[mc1].v[i];
		if(vertex != p[3] && vertex != p[1]) p[7] = vertex;
	}
	
	//conditions to abandon deleting.
	for(int i=0; i<3; ++i){
		for (int j=0; j<3; ++j){
			if( frontelements[m1].n[i] == frontelements[m2].n[j] && frontelements[m1].n[i] != m) {
				//i.e if m1 and m2 have a common neighbour other than m
				return false; 
			}
			if( frontelements[mc0].n[i] == frontelements[mc1].n[j] && frontelements[mc0].n[i] != mc) {
				//i.e if mc0 and mc1 have a common neighbour other than mc
				return false;
			}
		}
	}


	for(int i=0; i<3; ++i){
		frontpoints[p[0]].x[i] = 0.50 *( frontpoints[p[0]].x[i] + frontpoints[p[1]].x[i]) -
		                 0.125 *( frontpoints[p[2]].x[i] + frontpoints[p[3]].x[i]) +
		                 0.0625 *( frontpoints[p[4]].x[i] + frontpoints[p[5]].x[i] + frontpoints[p[6]].x[i] + frontpoints[p[7]].x[i]);
		
	}
	//delete elements m,mc and point p[1]
	delete_element(fr, m);
	delete_element(fr, mc);
	delete_point(fr, p[1]);

	//resetting corner
	mtemp = m1;
	while(1){
		for(int i=0; i<3; ++i){
			if(frontelements[mtemp].v[i] == p[1]) {
				frontelements[mtemp].v[i] = p[0];
				assert(frontelements[mtemp].p[i].pid == pid());
				frontelements[mtemp].p[i].alias = p[0];
				ntemp = i;
			}
		}
		if(frontelements[mtemp].n[ntemp] == mc) break;
		mtemp = frontelements[mtemp].n[ntemp]; 
	}

	//resetting neighbour elements around deleed elements
	for(int i=0; i<3; ++i){
		if(frontelements[m1].n[i] == m) {
			frontelements[m1].n[i] = m2;
			assert(frontelements[m1].e[i].pid == pid());
			frontelements[m1].e[i].alias = m2;
		}
		if(frontelements[m2].n[i] == m) {
			frontelements[m2].n[i] = m1;
			assert(frontelements[m2].e[i].pid == pid());
			frontelements[m2].e[i].alias = m1;
		}
		if(frontelements[mc0].n[i] == mc) {
			frontelements[mc0].n[i] = mc1;
			assert(frontelements[mc0].e[i].pid == pid());
			frontelements[mc0].e[i].alias = mc1;
		}
		if(frontelements[mc1].n[i] == mc) {
			frontelements[mc1].n[i] = mc0;
			assert(frontelements[mc1].e[i].pid == pid());
			frontelements[mc1].e[i].alias = mc0;
		}
	}

	return true;	//delete successful		

}
/**
Add an element.
*/
bool add_front_element(Front *fr, int elem_id, int vertex_id){

	Frontelement * frontelements = fr->frontelements;
	Frontpoint * frontpoints = fr->frontpoints;
	int m, n0, n1, n2, 
	    mc, nc0, nc1, nc2,
	    m2, mc0, mc1, p[8],
	    newpt, newm, newmc,
	    vertex;

	m = elem_id;
	n0 = vertex_id;

	local_valence(m, n0);
	local_valence(m, (n0+1)%3);

	nflocal(fr, &m, &n0, &n1, &n2,
	         &mc, &nc0, &nc1, &nc2, p);
	
	//m1 = Nbr[m][n1];
	m2  = frontelements[m].n[n2];
	mc0 = frontelements[mc].n[nc0];
	mc1 = frontelements[mc].n[nc1];
	for(int i=0; i<3; ++i){
		vertex = frontelements[mc0].v[i];
		if(vertex != p[0] && vertex != p[3]) p[6] = vertex;
		vertex = frontelements[mc1].v[i];
		if(vertex != p[3] && vertex != p[1]) p[7] = vertex;
	}

	//Create Elements newm, newmc and Point newpt
	add_local_element(fr, newm);
	add_local_element(fr, newmc);
	add_local_point(fr, newpt);
	//Coordinates of new point
	for(int i=0; i<3; ++i){
		frontpoints[newpt].x[i] = 0.50 *( frontpoints[p[0]].x[i] + frontpoints[p[1]].x[i]) -
		                 0.125 *( frontpoints[p[2]].x[i] + frontpoints[p[3]].x[i]) +
		                 0.0625 *( frontpoints[p[4]].x[i] + frontpoints[p[5]].x[i] + frontpoints[p[6]].x[i] + frontpoints[p[7]].x[i]);
	}

	//Resetting Corners and Nbrs
	for(int i=0;i<3; ++i) {
		frontelements[newm].e[i].pid = pid();
		frontelements[newm].p[i].pid = pid();
		frontelements[newmc].e[i].pid = pid();
		frontelements[newmc].p[i].pid = pid();
	}

	frontelements[m].v[n0] = newpt;
	frontelements[mc].v[nc0] = newpt;
	frontelements[m].p[n0].alias = newpt;
	frontelements[mc].p[nc0].alias = newpt;
	
	frontelements[newm].v[0] = newpt;
	frontelements[newm].v[1] = p[2];
	frontelements[newm].v[2] = p[0];
	frontelements[newm].p[0].alias = newpt;
	frontelements[newm].p[1].alias = p[2];
	frontelements[newm].p[2].alias = p[0];

	frontelements[newmc].v[0] = newpt;
	frontelements[newmc].v[1] = p[0];
	frontelements[newmc].v[2] = p[3];
	frontelements[newmc].p[0].alias = newpt;
	frontelements[newmc].p[1].alias = p[0];
	frontelements[newmc].p[2].alias = p[3];
	
	frontelements[m].n[n2] = newm;
	frontelements[mc].n[nc0] = newmc;
	frontelements[m].e[n2].alias = newm;
	frontelements[mc].e[nc0].alias = newmc;

	frontelements[newm].n[0] = m;
	frontelements[newm].n[1] = m2;
	frontelements[newm].n[2] = newmc;
	frontelements[newm].e[0].alias = m;
	frontelements[newm].e[1].alias = m2;
	frontelements[newm].e[2].alias = newmc;
		
	frontelements[newmc].n[0] = newm;
	frontelements[newmc].n[1] = mc0;
	frontelements[newmc].n[2] = mc;
	frontelements[newmc].e[0].alias = newm;
	frontelements[newmc].e[1].alias = mc0;
	frontelements[newmc].e[2].alias = mc;

	for (int i=0; i<3; ++i){
		if(frontelements[m2].n[i] == m)	{
			frontelements[m2].n[i] = newm;
			frontelements[m2].e[i].pid = pid();
			frontelements[m2].e[i].alias = newm;
		}
		if(frontelements[mc0].n[i] == mc) { 
			frontelements[mc0].n[i] = newmc;
			frontelements[mc0].e[i].alias = newmc;
			frontelements[mc0].e[i].pid = pid();
		}
	}

	return true;
	//assert(front_connectivity_valid(fr));

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

					if(add_front_element(fr, elem, n2)) {

					//if an element added, start the loop from the beginning
						ielem = 0;
						elem = *FirstElement;
						continue;					
					}
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
						ielem = 0;
						elem = *FirstElement;
						continue;
					}
					else {
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

#endif //_MPI