#if _MPI

#if dimension == 2
#define locate_elementrank() locate_rank(frontelement->c[0], frontelement->c[1])
#define locate_pointrank() locate_rank(frontpoint->x[0], frontpoint->x[1])
#define locate_coord(x) locate(x[0], x[1])
#elif dimension == 3
#define locate_elementrank() locate_rank(frontelement->c[0], frontelement->c[1], frontelement->c[2])
#define locate_pointrank() locate_rank(frontpoint->x[0], frontpoint->x[1], frontpoint->x[2])
#define locate_coord(x) locate(x[0], x[1], x[2])
#endif

typedef struct {
	MPI_Request sreq[2], rreq[2];
	bool reqstat;
	int pid;
	Array * q, * r;
	Array * temp;
}Fmpi_sr;

typedef struct{
	int npid;
	Fmpi_sr * sr;
	int * mappid;
}Fmpi;

@def foreach_fmpi_sr(f)
{
	Fmpi_sr * sr = f->sr; 	
	for(int p=0; p < f->npid; ++p, ++sr) {
@

@def end_foreach_fmpi_sr()
	}
}
@

Fmpi * fmpi_new() {
		
	MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
	RcvPid * rcvpid = (mpi->mpi_level).snd;

	int npid = rcvpid->npid;
	//assert(npid >0);
	Fmpi_sr * sr = (Fmpi_sr *)malloc(npid*sizeof(Fmpi_sr)), * s;
	s = sr;
	for(int p=0; p<npid; ++p, ++s) {
		s->pid = (rcvpid->rcv[p]).pid;
		s->q = array_new();
		s->r = array_new();
		s->temp = array_new();
	}

	//fixme : mappid slows down
	int * mappid = (int *)malloc(npe()*sizeof(int));
	for(int p=0; p<npe(); p++) 
		mappid[p] = -1;
	for(int p=0; p<npid; ++p)
		mappid[(rcvpid->rcv[p]).pid] = p;	

	Fmpi * f = (Fmpi *)malloc(sizeof(Fmpi));
	f->mappid = mappid;
	f->npid = npid;
	f->sr = sr;
	return f;
}

void fmpi_free(Fmpi * f){
	foreach_fmpi_sr(f){
		array_free(sr->q);
		array_free(sr->r);
		array_free(sr->temp);
	}
	free(f->sr);
	free(f->mappid);
	free(f);
}
/**
	Functions to send/rcv/wait for qbuff/rbuff
*/

#define FMPI_SR() 1003

typedef void (*fmpi_fn)(Array *, Array *);

void send_fmpi_general (Array * a, int to, MPI_Request * r)
{
  MPI_Isend (&a->len, 1, MPI_LONG, to, FMPI_SR(), MPI_COMM_WORLD, &r[0]);
  if (a->len > 0) 
    MPI_Isend (a->p, a->len, MPI_BYTE, to, FMPI_SR(), MPI_COMM_WORLD, &r[1]); 
}

long recv_only_fmpi_general (Array * a, int from, fmpi_fn fn)
{
	long n = 0;
  Array r;
  mpi_recv_check (&r.len, 1, MPI_LONG, from, FMPI_SR(),
		  MPI_COMM_WORLD, MPI_STATUS_IGNORE, "return_fmpi (len)");
  if (r.len > 0) {
    r.p = malloc (r.len);
    mpi_recv_check (r.p, r.len, MPI_BYTE, from, FMPI_SR(),
		    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "return_fmpi (p)");
      
		fn(&r, a); //a is updated. a = fn(&r)
		n += r.len;
    free (r.p);	
  }
	return n;
}

void recv_only_fmpi(Fmpi * f, fmpi_fn fn){
	foreach_fmpi_sr(f)
		recv_only_fmpi_general(sr->r, sr->pid, fn); 
}

long recv_return_fmpi_general (int from, MPI_Request *r, fmpi_fn fn)
{
  Array a, * s = array_new();
  mpi_recv_check (&a.len, 1, MPI_LONG, from, FMPI_SR(),
		  MPI_COMM_WORLD, MPI_STATUS_IGNORE, "recv_fmpi (len)");
  if (a.len > 0) {
    a.p = malloc (a.len);
    mpi_recv_check (a.p, a.len, MPI_BYTE, from, FMPI_SR(),
		    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "recv_fmpi (p)");
		fn(&a, s); //s is updated. s = fn(&a)	
    free (a.p);	
  }

  MPI_Isend (&s->len, 1, MPI_LONG, from, FMPI_SR(), MPI_COMM_WORLD, &r[0]);
  if (s->len > 0) 
    MPI_Isend (s->p, s->len, MPI_BYTE, from, FMPI_SR(), MPI_COMM_WORLD, &r[1]);
	long len = s->len;
	array_free(s);
	return len;
}

void recv_return_fmpi(Fmpi * f, fmpi_fn fn){
	foreach_fmpi_sr(f) {
		//sr->reqstat = false; 
		long len = recv_return_fmpi_general (sr->pid, sr->rreq, fn); 
		sr->reqstat = len > 0; 
	}
}

static void wait_fmpi_general (Array * a, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (a->len > 0)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

void wait_fmpi (Fmpi * f)
{
	foreach_fmpi_sr(f) 
		wait_fmpi_general (sr->q, sr->sreq);
}

static void wait2_fmpi_general (bool c, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (c)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

void wait2_fmpi (Fmpi * f)
{
	foreach_fmpi_sr(f) 
		wait2_fmpi_general (sr->reqstat, sr->sreq);
}

void reset_fmpi_arrays(Fmpi * f){
	foreach_fmpi_sr(f) {
		sr->q->len = 0; sr->r->len = 0;
		sr->temp->len = 0;
	}
}

void  fmpi(Fmpi * f, fmpi_fn fn1, fmpi_fn fn2)
{
	//send_fmpi (f);
	foreach_fmpi_sr(f)
		send_fmpi_general(sr->q, sr->pid, sr->sreq);

	//recv_return_fmpi(f, fn1);
	foreach_fmpi_sr(f) {
		long len = recv_return_fmpi_general (sr->pid, sr->rreq, fn1); 
		sr->reqstat = len > 0; 
	}
	
	//recv_only_fmpi(f, fn2);
	foreach_fmpi_sr(f)
		recv_only_fmpi_general(sr->r, sr->pid, fn2); 
}  

#endif