typedef struct {
 coord Center; 
 coord InwardNormal;
 double Delta;
} BoundaryFace;

typedef struct {
  double Kp;
  double Ti;
  double Td;
  double e; // current error, e = log(Q/Qtarget)
  double E; // integral of error
  double de_dt; // derivative 
  double fQ; // control variable, such that eta = eta_b(fQ*Q,inlet)
} PID_Controller;

void Reset_PID_Controller(PID_Controller * _C)
{
 _C->e = 0.;
 _C->E = 0.;
 _C->de_dt = 0.;
 _C->fQ = 1.0;
} 

struct Inlet {
  coord Segment[2]; // LEFT bank to RIGHT bank
  coord vec_n; // inward normal
  coord vec_t;
  // Rasterization on grid faces ("staircase" version of the boundary):
  int nf;
  BoundaryFace * BFaces;
  PID_Controller Controller;
//  char Name[100];
//  int nQ;
//  double dtQ;
//  double * Q;
};
 
int Init_Inlet(struct Inlet * _inlet, coord LeftBank, coord RightBank)
{
  scalar pscaln[], pscalt[];
  double xA = LeftBank.x, yA = LeftBank.y, xB = RightBank.x, yB = RightBank.y;
  (_inlet->Segment)[0].x = xA;
  (_inlet->Segment)[0].y = yA;
  (_inlet->Segment)[1].x = xB;
  (_inlet->Segment)[1].y = yB;

  _inlet->vec_t.x = xB-xA;
  _inlet->vec_t.y = yB-yA;
  _inlet->vec_n.x = yA-yB;
  _inlet->vec_n.y = xB-xA;

  double norm = sqrt(sq(xB-xA) + sq(yB-yA));
  assert (norm > 0.);
  _inlet->vec_n.x = _inlet->vec_n.x/norm + 1e-6, _inlet->vec_n.y = _inlet->vec_n.y/norm - 1.5e-6;
  _inlet->vec_t.x = _inlet->vec_t.x/norm, _inlet->vec_t.y = _inlet->vec_t.y/norm;

  Reset_PID_Controller(&(_inlet->Controller));

  ////////////////// Rasterize [AB] segment on grid faces //////////////////////////

  // Compute scalar products for each cell center
  foreach(){
    pscaln[] = (x-xA)*(_inlet->vec_n.x) + (y-yA)*(_inlet->vec_n.y);
    if(pscaln[]==0.) pscaln[] = -(1e-6); // "on the boundary" means "outside"
    pscalt[] = (x-xA)*(_inlet->vec_t.x) + (y-yA)*(_inlet->vec_t.y);
  }       
  boundary({pscaln,pscalt});    

  double xf,yf;
  int nf = 0.;
  int iig,jjg;
  BoundaryFace * sendFaces = NULL;

  foreach_face(x)
  {
    xf = x;
    yf = y;
    iig = pscaln[]>0. ? -ig : ig; // *inward* pointing normal
    double Metric = Delta;
    
    // check whether face-sharing cells are on either side of the segment
    if( (pscaln[]*pscaln[-1,0]<0.) && (pscalt[]>=0.) && (pscalt[-1,0]>=0) && (pscalt[]<=norm) && (pscalt[-1,0]<=norm) )
    {
      // locate "outer" cell
      Point point = locate(xf-iig*L0/N/2.,yf);

      if(point.level>0) // keep face only if "outer", pseudo-ghost cell is in the sub-domain (avoids duplicates when using MPI)
      {
        sendFaces = (BoundaryFace *) realloc(sendFaces,(nf+1)*sizeof(BoundaryFace));
	sendFaces[nf].Center.x = xf;
	sendFaces[nf].Center.y = yf;
	sendFaces[nf].InwardNormal.x = (double) iig;
	sendFaces[nf].InwardNormal.y = 0.;
	sendFaces[nf].Delta = Metric;
        nf++;
      }
    } // end if
  } // end foreach_face(x)

  foreach_face(y)
  {
    xf = x;
    yf = y;
    jjg = pscaln[]>0. ? -jg : jg; // *inward* pointing normal
    double Metric = Delta;

    // check whether face-sharing cells are on either side of the segment
    if( (pscaln[]*pscaln[0,-1]<=0.) && (pscalt[]>=0.) && (pscalt[0,-1]>=0) && (pscalt[]<=norm) && (pscalt[0,-1]<=norm) )
    {
      // locate "outer" cell
      Point point = locate(xf,yf-jjg*L0/N/2.);

      if(point.level>0) // keep face only if "outer", pseudo-ghost cell is in the sub-domain (avoids duplicates when using MPI)
      {
        sendFaces = (BoundaryFace *) realloc(sendFaces,(nf+1)*sizeof(BoundaryFace));
	sendFaces[nf].Center.x = xf;
	sendFaces[nf].Center.y = yf;
	sendFaces[nf].InwardNormal.x = 0.;
	sendFaces[nf].InwardNormal.y = (double) jjg;
	sendFaces[nf].Delta = Metric;
        nf++;
      }
    } // end if
  } // end foreach_face(y)

  // We want all the processes to share the full rasterized boundary, so we gather the results
  @if _MPI
    int nfa[npe()], nfat[npe()];
    int nft;
    nfat[0] = 0;
    // First count how many faces where found by each process, put the results in array nfa
    MPI_Allgather (&nf, 1, MPI_INT, &nfa[0], 1, MPI_INT, MPI_COMM_WORLD);

    // Compute displacements
    for(int j=1;j<npe();j++)
     nfat[j] = nfat[j-1]+nfa[j-1];

    // Compute total number of faces
    nft = nfat[npe()-1]+nfa[npe()-1];  
    // Allocate receiver buffer in the structure
    _inlet->BFaces = malloc(sizeof(BoundaryFace)*nft);
    // Scale size of counts (nfa) & displacements (nfat) to have them in bytes 
    for(int j=0;j<npe();j++){
      nfa[j]  *= sizeof(BoundaryFace);
      nfat[j] *= sizeof(BoundaryFace);
    }
    // Gather
    MPI_Allgatherv (&sendFaces[0], nfa[pid()], MPI_BYTE,
                    &(_inlet->BFaces[0]),nfa,nfat,MPI_BYTE,
		    MPI_COMM_WORLD); 
    _inlet->nf = nft;
  @else
    _inlet->BFaces = malloc(sizeof(BoundaryFace)*nf);
    for(int j=0;j<nf;j++){
      _inlet->BFaces[j].Center.x = sendFaces[j].Center.x;
      _inlet->BFaces[j].Center.y = sendFaces[j].Center.y;
      _inlet->BFaces[j].InwardNormal.x = sendFaces[j].InwardNormal.x;
      _inlet->BFaces[j].InwardNormal.y = sendFaces[j].InwardNormal.y;
      _inlet->BFaces[j].Delta = sendFaces[j].Delta;
    }
    _inlet->nf = nf;
  @endif

  free(sendFaces);      
  return(0);
}

int Deallocate_Inlet(struct Inlet * _inlet)
{
   free(_inlet->BFaces);
   _inlet->nf = 0.;
   return(0);
}

/******************************************************************************************************************************/

struct Eta_b {
  // compulsory arguments
  double Q_b;
  struct Inlet * ptr_inlet;
  // optional arguments
  double prec; // precision (default 0.1%)
};

static double bflux (struct Eta_b p, double eta_b)
{
  double Q = 0.;

  for(int k = 0;k<(p.ptr_inlet->nf);k++)
  {
    int iig = (int) p.ptr_inlet->BFaces[k].InwardNormal.x;
    int jjg = (int) p.ptr_inlet->BFaces[k].InwardNormal.y;
    // locate "outer" cell (is in subdomain by construction)
    double xint = (p.ptr_inlet->BFaces[k].Center.x) - iig*L0/N/2;
    double yint = (p.ptr_inlet->BFaces[k].Center.y) - jjg*L0/N/2;
    Point point = locate(xint,yint);

    if(point.level>0){
      double hp = max(h[iig,jjg],0.);  // depth in inner cell
      double hn = max(eta_b-zb[],0.);  // depth in current (outer) cell

      // Use the velocity field as it is
      double up = iig*u.x[iig,jjg]+jjg*u.y[iig,jjg]; // projection of inner velocity on face inward normal 
      double un = iig*u.x[]+jjg*u.y[] ;

      if (hp > dry || hn > dry) {
        double fh, fu, dtmax;
        kurganov (hn, hp, un, up, 0., &fh, &fu, &dtmax);
        Q += fh*(p.ptr_inlet->BFaces[k].Delta);
      }
    } // end if point.level>0
  } // end for

// Reduction

@if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &Q, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
@endif

//  printf("Q(eta=%g)=%g\n",eta_b,Q);
  return Q;

}


static double falsepos (struct Eta_b p,
			double binf, double qinf,
			double bsup, double qsup)
{
  int n = 0;
  double newb, newq;
  qinf -= p.Q_b;
  qsup -= p.Q_b;
  do {
    newb = (binf*qsup - bsup*qinf)/(qsup - qinf);
    newq = bflux (p, newb) - p.Q_b;
    if (newq > 0.)
      bsup = newb, qsup = newq;
    else
      binf = newb, qinf = newq;
    n++;
  } while (fabs(newq/p.Q_b) > p.prec && n < 100);

  if (n >= 100)
    fprintf (stderr, "WARNING: eta_b(): convergence not reached\n");
  
  return newb;
}

double eta_b (struct Eta_b p)
{
  double zmin = HUGE, etas = 0., hs = 0.;

  for(int k = 0;k<(p.ptr_inlet->nf);k++)
  {
    int iig = (int) p.ptr_inlet->BFaces[k].InwardNormal.x;
    int jjg = (int) p.ptr_inlet->BFaces[k].InwardNormal.y;
    // locate "outer" cell
    double xint = (p.ptr_inlet->BFaces[k].Center.x) - iig*L0/N/2;
    double yint = (p.ptr_inlet->BFaces[k].Center.y) - jjg*L0/N/2;
    Point point = locate(xint,yint);
    if(point.level>0){
      if (zb[] < zmin)
	zmin = zb[];
      etas += (p.ptr_inlet->BFaces[k].Delta)*h[]*eta[];
      hs += (p.ptr_inlet->BFaces[k].Delta)*h[];
    }
  }

// Reduction

  @if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &zmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &etas, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &hs  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  @endif

  if (p.Q_b <= 0.)
    return zmin - 1.;

  // We try to find good bounds on the solution.

  double etasup = hs > 0. ? etas/hs : zmin;
  double Qsup = bflux (p, etasup), etainf = zmin, Qinf = 0.;
  double h0 = etasup - zmin;
  if (h0 < dry)
    h0 = 1.;
  int n = 0;
  while (Qsup < p.Q_b && n++ < 100) {
    etainf = etasup, Qinf = Qsup;
    etasup += h0;
    Qsup = bflux (p, etasup);
  }
  if (n >= 100)
    fprintf (stderr, "WARNING: eta_b() not converged\n");
    
  if (!p.prec) p.prec = 0.001; // 0.1% by default
  return falsepos (p, etainf, Qinf, etasup, Qsup);
}

double segment_flux (coord segment[2], scalar h, vector u)  // vient de "layered/hydro.h"
{
  coord m = {segment[0].y - segment[1].y, segment[1].x - segment[0].x};
  normalize (&m);
  double tot = 0.;

  foreach_segment (segment, p) {
    double dl = 0.;
    foreach_dimension() {
      double dp = (p[1].x - p[0].x)*Delta/Delta_x*(fm.y[] + fm.y[0,1])/2.;
      dl += sq(dp);
    }
    dl = sqrt (dl);    
    for (int i = 0; i < 2; i++) {
      coord a = p[i];
      tot += dl/2.*
	interpolate_linear (point, (struct _interpolate)
			    {h, a.x, a.y, 0.})*
	(m.x*interpolate_linear (point, (struct _interpolate)
				 {u.x, a.x, a.y, 0.}) +
	 m.y*interpolate_linear (point, (struct _interpolate)
				 {u.y, a.x, a.y, 0.}));
    }
  }
  // reduction
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return tot;
}