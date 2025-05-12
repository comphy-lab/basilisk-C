/**
[Return to my homepage](http://basilisk.fr/sandbox/nlemoine/README)

# Toolbox for handling discharge inflows in a 2D domain without masks (MPI compatible)

## Structures

The inflow procedure relies on the definition of an inlet segment (left bank to right bank), which is then "rasterized" on grid faces. We first define a couple of useful structures:

- the *BoundaryFace* structure, which is the basic element constituting the rasterized inlet segment
*/

typedef struct {
 coord Center; 
 coord InwardNormal;
 double Delta;
} BoundaryFace;

/**
- the *PID_controller* structure which implements a [proportional-integral-derivative](https://en.wikipedia.org/wiki/PID_controller) control loop. It will help regulate inflow discharge according to the difference between the desired setpoint and the "measured" discharge.*/

typedef struct {
  double Kp;
  double Ti;
  double Td;
  double e; // current error, e = -log(Q/Qsetpoint)
  double E; // integral of error
  double de_dt; // derivative 
  double fQ; // control variable, such that eta = eta_b(fQ*Q,inlet)
} PID_Controller;

/**
- the *Inlet* structure itself, basically composed of an array of *BoundaryFace* structures, and fitted with a PID controller.
*/

struct Inlet {
  coord Segment[2]; // LEFT bank to RIGHT bank
  coord vec_n; // inward normal
  coord vec_t;
  // Rasterization on grid faces ("staircase" version of the boundary):
  int nf;
  BoundaryFace * BFaces;
  PID_Controller Controller;
//  char Name[100];
  int ID;
  bool with_velocity;
  int nQ;
  double * Time;
  double * Q;
};

/**

![Inlet geometry](https://dropsu.sorbonne-universite.fr/s/5Sd324owXNs7jxE/download){width="500px"}

## Initialization functions

### Function to reset the PID controller
*/

int Reset_PID_Controller(PID_Controller * _C)
{
 _C->e = 0.;
 _C->E = 0.;
 _C->de_dt = 0.;
 _C->fQ = 1.0;

 return(0);
} 

/**
### Function to initialize the inlet structure

Most importantly, this function performs the rasterization of the inlet segment on the grid. This rasterization is shared by all processes.
*/
int Init_Inlet(struct Inlet * _inlet, coord LeftBank, coord RightBank)
{
  scalar pscaln[], pscalt[];
  double xA = LeftBank.x, yA = LeftBank.y, xB = RightBank.x, yB = RightBank.y;
  (_inlet->Segment)[0] = (coord){xA,yA};
  (_inlet->Segment)[1] = (coord){xB,yB};

  double norm = sqrt(sq(xB-xA) + sq(yB-yA));
  assert (norm > 0.);
  
  // tangent and normal unit vector to the inflow section 
  _inlet->vec_t = (coord){ (xB-xA)/norm , (yB-yA)/norm };
  _inlet->vec_n = (coord){ (yA-yB)/norm + 1.e-6 , (xB-xA)/norm - 1.5e-6};

  (void) Reset_PID_Controller(&(_inlet->Controller));

  /* Rasterize [AB] segment on grid faces */

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
    double facesize = Delta;
    
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
	sendFaces[nf].Delta = facesize;
        nf++;
      }
    } // end if
  } // end foreach_face(x)

  foreach_face(y)
  {
    xf = x;
    yf = y;
    jjg = pscaln[]>0. ? -jg : jg; // *inward* pointing normal
    double facesize = Delta;

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
	sendFaces[nf].Delta = facesize;
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

  _inlet->nQ = 0;
  _inlet->Time = NULL;
  _inlet->Q = NULL;

  return(0);
}

int Deallocate_Inlet(struct Inlet * _inlet)
{
   free(_inlet->BFaces);
   _inlet->nf = 0.;
   return(0);
}

/**
## Discharge functions

The functions provided in this section are essentially a copy of the functions available in [discharge.h](http://basilisk.fr/src/discharge.h): we seek the value $\eta_b$ of the water surface elevation upstream of the inlet section, which yields the desired inflow discharge through the inlet section. The only difference is in the *bflux* function: because we do not have a true boundary condition on $z_b$ (typically Neumann), we have to perform the hydrostatic reconstruction of interface elevation for each boundary face composing the inlet section, just as in the Saint-Venant solver.
*/

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
  char fluxfile[200];
  bool write_fluxes = false ; // (p.ptr_inlet->ID)>0;
  FILE * fp ;
  if(write_fluxes){
   sprintf(fluxfile,"out/bflux/inlet_%d_t_%.3f_pid_%d.csv",p.ptr_inlet->ID,t,pid());
   fp = fopen(fluxfile,"w"); 
   fprintf(fp,"%s\n","x,y,ig,jg,eta_b,zb[out],h[in],dQ");
  }

// Center of inlet is taken as reference point for "lake-at-rest" correction with tilt
  double xcent = (p.ptr_inlet->Segment[0].x + p.ptr_inlet->Segment[1].x)/2.;
  double ycent = (p.ptr_inlet->Segment[0].y + p.ptr_inlet->Segment[1].y)/2.;

  for(int k = 0;k<(p.ptr_inlet->nf);k++)
  {
    int iig = (int) p.ptr_inlet->BFaces[k].InwardNormal.x;
    int jjg = (int) p.ptr_inlet->BFaces[k].InwardNormal.y;
    // locate "outer", pseudo-ghost cell (is in subdomain by construction)
    double xout = (p.ptr_inlet->BFaces[k].Center.x) - iig*L0/N/2;
    double yout = (p.ptr_inlet->BFaces[k].Center.y) - jjg*L0/N/2;
    Point point = locate(xout,yout);

    if(point.level>0){

 /** As the function can be used with [topography detrending](manning-tilt.h), we also correct $\eta_b$ in the inlet area as a function of pseudo-ghost cell coordinates $(x_\mathrm{out},y_\mathrm{out})$ so that the water surface elevation truly corresponds to a lake-at-rest condition:
   */   
      double corr_eta = tilt.x*(xout-xcent)+tilt.y*(yout-ycent);

/** Left / right states reconstruction (see [here](http://basilisk.fr/src/saint-venant.h#computing-fluxes) for details): */

      double dx = p.ptr_inlet->BFaces[k].Delta ;
      
      double up = iig*u.x[iig,0]+jjg*u.y[0,jjg]; // projection of inner velocity on face inward normal
      double upp = iig*u.x[2*iig,0]+jjg*u.y[0,2*jjg];
      double un = p.ptr_inlet->with_velocity ? iig*u.x[]+jjg*u.y[] : 0.;
      double unn = p.ptr_inlet->with_velocity ? iig*u.x[-iig,0]+jjg*u.y[0,-jjg] : 0.;

      double etap = eta[iig,jjg]; // inner cell
      double etapp = eta[2*iig,2*jjg];
      double etan = eta_b + corr_eta; // outer, pseudo-ghost cell
      double etann = eta_b + corr_eta - iig*dx*tilt.x - jjg*dx*tilt.y;

      double hp = h[iig,jjg];
      double hpp = h[2*iig,2*jjg];
      double hn = etan - zb[];
      double hnn = etann - zb[-iig,-jjg];

      // Estimate gradients

      double gr_eta_p = eta.gradient(etan,etap,etapp)/dx;
      double gr_eta_n = eta.gradient(etann,etan,etap)/dx;
      double gr_h_p = h.gradient(hn,hp,hpp)/dx;
      double gr_h_n = h.gradient(hnn,hn,hp)/dx;
      double gr_u_p = u.x.gradient(un,up,upp)/dx;
      double gr_u_n = u.x.gradient(unn,un,up)/dx;

      // Left/Right state reconstruction

      double zl = (etap-hp) - (dx/2.)*(gr_eta_p - gr_h_p);
      double zr = (etan-hn) + (dx/2.)*(gr_eta_n - gr_h_n);
      double zlr = max(zl, zr);
	
      double hl = hp - (dx/2.)*gr_h_p;
      double ul = up - (dx/2.)*gr_u_p;
      hl = max(0., hl + zl - zlr);
	
      double hr = hn + (dx/2.)*gr_h_n;
      double ur = un + (dx/2.)*gr_u_n;
      hr = max(0., hr + zr - zlr);

      hp = hl;
      up = ul;
      hn = hr;
      un = ur;

      if (hp > dry || hn > dry) {
        double fh, fu, dtmax;
        kurganov (hn, hp, un, up, 0., &fh, &fu, &dtmax);
        double dQ = fh*(p.ptr_inlet->BFaces[k].Delta); 
        Q += dQ;
        if(write_fluxes)
         fprintf(fp,"%g,%g,%d,%d,%g,%g,%g,%g\n", (p.ptr_inlet->BFaces[k].Center.x), (p.ptr_inlet->BFaces[k].Center.y),iig,jjg,eta_b,zb[],hp,dQ);
      }
    } // end if point.level>0
  } // end for

// Reduction

if(write_fluxes)
  fclose(fp);

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

int set_inlet_fields(struct Inlet * _inlet, double etain)
{
    double xcent = (_inlet->Segment[0].x + _inlet->Segment[1].x)/2.;
    double ycent = (_inlet->Segment[0].y + _inlet->Segment[1].y)/2.;

    foreach(){  // set eta to etain for all points in alcove
      double xrel = x-(_inlet->Segment)[0].x , yrel = y-(_inlet->Segment)[0].y;
      double norm = sqrt(sq(_inlet->Segment[0].x-_inlet->Segment[1].x)+sq(_inlet->Segment[0].y-_inlet->Segment[1].y)); 
      double psn = xrel*(_inlet->vec_n.x) + yrel*(_inlet->vec_n.y);
      double pst = xrel*(_inlet->vec_t.x) + yrel*(_inlet->vec_t.y);
      if(psn<0. && (psn>=-500.) && pst>=0 && pst<=norm)
      {
        double corr_eta = tilt.x*(x-xcent)+tilt.y*(y-ycent);
        h[] = fmax(0.,etain+corr_eta-zb[]);
        if( !(_inlet->with_velocity) | h[]<dry){ 
          u.x[] = 0.;
          u.y[] = 0.;
        }
      }
    } // end foreach()

    // Set pseudo-Neumann velocity on segment outer cells

    for(int k = 0;k<(_inlet->nf);k++)
    {
      int iig = (int) _inlet->BFaces[k].InwardNormal.x;
      int jjg = (int) _inlet->BFaces[k].InwardNormal.y;
      // locate "outer" cell (is in subdomain by construction)
      double xint = (_inlet->BFaces[k].Center.x) - iig*L0/N/2;
      double yint = (_inlet->BFaces[k].Center.y) - jjg*L0/N/2;
      Point point = locate(xint,yint);
      if(point.level>0){
        double corr_eta = tilt.x*(xint-xcent)+tilt.y*(yint-ycent);
        h[] = fmax(0.,etain+corr_eta-zb[]);
         if( !(_inlet->with_velocity) | h[]<dry){
          u.x[] = 0.;
          u.y[] = 0.;
         }
      }
    } // end for(int k...

    // end

    return(0);
}

int update_PID_error(double Qm, double Qsp, PID_Controller * _C)
{
   double newe = -log(Qm/Qsp); // setpoint-measurement, in log
   double de_dt = (newe-(_C->e))/dt; 
   double rho = exp(-dt/(_C->Ti));
   // Integral error with exp(t-t0) weight function (~Gauss-Laguerre quadrature)   
   _C->E = rho*(_C->E) + de_dt*dt + ((_C->e)-de_dt*(_C->Ti))*(1.-rho);

   double newcontrol =  (_C->Kp)*(newe
                                  + (_C->E)
                                  + de_dt*(_C->Td)
                                 );
   _C->e = newe;
   (_C->fQ) *= exp(newcontrol);
   
   return(0);
}

/**
## Managing inflow hydrographs

*/

/**
### Read hydrograph from file

This function has 3 arguments: a pointer to an *Inlet* structure, the path to the text file containing (time,Q) pairs, and an offset expressed in fractional days. In this implementation, times in the file are assumed to be expressed in fractional days as returned by the function *datenum()* in Matlab, Scilab, or Octave (e.g. 2022/01/01 14:30:00 is 738522.604167). As the file is parsed, times are converted to seconds elapsed from the offset *datenum0*.   

*/

int load_hydrograph(struct Inlet * _inlet, char * hydrofile, double datenum0)
{
  char buffer[200]; 
  int nQ;
  double time,Q;
  FILE * fp;

  if(!(fp = fopen (hydrofile, "r")))
  {
    printf("Failed to open hydrograph file!\n");
    return -1;
  }
  
  nQ = 0;

  // Read (time,Q) pairs
  while( fgets(buffer, 200, fp)!= NULL ) // parse until end-of-file
  {
  	sscanf(buffer,"%lf %lf",&time, &Q);
        double t_sec = (time-datenum0)*86400.;
        if(t_sec>=0)
	{
	  nQ++;
	  _inlet->Time = (double *) realloc( _inlet->Time, nQ*sizeof(double));
	  _inlet->Q    = (double *) realloc( _inlet->Q   , nQ*sizeof(double));
	  _inlet->Time[nQ-1] = t_sec;
	  _inlet->Q[nQ-1] = Q;
	}
  }

  _inlet->nQ = nQ;
  fclose(fp);
  printf("%d records loaded from %s\n",nQ,hydrofile);

  return(0);   
}

/**
### Linear interpolation in hydrograph
*/

double interpolate_discharge(struct Inlet * _inlet, double time)
{
  if(time <= _inlet->Time[0])
    return(_inlet->Q[0]);

  if(time >= _inlet->Time[(_inlet->nQ)-1])
    return(_inlet->Q[(_inlet->nQ)-1]);

  int pos = 0;
  while(time > (_inlet->Time[pos]) )
   pos++;

  double t0 = _inlet->Time[pos-1];
  double t1 = _inlet->Time[pos];
  double Q0 = _inlet->Q[pos-1];
  double Q1 = _inlet->Q[pos];

  return (Q0 + (time-t0)*(Q1-Q0)/(t1-t0));
}