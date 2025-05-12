/**
[Return to my homepage](http://basilisk.fr/sandbox/nlemoine/README)


# Basilisk interface for computing the *Basal Envelope Surface of Talwegs* using the Analytical Element Method

This is the source code for the study published in [Le Moine, 2023](https://doi.org/10.5802/crgeos.164).

This code mixes Basilisk C code with a [FORTRAN library](./libaem.h) for computing [slit complex potential](./slit.f) and [areal sink potential](./areal.f), and finally [linear least-square LAPACK subroutines](./lsq_lapack.f). In order to run the code, you will need to [install Basilisk](http://basilisk.fr/src/INSTALL), get the content of the [envelope](./) directory, and then run:

make envelope.tst

*/

#include "grid/quadtree.h"
#include "libaem.h"
#include "run.h"
#include "nlemoine/esri-binary-grids.h"
#include "envelope.h"
#pragma autolink -L. -laem -lgfortran -lblas -llapack

#define M_PI acos(-1.0)
#define LEVEL 10

scalar Phi[];

/**
## The main initializes the grid (Evel catchment, Lambert 93 coord.)
*/

int main (int argc, char * argv[])
{  
  N = 1 << LEVEL;
  L0 = 25600.;

  size (L0);
  X0 = 267000.-L0/2.;
  Y0 = 6782500.-L0/2.;
  origin (X0,Y0);
  init_grid (1 << LEVEL);

  run();
  return(0);
}

event init (i=0)
{
 	double * ReZstart, * ImZstart, * Phistart, * ReZend, * ImZend, * Phiend , * PHI_topo, * Aup, * dA;
        int * MinUpID, * DownID;
	double * ReCn, * ImCn, * Q;
	double v0x,v0y,Phi0;
	double * ReZ, * ImZ, * PHI, * PSI;
	double * J_PHI, * J_PSI, * J_Vx, * J_Vy;
        double * SlitLength;

	int j,k;

	char buffer[200]; 
        char datafile[200];
	FILE * fp;

        int one = 1;
        double PHI1, PSI1, Vx1, Vy1;
        double Xp, Yp;

	int ncoef = 2;
	int ntheta = 4*ncoef;
	int nslit = 0;
        int nsink,nbound;

	ReZstart = NULL;
	ImZstart = NULL;
	Phistart = NULL;
	ReZend = NULL;
	ImZend = NULL;
	Phiend = NULL;
	MinUpID = NULL;
	DownID = NULL;
        SlitLength = NULL;
        Aup = NULL;
        dA = NULL;
/**
## Read data for stream slit elements
*/
	sprintf(datafile,"../sinks.dat");

	if(!(fp = fopen (datafile, "r")))
	{
	  printf("Failed to open sink data file!\n");
	  return -1;
	}
  
	while( fgets(buffer, 200, fp)!= NULL ) // parse until end-of-file
	{
	   ReZstart = (double *) realloc(ReZstart,(nslit+1)*sizeof(double));
	   ImZstart = (double *) realloc(ImZstart,(nslit+1)*sizeof(double));
	   Phistart = (double *) realloc(Phistart,(nslit+1)*sizeof(double));

	   ReZend = (double *) realloc(ReZend,(nslit+1)*sizeof(double));
	   ImZend = (double *) realloc(ImZend,(nslit+1)*sizeof(double));
	   Phiend = (double *) realloc(Phiend,(nslit+1)*sizeof(double));

	   Aup 	   = (double *) realloc(Aup,(nslit+1)*sizeof(double));
	   dA 	   = (double *) realloc(dA,(nslit+1)*sizeof(double));
	   MinUpID    = (int *) realloc(MinUpID,(nslit+1)*sizeof(int));
	   DownID     = (int *) realloc(DownID,(nslit+1)*sizeof(int));

	   sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %lf %d %d",ReZstart+nslit,ImZstart+nslit,Phistart+nslit,ReZend+nslit,ImZend+nslit,Phiend+nslit,
                                                         Aup+nslit,dA+nslit,MinUpID+nslit,DownID+nslit);
	   	
	   nslit++;
	}

        fclose(fp);

        fprintf(stderr,"%d sinks read\n",nslit);
	fflush(stderr);

/**
## Prune tree at higher support area if desired
*/
	// Prune tree

	double Aup_min = 0.0e6;

	(void) prune_tree ( 	& ReZstart, & ImZstart, & Phistart, 
		 		& ReZend, & ImZend, & Phiend,
		 		& Aup, & dA, & DownID, & nslit, Aup_min);

        fprintf(stderr,"Tree pruned down to %d sinks\n",nslit);
	fflush(stderr);

	SlitLength = (double *) malloc(nslit*sizeof(double));

	for(int s=0;s<nslit;s++)
	   SlitLength[s] = sqrt( sq(ReZend[s]-ReZstart[s]) + sq(ImZend[s]-ImZstart[s]) );	

        nsink = nslit;

/**
## Read boundary data
*/  
	sprintf(datafile,"../boundary.dat");

	if(!(fp = fopen (datafile, "r")))
	{
	  printf("Failed to open boundary data file!\n");
	  return -1;
	}

	while( fgets(buffer, 200, fp)!= NULL ) // parse until end-of-file
	{
	   ReZstart = (double *) realloc(ReZstart,(nslit+1)*sizeof(double));
	   ImZstart = (double *) realloc(ImZstart,(nslit+1)*sizeof(double));

	   ReZend = (double *) realloc(ReZend,(nslit+1)*sizeof(double));
	   ImZend = (double *) realloc(ImZend,(nslit+1)*sizeof(double));

	   sscanf(buffer,"%lf %lf %lf %lf",ReZstart+nslit,ImZstart+nslit,ReZend+nslit,ImZend+nslit);

	   SlitLength = (double *) realloc(SlitLength,(nslit+1)*sizeof(double));
	   SlitLength[nslit] = sqrt( sq(ReZend[nslit]-ReZstart[nslit]) + sq(ImZend[nslit]-ImZstart[nslit]) );	
           
	   nslit++;
	}

        fclose(fp);
        
        nbound = nslit-nsink; 
        fprintf(stderr,"%d boundary segments read\n",nbound);
	fflush(stderr);

        double * THETA = (double *) malloc( ntheta*sizeof(double) );

	for(int k=0;k<ntheta;k++)
	  THETA[k] = M_PI*(k+0.5)/ntheta;

/**
## Get Jacobian blocks for regression
*/
        double * UL = NULL, * UC = NULL, * UR = NULL,
 	       * LL = NULL, * LC = NULL, * LR = NULL;

        (void) get_LHS_blocks( nsink, nbound, ncoef, ntheta, & THETA, 
                               & ReZstart, & ImZstart, & ReZend, & ImZend,
                               & UL, & UC, & UR, & LL, & LC, & LR,
                               & ReZ, & ImZ);

        // Get potential & velocity created by areal sink of unit infitration rate (γ0 = -1)

        double gamma0 = -1.0;
  
        double * PHI_areal1 = (double *) malloc(nslit*ntheta*sizeof(double));
        double * Vx_areal1 = (double *) malloc(nslit*ntheta*sizeof(double));
        double * Vy_areal1 = (double *) malloc(nslit*ntheta*sizeof(double));
        double A;
        int nn = nslit*ntheta;

        (void) arealsink_ ( &ReZstart[nsink],&ImZstart[nsink],&nbound,   
                            ReZ,ImZ,&nn,&gamma0,
                            PHI_areal1,Vx_areal1,Vy_areal1,&A);

        // Compute array of normal velocity created by unit-rate areal sink

        double * Vn_areal1 = (double *) malloc(nbound*ntheta*sizeof(double)); 

        for(j=0;j<nbound;j++)
        {
           // inward normal
           double complex vv = I*(ReZend[nsink+j]-ReZstart[nsink+j] + I*ImZend[nsink+j]-I*ImZstart[nsink+j]);
           vv /= cabs(vv);
           double nx = creal(vv);
           double ny = cimag(vv);
	   for(k=0;k<ntheta;k++)
             Vn_areal1[ntheta*j+k] = nx*Vx_areal1[ntheta*(nsink+j)+k] + ny*Vy_areal1[ntheta*(nsink+j)+k];
        }     
 
        free(Vx_areal1);
        free(Vy_areal1); 

	PHI_topo = (double *) malloc(nsink*ntheta*sizeof(double));
	PHI      = (double *) malloc(nsink*ntheta*sizeof(double));

	for(j=0;j<nsink;j++){
	   for(k=0;k<ntheta;k++){
             int row = ntheta*j+k;
	     PHI_topo[row] = cos(THETA[k])*(Phiend[j]-Phistart[j])/2.+(Phiend[j]+Phistart[j])/2.;
 	   }
	}

        Q = (double *) malloc(nslit*sizeof(double));
        ReCn = (double *) malloc(nslit*ncoef*sizeof(double));
        ImCn = (double *) malloc(nslit*ncoef*sizeof(double));

/**
## Read  groundwater data
*/

	double * Xgw = NULL, * Ygw = NULL, * PHI_gw = NULL;
        int ngw = 0;

	sprintf(datafile,"../groundwater_feb.dat");

	if(!(fp = fopen (datafile, "r")))
	{
	  printf("Failed to open groundwater data file!\n");
	  return -1;
	}

	while( fgets(buffer, 200, fp)!= NULL ) // parse until end-of-file
	{
	   Xgw = (double *) realloc(Xgw,(ngw+1)*sizeof(double));
	   Ygw = (double *) realloc(Ygw,(ngw+1)*sizeof(double));
	   PHI_gw = (double *) realloc(PHI_gw,(ngw+1)*sizeof(double));

	   sscanf(buffer,"%lf %lf %lf",Xgw+ngw,Ygw+ngw,PHI_gw+ngw);
	   ngw++;
	}
	fclose(fp);     
        fprintf(stderr,"Groundwater data read for %d wells.\n",ngw); 

	double * PL, * PC, * PR, * PHI_areal1_gw;

	(void) get_LHS_blocks_p ( ngw, & Xgw, & Ygw, 
	                	  nsink, nbound, ncoef,
				  & ReZstart, & ImZstart, & ReZend, & ImZend,
				  & ReCn, & ImCn, & Q,
				  & PL, & PC, & PR, & PHI_areal1_gw );
/**
## If sink terms $Q_j$'s are enabled, allow estimate of areal exfiltration rate $\gamma_0$ with groundwater data
*/
        bool with_sinks = true;
  
	if(with_sinks)
	{
	    (void) solve_recharge_p ( nsink, nbound, ncoef, ntheta, ngw,
                       		      & UL, & UC, & UR,
                                      & LL, & LC, & LR,
                                      & PL, & PC, & PR,
                                      & PHI_areal1, & Vn_areal1, & PHI_areal1_gw,  
				      & ReZ, & ImZ, & Xgw, & Ygw,
				      & SlitLength, & PHI_topo, & PHI_gw,
                       		      & dA, & Aup,
				      & gamma0, & Q, & ReCn, & ImCn, & v0x, & v0y, & Phi0,
                       		      & PHI );  
	}
	else
/**
## Otherwise only estimate slit degrees of freedom for computing the Basal Envelope Surface of Talweg $\Phi_\mathrm{BEST}$
*/
        {
	  gamma0 = 0.;

	  for(int s=0;s<nslit;s++)
	    Q[s] = 0.;

	  (void) solve_correcting_dipoles(nsink, nbound, ncoef, ntheta,
 					  & UC, & UR, & LC, & LR,
					  & ReZ, & ImZ,
					  & SlitLength, & PHI_topo,
					  & PHI_areal1, & Vn_areal1,
					  & UL, & LL, gamma0, & dA,
					  & Q, & ReCn, & ImCn, & v0x, & v0y, & Phi0,
					  & PHI );
        }

	Phi0 += (v0x*X0 + v0y*Y0);

	int varout_grid[3] = {0,0,0};
        J_PHI = NULL; J_PSI = NULL;
        J_Vx = NULL; J_Vy = NULL;

/**
## Compute map of potential and write as .FLT
*/ 
	foreach()
	{
	   Xp = x;
	   Yp = y;
           double res;

           (void) inpolygon_ ( &ReZstart[nsink],&ImZstart[nsink],&nbound,   
                               &Xp,&Yp,&one,&res);

           if(res>0.5)
           {
	     (void) slit_ ( &Xp,&Yp,&one,
		   	    ReZstart,ImZstart,ReZend,ImZend,&nslit,   
			    ReCn,ImCn,&ncoef,Q,&v0x,&v0y,&Phi0,
			    varout_grid,
			    &PHI1,&PSI1,&Vx1,&Vy1,J_PHI,J_PSI,J_Vx,J_Vy);
             Phi[] = PHI1; 

             if(gamma0!=0.)
             {
                (void) arealsink_ ( &ReZstart[nsink],&ImZstart[nsink],&nbound,   
                                    &Xp,&Yp,&one,&gamma0,
                                    &PHI1,&Vx1,&Vy1,&A);
                Phi[] += PHI1;
             }
	   }
           else
             Phi[] = -9999.;
	}

        fprintf(stderr,"Done computing Phi[] map. Writing to file...\n");
        fflush(stderr);

	char fltfile[200];
  	sprintf(fltfile,"Phi.flt");
  	output_flt(Phi,fltfile);

	// Deallocate other pointers

	free(ReZstart);
	free(ImZstart);
	free(Phistart);
	free(ReZend);
	free(ImZend);
	free(Phiend);
	free(MinUpID);
	free(ReCn);
	free(ImCn);
	free(Q);
	free(ReZ);		
	free(ImZ);		
	free(PHI);
        free(PHI_areal1);
        free(Vn_areal1);
        free(UL);free(UC);free(UR);
        free(LL);free(LC);free(LR);
        free(PL);free(PC);free(PR);
        free(PHI_areal1_gw);
}
/**
![](https://comptes-rendus.academie-sciences.fr/geoscience/article/CRGEOS_2023__355_S1_79_0/jpg/src/tex/figures/fig14.jpg){width="800px"}
*/