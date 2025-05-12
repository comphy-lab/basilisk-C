#include "grid/octree.h" 
#include "utils.h"
#include "fractions.h"
#include "tag.h"
#include "output.h"

#define d0  8.0e-4
#define Nbin 25
#define lmax  3.0e-4 
#define l0 0.0 
double init= 0.0069;
double end=0.0069;
#define x_b 0.0123
#define y_b 0.00276
#define x_r 1.0e2    // HUGE
#define y_t 1.0e2    // HUGE
#define dis 19.5
#define L  0.03125
#define cut_r 2

int maxlevel=11;
double cell_size, l_cut;

int * o=0;


double l[Nbin], pdf_domain_avg[Nbin], pdf_far_avg[Nbin], pdf_near_avg[Nbin];
double sum_domain,  sum_plate,  sum_box;

char name[80];
char name2[80];
char name3[80];
char name4[80];

double bin=0;    
double range=0;  
         

coord p1={x_r, y_t, 0};   // main box {jet_dir (x_r),  cross_dir(y_r),  width}
coord p2={x_b, y_b, 0};   // near field {height_breakup,  horizental_breakup, width}
coord p3={0, dis*d0, 0};   // far-field plate { , 19.5*d  , }

//int v_cut=4/3*pi*pow((0.5*d_cut),3);
//int v_max=4/3*pi*pow((0.5*lmax),3);


    
int main (int argc, char * argv[])
 {
if (argc > 1)
    p1.x = atof (argv[1]);
if (argc > 2)
    p1.y = atof (argv[2]);
if (argc > 3)
    p2.x = atof (argv[3]);
if (argc > 4)
    p2.y = atof (argv[4]);
if (argc > 5)
    init = atof (argv[5]);
if (argc > 6)
    end = atof (argv[6]);
if (argc > 7)
    maxlevel = atoi (argv[7]);


range=lmax-l0;
bin=range/Nbin;

cell_size=L/pow(2, maxlevel);
l_cut= cut_r*cell_size;



fprintf(fout, "t n domain box plate V_tot V_domain_filtered V_box_filtered V_plate_filtered\n");
/**
fprintf (fout, "%g %g %g \n",  p1.x, p1.y, p1.z);
fprintf (fout, "%g %g %g \n",  p2.x, p2.y, p2.z);
fprintf (fout, "%g %g %g \n",  p3.x, p3.y, p3.z);
fprintf (fout, "%g %g %g %d \n", range, bin, l_cut, maxlevel);
//scanf("%d\n", o);
*/
  
for (int k = 0; k < Nbin; k++) {
 
 pdf_domain_avg[k]=pdf_far_avg[k]=pdf_near_avg[k]=0;
 
 }

for (double t=init; t <= end; t+=2e-5){

    
  sprintf (name, "snapshot-%g", t);
  sprintf (name2, "reduced_size_dump-%g", t);   
  sprintf (name3, "jet_column-new%g", t);   
  sprintf (name4, "gfsfield-%g", t);
  init_grid(2);
  vector u[];
  scalar f[];
  f.prolongation = fraction_refine; //for volume fraction field

  scalar * interesting_scalars =  {f} ; 
  restore (name2, interesting_scalars);
  
  boundary(interesting_scalars);


   scalar m[];
   foreach()
    m[] = f[] > 1e-3;
    int n = tag (m);

   /**
   Once each cell is tagged with a unique droplet index, we can easily
   compute the volume *v* and position *b* of each droplet. Note that
   we use *foreach_leaf()* rather than *foreach()* to avoid doing a
   parallel traversal when using OpenMP. This is because we don't have
   reduction operations for the *v* and *b* arrays (yet). */


   double v[n];
   coord b[n];    
   for (int j = 0; j < n; j++)
     v[j] = b[j].x = b[j].y = b[j].z = 0.;
   foreach_leaf()
     if (m[] > 0) {
       int j = m[] - 1;
       v[j] += dv()*f[];
       coord p = {x,y,z};
       foreach_dimension()
         b[j].x += dv()*f[]*p.x;   
     }

  /**
  When using MPI we need to perform a global reduction to get the
  volumes and positions of droplets which span multiple processes. */


 #if _MPI
   MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // We need to add the samething for velocity as well, Vp 
 #endif


/**

for (int j = 0; j < n; j++)
    fprintf (fout, "%g %d %g %g %g\n", t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j]);
  fflush (fout);
*/
 

 if (mpi_rank==0) {
     
       coord cenoid[n];
       double d[n];
  
       for (int j=0; j< n; j++){
        
	         d[j]= 2*cbrt(0.75*v[j]/pi);    
           foreach_dimension()
           cenoid[j].x=b[j].x/v[j];            
       }
      
// fprintf(fout, "time index volume cen_X cen_Y size\n");

/**
for (int j = 0; j < n; j++)
    fprintf (fout, "%g %d %g %g %g %g \n", t,
             j, v[j], cenoid[j].x, cenoid[j].y, d[j]);
  fflush (fout);
*/



//droplet PDF:
    
    double l[Nbin], pdf[Nbin], pdf_plate[Nbin], pdf_box[Nbin], V_tot, V_domain_filtered, V_box_filtered, V_plate_filtered;  

/** It is being done in one processor,
     so we need to store everything in one processor*/

    int domain=0;  // whole domain
    int box=0;     //near_field 
    int plate=0;   //plate

    V_tot=V_domain_filtered=V_box_filtered=V_plate_filtered=0;
   
    for ( int k=0; k<Nbin; k++){

       l[k]=l0+k*bin;
       pdf[k]=pdf_plate[k]=pdf_box[k]=0;

      }

    for ( int j=0; j<n; j++) {  
         
	 V_tot += v[j];

        if ( cenoid[j].x <= p1.x  &&          // jet direction
             cenoid[j].y <= p1.y  &&          // crossflow direction
	            // z_min <= cenoid[j].z <= z_max  && 
              l_cut  <= d[j] && d[j] <= lmax ) {      // treshhold for poor-resolved droplets && volume fraction in unrefined cells && large              

            V_domain_filtered += v[j];  
            domain++;  
            int k2 = floor((d[j]-l0)/bin);
            pdf[k2]++;

           }
  
        if ( cenoid[j].x <= p2.x  &&
             cenoid[j].y <= p2.y  &&
             //z_min <= cenoid[j].z <= z_max &&
      	     l_cut  <= d[j] && d[j] <= lmax ) {

             V_box_filtered += v[j];
            box++;
            int k2 = floor((d[j]-l0)/bin);
            pdf_box[k2]++;

           }

        if ( // p3.x <= cenoid[j].x <= p3.x+d0  &&
             p3.y <= cenoid[j].y && cenoid[j].y <= p3.y+d0  &&
             // z_min <= cenoid[j].z <= z_max  &&  
	           l_cut  <= d[j] && d[j] <= lmax ) {

            V_plate_filtered += v[j];
            plate++;
            int k2 = floor((d[j]-l0)/bin);
            pdf_plate[k2]++;

            }

  }


//  printf("******postprocess is done for t=%g\n", t);
//  fprintf(fout, "t n domain box plate V_tot V_domain_filtered V_box_filtered V_plate_filtered\n");
  fprintf (fout, "%g %d %d %d %d %g %g %g %g \n", t, n, domain, box, plate,
          V_tot, V_domain_filtered, V_box_filtered, V_plate_filtered);
 }

 delete (interesting_scalars);

}
        

}




