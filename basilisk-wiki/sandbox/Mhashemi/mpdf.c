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
double mpdf_domain_avg[Nbin], mpdf_far_avg[Nbin], mpdf_near_avg[Nbin];
double sum_domain,  sum_plate,  sum_box,  mass_domain,  mass_plate,  mass_box;

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


fprintf (fout, "%g %g %g \n",  p1.x, p1.y, p1.z);
fprintf (fout, "%g %g %g \n",  p2.x, p2.y, p2.z);
fprintf (fout, "%g %g %g \n",  p3.x, p3.y, p3.z);
fprintf (fout, "%g %g %g %d \n", range, bin, l_cut, maxlevel);


  
for (int k = 0; k < Nbin; k++) {
 
pdf_domain_avg[k]=pdf_far_avg[k]=pdf_near_avg[k]=0;
mpdf_domain_avg[k]=mpdf_far_avg[k]=mpdf_near_avg[k]=0;

 }

sum_domain=sum_plate=sum_box=mass_domain=mass_plate=mass_box=0;

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
   MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
 #endif

 

 if (mpi_rank==0) {
     
       coord cenoid[n];
       double d[n];
  
       for (int j=0; j< n; j++){
        
	         d[j]= 2*cbrt(0.75*v[j]/pi);    
           foreach_dimension()
           cenoid[j].x=b[j].x/v[j];            
       }
      
 fprintf(fout, "time index volume cen_X cen_Y size\n");


for (int j = 0; j < n; j++)
    fprintf (fout, "%g %d %g %g %g %g \n", t,
             j, v[j], cenoid[j].x, cenoid[j].y, d[j]);
  fflush (fout);




//droplet PDF:
    
    double l[Nbin], pdf[Nbin], pdf_plate[Nbin], pdf_box[Nbin], mdomain[Nbin], mbox[Nbin], mplate[Nbin];  

/** It is being done in one processor,
     so we need to store everything in one processor*/

    int domain=0;  // whole domain
    int box=0;     //near_field 
    int plate=0;   //plate
    double domain_mass=0;
    double box_mass=0;
    double plate_mass=0;
        
    for (int k=0; k<Nbin; k++){

       l[k]=l0+k*bin;
       pdf[k]=pdf_plate[k]=pdf_box[k]=0;
       mdomain[k]=mbox[k]=mplate[k]=0;

      }

    for ( int j=0; j<n; j++) {  

        if ( cenoid[j].x <= p1.x  &&          // jet direction
             cenoid[j].y <= p1.y  &&          // crossflow direction
	            // z_min <= cenoid[j].z <= z_max  && 
              l_cut  <= d[j] && d[j] <= lmax ) {      // treshhold for poor-resolved droplets && volume fraction in unrefined cells && large              
	  
            domain++;  
	    domain_mass += v[j] ;

            int k2 = floor((d[j]-l0)/bin);
            pdf[k2]++;
	    mdomain[k2] += v[j]; 
           
	   /** if (k2==0) {

            fprintf (fout, "this is index and size %d %g\n",
             j, d[j]);
            fflush (fout);


	    } */

           }
  
        if ( cenoid[j].x <= p2.x  &&
             cenoid[j].y <= p2.y  &&
             //z_min <= cenoid[j].z <= z_max &&
      	     l_cut  <= d[j] && d[j] <= lmax ) {

            box++;
	    box_mass += v[j];

            int k2 = floor((d[j]-l0)/bin);
            pdf_box[k2]++;
	    mbox[k2] += v[j];

           }

        if ( // p3.x <= cenoid[j].x <= p3.x+d0  &&
             p3.y <= cenoid[j].y && cenoid[j].y <= p3.y+d0  &&
             // z_min <= cenoid[j].z <= z_max  &&  
	           l_cut  <= d[j] && d[j] <= lmax ) {

            plate++;
	    plate_mass += v[j] ;

            int k2 = floor((d[j]-l0)/bin);            
	    pdf_plate[k2]++;
	    mplate[k2] += v[j];

            }

  }

//print the Frequency
   fprintf(fout, "bin_center Frequency_Domain Frequency_Near Frequency_Far mass_Domain mass_box mass_plate\n");
   for (int k = 0; k < Nbin; k++)
        fprintf (fout, "%g %g %g %g %g %g %g \n",  l[k]+bin/2,
                 pdf[k], pdf_box[k], pdf_plate[k], mdomain[k], mbox[k], mplate[k]);
        fflush (fout);
//average the trajectory based on time 
  for (int k = 0; k < Nbin; k++){

          pdf_domain_avg[k]+=pdf[k]; 
          pdf_near_avg[k]+=pdf_box[k]; 
      	  pdf_far_avg[k]+=pdf_plate[k]; 
          
	  mpdf_domain_avg[k] += mdomain[k];
          mpdf_near_avg[k] += mbox[k];
          mpdf_far_avg[k] += mplate[k];

	}

  sum_domain +=domain;
  sum_box += box;
  sum_plate += plate;
  mass_domain += domain_mass;
  mass_box  += box_mass;
  mass_plate += plate_mass;
 

  printf("******postprocess is done for t=%g\n", t);

  fprintf (fout, "%d  \n", domain);
  fprintf (fout, "%d  \n", box);
  fprintf  (fout, "%d  \n", plate);
  fprintf  (fout, "%g  \n", domain_mass);
  fprintf  (fout, "%g  \n", box_mass);
  fprintf  (fout, "%g  \n", plate_mass);




 // scanf("%d\n", o);

 }

 delete (interesting_scalars);

}


if (mpi_rank==0) {

  fprintf(fout, "bin_center PDF_domain_avg PDF_Near_avg  PDF_Far_avg Mass_domain_avg Mass_PDF_Near_avg  Mass_PDF_Far_avg\n");
   for (int k = 0; k < Nbin; k++){
        
	 l[k]=l0+k*bin;
        fprintf (fout, "%g %g %g %g %g %g %g\n", (l[k]+bin/2)*1e+6,
                 pdf_domain_avg[k]/(bin*1e+6*sum_domain), pdf_near_avg[k]/(bin*1e+6*sum_box), pdf_far_avg[k]/(bin*1e+6*sum_plate),
		 mpdf_domain_avg[k]/(bin*1e+6*mass_domain), mpdf_near_avg[k]/(bin*1e+6*mass_box), mpdf_far_avg[k]/(bin*1e+6*mass_plate)
		 );   
     }

     fflush (fout);     
 }        

}












