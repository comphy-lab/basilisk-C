#define L_domain 2.0           // the length of computation domian
#define num_initgrid 256    // initial grid 
#define numofcell num_initgrid // 2^Max_Level
#define dxmin (L0/numofcell)   // size of the finest mesh 
#define Max_Level 6           // the finest level
#define DT_Max  1.e-6          // dt
#define rd  (1.0)             // ratio of Density
#define rv  (1.0)             // ratio of Viscosity
#define rk  (1.0)             // ratio of thermal conductivity k_s/k_l
#define rCp (1.0)             // ratio of Heat capacity Cp_s/Cp_l


#define dtoutput 0.1          // the time interval to output
#define endoutput (dtoutput*10.00001) 


#define A_AC (1.0)

#define Re (100.0)  
#define Pe (1.0)
#define We (10.)             // Weber number    We = rho*u*u*a/sigma
#define St (1.0)             // Stefan number    St=L/Cpl △T
#define Tm (1.0)                //Dimensionless melting temperature

#define two_phase
#define solidification


#define TmpE_on               // define TmpE_on mean ACE is used 

//#define NS_ON  
//#define DT_adapt_on           
//----------------------------------------------------------------------------------------

double Cn,D_AC,xdomain=L_domain,ydomain=(L_domain*1.0/num_initgrid),zdomain=L_domain;   // the thickness of interface;
double eta=2.0*DT_Max;
double totalphi1=0.,totalphi10 = 0.;
scalar phi1[],convprev[],Tm_eff[];
scalar phi1_ns[],phi1new[];
scalar phi1old[];


scalar T[],convprev_TE[];  //define temperature equation variables
FILE *fp1,*fp2,*fp3;

#ifdef  TmpE_on
	scalar * tracers = {T,phi1};
#else
	scalar * tracers = {phi1};
#endif 




//-----------------------------------------------------------
void WENO();
void wirte_ASCII_of_str();
double  p_interpolation(double x);
void output_a_field( scalar phi1,FILE * fp);
void output_a_field_imagestec(scalar scalar_f,FILE * fp);