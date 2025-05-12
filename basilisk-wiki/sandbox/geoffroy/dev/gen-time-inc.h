/**
# Generic time-integration scheme
*/
#include "utils.h"

// Required from solver
// fields updated by time-integration
extern scalar * evolving;

///////////////////////////////////////////////////////////////////////////
// Butcher tables, simplified by Kopal and re-organized by Williamson
// See "Low-storage runge-kutta schemes", williamson,1980
// & "Fourth-order 2N-storage rk schemes", Carpenter & Kennedy, 1994, NASA
static double ** butcherd(order){
  // butcher[0][i] = ai
  // butcher[1][i] = bi
  double ** butcher;
  butcher = malloc(2*order*sizeof(double));
  butcher[0] = malloc(order*sizeof(double));
  butcher[1] = malloc(order*sizeof(double));

  butcher[0][0] = 0; //a1
  if(order == 1)
    butcher[1][0] = 1;
  else if(order == 2){
//      a[]                         b[]  
                           butcher[1][0] = 0.5;//1
  butcher[0][1] = -0.5;    butcher[1][1] = 1;//2
  }
  else if(order == 3){
//      a[]                         b[]   
                             butcher[1][0] = 1/3.; //1
  butcher[0][1] = -5/9.;     butcher[1][1] = 15/16.; //2
  butcher[0][2] = -153/128.; butcher[1][2] = 8/15.; //3
  }
  else if(order == 5){
    // 5 steps for the 4th-order scheme
//      a[]                        
  butcher[0][1] = -1 ;     
  butcher[0][2] = -1/3.+pow(2,2/3.)/6. - 2*cbrt(2)/3.;     
  butcher[0][3] = -cbrt(2) - pow(2,2/3.) - 2;      
  butcher[0][4] = -1 + cbrt(2);      

  //   b[]
  butcher[1][0] = 2/3. + cbrt(2)/3. + pow(2,2/3.)/6.; //1
  butcher[1][1] = -pow(2,2/3.)/6. + 1/6.;    //2
  butcher[1][2] = -1/3. - 2*cbrt(2)/3. - pow(2,2/3.)/3.; //3
  butcher[1][3] = 1/3. - cbrt(2)/3. - pow(2,2/3.)/6.; //3
  butcher[1][4] = 1/3. + cbrt(2)/6. + pow(2,2/3.)/12.; //3
  }
  else assert(false);
 return butcher;
}

//////////////////////////////////////////////////////
// User-provided functions
// gradient
double (* gradient)  (double, double, double) = minmod2;

// how to compute updates
double (* update) (scalar * evolving, scalar * updates, double dtmax) = NULL;

double t = 0., dt = 0.;
static void advance_generic (scalar * output, scalar * input, scalar * updates,
			     double dt){
  trash (output);  
  foreach() {
    scalar o,i,u;
    for(o,i,u in output,input,updates)
       o[] = i[]+dt*u[];
  }
  boundary (output);
}

static void (* advance) (scalar * output, scalar * input, scalar * updates,
			 double dt) = advance_generic;
event defaults (i = 0)
{
  // limiting
  for (scalar s in all)
    s.gradient = gradient;

  // default values
  foreach()
    for (scalar s in all)
      s[] = 0.;
  boundary (all);
}


//////////////////////////////////////////////////////
// Run function
void run()
{
  t = 0.;
  init_grid (N);
  timer start = timer_start();
  int i = 0; long tnc = 0;

  // Take care of the order 4 which is in 5 steps
  int order_k = order_t < 4 ? order_t : 5;
  // williamson/carpenter tables from butcher tables
  double ** butcher = butcherd(order_k);
      
  while (events (i, t, true)) {
    
    // list of updates and intermediary fields
    scalar * updates = list_clone (evolving);
    scalar * interm =  list_clone (evolving); 
    // interm = 0
    foreach(){
      for( scalar s in interm)
	s[] = 0;
    }
    
    dt = dtnext (t, update (evolving, updates, DT));
    for( int order_ = 0; order_ < order_k ; order_ ++ ){
      /* Nnd-order time-integration */
      foreach() {
	for(scalar s in interm)
	  s[] *= butcher[0][order_];
      }
      //interm = a*qj-1
      
      update(evolving,updates,dt); // update = F(Yj-1)
      advance (interm, interm, updates, dt);  // interm = qj = a*qj-1 + dt F(Yj-1) 
     
      foreach(){
	scalar s1,s2;
	for(s1,s2 in evolving,interm)
	  s1[] += s2[]*butcher[1][order_];

      }
      // Yj = Yj-1 + bj*qj
    }
    
    delete(interm);
    free(interm);
    delete (updates);
    free (updates);
    
    long nc = 0;
    foreach (reduction(+:nc)) 
      nc++;
    tnc += nc;
    i++; t = tnext;
  }
 
  timer_print (start, i, tnc);
  free_grid();
  free(butcher);
}
