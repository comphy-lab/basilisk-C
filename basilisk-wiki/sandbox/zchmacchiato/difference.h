// foreach_dimension()
// inline double iso_x(Point point, scalar phi){
//     return(((phi[1,0]-phi[-1,0])/3.0+(phi[1,1]-phi[-1,-1])/12.0+(phi[1,-1]-phi[-1,1])/12.0)/Delta);

// } 

// foreach_dimension()
// inline double bias_x(Point point, scalar phi){
//     return((-phi[2,0]+4.0*phi[1,0]+phi[-2,0]-4.0*phi[-1,0])/3.0+(-phi[2,2]+4.0*phi[1,1]+phi[-2,2]-4.0*phi[-1,1])/12.0+(-phi[2,-2]+4.0*phi[1,-1]+phi[-2,-2]-4.0*phi[-1,-1])/12.0);

// } 

// foreach_dimension()
// inline double mix_x(Point point, scalar phi){
//     double mix;//
//     mix=0.5*iso_x(point,phi)+0.5*bias_x(point,phi);
//     return(mix);

// } 
// foreach_dimension()
// inline double isof_x(Point point, scalar phi){
//     return((phi[1,0]-phi[-1,0])/3.0+(phi[1,1]-phi[-1,-1])/12.0+(phi[1,-1]-phi[-1,1])/12.0);

// } 

// foreach_dimension()
// inline double isob_x(Point point, scalar phi){
//     return((phi[0,0]-phi[-2,0])/3.0+(phi[0,1]-phi[-2,-1])/12.0+(phi[0,-1]-phi[-2,1])/12.0);

// }
// inline double iso_lap(Point point, scalar phi ){
//     return((phi[1,1]+phi[-1,1]+phi[-1,-1]+phi[1,-1]+4.0*phi[1,0]+4.0*phi[-1,0]+4.0*phi[0,1]+4.0*phi[0,-1]-20*phi[0,0])/6.0/sq(Delta));

// } 
foreach_dimension()
inline double forward_x(Point point, scalar phi){
        
    return((phi[1,0]-phi[0,0])/Delta);

} 

foreach_dimension()
inline double backward_x(Point point, scalar phi){
        
    return((phi[0,0]-phi[-1,0])/Delta);

} 



foreach_dimension()
inline double center_x(Point point, scalar phi){
     return((phi[1,0]-phi[-1,0])/2.0/Delta);

} 


inline double laplace(Point point, scalar phi){
     return((phi[1,0]+phi[-1,0]+phi[0,1]+phi[0,-1]-4.0*phi[])/sq(Delta));

} 

