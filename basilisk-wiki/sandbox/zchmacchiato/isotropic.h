inline double isotropic_x(Point point, scalar phi){
    return(((phi[1,0]-phi[-1,0])/3.0+(phi[1,1]-phi[-1,-1])/12.0+(phi[1,-1]-phi[-1,1])/12.0)/Delta);

} 
inline double isotropic_y(Point point, scalar phi){
    return(((phi[0,1]-phi[0,-1])/3.0+(phi[1,1]-phi[-1,-1])/12.0+(phi[-1,1]-phi[1,-1])/12.0)/Delta);

} 

inline double isotropic_laplace(Point point, scalar phi ){
    return((phi[1,1]+phi[-1,1]+phi[-1,-1]+phi[1,-1]+4.0*phi[1,0]+4.0*phi[-1,0]+4.0*phi[0,1]+4.0*phi[0,-1]-20*phi[0,0])/6.0/sq(Delta));

} 