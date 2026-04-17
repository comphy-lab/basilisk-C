/**
 
 # Example of "Point point locate()" not working on GPU.

 the code with #if 1 : correct version (foreach_point())
 the code with #if 0 : incorrect version (Point point locate())

Using Point point locate() should raise an error (see [here](https://basilisk.fr/sandbox/bugs/pointpoint.c)).

 */
scalar * list_out;
float data;
scalar * nc_scalar_list;
scalar eta;
vector u;

int main(int argc, char *argv[])  
{

  init_grid (1);
  eta = new scalar;
  u = new vector[1];

  list_out = {eta, u};
  nc_scalar_list = list_copy(list_out);

  float xp=0.5;
  float yp=0.5;
  
  //set value at [0,0] to be 0.
  foreach(){
    u.x[0,0] = 0.0;
    u.y[0,0] = 0.0;
    eta[0,0] = 0.0;
  }

  // print value at [0,0]
  for (scalar s in nc_scalar_list){
    fprintf(stderr, "%s ", s.name); 
  #if 1
    foreach_point(xp,yp) 
      fprintf(stderr,"i=%d, j=%d, %f\n", 0, 0, val(s));
  #else
    Point point = locate (xp, yp);
    data = point.level >= 0 ? val(s) : nodata;
    fprintf(stderr,"i=%d, j=%d, %f\n", 0, 0, data);
  #endif
  }

  //update value at [0,0]
  foreach(){
    u.x[0,0] = 1.0;
    u.y[0,0] = 2.0;
    eta[0,0] = -1.0;
  }

  // print again value at [0,0], should be 1 for u.x, 0 for eta
  for (scalar s in nc_scalar_list){
    fprintf(stderr, "%s ", s.name); 
  #if 1
    foreach_point(xp,yp) 
      fprintf(stderr,"i=%d, j=%d, %f\n", 0, 0, val(s));
  #else
    Point point = locate (xp, yp);
    data = point.level >= 0 ? val(s) : nodata;
    fprintf(stderr,"i=%d, j=%d, %f\n", 0, 0, data);
  #endif
  }

} // end main
