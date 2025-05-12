/**
# Output functions

## *output_field()*: Multiple fields interpolated on a regular grid (text format)

This function interpolates a list of fields onto a regular N\*N\*N Cartesian grid. It outputs each processor's domain data individually and you have to catch the data together for post-processing.

The arguments and their default values are:

*list*
: list of fields to output. Default is *all*.

*n*
: number of points along each dimension. Default is *N*.

*linear*
: whether to use linear interpolation or not. 

*box*
: the lower-left-back and upper-right-front coordinates of the domain to consider.
 Default is the entire domain. */

struct OutputField_3d {
  scalar * list;
  int n;
  bool linear;
  double box[2][3];
};

trace
void output_field_3d (struct OutputField_3d p)
{
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && p.box[0][2]==0.0 &&
      p.box[1][0] == 0. && p.box[1][1] == 0. && p.box[1][2]==0.0) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;      p.box[0][2]=Z0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0; p.box[1][2]=Z0+L0;
  }
  
  double Delta = 0.999999*(p.box[1][0] - p.box[0][0])/(p.n);
  int nx = p.n;
  int ny = p.n;
  int nz = p.n;
  double field;
  char name[80];
  sprintf(name, "outputfield_%g_%d.csv", t, pid());
  FILE * fp = fopen(name, "a");
  for (int i = 0; i < nx; i++) {
    double xp = Delta*(i+0.5) + p.box[0][0];
    for (int j = 0; j < ny; j++) {
      double yp = Delta*(j+0.5) + p.box[0][1];
      for (int k=0; k<nz; k++) {
        double zp=Delta*(k+0.5)+p.box[0][2];
        if (p.linear) {
	  //2023.11.25: judge if the point is on the current processor.
	  if(locate(xp,yp,zp).level>0) {
	    fprintf(fp, "%g %g %g", xp, yp, zp);
	    for (scalar s in p.list) {
	      field = interpolate (s, xp, yp, zp);
	      fprintf(fp, " %g", field);
	    }
	    fputc('\n', fp);  //not ""
	  }
        }
        else {
          if(locate(xp,yp,zp).level>0) {
            Point point = locate (xp, yp, zp);
            fprintf(fp, "%g %g %g", xp, yp, zp);
            for (scalar s in p.list) {
              field = s[];
              fprintf(fp, " %g", field);
            }
            fputc('\n', fp);  //not ""
          }
        }
      }
    }
  }
  fclose(fp);
}
/** After the calculation, you may integrate the individual data files into one using the command  

~~~bash
cat outputfield_* > outputfield_total
~~~

*/