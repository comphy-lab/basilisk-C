#include <stdio.h>

#include "interpolation_data.h"

#define verbose 1
#define verbose_Pascal 0



double lerp (double a, double b, double t) { // Linear interpolation
  double x = (1. - t)*a + t*b ;
  return x;
}


double Lagrange (int n, double Ptsx[n], double Ptsy[n], double x) { // Lagrange interpolation of order n
  double L = 0., l;
  int i, j;
  for (i = 0; i < n; i++) {
    l = 1.;
    for (j = 0; j < n; j++) {
      if (i != j)
	l *= (x - Ptsx[j])/(Ptsx[i] - Ptsx[j]);
    }
    L += l * Ptsy[i];
  }
  return L;  
}


double Bezier (int n, double Ctrl_Pts[n+1], double t) { // Bézier approximation of order n (n+1 control points)
  double B = 0., Bernstein;
  long long P[n+1][n+1]; // TODO : compute Pascal's triangle only once, when necessary (ideas : compute it outside, or make a test every time ?)
  int i, j;
  // Pascal's triangle (in order to "compute" binomial coefficient easily just after)
  for (i = 0; i < n+1; i++) { // We could optimize the process by only computing the last 2 lines, as they're the only ones used here
    for (j = 0; j < n+1; j++) {
      if ((j == 0) || (i == j))
	P[i][j] = 1;
      else if (j < i) {
	P[i][j] = P[i-1][j-1] + P[i-1][j];
      } else
	P[i][j] = 0;
    }
  }
  #if verbose_Pascal
  fprintf (stdout, "Pascal's triangle :\n");
  for (i = 0; i < n+1; i++) {
    for (j = 0; j < n+1; j++) {
      fprintf (stdout, "%lli ", P[i][j]);
    }
    fprintf (stdout, "\n");
  }
  #endif
  // Proper Bézier curve computation
  for (i = 0; i <= n; i++) {
    Bernstein = P[n][i] * pow(t, i) * pow(1. - t, n - i);
    B += Bernstein * Ctrl_Pts[i];
  }
  return B;      
}





void interpolate_bathymetry_from_datafile (scalar zb) {
  // Read the datafile containing the topography points to interpolate
  const int topo_size = min (N, interpolation_data_N);
  const int interpolation_level = interpolation_order;
  const int LorB = interpolation_type;
  int ii, jj;
  double tmpx[interpolation_level+1], tmpy[interpolation_level+1];
  
#if verbose
  fprintf (stdout, "TOPOGRAPHY DATA\n");
  fprintf (stdout, "Number of data points = %i\n", topo_size);
  fprintf (stdout, "Data points :\n");
#endif
  double topo_data[topo_size+1]; // Stores the values in the datafile "topo_data"
  for (ii = 0; ii < topo_size; ii++) {
    topo_data[ii] = interpolation_data_values[ii];
#if verbose
    fprintf (stdout, "%g ", topo_data[ii]);
#endif
  }
  topo_data[topo_size] = topo_data[topo_size-1]; // We add a "ghost" datapoint
#if verbose
  fprintf (stdout, "(%g)\n", topo_data[topo_size]);
  fprintf (stdout, "Indices :\n");
#endif

  
  if (topo_size == 1) { // No interpolation, flat batyhmetry at the desired elevation
    foreach()
      zb[] = topo_data[0];
    
  } else {
    // Create a table containing indices of the mesh at which the bathymetry will be interpolated
    int indices[topo_size+1];
    for (ii = 0; ii < topo_size; ii++) {
      indices[ii] = (int) ((double)ii * ((double)N / (double)(topo_size-1)));
#if verbose
      fprintf (stdout, "%i ", indices[ii]);
#endif
    }
    indices[topo_size] = indices[topo_size-1] + 1; // We add a "ghost" index here too
#if verbose
    fprintf (stdout, "(%i)\n", indices[topo_size]);
    fprintf (stdout, "Positions of indices :\n");
#endif
    // Idem but for the positions x of the interpolated points (corresponding to the indices calculated just above)
    double positions_indices[topo_size];
    ii = 0; jj = 0;
    foreach(serial) {
      jj += 1;
      if (jj >= indices[ii]) {
	positions_indices[ii] = x;
	ii += 1;
      }
    }
#if verbose
    for (ii = 0; ii < topo_size; ii++)
      fprintf (stdout, "%.4g ", positions_indices[ii]);
    fprintf (stdout, "\n\n");
#endif
    // Properly initialize bathymetry zb by interpolating the data
    ii = -1; jj = 0;
    foreach(serial) {
      jj += 1; // Index of the cell (from the global mesh) we're looking at
      if (jj >= indices[ii+1]) {
	ii++; // Index of the segment containing cell jj
      }
    
      if (interpolation_level == 0) // Nearest neighbor interpolation
	zb[] = (jj-indices[ii] < indices[ii+1]-jj ? topo_data[ii] : topo_data[ii+1]);
    
      else if (interpolation_level == 1) { // Linear interpolation
	double tt = (double)(jj - indices[ii])/(double)(indices[ii+1] - indices[ii]); // "part" of the segment where we're at (in [0;1])
	zb[] = lerp (topo_data[ii], topo_data[ii+1], tt);
      
      } else if (interpolation_level > 0) { // By part interpolation
	int ij;
	static int interpol_order = 1; // "True" order of the interpolation (especially useful for the last points)
	static int previous_part_ii = -1; // To only define interpolation points once, when "entering" the segment(s) to be interpolated
	if ((ii % interpolation_level == 0) && (previous_part_ii != ii)) { // cf just above
	  previous_part_ii = ii;
	  for (ij = 0; ij < interpolation_level+1; ij++) { // We set the interpolation points needed
	    tmpx[ij] = positions_indices[ii + min(topo_size, ij)];
	    tmpy[ij] = topo_data[ii + min(topo_size+1, ij)];
	  }
	  interpol_order = (topo_size > ii + interpolation_level ? interpolation_level : max(1, topo_size-ii-1)); // Order of the last interpolation if there are not enough points left for the desired interpolation order
	}
	if (LorB) {
	  double tt = (x - tmpx[0])/(tmpx[interpol_order] - tmpx[0]);
	  zb[] = Bezier (interpol_order, tmpy, tt);
	} else
	  zb[] = Lagrange (interpol_order+1, tmpx, tmpy, x);
      
      } else // Interpolation over all points
	if (LorB) {
	  double tt = (x - X0)/L0;
	  zb[] = Bezier (topo_size-1, topo_data, tt);
	} else
	  zb[] = Lagrange (topo_size, positions_indices, topo_data, x);
    }
  
    // Write a file describing the interpolation
    FILE * fp = fopen ("interpolation_points", "w");
    for (ii = 0; ii < topo_size; ii++)
      fprintf (fp, "%g %g\n", positions_indices[ii], topo_data[ii]);
    fclose (fp);
  }
}
