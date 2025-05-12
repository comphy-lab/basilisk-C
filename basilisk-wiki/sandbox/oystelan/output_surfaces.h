/**
Various output functions for extracting a surface from field data using the build in VOF PLIC surface or the isosurface reconstruction in bview
Three different formats are supported:
- .ply
- .vtu
- .vtp
all of which can be loaded into Paraview. for the two later alternatives, field data may be added by using the functions with "_w_fielddata".
This is an updated version of the previous fractions_output.h

Author: Oystein Lande
Date: July 2019
*/

#include "geometry.h"
#include "fractions.h"

#if dimension == 1
coord mycs (Point point, scalar c) {
    coord n = {1.};
    return n;
}
#elif dimension == 2
# include "myc2d.h"
#else // dimension == 3
# include "myc.h"
#endif


/**
 This function outputs the VOF surface to a .ply file
 */
void output_ply (struct OutputFacets p)
{

    #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
    #endif
    scalar c = p.c;
    face vector s = p.s;
    if (!p.fp) p.fp = stdout;

    // print header text
    fputs ("ply\n", p.fp);
    fputs ("format ascii 1.0\n", p.fp);

    int nverts = 0;
    int nfacets = 0;

    foreach()
        if (c[] > 1e-6 && c[] < 1. - 1e-6) {
            coord n;
            if (!s.x.i)
                // compute normal from volume fraction
                n = mycs (point, c);
            else {
                // compute normal from face fractions
                double nn = 0.;
                foreach_dimension() {
                    n.x = s.x[] - s.x[1];
                    nn += fabs(n.x);
                }
                assert (nn > 0.);
                foreach_dimension()
                    n.x /= nn;
            }
            double alpha = plane_alpha (c[], n);

            coord v[12];
            int m = facets (n, alpha, v, 1.);
            for (int i = 0; i < m; i++) {
                nverts ++;
            }
            if (m > 0) {
                nfacets ++;
            }
        }

        fprintf (p.fp, "element vertex %i\n", nverts);
        fputs ("property float x\n", p.fp);
        fputs ("property float y\n", p.fp);
        fputs ("property float z\n", p.fp);
        fprintf (p.fp, "element face %i\n", nfacets);
        fputs ("property list uchar int vertex_index\n", p.fp);
        fputs ("end_header\n", p.fp);

        int facet_num[nfacets];

        int ifacet = 0;
        int ivert = 0;

        foreach()
            if (c[] > 1e-6 && c[] < 1. - 1e-6) {
                coord n;
                if (!s.x.i)
                    // compute normal from volume fraction
                    n = mycs (point, c);
                else {
                    // compute normal from face fractions
                    double nn = 0.;
                    foreach_dimension() {
                        n.x = s.x[] - s.x[1];
                        nn += fabs(n.x);
                    }
                    assert (nn > 0.);
                    foreach_dimension()
                        n.x /= nn;
                }
                double alpha = plane_alpha (c[], n);

                coord v[12];
                int m = facets (n, alpha, v, 1.);
                for (int i = 0; i < m; i++) {
                    fprintf (p.fp, "%g %g %g\n",
                             x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
                }
                if (m > 0) {
                    facet_num[ifacet] = m;
                    ifacet ++;
                }
            }

            // print face list
            for (ifacet = 0; ifacet < nfacets; ifacet++) {
                fprintf (p.fp, "%i ", facet_num[ifacet]);
                for (int iv = 0; iv < facet_num[ifacet]; iv ++) {
                    fprintf (p.fp, "%i ", ivert);
                    ivert ++;
                }
                fputc ('\n', p.fp);
            }

            fflush (p.fp);
            #if defined(_OPENMP)
            omp_set_num_threads(num_omp);
            #endif
}

/**
 This function outputs the iso surface to a .ply file
 */
void output_ply_iso (struct OutputFacets p)
{

    #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
    #endif
    scalar c = p.c;
    face vector s = p.s;
    if (!p.fp) p.fp = stdout;

    // print header text
    fputs ("ply\n", p.fp);
    fputs ("format ascii 1.0\n", p.fp);

    // Start by creating the vertex and normal field
    vertex scalar v[];
    foreach_vertex()
        v[] = (f[] + f[-1] + f[0,-1] + f[-1,-1] +
        f[0,0,-1] + f[-1,0,-1] + f[0,-1,-1] + f[-1,-1,-1])/8.;

    vector n[];
    foreach()
        foreach_dimension()
            n.x[] = (f[1] - f[-1])/(2.*Delta);
        boundary ((scalar *){n});


    // Loop through all surface cells
    // The point of this first round is to count the number of isosurface triangles

    int nverts = 0;
    int nfacets = 0;
    foreach(){
        //if (c[] > 1e-7 && c[] < 1. - 1e-7) {
        double val[8] = {
            v[0,0,0], v[1,0,0], v[1,0,1], v[0,0,1],
            v[0,1,0], v[1,1,0], v[1,1,1], v[0,1,1]
        };
        double t[5][3][3];
        int nt = polygonize (val, 0.5, t);
        nfacets += nt;
        nverts += nt*3;
    }

    fprintf (p.fp, "element vertex %i\n", nverts);
    fputs ("property float x\n", p.fp);
    fputs ("property float y\n", p.fp);
    fputs ("property float z\n", p.fp);
    fprintf (p.fp, "element face %i\n", nfacets);
    fputs ("property list uchar int vertex_index\n", p.fp);
    fputs ("end_header\n", p.fp);


    foreach(){
        //if (c[] > 1e-7 && c[] < 1. - 1e-7) {

        double val[8] = {
            v[0,0,0], v[1,0,0], v[1,0,1], v[0,0,1],
            v[0,1,0], v[1,1,0], v[1,1,1], v[0,1,1]
        };
        double t[5][3][3];
        int nt = polygonize (val, 0.5, t);
        for (int i = 0; i < nt; i++) {
            for (int j = 0; j < 3; j++) {
                coord v = {t[i][j][0], t[i][j][1], t[i][j][2]}, np;
                //foreach_dimension()
                //np.x = interp (point, v, n.x);
                //glNormal3d (np.x, np.y, np.z);
                //color_vertex (p, interp (point, v, col));
                //glvertex3d (view, x + v.x*Delta, y + v.y*Delta, z + v.z*Delta);
                fprintf (p.fp, "%g %g %g\n",
                         x + v.x*Delta, y + v.y*Delta, z + v.z*Delta);
            }
        }
    }

    int ifacet = 0;
    int ivert = 0;
    // print face list
    for (ifacet = 0; ifacet < nfacets; ifacet++) {
        fprintf (p.fp, "%i ", 3);
        for (int iv = 0; iv < 3; iv ++) {
            fprintf (p.fp, "%i ", ivert);
            ivert ++;
        }
        fputc ('\n', p.fp);
    }

    fflush (p.fp);
    #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
    #endif
}



/**
 This function outputs the VOF surface to a .vtu format file
 */
void output_vtu (struct OutputFacets p)
{
    #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
    #endif
    scalar c = p.c;
    face vector s = p.s;
    if (!p.fp) p.fp = stdout;

    // print header text
    fputs ("<?xml version=\"1.0\"?>\n", p.fp);
    fputs ("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n", p.fp);
    fputs ("  <UnstructuredGrid>\n", p.fp);

    int nverts = 0;
    int nfacets = 0;

    foreach()
        if (c[] > 1e-6 && c[] < 1. - 1e-6) {
            coord n;
            if (!s.x.i)
                // compute normal from volume fraction
                n = mycs (point, c);
            else {
                // compute normal from face fractions
                double nn = 0.;
                foreach_dimension() {
                    n.x = s.x[] - s.x[1];
                    nn += fabs(n.x);
                }
                assert (nn > 0.);
                foreach_dimension()
                    n.x /= nn;
            }
            double alpha = plane_alpha (c[], n);

            coord v[12];
            int m = facets (n, alpha, v, 1.);
            for (int i = 0; i < m; i++) {
                nverts ++;
            }
            if (m > 0) {
                nfacets ++;
            }
        }

        fprintf (p.fp, "    <Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n", nverts, nfacets);
        fputs ("      <Points>\n", p.fp);
        fputs ("        <DataArray type=\"Float32\" Name=\"vertices\" NumberOfComponents=\"3\" format=\"ascii\">\n", p.fp);

        int offsets[nfacets];

        int ifacet = 0;
        int offset = 0;

        foreach()
            if (c[] > 1e-6 && c[] < 1. - 1e-6) {
                coord n;
                if (!s.x.i)
                    // compute normal from volume fraction
                    n = mycs (point, c);
                else {
                    // compute normal from face fractions
                    double nn = 0.;
                    foreach_dimension() {
                        n.x = s.x[] - s.x[1];
                        nn += fabs(n.x);
                    }
                    assert (nn > 0.);
                    foreach_dimension()
                        n.x /= nn;
                }
                double alpha = plane_alpha (c[], n);

                coord v[12];
                int m = facets (n, alpha, v, 1.);
                for (int i = 0; i < m; i++) {
                    fprintf (p.fp, "%g %g %g ",
                             x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
                }
                if (m > 0) {
                    offset += m;
                    offsets[ifacet] = offset;
                    ifacet ++;
                }
            }


            fputs ("        </DataArray>\n", p.fp);
            fputs ("      </Points>\n", p.fp);
            fputs ("      <Cells>\n", p.fp);

            fputs ("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n", p.fp);

            // print vert numbers
            for (int ivert = 0; ivert < nverts; ivert++)
                fprintf (p.fp, "%i ", ivert);

            fputs ("        </DataArray>\n", p.fp);
            fputs ("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n", p.fp);

            // print offsets
            for (ifacet = 0; ifacet < nfacets; ifacet++)
                fprintf (p.fp, "%i ", offsets[ifacet]);

            fputs ("        </DataArray>\n", p.fp);
            fputs ("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", p.fp);

            // print cell type list
            for (ifacet = 0; ifacet < nfacets; ifacet++)
                fprintf (p.fp, "7 ");

            fputs ("        </DataArray>\n", p.fp);
            fputs ("      </Cells>\n", p.fp);
            fputs ("      <PointData>\n", p.fp);
            fputs ("      </PointData>\n", p.fp);
            fputs ("      <CellData>\n", p.fp);
            fputs ("      </CellData>\n", p.fp);
            fputs ("    </Piece>\n", p.fp);
            fputs ("  </UnstructuredGrid>\n", p.fp);
            fputs ("</VTKFile>\n", p.fp);

            fflush (p.fp);
            #if defined(_OPENMP)
            omp_set_num_threads(num_omp);
            #endif
}


struct OutputFacets_scalar {
    scalar c;
    FILE * fp;     // optional: default is stdout
    scalar * list;  // List of scalar fields to include when writing vtu surface to file
    vector * vlist; // List of vector fields to include.
    face vector s; // optional: default is none
};


/**
 Outputs VOF surface with fielddata in .vtu format 
 */
void output_vtu_w_fielddata (struct OutputFacets_scalar p)
{
    #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
    #endif
    scalar c = p.c;
    face vector s = p.s;
    if (!p.fp) p.fp = stdout;

    // print header text
    fputs ("<?xml version=\"1.0\"?>\n", p.fp);
    fputs ("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n", p.fp);
    fputs ("\t<UnstructuredGrid>\n", p.fp);

    int nverts = 0;
    int nfacets = 0;

    foreach()
        if (c[] > 1e-6 && c[] < 1. - 1e-6) {
            coord n;
            if (!s.x.i)
                // compute normal from volume fraction
                n = mycs (point, c);
            else {
                // compute normal from face fractions
                double nn = 0.;
                foreach_dimension() {
                    n.x = s.x[] - s.x[1];
                    nn += fabs(n.x);
                }
                assert (nn > 0.);
                foreach_dimension()
                    n.x /= nn;
            }
            double alpha = plane_alpha (c[], n);

            coord v[12];
            int m = facets (n, alpha, v, 1.);
            for (int i = 0; i < m; i++) {
                nverts ++;
            }
            if (m > 0) {
                nfacets ++;
            }
        }

        fprintf (p.fp, "\t\t<Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n", nverts, nfacets);

        // Write list of scalar field values to file
        fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", p.fp);
        for (scalar s in p.list) {
            fprintf (p.fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
            foreach(){
                if (c[] > 1e-6 && c[] < 1. - 1e-6) {
                    fprintf (p.fp, "%g\n", val(s));
                }
            }
            fputs ("\t\t\t\t </DataArray>\n", p.fp);
        }
        // Write list of vector field values to file
        for (vector v in p.vlist) {
            fprintf (p.fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
            foreach(){
                if (c[] > 1e-6 && c[] < 1. - 1e-6) {
                    #if dimension == 2
                    fprintf (p.fp, "%g %g 0.\n", val(v.x), val(v.y));
                    #endif
                    #if dimension == 3
                    fprintf (p.fp, "%g %g %g\n", val(v.x), val(v.y), val(v.z));
                    #endif
                }
            }
            fputs ("\t\t\t\t </DataArray>\n", p.fp);
        }
        fputs ("\t\t\t </CellData>\n", p.fp);


        // Write points to file
        fputs ("      <Points>\n", p.fp);
        fputs ("        <DataArray type=\"Float32\" Name=\"vertices\" NumberOfComponents=\"3\" format=\"ascii\">\n", p.fp);

        int offsets[nfacets];

        int ifacet = 0;
        int offset = 0;

        foreach()
            if (c[] > 1e-6 && c[] < 1. - 1e-6) {
                coord n;
                if (!s.x.i)
                    // compute normal from volume fraction
                    n = mycs (point, c);
                else {
                    // compute normal from face fractions
                    double nn = 0.;
                    foreach_dimension() {
                        n.x = s.x[] - s.x[1];
                        nn += fabs(n.x);
                    }
                    assert (nn > 0.);
                    foreach_dimension()
                        n.x /= nn;
                }
                double alpha = plane_alpha (c[], n);

                coord v[12];
                int m = facets (n, alpha, v, 1.);
                for (int i = 0; i < m; i++) {
                    fprintf (p.fp, "%g %g %g ",
                             x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
                }
                if (m > 0) {
                    offset += m;
                    offsets[ifacet] = offset;
                    ifacet ++;
                }
            }


            fputs ("        </DataArray>\n", p.fp);
            fputs ("      </Points>\n", p.fp);
            fputs ("      <Cells>\n", p.fp);

            fputs ("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n", p.fp);

            // print vert numbers
            for (int ivert = 0; ivert < nverts; ivert++)
                fprintf (p.fp, "%i ", ivert);

            fputs ("        </DataArray>\n", p.fp);
            fputs ("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n", p.fp);

            // print offsets
            for (ifacet = 0; ifacet < nfacets; ifacet++)
                fprintf (p.fp, "%i ", offsets[ifacet]);

            fputs ("        </DataArray>\n", p.fp);
            fputs ("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", p.fp);

            // print cell type list
            for (ifacet = 0; ifacet < nfacets; ifacet++)
                fprintf (p.fp, "7 ");

            fputs ("        </DataArray>\n", p.fp);
            fputs ("      </Cells>\n", p.fp);
            fputs ("      <PointData>\n", p.fp);
            fputs ("      </PointData>\n", p.fp);
            //fputs ("      <CellData>\n", p.fp);
            //fputs ("      </CellData>\n", p.fp);
            fputs ("    </Piece>\n", p.fp);
            fputs ("  </UnstructuredGrid>\n", p.fp);
            fputs ("</VTKFile>\n", p.fp);

            fflush (p.fp);
            #if defined(_OPENMP)
            omp_set_num_threads(num_omp);
            #endif
}

struct _interpolate_weighted {
  scalar v;
  scalar f;
  double x, y, z;
};

/** Experimental function for interpolating values with a "skewed" weight to one of the two phases */

static double interpolate_momentum_weighted (Point point, struct _interpolate_weighted p)
{
  scalar v = p.v;
  scalar wmat = p.f;
#if dimension == 1
  x = (p.x - x)/Delta - v.d.x/2.;
  int i = sign(x);
  x = fabs(x);
  double fsum = p.f[] + p.f[i];
  /* linear interpolation */
  return (v[]*(1. - x)*p.f[] + v[i]*x*p.f[i])/fsum;
#elif dimension == 2
  x = (p.x - x)/Delta - v.d.x/2.;
  y = (p.y - y)/Delta - v.d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);

  double fsum = p.f[] + p.f[i]+ p.f[0,j] + p.f[i,j];
  /* bilinear interpolation */
  return(((v[]*(1. - x)*p.f[] + v[i]*x*p.f[i])*(1. - y) +
	  (v[0,j]*(1. - x)*p.f[0,j] + v[i,j]*x*p.f[i,j])*y))/fsum;
#else // dimension == 3
  x = (p.x - x)/Delta - v.d.x/2.;
  y = (p.y - y)/Delta - v.d.y/2.;
  z = (p.z - z)/Delta - v.d.z/2.;
  int i = sign(x), j = sign(y), k = sign(z);
  x = fabs(x); y = fabs(y); z = fabs(z);

  double fsum = wmat[] + wmat[i] + wmat[0,j] + wmat[i,j] + wmat[0,0,k] + wmat[i,0,k] + wmat[0,j,k] + wmat[i,j,k];
  /* trilinear interpolation */
  return (((v[]*(1. - x)*wmat[] + v[i]*x*wmat[i])*(1. - y) +
	   (v[0,j]*(1. - x)*wmat[0,j] + v[i,j]*x*wmat[i,j])*y)*(1. - z) +
	  ((v[0,0,k]*(1. - x)*wmat[0,0,k] + v[i,0,k]*x*wmat[i,0,k])*(1. - y) +
	   (v[0,j,k]*(1. - x)*wmat[0,j,k] + v[i,j,k]*x*wmat[i,j,k])*y)*z)/fsum;
#endif
}


static inline double interp3_skewed (Point point, coord p, scalar col, scalar weights) {
  struct _interpolate_weighted _r = { col, weights, x + p.x*Delta, y + p.y*Delta, z + p.z*Delta };
  return interpolate_momentum_weighted (point, _r);
}

static inline double interp3 (Point point, coord p, scalar col) {
  struct _interpolate _r = { col, x + p.x*Delta, y + p.y*Delta, z + p.z*Delta };
  return interpolate_linear (point, _r);
}


/**
 Outputs ISO surface with fielddata in .vtp format 
 */
void output_vtp_iso_w_fielddata (struct OutputFacets_scalar p)
{
    #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
    #endif
    scalar c = p.c;
    //face vector s = p.s;
    if (!p.fp) p.fp = stdout;

    // print header text
    fputs ("<?xml version=\"1.0\"?>\n", p.fp);
    fputs ("<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n", p.fp);
    fputs ("\t<PolyData>\n", p.fp);

    // Start by creating the vertex and smoothed normal field
    vertex scalar v[];
    foreach_vertex()
        v[] = (c[] + c[-1] + c[0,-1] + c[-1,-1] +
        c[0,0,-1] + c[-1,0,-1] + c[0,-1,-1] + c[-1,-1,-1])/8.;


    /** Loop through all surface cells
     *       The point of this first round is to count the number of isosurface triangles. Should be improved...
     */
    int nverts = 0;
    int nfacets = 0;
    foreach(){
        //if (c[] > 1e-7 && c[] < 1. - 1e-7) {
        double val[8] = {
            v[0,0,0], v[1,0,0], v[1,0,1], v[0,0,1],
            v[0,1,0], v[1,1,0], v[1,1,1], v[0,1,1]
        };
        double t[5][3][3];
        int nt = polygonize (val, 0.5, t);
        nfacets += nt;
        nverts += nt*3;
    }

    fprintf (p.fp, "\t\t<Piece NumberOfPoints=\"%i\" NumberOfPolys=\"%i\">\n", nverts, nfacets);

    fputs ("      <CellData>\n", p.fp);
    fputs ("      </CellData>\n", p.fp);

    // Write list of scalar field values to file
    fputs ("\t\t\t <PointData Normals=\"Normals\">\n", p.fp);
    for (scalar s in p.list) {
        fprintf (p.fp,"\t\t\t\t <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", s.name);
        foreach() {
	      	// Rearranging v[]
	      	double val[8] = {
	      		v[0,0,0], v[1,0,0], v[1,0,1], v[0,0,1],
	      		v[0,1,0], v[1,1,0], v[1,1,1], v[0,1,1]
	     	 };
	      	double t[5][3][3];
	      	int nt = polygonize (val, 0.5, t);
	      	for (int i = 0; i < nt; i++) {
		      	for (int j = 0; j < 3; j++) {
	      	  		coord v = {t[i][j][0], t[i][j][1], t[i][j][2]};
	      	  		fprintf (p.fp, "%g\n", interp3 (point, v, s));
	      		}
	      	}
    	}
        fputs ("\t\t\t\t </DataArray>\n", p.fp);
    }
    // Write list of vector field values to file
    for (vector ve in p.vlist) {
        fprintf (p.fp,"\t\t\t\t <DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", ve.x.name);
        foreach() {
	      	// Rearranging v[]
	      	double val[8] = {
	      		v[0,0,0], v[1,0,0], v[1,0,1], v[0,0,1],
	      		v[0,1,0], v[1,1,0], v[1,1,1], v[0,1,1]
	     	 };
	      	double t[5][3][3];
	      	int nt = polygonize (val, 0.5, t);
	      	for (int i = 0; i < nt; i++) {
		      	for (int j = 0; j < 3; j++) {
	      	  		coord v = {t[i][j][0], t[i][j][1], t[i][j][2]};
	      	  		#if dimension == 2
	                fprintf (p.fp, "%g %g 0.\n", interp3 (point, v, ve.x), interp3 (point, v, ve.y));
	                #endif
	                #if dimension == 3
	                fprintf (p.fp, "%g %g %g\n", interp3 (point, v, ve.x), interp3 (point, v, ve.y), interp3 (point, v, ve.z));
	                #endif
	      	  		//fprintf (p.fp, "%g\n", interp3 (point, v, s));
	      		}
	      	}
    	}
        fputs ("\t\t\t\t </DataArray>\n", p.fp);
    }
    fputs ("\t\t\t </PointData>\n", p.fp);


    // Write points to file
    fputs ("      <Points>\n", p.fp);
    fputs ("        <DataArray type=\"Float32\" Name=\"vertices\" NumberOfComponents=\"3\" format=\"ascii\">\n", p.fp);


    foreach(){
        //if (c[] > 1e-7 && c[] < 1. - 1e-7) {

        double val[8] = {
            v[0,0,0], v[1,0,0], v[1,0,1], v[0,0,1],
            v[0,1,0], v[1,1,0], v[1,1,1], v[0,1,1]
        };
        double t[5][3][3];
        int nt = polygonize (val, 0.5, t);
        for (int i = 0; i < nt; i++) {
            for (int j = 0; j < 3; j++) {
                coord v = {t[i][j][0], t[i][j][1], t[i][j][2]}, np;

                fprintf (p.fp, "%g %g %g\n",
                         x + v.x*Delta, y + v.y*Delta, z + v.z*Delta);
            }
        }
    }

    fputs ("        </DataArray>\n", p.fp);
    fputs ("      </Points>\n", p.fp);
    fputs ("      <Polys>\n", p.fp);

    fputs ("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n", p.fp);

    // print vert numbers
    for (int ivert = 0; ivert < nverts; ivert++)
        fprintf (p.fp, "%i ", ivert);

    fputs ("        </DataArray>\n", p.fp);
    fputs ("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n", p.fp);

    // print offsets
    for (int ifacet = 0; ifacet < nfacets; ifacet++)
        fprintf (p.fp, "%i ", ifacet*3+3);


    fputs ("        </DataArray>\n", p.fp);
    fputs ("      </Polys>\n", p.fp);
    fputs ("    </Piece>\n", p.fp);
    fputs ("  </PolyData>\n", p.fp);
    fputs ("</VTKFile>\n", p.fp);

    fflush (p.fp);
    #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
    #endif
}