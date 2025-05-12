/**
HyperTreeGrid (.htg)
This function writes basilisk scalar and vector field data to a hyperTreegrid file (.htg), which can be
loaded into the later versions of Paraview (>v5.6.0). Hypertree is an octree storage format similar to 
Basilisks, using x-ordering in oppose to Basilisks z-ordering. This means that the exported x and z axis
will swap when exporting to this format.
Using hypertree instead of vtu allows for the use of the latest hypertree filters, which support
octree meshes and among other things support the generation of seamless contour meshes.
The function does currently NOT support mpi.

Author: Arne Bockmann, Oystein Lande
Date: July, 2019

*/

int nlcount;

/** Some recursive subfunctions necessary for grid traversal */
void traverse_children_to_leaf_scalar(FILE * fp, Point point, scalar f, int lvl, int ii, int jj, int kk) {

  point.i = ii;
  point.j = jj;
  point.k = kk;

  point.level = lvl + 1;
  for (int di = 0; di < 2; di++) {
    for (int dj = 0; dj < 2; dj++) {
      for (int dk = 0; dk < 2; dk++) {
        point.i = 2*ii - GHOSTS + di;
        point.j = 2*jj - GHOSTS + dj;
        point.k = 2*kk - GHOSTS + dk;
        if is_leaf(cell) {
          fprintf (fp, "%g ", f[]); // print real value
          //fprintf (fp, "%d ", point.level); // print real value
        } else {
          fprintf (fp, "0 "); // print dummy value
        }
        nlcount++;
        if (nlcount%6 == 0)
          fprintf(fp,"\n");           
      }
    }
  }
  
  for (int di = 0; di < 2; di++) {
    for (int dj = 0; dj < 2; dj++) {
      for (int dk = 0; dk < 2; dk++) {

        // Her er de manglende tre linjene med kode som var kritisk for at traverseringen skulle bli riktig
        point.i = 2*ii - GHOSTS + di;
        point.j = 2*jj - GHOSTS + dj;
        point.k = 2*kk - GHOSTS + dk;

        if is_leaf(cell) {
        } else {
          traverse_children_to_leaf_scalar(fp, point, f, lvl+1, 2*ii-GHOSTS+di, 2*jj-GHOSTS+dj, 2*kk-GHOSTS+dk);
        }
      }
    }
  }
}

void traverse_children_to_leaf_vector(FILE * fp, Point point, vector v, int lvl, int ii, int jj, int kk) {

  point.i = ii;
  point.j = jj;
  point.k = kk;

  point.level = lvl + 1;
  for (int di = 0; di < 2; di++) {
    for (int dj = 0; dj < 2; dj++) {
      for (int dk = 0; dk < 2; dk++) {
        point.i = 2*ii - GHOSTS + di;
        point.j = 2*jj - GHOSTS + dj;
        point.k = 2*kk - GHOSTS + dk;
        if is_leaf(cell) {
          fprintf (fp, "%g %g %g ", v.x[], v.y[], v.z[]); // print real value
          //fprintf (fp, "%d ", point.level); // print real value
        } else {
          fprintf (fp, "0 0 0 "); // print dummy value
        }
        nlcount++;
        if (nlcount%6 == 0)
          fprintf(fp,"\n");           
      }
    }
  }
  
  for (int di = 0; di < 2; di++) {
    for (int dj = 0; dj < 2; dj++) {
      for (int dk = 0; dk < 2; dk++) {

        // Her er de manglende tre linjene med kode som var kritisk for at traverseringen skulle bli riktig
        point.i = 2*ii - GHOSTS + di;
        point.j = 2*jj - GHOSTS + dj;
        point.k = 2*kk - GHOSTS + dk;

        if is_leaf(cell) {
        } else {
          traverse_children_to_leaf_vector(fp, point, v, lvl+1, 2*ii-GHOSTS+di, 2*jj-GHOSTS+dj, 2*kk-GHOSTS+dk);
        }
      }
    }
  }
}
/** This is the function you want to call. Writes scalar and vector data to .htg format */
void output_htg (scalar * list, vector * vlist, FILE * fp){

  #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
  #endif

  Point point;
  
  // Define the xml header of vtkHyperTreeGrid
  fputs ("<VTKFile type=\"HyperTreeGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n",fp);
  fprintf (fp, "\t<HyperTreeGrid Dimension=\"3\" BranchFactor=\"2\" TransposedRootIndexing=\"0\" GridSize=\"1 1 1\" GridOrigin=\"0 0 0\" GridScale=\"1 1 1\">\n");
  fputs ("\t\t<Coordinates>\n",fp);
  fputs ("\t\t\t<DataArray type=\"Float64\" Name=\"XCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"1\">\n",fp);
  fprintf (fp, "\t\t\t%g %g\n",X0,L0+X0);
  fputs ("\t\t\t</DataArray>\n",fp);
  fputs ("\t\t\t<DataArray type=\"Float64\" Name=\"YCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"1\">\n",fp);
  fprintf (fp, "\t\t\t%g %g\n",Y0,L0+Y0);
  fputs ("\t\t\t</DataArray>\n",fp);
  fputs ("\t\t\t<DataArray type=\"Float64\" Name=\"ZCoordinates\" NumberOfTuples=\"2\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"1\">\n",fp);
  fprintf (fp, "\t\t\t%g %g\n",Z0,L0+Z0);
  fputs ("\t\t\t</DataArray>\n",fp);
  fputs ("\t\t</Coordinates>\n",fp);
  fputs ("\t\t<Topology>\n",fp);
  fputs ("\t\t\t<DataArray type=\"Int64\" Name=\"MaterialMaskIndex\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n",fp);
  fputs ("\t\t\t0\n",fp);
  fputs ("\t\t\t</DataArray>\n",fp);
  // Find number of cells (excluding halo cells)
  nlcount = 0;
  for (int lvl = 0; lvl <= depth(); lvl++) {
      foreach_level(lvl) {
        nlcount ++;
      }
  }
  fprintf(fp, "\t\t\t<DataArray type=\"Bit\" Name=\"Descriptor\" NumberOfTuples=\"%d\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"1\">\n",nlcount);
   
  // Set bit
  nlcount = 0;
  for (int lvl = 0; lvl <= depth(); lvl++) {
      foreach_level(lvl) {
          if (is_leaf(cell)) {
              fputs("0 ", fp);
          }
          else {
              fputs("1 ", fp);
          }
          nlcount++;
          if (nlcount%6 == 0)
              fputs("\n", fp); // line break every 6th number
      }
  }
  if (nlcount%6 != 0)
      fputs("\n", fp); // line break every 6th number

  fputs ("\t\t\t</DataArray>\n",fp);
  fputs ("\t\t</Topology>\n",fp);

  // insert pointdata arrays
  fputs ("\t\t<PointData>\n",fp);
  nlcount = 1;
  for (scalar s in list) {
    fprintf(fp,"\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n",s.name);
    fprintf(fp, "-1 ");
    traverse_children_to_leaf_scalar(fp, point, s, 0, GHOSTS, GHOSTS, GHOSTS);
    fputs ("\n\t\t\t</DataArray>\n",fp);
  }
  for (vector v in vlist) {
    char *vname = strtok(v.x.name, ".");
    fprintf (fp, "\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n",vname);
    fprintf (fp, "-1 -1 -1 ");
    traverse_children_to_leaf_vector(fp, point, v, 0, GHOSTS, GHOSTS, GHOSTS);
    fputs ("\n\t\t\t</DataArray>\n",fp);
  }
  fputs ("\t\t</PointData>\n",fp);
  fputs ("\t</HyperTreeGrid>\n",fp);
  fputs ("</VTKFile>\n",fp);

  #if defined(_OPENMP)
    omp_set_num_threads(num_omp);
  #endif
}



