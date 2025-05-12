/**
# Helper functions for output_vtu.h
*/

/** ## Functions to write light data */
/** ### Write the VTU file header */
void write_vtu_header(FILE *fp, long no_points, long no_cells) {
  fputs("<?xml version=\"1.0\"?>\n", fp);
  fputs("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs("\t<UnstructuredGrid>\n", fp);
  fprintf(fp, "\t\t<Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\">\n", no_points, no_cells);
}

void write_pvtu_header(FILE *fp) {
  fputs("<?xml version=\"1.0\"?>\n", fp);
  fputs("<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs("\t<PUnstructuredGrid GhostLevel=\"0\">\n", fp);
}

/** ### Write scalar data arrays */
void write_scalar_light_data(FILE *fp, scalar *list, vector *vlist, long *count, long no_cells) {
  fputs("\t\t\t<CellData Scalars=\"scalars\">\n", fp);
  for (scalar s in list){
    fprintf(fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%ld\"/>\n", s.name, *count);
    *count += (no_cells * sizeof(double)) + sizeof(long);
  }

  for (vector v in vlist){
    fprintf(fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\"/>\n", v.x.name, *count);
    *count += (3 * no_cells * sizeof(double)) + sizeof(long);
  }
  fputs("\t\t\t</CellData>\n", fp);
}

void write_scalar_light_pdata(FILE *fp, scalar *list, vector *vlist) {
  fputs("\t\t<PCellData Scalars=\"scalars\">\n", fp);

  for (scalar s in list){
    fprintf(fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"%s\" format=\"appended\"/>\n", s.name);
  }
  for (vector v in vlist){
    fprintf(fp, "\t\t\t<PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"appended\"/>\n", v.x.name);
  }

  fputs("\t\t</PCellData>\n", fp);
}

/** ### Write points data array  */
void write_points_light_data(FILE *fp, long *count, long no_points) {
  fputs("\t\t\t<Points>\n", fp);
  fprintf(fp, "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\"/>\n", *count);
  fputs("\t\t\t</Points>\n", fp);

  *count += (3 * no_points * sizeof(double)) + sizeof(long);
}

void write_points_light_pdata(FILE *fp) {
  fputs("\t\t<PPoints>\n", fp);
  fputs("\t\t\t<PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\"/>\n", fp);
  fputs("\t\t</PPoints>\n", fp);
}


/** ### Write cells data arrays */
void write_cells_light_data(FILE *fp, long *count, long no_cells, long no_cells_offset) {
  fputs("\t\t\t<Cells>\n", fp);

  fprintf(fp, "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"%ld\"/>\n", *count);
  *count += (no_cells * sizeof(long)) + sizeof(long);

  fprintf(fp, "\t\t\t\t<DataArray type=\"Int8\" Name=\"types\" format=\"appended\" offset=\"%ld\"/>\n", *count);
  *count += (no_cells * sizeof(char)) + sizeof(long);

  fprintf(fp, "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"%ld\"/>\n", *count);
  *count += (no_cells_offset * sizeof(long)) + sizeof(long);

  fputs("\t\t\t</Cells>\n", fp);
}

/** ### Write appended data section */
void write_vtu_appended(FILE *fp) {
  fputs("\t\t</Piece>\n", fp);
  fputs("\t</UnstructuredGrid>\n", fp);
  fputs("\t<AppendedData encoding=\"raw\">\n", fp);
  fputs("_", fp);
}

/** ### Write piece references for each process */
void write_pieces_light_pdata(FILE *fp, char *subname) {
  for (int i = 0; i < npe(); i++){
    fprintf(fp, "\t\t<Piece Source=\"%s_n%3.3d.vtu\"/>\n", subname, i);
  }
  fputs("\t</PUnstructuredGrid>\n", fp);
  fputs("</VTKFile>\n", fp);
}


/** ## Functions to write heavy data */
/** ### Write scalar field data */
void write_scalar_heavy_data(FILE *fp, scalar *list, scalar per_mask, long no_cells){
  long block_len = no_cells * sizeof(double);
  for (scalar s in list){
    fwrite(&block_len, sizeof(long), 1, fp);
    foreach (){
      if (per_mask[]){
        fwrite(&val(s), sizeof(double), 1, fp);
      }
    }
  }
}

void write_scalar_heavy_data_slice(FILE *fp, scalar *list, scalar per_mask, long no_cells, coord n = {0, 0, 1}, double _alpha = 0){
  long block_len = no_cells * sizeof(double);
  for (scalar s in list){
    fwrite(&block_len, sizeof(long), 1, fp);
    foreach (){
      if (per_mask[]){
        double sval;
        if (n.x == 1)
          sval = 0.5 * (val(s) + val(s, 1, 0, 0));
        else if (n.y == 1)
          sval = 0.5 * (val(s) + val(s, 0, 1, 0));
        else
          sval = 0.5 * (val(s) + val(s, 0, 0, 1));
        fwrite(&sval, sizeof(double), 1, fp);
      }
    }
  }
}

void write_scalar_heavy_data_array(FILE *fp, long no_cells, double *pt_array_s){
  long block_len = no_cells * sizeof(double);
  fwrite(&block_len, sizeof(long), 1, fp);
  for (int i = 0; i < no_cells; i++){
    fwrite(&pt_array_s[i], sizeof(double), 1, fp);
  }
}

/** ### Write vector field data */
void write_vector_heavy_data(FILE *fp, vector *vlist, scalar per_mask, long no_cells){
  long block_len = no_cells * 3 * sizeof(double);
  for (vector v in vlist){
    fwrite(&block_len, sizeof(long), 1, fp);
    foreach (){
      if (per_mask[]){
        fwrite(&val(v.x), sizeof(double), 1, fp);
        fwrite(&val(v.y), sizeof(double), 1, fp);
        #if dimension == 2
          double vz = 0;
          fwrite(&vz, sizeof(double), 1, fp);
        #elif dimension == 3
          fwrite(&val(v.z), sizeof(double), 1, fp);
        #endif
      }
    }
  }
}

void write_vector_heavy_data_slice(FILE *fp, vector *vlist, scalar per_mask, long no_cells, coord n = {0, 0, 1}, double _alpha = 0){
  long block_len = no_cells * 3 * sizeof(double);
  for (vector v in vlist){
    fwrite(&block_len, sizeof(long), 1, fp);
    foreach (){
      if (per_mask[]){
        double xval, yval, zval;
        if (n.x == 1){
          xval = 0.5 * (val(v.x) + val(v.x, 1, 0, 0));
          yval = 0.5 * (val(v.y) + val(v.y, 1, 0, 0));
          #if dimension == 3
            zval = 0.5 * (val(v.z) + val(v.z, 1, 0, 0));
          #else
            zval = 0;
          #endif
        }
        else if (n.y == 1){
          xval = 0.5 * (val(v.x) + val(v.x, 0, 1, 0));
          yval = 0.5 * (val(v.y) + val(v.y, 0, 1, 0));
          #if dimension == 3
            zval = 0.5 * (val(v.z) + val(v.z, 0, 1, 0));
          #else
            zval = 0;
          #endif
        }
        else {
          xval = 0.5 * (val(v.x) + val(v.x, 0, 0, 1));
          yval = 0.5 * (val(v.y) + val(v.y, 0, 0, 1));
          #if dimension == 3
            zval = 0.5 * (val(v.z) + val(v.z, 0, 0, 1));
          #else
            zval = 0;
          #endif
        }
        fwrite(&xval, sizeof(double), 1, fp);
        fwrite(&yval, sizeof(double), 1, fp);
        fwrite(&zval, sizeof(double), 1, fp);
      }
    }
  }
}

/** ### Write points data */
void write_points_heavy_data(FILE *fp, long no_points) {
  long block_len = no_points * 3 * sizeof(double);
  fwrite(&block_len, sizeof(long), 1, fp);
  foreach_vertex(){
    fwrite(&x, sizeof(double), 1, fp);
    fwrite(&y, sizeof(double), 1, fp);
    fwrite(&z, sizeof(double), 1, fp);
  }
}

void write_points_heavy_data_slice(FILE *fp, long no_points, coord n = {0,0,1}, double _alpha = 0) {
  long block_len = no_points * 3 * sizeof(double);
  fwrite(&block_len, sizeof(long), 1, fp);
  foreach_vertex(){
    shortcut_slice(n, _alpha);
    fwrite(&x, sizeof(double), 1, fp);
    fwrite(&y, sizeof(double), 1, fp);
    fwrite(&z, sizeof(double), 1, fp);
  }
}

void write_points_heavy_data_array(FILE *fp, long no_points, double *pt_array_x, double *pt_array_y, double *pt_array_z) {
  long block_len = no_points * 3 * sizeof(double);
  fwrite(&block_len, sizeof(long), 1, fp);
  for (int i = 0; i < no_points; i++){
    fwrite(&pt_array_x[i], sizeof(double), 1, fp);
    fwrite(&pt_array_y[i], sizeof(double), 1, fp);
    #if dimension == 2
      double vz = 0;
      fwrite(&vz, sizeof(double), 1, fp);
    #elif dimension == 3
      fwrite(&pt_array_z[i], sizeof(double), 1, fp);
    #endif
  }
}

/** ### Write cell offsets */
void write_cell_offsets(FILE *fp, long no_cells, char noffset) {
  long block_len = no_cells * sizeof(long);
  fwrite(&block_len, sizeof(long), 1, fp);
  for (int i = 0; i < no_cells; i++){
    long offset = (i + 1) * noffset;
    fwrite(&offset, sizeof(int64_t), 1, fp);
  }
}

void write_cell_offsets2(FILE *fp, long nfacets, long *offsets) {
  long block_len = nfacets * sizeof(long);
  fwrite(&block_len, sizeof(long), 1, fp);
  for (int ii = 0; ii < nfacets; ii++)
    fwrite(&offsets[ii], sizeof(long), 1, fp);
}

void write_cell_offsets3(FILE *fp, long no_cells) {
  long block_len = no_cells * sizeof(long);
  fwrite(&block_len, sizeof(long), 1, fp);
  for (long i = 0; i < no_cells; i++){
    fwrite(&i, sizeof(int64_t), 1, fp);
  }
}

/** ### Write cell types */
void write_cell_types(FILE *fp, long no_cells, char type) {
  long block_len = no_cells * sizeof(char);
  fwrite(&block_len, sizeof(long), 1, fp);
  for (int i = 0; i < no_cells; i++){
    fwrite(&type, sizeof(char), 1, fp);
  }
}

/** ### Write cell connectivity */
void write_cell_connectivity(FILE *fp, vertex scalar marker, scalar per_mask, long no_cells, char noffset) {
  long block_len = no_cells * noffset * sizeof(long);
  fwrite(&block_len, sizeof(long), 1, fp);
  foreach (serial, noauto){
    if (per_mask[]){
      long connectivity[noffset];
      connectivity[0] = (long)marker[];
      connectivity[1] = (long)marker[1];
      connectivity[2] = (long)marker[1, 1];
      connectivity[3] = (long)marker[0, 1];
#if dimension == 3
      connectivity[4] = (long)marker[0, 0, 1];
      connectivity[5] = (long)marker[1, 0, 1];
      connectivity[6] = (long)marker[1, 1, 1];
      connectivity[7] = (long)marker[0, 1, 1];
#endif
      fwrite(connectivity, sizeof(long), noffset, fp);
    }
  }
}

void write_cell_connectivity_slice(FILE *fp, vertex scalar marker, scalar per_mask, long no_cells, char noffset, coord n = {0,0,1} ) {
  long block_len = no_cells * noffset * sizeof(long);
  fwrite(&block_len, sizeof(long), 1, fp);
  foreach (serial, noauto){
    if (per_mask[]){
      long connectivity[noffset];
      if (n.x == 1){
        connectivity[0] = (long)marker[1, 0, 0];
        connectivity[1] = (long)marker[1, 1, 0];
        connectivity[2] = (long)marker[1, 1, 1];
        connectivity[3] = (long)marker[1, 0, 1];
      }
      else if (n.y == 1){
        connectivity[0] = (long)marker[0, 1, 0];
        connectivity[1] = (long)marker[1, 1, 0];
        connectivity[2] = (long)marker[1, 1, 1];
        connectivity[3] = (long)marker[0, 1, 1];
      }
      else{
        connectivity[0] = (long)marker[0, 0, 1];
        connectivity[1] = (long)marker[1, 0, 1];
        connectivity[2] = (long)marker[1, 1, 1];
        connectivity[3] = (long)marker[0, 1, 1];
      }
      fwrite(connectivity, sizeof(long), noffset, fp);
    }
  }
}

/**
## Other functions

### Count the number of vertices and facets */

void count_vertices_and_facets(scalar c, long *nverts, long *nfacets) {
  foreach (serial, noauto){
    #if EMBED
    if ((c[] > 1e-6 && c[] < 1. - 1e-6) && cs[] == 1)
    #else
    if (c[] > 1e-6 && c[] < 1. - 1e-6)
    #endif
    {
      shortcut_facets
      for (int i = 0; i < m; i++)(*nverts)++;
      if (m > 0)
        (*nfacets)++;
    }
  }
}

/** ### Populate vertex coordinates and facet offsets */
void populate_vertex_and_offset_arrays(scalar c, long nverts, long nfacets, double *pt_array_x, double *pt_array_y, double *pt_array_z, long *offsets) {
  long iverts = 0, ifacet = 0, offset = 0;
  foreach (serial, noauto){
    #if EMBED
    if ((c[] > 1e-6 && c[] < 1. - 1e-6) && cs[] == 1)
    #else
    if (c[] > 1e-6 && c[] < 1. - 1e-6)
    #endif
    {
      shortcut_facets // we cycle if cell is not at the interface

      // Calculate and store vertex coordinates
      coord _p = {x, y, z};
      for (int i = 0; i < m; i++){
        pt_array_x[iverts] = _p.x + v[i].x * Delta;
        pt_array_y[iverts] = _p.y + v[i].y * Delta;
        #if dimension == 3
          pt_array_z[iverts] = _p.z + v[i].z * Delta;
        #endif
        iverts++;
      }
      // Store facet offset if there are vertices in the facet
      if (m > 0){
        offset += m;
        offsets[ifacet] = offset;
        ifacet++;
      }
    }
  }
}

/** ### Populate arrays with values at the facets */
void populate_facet_arrays(scalar c, scalar s, long nverts, long nfacets, double *val_array_s) {
  long ifacet = 0;
  foreach (serial, noauto){
    #if EMBED
    if ((c[] > 1e-6 && c[] < 1. - 1e-6) && cs[] == 1)
    #else
    if (c[] > 1e-6 && c[] < 1. - 1e-6)
    #endif
    {
      shortcut_facets // we cycle if cell is not at the interface

      if (m > 0){
        val_array_s[ifacet] = s[];
        ifacet++;
      }
    }
  }
}
