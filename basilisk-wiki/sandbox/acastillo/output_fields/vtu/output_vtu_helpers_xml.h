#ifndef OUTPUT_VTU_HELPERS_XML_H
#define OUTPUT_VTU_HELPERS_XML_H

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

#endif
