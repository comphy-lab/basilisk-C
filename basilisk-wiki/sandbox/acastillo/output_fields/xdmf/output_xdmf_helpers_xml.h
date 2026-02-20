#ifndef OUTPUT_XDMF_HELPERS_XML_H
#define OUTPUT_XDMF_HELPERS_XML_H

/** 
## Functions to write light data 
*/

/** ### Write the XDMF header elements to the file */ 
void write_xdmf_header(FILE *fp, const char *file_name){
  fputs("<?xml version=\"1.0\"?>\n", fp);
  fputs("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n", fp);
  fprintf(fp, "<!ENTITY HeavyData \"%s:\">\n", file_name);
  fputs("]>\n", fp);
  fputs("<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"3.0\">\n", fp);
}


/** ### Write points data array */ 
void write_xdmf_topology(FILE *fp, int dim, long num_cells, long num_points, double t) {

  fputs("\t<Domain>\n", fp);
  fputs("\t\t<Grid Name=\"Unstructured Grid\" GridType=\"Uniform\">\n", fp);
  fprintf(fp, "\t\t\t<Time Type=\"Single\" Value=\"%g\" />\n", t);

  // Write topology based on the dimension
  if (dim == 2){
    // Write 2D topology (Quadrilateral)
    fprintf(fp, "\t\t\t<Topology TopologyType=\"Quadrilateral\" NumberOfElements=\"%ld\">\n", num_cells);
    fprintf(fp, "\t\t\t\t<DataItem Format=\"HDF\" Dimensions=\"%ld 4\" DataType=\"Int\" Precision=\"8\" >\n", num_cells);
  }
  else if (dim == 3){
    // Write 3D topology (Hexahedron)
    fprintf(fp, "\t\t\t<Topology TopologyType=\"Hexahedron\" NumberOfElements=\"%ld\">\n", num_cells);
    fprintf(fp, "\t\t\t\t<DataItem Format=\"HDF\" Dimensions=\"%ld 8\" DataType=\"Int\" Precision=\"8\" >\n", num_cells);
  }

  // Write data item and close tags
  fputs("\t\t\t\t\t&HeavyData;/Topology\n", fp);
  fputs("\t\t\t\t</DataItem>\n", fp);
  fputs("\t\t\t</Topology>\n", fp);

  // Write geometry information
  fputs("\t\t\t<Geometry GeometryType=\"XYZ\">\n", fp);
  fprintf(fp, "\t\t\t\t<DataItem Format=\"HDF\" NumberType=\"Float\" Dimensions=\"%ld 3\" Precision=\"8\" >\n", num_points);
  fputs("\t\t\t\t\t&HeavyData;/Geometry/Points\n", fp);
  fputs("\t\t\t\t</DataItem>\n", fp);
  fputs("\t\t\t</Geometry>\n", fp);
}


/** ### Write attributes for scalars and vectors */ 
void write_xdmf_attributes(FILE *fp, long num_cells, scalar *slist, vector *vlist) {

  // Loop over scalars in list and write attributes
  for (scalar s in slist){
    fprintf(fp, "\t\t\t<Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", s.name);
    fprintf(fp, "\t\t\t\t<DataItem Dimensions=\"%ld\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", num_cells);
    fprintf(fp, "\t\t\t\t\t&HeavyData;/Cells/%s\n", s.name);
    fputs("\t\t\t\t</DataItem>\n", fp);
    fputs("\t\t\t</Attribute>\n", fp);
  }

  // Loop over vectors in list and write attributes
  for (vector v in vlist){
    fprintf(fp, "\t\t\t<Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Cell\">\n", v.x.name);
    fprintf(fp, "\t\t\t\t<DataItem Dimensions=\"%ld 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", num_cells);
    fprintf(fp, "\t\t\t\t\t&HeavyData;/Cells/%s\n", v.x.name);
    fputs("\t\t\t\t</DataItem>\n", fp);
    fputs("\t\t\t</Attribute>\n", fp);
  }
}

/** ### Write the XDMF footer elements to the file */ 
void write_xdmf_footer(FILE *fp) {
  // Write the closing tags for the XDMF file
  fputs("\t\t</Grid>\n", fp);
  fputs("\t</Domain>\n", fp);
  fputs("</Xdmf>\n", fp);
}

/** 
## write_xdmf_light_data(): write an `.xdmf` file containing the light data
*/

void write_xdmf_light_data(scalar *slist, vector *vlist, char *file_name, char *subname, long num_cells = 0, long num_points = 0, int dim = dimension){

#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  // Open a file for writing
  char name[111];
  sprintf(name, "%s.xmf", subname);
  FILE *fp = fopen(name, "w");

  // Write the header of the xdmf file.
  write_xdmf_header(fp, file_name);
  write_xdmf_topology(fp, dim, num_cells, num_points, t);
  write_xdmf_attributes(fp, num_cells, slist, vlist);
  write_xdmf_footer(fp);

  // Close the file
  fflush(fp);
  fclose(fp);

#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}

#endif // OUTPUT_XDMF_HELPERS_XML_H
