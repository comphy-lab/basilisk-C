/**
# Parallel output of skeleton in a paraview compatible format


*CAREFUL* This assumes the existence of a global variable n_part.
*/

void output_skeleton (coord * loc, double * myfield, char * subname)
{
  char name[80];
  // uint64_t count = 0; => could be used for appended datas, but I didn't manage to make it work

  sprintf(name, "%s.vtp", subname);
  FILE * fp = fopen(name, "w");
  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <PolyData>\n", fp);
  fprintf (fp,"\t\t<Piece NumberOfPoints=\"%ld\" NumberOfVerts=\"1\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n",n_part);
  fputs("\t\t<PointData Scalars=\"scalar\">\n",fp);

/**
 * Here I should code something that uses a list of pointers to all the fields that I want to ouput.
 * This is not the case for now.
 *  */
    fprintf(fp,"\t\t\t<DataArray type =\"Float64\" Name =\"PID\" format=\"ascii\">\n");
      for (int i = 0; i < n_part; i++){
        fprintf (fp, "%g ", myfield[i]);
      }  
      fputs("\n",fp);
  fputs("\t\t\t</DataArray>\n",fp);
  fputs ("\t\t</PointData>\n", fp);
  fputs("\t\t<CellData/>\n",fp);
  fputs ("\t\t<Points>\n", fp);
    fprintf (fp,"\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"ascii\">\n");
    for (int i = 0; i < n_part; i++){
      fprintf (fp,"%g %g %g ",loc[i].x,loc[i].y,loc[i].z);
    }
    fputs("\n",fp);
    fputs("\t\t\t</DataArray>\n",fp);
    fputs ("\t\t</Points>\n", fp);
  fputs("\t\t<Verts>\n",fp);
  fprintf(fp,"\t\t\t<DataArray type=\"Int64\" Name= \"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%ld\">\n",n_part-1);
  for (long unsigned int i = 0; i < n_part; i++){
      fprintf (fp,"%ld ",i);
  }
  fputs("\n",fp);
  fputs("\t\t\t</DataArray>\n",fp);
  fprintf(fp,"\t\t\t<DataArray type=\"Int64\" Name= \"offsets\" format=\"ascii\" RangeMin=\"%ld\" RangeMax=\"%ld\">\n",n_part,n_part);
  fprintf(fp,"\t\t\t\t %ld\n",n_part);
  fputs("\t\t\t</DataArray>\n",fp);
  fputs("\t\t</Verts>\n",fp);
  fputs("\t</Piece>\n",fp);
  fputs ("\t </PolyData>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
  fclose(fp);
}

void output_pskeleton (coord * loc, double * myfield, char * subname){
  char name[80];
  FILE * fp;

  // Indentation should the indentation in the output file...
  if (pid() == 0) {
    sprintf(name, "%s.pvtp", subname);
    fp = fopen(name, "w");
    fputs ("<?xml version=\"1.0\"?>\n"
           "<VTKFile type=\"PPolyData\">\n",fp);
    fputs ("\t <PPolyData GhostLevel=\"0\">\n", fp);
        fputs ("\t\t\t <PCellData></PCellData>\n", fp);
        fputs ("\t\t\t <PPointData Scalars=\"scalars\">\n",fp);
          fputs("\t\t\t\t <PDataArray type =\"Float64\" Name =\"PID\" format=\"ascii\"/>\n",fp);
        fputs("\t\t\t </PPointData>\n",fp);
        fputs ("\t\t\t <PPoints>\n", fp);
          fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>\n", fp);
        fputs ("\t\t\t </PPoints>\n", fp);
      for (int i = 0; i < npe(); i++)
        fprintf (fp, "<Piece Source=\"%s_n%3.3d.vtp\"/> \n", subname, i);
      fputs ("\t </PPolyData>\n", fp);
    fputs ("</VTKFile>\n", fp);
    fflush (fp);
    fclose (fp);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  sprintf(name, "%s_n%3.3d", subname, pid());
  output_skeleton (loc, myfield, name);
}




/*
TODO: Switch to binary data with appended fields ...
*/