/**
# Wrapper for output functions with Paraview

File copy-pasted (with minor modifications by Can Celçuk and Guodong Gai) from [Oystein Lande's sandbox](http://basilisk.fr/sandbox/oystelan/output_vtu_foreach.h). All credit goes to Oystein Lande, Can Selçuk and Guodong Gai.

*/
/*Here we defined the directory of the results*/
# define result_dir "Res" 
# define result_fluid_rootfilename "fluid_basilisk"



/**
# output_pvtu_bin
This function writes one XML file which allows to read the *.vtu files generated
by output_vtu_bin_foreach() when used in MPI. Tested in (quad- and oct-)trees
using MPI.
*/
void output_pvtu_bin (scalar * list, vector * vlist,  FILE * fp, char * subname)
{
  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
  fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"appended\">\n", s.name);
    fputs ("\t\t\t\t </PDataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"appended\">\n", v.x.name);
    fputs ("\t\t\t\t </PDataArray>\n", fp);
  }
  fputs ("\t\t\t </PCellData>\n", fp);
  fputs ("\t\t\t <PPoints>\n", fp);
  fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
  fputs ("\t\t\t\t </PDataArray>\n", fp);
  fputs ("\t\t\t </PPoints>\n", fp);

  for (int i = 0; i < npe(); i++)
    fprintf (fp, "<Piece Source=\"%s_%d.vtu\"/> \n", subname, i);

  fputs ("\t </PUnstructuredGrid>\n", fp);
  fputs ("</VTKFile>\n", fp);
}

/**
# output_vtu_bin_foreach
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on binary format. If one writes
one *.vtu file per PID process this function may be combined with
output_pvtu_bin() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).
Bug correction: %g turns into scientific notation for high integer values. This is 
not supported by paraview. Hence a fix was needed.
Oystein Lande 2017
*/
void output_vtu_bin_foreach (scalar * list, vector * vlist, FILE * fp, bool linear)
{
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  vertex scalar marker[];
  int no_cells = 0;
  int no_points = 0;
 
  foreach_cache(tree->vertices){
    marker[] = _k;
    no_points += 1;
  }
  foreach_cache(tree->leaves){
   no_cells += 1;
  }

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);

  int count = 0;
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%d\">\n", s.name,count);
    count += ((no_cells)+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"Vect-%s\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", v.x.name,count);
    count += ((no_cells*3)+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n",count);
  count += ((no_points*3)+1)*8;
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
  foreach(){
#if dimension == 2
// Edit OLAND:
    // %g will turn into scientific notation for high integer values. this does not work with paraview
    //fprintf (fp, "%g %g %g %g \n", marker[], marker[1,0], marker[1,1], marker[0,1]);
    int ape1 = marker[];
    int ape2 = marker[1,0];
    int ape3 = marker[1,1];
    int ape4 = marker[0,1];
    fprintf (fp, "%u %u %u %u \n", ape1, ape2, ape3, ape4);
#endif
#if dimension == 3
    // Edit OLAND:
    // %g will turn into scientific notation for high integer values. this does not work with paraview
    //fprintf (fp, "%g %g %g %g %g %g %g %g\n", marker[], marker[1,0,0], marker[1,1,0], marker[0,1,0],marker[0,0,1], marker[1,0,1], marker[1,1,1], marker[0,1,1]);
    int ape1 = marker[];
    int ape2 = marker[1,0,0];
    int ape3 = marker[1,1,0];
    int ape4 = marker[0,1,0];
    int ape5 = marker[0,0,1];
    int ape6 = marker[1,0,1];
    int ape7 = marker[1,1,1];
    int ape8 = marker[0,1,1];
    fprintf (fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, ape6, ape7, ape8);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);
  for (int i = 1; i < no_cells+1; i++){
#if dimension == 2
    fprintf (fp, "%d \n", i*4);
#endif
#if dimension == 3
    fprintf (fp, "%d \n", i*8);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
  foreach(){
#if dimension == 2
    fputs ("9 \n", fp);
#endif
#if dimension == 3
    fputs ("12 \n", fp);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
  fputs ("_", fp);
  unsigned long long block_len=no_cells*8;
#if dimension == 2
  double vz=0;
#endif
  for (scalar s in list) {
    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
    foreach()
      fwrite (&val(s), sizeof (double), 1, fp);
  }
  block_len=no_cells*8*3;
  for (vector v in vlist) {
    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
    foreach(){
      fwrite (&val(v.x), sizeof (double), 1, fp);
      fwrite (&val(v.y), sizeof (double), 1, fp);
#if dimension == 2
      fwrite (&vz, sizeof (double), 1, fp);
#endif
#if dimension == 3
      fwrite (&val(v.z), sizeof (double), 1, fp);
#endif
    }
  }
  block_len=no_points*8*3;
  fwrite (&block_len, sizeof (unsigned long long), 1, fp);
  foreach_vertex(){
    fwrite (&x, sizeof (double), 1, fp);
    fwrite (&y, sizeof (double), 1, fp);
    fwrite (&z, sizeof (double), 1, fp);
  }
  fputs ("\t\n", fp);
  fputs ("\t </AppendedData>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}




//----------------------------------------------------------------------------
void output_pvd( FILE * fp, char const* times_series )
//----------------------------------------------------------------------------
{
  fputs( "<?xml version=\"1.0\"?>\n"
  	"<VTKFile type=\"Collection\" version=\"1.0\" "
	"byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp );
  fputs( "<Collection>\n", fp );
  fputs( times_series, fp );
  fputs( "</Collection>\n", fp );
  fputs( "</VTKFile>\n", fp );
}




//----------------------------------------------------------------------------
void save_data( scalar * list, vector * vlist, double const time )
//----------------------------------------------------------------------------
{

  char vtk_times_series[100000] = "";
  static int cycle_number = 0;
  //if ( !cycle_number ) cycle_number = init_cycle_number;//ggd

  FILE * fpvtk;
  char filename_vtu[80] = "";
  char filename_pvtu[80] = "";
  char suffix[80] = "";

  // Write the VTU file
  sprintf( filename_vtu, "%s", result_dir );
  strcat( filename_vtu, "/" );
  strcat( filename_vtu, result_fluid_rootfilename );

  sprintf( suffix, "_T%d_%d.vtu", cycle_number, pid());
  strcat( filename_vtu, suffix );

  fpvtk = fopen( filename_vtu, "w" );
  output_vtu_bin_foreach( list, vlist, fpvtk, false );
  fclose( fpvtk );

  // Write the PVTU file
  if ( pid() == 0 )
  {
    sprintf( filename_pvtu, "%s", result_dir );
    strcat( filename_pvtu, "/" );
    strcat( filename_pvtu, result_fluid_rootfilename );
    sprintf( suffix, "_T%d.pvtu", cycle_number );
    strcat( filename_pvtu, suffix );

    fpvtk = fopen( filename_pvtu, "w" );

    sprintf ( filename_vtu, "%s", result_fluid_rootfilename );
    sprintf( suffix, "_T%d", cycle_number );
    strcat( filename_vtu, suffix );
    output_pvtu_bin ( list, vlist, fpvtk, filename_vtu );

    fclose( fpvtk );
  }

  // Write the PVD file
  if ( pid() == 0 )
  {
    char filename_pvd[80] = "";
    sprintf( filename_pvd, "%s", result_dir );
    strcat( filename_pvd, "/" );
    strcat( filename_pvd, result_fluid_rootfilename );
    strcat( filename_pvd, ".pvd" );

    fpvtk = fopen( filename_pvd, "w" );

    char time_line[200] = "";
    strcpy( time_line, "<DataSet timestep=" );
    sprintf( suffix, "\"%.4e\"", time );
    strcat( time_line, suffix );
    strcat( time_line, " group=\"\" part=\"0\" file=\"" );
    strcpy( filename_pvtu, result_fluid_rootfilename );
    sprintf( suffix, "_T%d.pvtu", cycle_number );
    strcat( filename_pvtu, suffix );
    strcat( time_line, filename_pvtu );
    strcat( time_line, "\"/>\n" );
    strcat( vtk_times_series, time_line );
    output_pvd( fpvtk, vtk_times_series );

    fclose( fpvtk );
  }

  // Write the last cycle number in a file for restart
  if ( pid() == 0 )
  {
    char filename_lcn[256] = "";
    sprintf( filename_lcn, "%s", result_dir );
    strcat( filename_lcn, "/" );
    strcat( filename_lcn, result_fluid_rootfilename );
    strcat( filename_lcn, "_lcn_vtk.txt" );

    fpvtk = fopen( filename_lcn, "w" );

    fprintf( fpvtk, "%d\n", cycle_number );

    fclose( fpvtk );
  }

  ++cycle_number;
}

/*
//----------------------------------------------------------------------------
void reinitialize_vtk_restart( void )
//----------------------------------------------------------------------------
{
  // Get the last cycle cumber from previous simulation
  char filename_lcn[80] = "";
  sprintf( filename_lcn, "%s", result_dir );
  strcat( filename_lcn, "/" );
  strcat( filename_lcn, result_fluid_rootfilename );
  strcat( filename_lcn, "_lcn_vtk.txt" );

  FILE * fpvtk = fopen( filename_lcn, "r" );

  fscanf ( fpvtk, "%d", &init_cycle_number );
  ++init_cycle_number;

  fclose( fpvtk );

  // Re-initialize the time output series string
  if ( pid() == 0 )
  {
    char filename_pvd[80] = "";
    char time_line[256] = "";
    char start[20] = "";
    char start_ref[20] = "<DataSet";
    sprintf( filename_pvd, "%s", result_dir );
    strcat( filename_pvd, "/" );
    strcat( filename_pvd, result_fluid_rootfilename );
    strcat( filename_pvd, ".pvd" );

    fpvtk = fopen( filename_pvd, "r" );

    while ( fgets( time_line, sizeof(time_line), fpvtk ) )
    {
      // Extract 8 first characters
      strncpy( start, time_line, 8 );

      // If 8 first characters equal "<DataSet", it is an output time line
      // We add to the vtk time series string
      if ( ! strcmp( start, start_ref ) )
        strcat( vtk_times_series, time_line );
    }

    fclose( fpvtk );
  }
}
*/


