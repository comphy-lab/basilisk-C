char *folder_name = NULL;
FILE *log_file = NULL;

/** 
Checks a condition and, if false, prints an error message and ends execution if requested
(Similar to the standard "assert" macro). */
void ErrorMessage(int Condition, const char *Message, const char *FunctionName, char EndExecution)
{   
  // If the condition is satisfied, we do nothing
  if( Condition )
    return;

  printf("  ===== ERROR in function %s ...\n", (FunctionName) ? FunctionName : "NOT-SPECIFIED");
  printf("  ===== %s\n", Message ? Message : "");
  if( EndExecution ) {
    printf("  ===== Ending execution...\n\n");
    exit(0);
  }

  return;
}

/**
Returns a command line argument given by the user when calling the program
Input Parameter: the index (starting from zero) of the argument to be returned
Output Parameter: string location to store the argument. */
void GetCommandLineArgument(int ArgumentIndex, char *ReturnArgument, int argc, char *argv[])
{   
  /// Checking if the user has provided the argument with the requested index
  if( argc<=ArgumentIndex+1 ) {
    printf("  ===== ERROR in function GetCommandLineArgument ...\n");
    printf("  ===== The requested parameter %d was not provided when calling the program.\n", ArgumentIndex);
    printf("  ===== Ending execution...\n\n");
    exit(0);
  }

  /// Copying the argument into the return string
  strcpy(ReturnArgument, argv[ArgumentIndex + 1]);
  return;
}

/**
Creates a folder for this simulation
All the outputs from functions "PrintLog", "PrintMeshVTK" and 
"PrintInterfaceVTK" will go into this folder

Parameter: the name of the folder to be created. You can use the same syntax as "printf" 
to format this folder name */
void OpenSimulationFolder(const char* format, ...)
{
  // === Only root process does this
  if( pid() )
    return;

  char *folder_temp = (char *)malloc(1000*sizeof(char));
  folder_name = (char *)malloc(1100*sizeof(char));

  // === Formatting the string with the function parameters
  va_list argptr;
  va_start(argptr, format);
  vsprintf(folder_temp, format, argptr);
  va_end(argptr);

  // === Creating the base "outputs" folder
  system("mkdir outputs/");

  // === Creating the simulation folder
  char command[1000];
  sprintf(folder_name, "outputs/%s/", folder_temp);
  sprintf(command, "mkdir %s", folder_name);
  system(command);

  /// === Opening the new log file
  if( log_file ) 
    fclose(log_file);
  sprintf(command, "%s/log_file.txt", folder_name);
  log_file = fopen(command, "wt");
  if( !log_file ) {
    printf("  ===== ERROR in function OpenSimulationFolder ...\n");
    printf("  ===== Unable to create the folder and log_file. Do you have permission for creating this folder?\n");
    printf("  ===== Ending execution...\n\n");
    exit(0);
  }
  
  free(folder_temp);
  return;
}

/**
Closes the current simulation folder and output files currently open. */
void CloseSimulationFolder() {

  // === Only root process does this
  if( pid() )
    return;

  if( log_file )
    fclose(log_file);
  log_file = NULL;

  return;
}


/**
Read parameters from a file
Parameters: 
1. The name of the file containing the parameters.
2. A list of parameter names to look for in the file. Separated by semi-colon.
3. Pointers to variables where each of the parameters will be stored
NOTE: currently I'm reading all values as "double" variables. 
I can change this in future if there is interest for it */
// void ReadParametersFromFile(const char *file_name, const char *parameters_names, ...)
// {
//   char names_copy[1000], *token;
//   int number_of_parameters = 0;
//   FILE *file;

//   // == Counting how many parameters the user wants to read
//   strcpy(names_copy, parameters_names);
//   token = strtok(names_copy, ";");
//   while( token ) {
//     number_of_parameters++;
//     token = strtok(NULL, ";");
//   }

//   // == Looking for the parameters in the text file
//   if( !(file=fopen(file_name, "rt")) ) {
//     printf("  ===== ERROR in function ReadParametersFromFile ...\n");
//     printf("  ===== Unable to open the input text file \"%s\".\n", file_name);
//     printf("  ===== Ending execution...\n\n");
//     exit(0);
//   }
//   strcpy(names_copy, parameters_names);
//   token = strtok(names_copy, ";");
//   int i_token = 0;
//   while( token ) {
//     char name_read[1000];
//     double value_read;
    
//     // Scanning the entire file from the beginning loking for this token
//     rewind(file);
//     int found_parameter = 0;
//     while( 2==fscanf(file, "%s %lf\n", name_read, &value_read) ) {
//       if( !strcmp(name_read, token) ) {
//         found_parameter = 1;
//         break;
//       }
//     }
//     if( !found_parameter ) {
//       printf("  ===== ERROR in function ReadParametrsFromFile ...");
//       printf("  ===== Could not find the parameter [%s] in the input file.\n", token);
//       printf("  ===== Ending execution...\n\n");
//       exit(0);
//     }

//     // Finding the pointer in the list of pointers given by the user in the function parameters
//     va_list ap;
//     va_start(ap, parameters_names);
//     int i;
//     double *pointer;
//     for( i=0; i<=i_token; i++ )
//       pointer = va_arg(ap, double*);
//     va_end(ap);

//     // Setting the value read in the pointer
//     *pointer = value_read;

//     token = strtok(NULL, ";");
//     i_token++;
//   }
//   fclose( file );

//   return;
// }

/**
Prints something to the log file of this simulation */
void PrintLog(const char* format, ...)
{
  // === Only root process will print
  if( pid() )
    return;
  
  if( log_file==NULL )
    return;

  // Formatted printing
  va_list argptr;
  va_start(argptr, format);
  vfprintf(log_file, format, argptr);
  va_end(argptr);
  fflush(log_file);
  
  return;
}



/**
This is basically a wrapper for the MPI_Gatherv function.
It automatically counts how much data will be sent by each processor */
#if _MPI
void MPI_Gather_Uneven(const void *sendbuf, int sendcount, MPI_Datatype sendreceivetype, void *recvbuf, int root)
{
  int array_counts[npe()];
  int displs[npe()];

  if( pid()==root ) {
    int i;
    for( i=0; i<npe(); i++ ) {
      if( i==root )
        array_counts[i] = sendcount;
      else
        MPI_Recv(array_counts+i, 1, MPI_INT, i, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      displs[i] = (i==0) ? 0 : displs[i-1] + array_counts[i-1];
    }
  }
  else
    MPI_Send(&sendcount, 1, MPI_INT, 0, 1000, MPI_COMM_WORLD);

  MPI_Gatherv(sendbuf, sendcount, sendreceivetype, recvbuf, array_counts, displs, sendreceivetype, root, MPI_COMM_WORLD);
}
#endif

// #define PRINT_CELL_CRITERIA 1 // Printing all the cells in the mesh
#define PRINT_CELL_CRITERIA f[]>0.05 // Only printing cells with some volume in it

/**
Prints the mesh and (optionally) scalar data to an ASCII VTK file. */
void PrintMeshVTK_ASCII(int n, double time, scalar *list_scalar_data, const char **list_scalar_names)
{
    FILE *arq;
    char nomeArq[900];

    #if dimension==2
        // Each cell is a square (4 vertices)
        int vertices_per_cell = 4; 

        // VTK cell code that represents quads
        // (check Figure 2 here to understand: https://kitware.github.io/vtk-examples/site/VTKFileFormats/)
        int vtk_cell_type = 9; 
    #else
        int vertices_per_cell = 8; // Each cell is a cube

        // VTK cell code that represents voxels
        // (check Figure 2 here to understand: https://kitware.github.io/vtk-examples/site/VTKFileFormats/)
        int vtk_cell_type = 11; 
    #endif

    // Counting how many local cells we have in the mesh
    // NOTE: whenever I say "local", i mean things that are locally in this processor (in case of parallel simulation)
    int local_num_cells = 0;
    foreach(serial) {
        // if( PRINT_CELL_CRITERIA )
            local_num_cells++;
    }

    /// === Allocatting memory for the local vertex arrays
    double *local_vertices_x = (double *)malloc( vertices_per_cell*local_num_cells*sizeof(double) );
    double *local_vertices_y = (double *)malloc( vertices_per_cell*local_num_cells*sizeof(double) );
    double *local_vertices_z = (double *)malloc( vertices_per_cell*local_num_cells*sizeof(double) );

    /// === Allocating memory for ALL the local scalar data arrays
    int number_of_scalar_fields = list_len(list_scalar_data);
    typedef double** doublepp; // Hiding the double pointer with a typedef because qcc gets really annoying if i dont do this (why??)
    doublepp local_data = number_of_scalar_fields ? (double **)malloc( number_of_scalar_fields*sizeof(double *) ) : NULL;
    for( int k=0; k<number_of_scalar_fields; k++ )
        local_data[k] = (double *)malloc( local_num_cells*sizeof(double) );

    /// === Storing all the vertices coordinates in the arrays and also the scalar data to be printed data
    int i = 0, i_cell = 0;
    foreach(serial) {
        if( PRINT_CELL_CRITERIA ) {
            #if dimension==2
                local_vertices_x[i] = x-0.5*Delta; local_vertices_y[i] = y-0.5*Delta; local_vertices_z[i++] = 0.0;
                local_vertices_x[i] = x+0.5*Delta; local_vertices_y[i] = y-0.5*Delta; local_vertices_z[i++] = 0.0;
                local_vertices_x[i] = x+0.5*Delta; local_vertices_y[i] = y+0.5*Delta; local_vertices_z[i++] = 0.0;
                local_vertices_x[i] = x-0.5*Delta; local_vertices_y[i] = y+0.5*Delta; local_vertices_z[i++] = 0.0;
            #else // dimension==3
                local_vertices_x[i] = x-0.5*Delta; local_vertices_y[i] = y-0.5*Delta; local_vertices_z[i++] = z-0.5*Delta;
                local_vertices_x[i] = x+0.5*Delta; local_vertices_y[i] = y-0.5*Delta; local_vertices_z[i++] = z-0.5*Delta;
                local_vertices_x[i] = x-0.5*Delta; local_vertices_y[i] = y+0.5*Delta; local_vertices_z[i++] = z-0.5*Delta;
                local_vertices_x[i] = x+0.5*Delta; local_vertices_y[i] = y+0.5*Delta; local_vertices_z[i++] = z-0.5*Delta;

                local_vertices_x[i] = x-0.5*Delta; local_vertices_y[i] = y-0.5*Delta; local_vertices_z[i++] = z+0.5*Delta;
                local_vertices_x[i] = x+0.5*Delta; local_vertices_y[i] = y-0.5*Delta; local_vertices_z[i++] = z+0.5*Delta;
                local_vertices_x[i] = x-0.5*Delta; local_vertices_y[i] = y+0.5*Delta; local_vertices_z[i++] = z+0.5*Delta;
                local_vertices_x[i] = x+0.5*Delta; local_vertices_y[i] = y+0.5*Delta; local_vertices_z[i++] = z+0.5*Delta;
            #endif

            int list_index = 0;
            for(scalar s in list_scalar_data)
                local_data[list_index++][i_cell] = s[];

            i_cell++;
        }
    }

    /// === Getting how many cells we have in total between all processes
    /// === And then gathering all the local vertices into the root process
    int total_num_cells;
    double *vertices_x = NULL, *vertices_y = NULL, *vertices_z = NULL;
    doublepp data = NULL;
    #if _MPI
        /// Total number of cels
        MPI_Allreduce(&local_num_cells, &total_num_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        /// === Gathering all the local vertices into the root process
        if( pid()==0 ) {
            vertices_x = (double *)malloc( total_num_cells*vertices_per_cell*sizeof(double) );
            vertices_y = (double *)malloc( total_num_cells*vertices_per_cell*sizeof(double) );
            vertices_z = (double *)malloc( total_num_cells*vertices_per_cell*sizeof(double) );
            
            data = number_of_scalar_fields ? (double **)malloc( number_of_scalar_fields*sizeof(double*) ) : NULL;
            for(int k=0; k<number_of_scalar_fields; k++)
                data[k] = (double *)malloc( total_num_cells*sizeof(double) );
        }
        MPI_Gather_Uneven(local_vertices_x, vertices_per_cell*local_num_cells, MPI_DOUBLE, vertices_x, 0);
        MPI_Gather_Uneven(local_vertices_y, vertices_per_cell*local_num_cells, MPI_DOUBLE, vertices_y, 0);
        MPI_Gather_Uneven(local_vertices_z, vertices_per_cell*local_num_cells, MPI_DOUBLE, vertices_z, 0);

        int list_index = 0;
        for( list_index=0; list_index<number_of_scalar_fields; list_index++ )
            MPI_Gather_Uneven(local_data[list_index], local_num_cells, MPI_DOUBLE, (pid()==0) ? data[list_index] : NULL, 0);

        /// === Releasing local memory
        free(local_vertices_x);
        free(local_vertices_y);
        free(local_vertices_z);

        for(int k=0; k<number_of_scalar_fields; k++)
            free(local_data[k]);
        if( local_data )
            free(local_data);

    #else
        total_num_cells = local_num_cells;
        vertices_x = local_vertices_x;
        vertices_y = local_vertices_y;
        vertices_z = local_vertices_z;
        data = local_data;
    #endif

    
    /// === From now on, we only do file-printing stuff. Only the root process will do it
    if( pid()==0 ) {

        /// === Opening the file to print the mesh
        sprintf(nomeArq, "%s/Mesh-N%d.vtk", folder_name, n);
        arq = fopen(nomeArq, "wt");


        /// === Printing the VTK header information
        fprintf(arq, "# vtk DataFile Version 2.0\n");
        fprintf(arq, "MESH. step %d time %lf\n", n, time);
        fprintf(arq, "ASCII\n");
        fprintf(arq, "DATASET UNSTRUCTURED_GRID\n");

        /// === Printing all the vertices coordinates
        fprintf(arq, "POINTS %d float\n", total_num_cells*vertices_per_cell);
        int total_num_vertices = vertices_per_cell*total_num_cells;
        for( i=0; i<total_num_vertices; i++ )
            fprintf(arq, "%lf %lf %lf\n", vertices_x[i], vertices_y[i], vertices_z[i]);

        /// === Printing all the cells (each cell contains 4 (or 8) indices referring to the vertices printed above)
        fprintf(arq, "CELLS %d %d\n", total_num_cells, (vertices_per_cell + 1)*total_num_cells);
        int index = 0;
        for( i=0; i<total_num_cells; i++ ) {
            #if dimension==2
                fprintf(arq, "%d %d %d %d %d\n", vertices_per_cell, index, index+1, index+2, index+3);
            #else        
                fprintf(arq, "%d %d %d %d %d %d %d %d %d\n", vertices_per_cell, index, index+1, index+2, index+3, index+4, index+5, index+6, index+7);
            #endif
            index += vertices_per_cell;
        }

        /// === Printing the VTK cell_types (quads or voxels)
        fprintf(arq, "CELL_TYPES %d\n", total_num_cells);
        for( i=0; i<total_num_cells; i++ )
            fprintf(arq, "%d\n", vtk_cell_type);

        

        /// === Printing the actual simulation data that is stored in the cells
        if(data) {
            int list_index = 0;
            fprintf(arq, "CELL_DATA %d\n", total_num_cells);
            for( scalar s in list_scalar_data ) {
                fprintf(arq, "SCALARS %s float 1\n", list_scalar_names[list_index]);
                fprintf(arq, "LOOKUP_TABLE default\n");
                for( i=0; i<total_num_cells; i++ )
                    fprintf(arq, "%lf\n", data[list_index][i]);
                list_index++;
            }
        }


        /// === Closing the file
        fclose(arq);

        /// === Releasing memory
        free(vertices_x);
        free(vertices_y);
        free(vertices_z);
        for(int k=0; k<number_of_scalar_fields; k++)
            free(data[k]);
        if(data)
            free(data);
    }

    /// === All processes wait for the file printing to finish before continuing with the simulation
    #if _MPI
        MPI_Barrier(MPI_COMM_WORLD);
    #endif

    return;
}

// This function transforms the little-endian entries of an array into big-endian (or vice-versa). This will be used in the BINARY vtk printing functions
// IMPORTANT: currently, this function will also swap the order of the array entries as well as the byte ordering of each individual entry
void SwapArrayBytes(void *Array, size_t number_of_bytes)
{
    int i;
    size_t half_of_bytes = number_of_bytes/2;

    unsigned char *byte_array = (unsigned char *)Array;
    for( i=0; i<half_of_bytes; i++ ) {
        unsigned char temporary = byte_array[i];
        byte_array[i] = byte_array[number_of_bytes - i - 1];
        byte_array[number_of_bytes - i - 1] = temporary;
    }

    return;
}

/**
Prints mesh and (optionally) scalars to a double-precision BINARY VTK file. */
void PrintMeshVTK_Binary_Double(int n, double time, scalar *list_scalar_data, const char **list_scalar_names)
{
  FILE *arq;
  char nomeArq[900];

  // === Cell is either a square (2D) or cube (3D), meaning either 4 or 8 vertices per cell
  int vertices_per_cell = (dimension==2) ? 4 : 8;

  // === VTK cell code that represents voxels
  // === (check Figure 2 here to understand: https://kitware.github.io/vtk-examples/site/VTKFileFormats/)
  int vtk_cell_type = (dimension==2) ? 9 : 11;

  // === Counting how many local cells we have in the mesh
  // === NOTE: whenever I say "local", i mean things that are locally
  // === in this processor (in case of parallel simulation)
  int local_num_cells = 0;
  foreach(serial) {
    if( PRINT_CELL_CRITERIA )
      local_num_cells++;
  }

  // === Allocatting memory for the local vertex arrays
  // === Note: the x,y,z coordinates will all go into the same array, so the array will be [x1, y1, z1, x2, y2, z2, ...]
  double *local_vertices = (double *)malloc( 3*vertices_per_cell*local_num_cells*sizeof(double) );

  // === Allocating memory for ALL the local scalar data arrays
  int number_of_scalar_fields = list_len(list_scalar_data);
  typedef double** doublepp; // Hiding the double pointer with a typedef because qcc gets really annoying if i dont do this (why??)
  doublepp local_data = number_of_scalar_fields ? (double **)malloc( number_of_scalar_fields*sizeof(double *) ) : NULL;
  for( int k=0; k<number_of_scalar_fields; k++ )
    local_data[k] = (double *)malloc( local_num_cells*sizeof(double) );

  // === Storing all the vertices coordinates in the arrays
  int i = 0, i_cell = 0;
  foreach(serial) {
    if( PRINT_CELL_CRITERIA ) {
      // Using a macro conditional to avoid checking dimension==2 every loop during execution time...
      #if dimension==2
        local_vertices[i++] = 0.0; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = 0.0; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = 0.0; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = 0.0; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x-0.5*Delta;
      #else // dimension
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x+0.5*Delta;
      #endif

      int list_index = 0;
      for(scalar s in list_scalar_data) {
        local_data[list_index][i_cell] = s[];
        SwapArrayBytes(&local_data[list_index][i_cell], sizeof(double));
        list_index++;
      }

      i_cell++;
    }
  }

  // === Getting how many cells we have in total between all processes
  // === And then gathering all the local vertices into the root process
  int total_num_cells;
  double *vertices = NULL;
  doublepp data = NULL;
  #if _MPI
    /// === Total number of cells
    MPI_Allreduce(&local_num_cells, &total_num_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    /// === Gathering all the local vertices and the scalar fields data into the root process
    if( pid()==0 ) {
      vertices = (double *)malloc( 3*total_num_cells*vertices_per_cell*sizeof(double) );
      data = (double **)malloc( number_of_scalar_fields*sizeof(double*) );
      for(int k=0; k<number_of_scalar_fields; k++)
        data[k] = (double *)malloc( total_num_cells*sizeof(double) );
    }
    MPI_Gather_Uneven(local_vertices, 3*vertices_per_cell*local_num_cells, MPI_DOUBLE, vertices, 0);
    int list_index = 0;
    for( list_index=0; list_index<number_of_scalar_fields; list_index++ )
      MPI_Gather_Uneven(local_data[list_index], local_num_cells, MPI_DOUBLE, (pid()==0) ? data[list_index] : NULL, 0);

    /// === Releasing local memory
    free(local_vertices);
    for(int k=0; k<number_of_scalar_fields; k++)
      free(local_data[k]);
    if( local_data )
      free(local_data);
  #else
    total_num_cells = local_num_cells;
    vertices = local_vertices;
    data = local_data;
  #endif

  /**
    From now on, we only do file-printing stuff. Only the root process will do it. */
  if( pid()==0 ) {

    // === Opening the file to print the mesh
    sprintf(nomeArq, "%s/Mesh-N%d.vtk", folder_name, n);
    arq = fopen(nomeArq, "wt");

    // === Printing the VTK header information (printed as ASCII text)
    fprintf(arq, "# vtk DataFile Version 2.0\n");
    fprintf(arq, "MESH. step %d time %lf\n", 0, 0.0);
    fprintf(arq, "BINARY\n");
    fprintf(arq, "DATASET UNSTRUCTURED_GRID\n");

    // === Printing all the vertices coordinates (as BINARY)
    fprintf(arq, "POINTS %d double\n", total_num_cells*vertices_per_cell);
    SwapArrayBytes(vertices, 3*vertices_per_cell*total_num_cells*sizeof(double));
    fwrite(vertices, sizeof(double), 3*vertices_per_cell*total_num_cells, arq);
    fprintf(arq, "\n");

    // === Printing all the cells 
    // === Each cell contains 4 (or 8) indices referring to the vertices above
    fprintf(arq, "CELLS %d %d\n", total_num_cells, (vertices_per_cell + 1)*total_num_cells);
    int *array_cell_indices = malloc( (vertices_per_cell + 1)*total_num_cells*sizeof(int) );
    int offset = 0, vertex_index = 0;
    for( i=0; i<total_num_cells; i++ ) {
      array_cell_indices[offset] = vertex_index;
      array_cell_indices[offset + 1] = vertex_index + 1;
      array_cell_indices[offset + 2] = vertex_index + 2;
      array_cell_indices[offset + 3] = vertex_index + 3;
      #if dimension==3
        array_cell_indices[offset + 4] = vertex_index + 4;
        array_cell_indices[offset + 5] = vertex_index + 5;
        array_cell_indices[offset + 6] = vertex_index + 6;
        array_cell_indices[offset + 7] = vertex_index + 7;
      #endif
      array_cell_indices[offset + vertices_per_cell] = vertices_per_cell;
      offset += vertices_per_cell + 1;
      vertex_index += vertices_per_cell;
    }
    SwapArrayBytes(array_cell_indices, (vertices_per_cell + 1)*total_num_cells*sizeof(int));
    fwrite(array_cell_indices, sizeof(int), (vertices_per_cell + 1)*total_num_cells, arq);
    fprintf(arq, "\n");
        
    // === Printing cell types (squares or cubes)
    fprintf(arq, "CELL_TYPES %d\n", total_num_cells);
    SwapArrayBytes(&vtk_cell_type, sizeof(int));
    for( i=0; i<total_num_cells; i++ )
      array_cell_indices[i] = vtk_cell_type;
    fwrite(array_cell_indices, sizeof(int), total_num_cells, arq);
    fprintf(arq, "\n");


    // === Printing the actual simulation data that is stored in the cells
    int list_index = 0;
    fprintf(arq, "CELL_DATA %d\n", total_num_cells);
    for( scalar s in list_scalar_data ) {
      fprintf(arq, "SCALARS %s double 1\n", list_scalar_names[list_index]);
      fprintf(arq, "LOOKUP_TABLE default\n");
      fwrite(data[list_index], sizeof(double), total_num_cells, arq);
      fprintf(arq, "\n");
      list_index++;
    }

    // === Releasing memory
    free(array_cell_indices);
    free(vertices);
    for( list_index=0; list_index<number_of_scalar_fields; list_index++ )
      free(data[list_index]);
    free(data);
    fclose(arq);
  }

  return;
}

/**
Prints mesh and (optionally) scalars to a single-precision BINARY VTK file. */
void PrintMeshVTK_Binary_Float(int n, double time, scalar *list_scalar_data, const char **list_scalar_names)
{
  FILE *arq;
  char nomeArq[900];

  // ===  Cell is either a square (2D) or cube (3D), meaning either 4 or 8 vertices per cell
  int vertices_per_cell = (dimension==2) ? 4 : 8;

  // === VTK cell code that represents voxels
  // === (check Figure 2 here to understand: https://kitware.github.io/vtk-examples/site/VTKFileFormats/)
  int vtk_cell_type = (dimension==2) ? 9 : 11;

  // === Counting how many local cells we have in the mesh
  // === NOTE: whenever I say "local", i mean things that are locally in this processor (in case of parallel simulation)
  int local_num_cells = 0;
  foreach(serial) {
    if( PRINT_CELL_CRITERIA )
      local_num_cells++;
  }

  // === Allocatting memory for the local vertex arrays
  // === Note: the x,y,z coordinates will all go into the same array, 
  // === so the array will be [x1, y1, z1, x2, y2, z2, ...]
  float *local_vertices = (float *)malloc( 3*vertices_per_cell*local_num_cells*sizeof(float) );

  // === Allocating memory for ALL the local scalar data arrays
  int number_of_scalar_fields = list_len(list_scalar_data);
  typedef float** floatpp; // qcc gets uncomfortable if i dont hide the double pointer (?????)
  floatpp local_data = number_of_scalar_fields ? (float **)malloc( number_of_scalar_fields*sizeof(float *) ) : NULL;
  for( int k=0; k<number_of_scalar_fields; k++ )
    local_data[k] = (float *)malloc( local_num_cells*sizeof(float) );

  // === Storing all the vertices coordinates in the arrays
  int i = 0, i_cell = 0;
  foreach(serial) {
    if( PRINT_CELL_CRITERIA ) {

      /// Using a macro conditional to avoid checking dimension==2 every loop during execution time...
      #if dimension==2
        local_vertices[i++] = 0.0; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = 0.0; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = 0.0; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = 0.0; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x-0.5*Delta;
      #else // dimension
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x+0.5*Delta;
      #endif

      int list_index = 0;
      for(scalar s in list_scalar_data) {
        local_data[list_index][i_cell] = (float) s[];
        SwapArrayBytes(&local_data[list_index][i_cell], sizeof(float));
        list_index++;
      }

      i_cell++;
    }
  }

  // === Getting how many cells we have in total between all processes
  // === And then gathering all the local vertices into the root process
  int total_num_cells;
  float *vertices = NULL;
  floatpp data = NULL;
  #if _MPI
    // === Total number of cells
    MPI_Allreduce(&local_num_cells, &total_num_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // === Gathering all the local vertices and the scalar fields data into the root process
    if( pid()==0 ) {
      vertices = (float *)malloc( 3*total_num_cells*vertices_per_cell*sizeof(float) );
      data = (float **)malloc( number_of_scalar_fields*sizeof(float*) );
      for(int k=0; k<number_of_scalar_fields; k++)
        data[k] = (float *)malloc( total_num_cells*sizeof(float) );
    }
    MPI_Gather_Uneven(local_vertices, 3*vertices_per_cell*local_num_cells, MPI_FLOAT, vertices, 0);
    int list_index = 0;
    for( list_index=0; list_index<number_of_scalar_fields; list_index++ )
      MPI_Gather_Uneven(local_data[list_index], local_num_cells, MPI_FLOAT, (pid()==0) ? data[list_index] : NULL, 0);

    // === Releasing local memory
    free(local_vertices);
    for(int k=0; k<number_of_scalar_fields; k++)
      free(local_data[k]);
    if( local_data )
      free(local_data);
  #else
    total_num_cells = local_num_cells;
    vertices = local_vertices;
    data = local_data;
  #endif

  
  
  // === From now on, we only do file-printing stuff. Only the root process will do it
  if( pid()==0 ) {
    // === Opening the file to print the mesh
    sprintf(nomeArq, "%s/Mesh-N%d.vtk", folder_name, n);
    arq = fopen(nomeArq, "wt");

    // === Printing the VTK header information (printed as ASCII text)
    fprintf(arq, "# vtk DataFile Version 2.0\n");
    fprintf(arq, "MESH. step %d time %lf\n", 0, 0.0);
    fprintf(arq, "BINARY\n");
    fprintf(arq, "DATASET UNSTRUCTURED_GRID\n");

    // === Printing all the vertices coordinates (as BINARY)
    fprintf(arq, "POINTS %d float\n", total_num_cells*vertices_per_cell);
    SwapArrayBytes(vertices, 3*vertices_per_cell*total_num_cells*sizeof(float));
    fwrite(vertices, sizeof(float), 3*vertices_per_cell*total_num_cells, arq);
    fprintf(arq, "\n");

    // === Printing all the cells (each cell contains 4 (or 8) indices referring to the vertices printed above)
    fprintf(arq, "CELLS %d %d\n", total_num_cells, (vertices_per_cell + 1)*total_num_cells);
    int *array_cell_indices = malloc( (vertices_per_cell + 1)*total_num_cells*sizeof(int) );
    int offset = 0, vertex_index = 0;
    for( i=0; i<total_num_cells; i++ ) {
      array_cell_indices[offset] = vertex_index;
      array_cell_indices[offset + 1] = vertex_index + 1;
      array_cell_indices[offset + 2] = vertex_index + 2;
      array_cell_indices[offset + 3] = vertex_index + 3;
      #if dimension==3
        array_cell_indices[offset + 4] = vertex_index + 4;
        array_cell_indices[offset + 5] = vertex_index + 5;
        array_cell_indices[offset + 6] = vertex_index + 6;
        array_cell_indices[offset + 7] = vertex_index + 7;
      #endif
      array_cell_indices[offset + vertices_per_cell] = vertices_per_cell;
      offset += vertices_per_cell + 1;
      vertex_index += vertices_per_cell;
    }
    SwapArrayBytes(array_cell_indices, (vertices_per_cell + 1)*total_num_cells*sizeof(int));
    fwrite(array_cell_indices, sizeof(int), (vertices_per_cell + 1)*total_num_cells, arq);
    fprintf(arq, "\n");
        

    fprintf(arq, "CELL_TYPES %d\n", total_num_cells);
    SwapArrayBytes(&vtk_cell_type, sizeof(int));
    for( i=0; i<total_num_cells; i++ )
      array_cell_indices[i] = vtk_cell_type;
    fwrite(array_cell_indices, sizeof(int), total_num_cells, arq);
    fprintf(arq, "\n");


    // === Printing the actual simulation data that is stored in the cells
    int list_index = 0;
    fprintf(arq, "CELL_DATA %d\n", total_num_cells);
    for( scalar s in list_scalar_data ) {
      fprintf(arq, "SCALARS %s float 1\n", list_scalar_names[list_index]);
      fprintf(arq, "LOOKUP_TABLE default\n");
      fwrite(data[list_index], sizeof(float), total_num_cells, arq);
      fprintf(arq, "\n");
      list_index++;
    }

    free(array_cell_indices);
    free(vertices);
    for( list_index=0; list_index<number_of_scalar_fields; list_index++ )
      free(data[list_index]);
    free(data);
    fclose(arq);
  }

  return;
}

/**
Structure used to pass optional parameters to the PrintMeshVTK and PrintMeshVTK_Binary. */
typedef enum {VTK_TYPE_ASCII, VTK_TYPE_BINARY} VTK_FILE_TYPE;
typedef enum {VTK_PRECISION_FLOAT, VTK_PRECISION_DOUBLE} VTK_FILE_PRECISION;
struct StructPrintMesh {
  int n;
  double time;
  scalar *list_scalar_data;
  const char **list_scalar_names;
  VTK_FILE_TYPE vtk_type;
  VTK_FILE_PRECISION vtk_precision; // only relevant if vtk_type==binary
};

/**
  Prints mesh and (optionally) scalar data to a VTK file. */
void PrintMeshVTK(struct StructPrintMesh spm)
{
  // if( spm.vtk_type==VTK_TYPE_ASCII )
  //   PrintMeshVTK_ASCII(spm.n, spm.time, spm.list_scalar_data, spm.list_scalar_names);
  // else if( spm.vtk_type==VTK_TYPE_BINARY ) {        
  //   if( spm.vtk_precision==VTK_PRECISION_FLOAT )
  //     PrintMeshVTK_Binary_Float(spm.n, spm.time, spm.list_scalar_data, spm.list_scalar_names);
  //   else if( spm.vtk_precision==VTK_PRECISION_DOUBLE )
  //     PrintMeshVTK_Binary_Double(spm.n, spm.time, spm.list_scalar_data, spm.list_scalar_names);
  //   else
  //     ErrorMessage( 0, "The precision parameter for VTK Binary files should be either VTK_PRECISION_FLOAT or VTK_PRECISION_DOUBLE.", "PrintMeshVTK", 1 );
  // }
  // else
  //   ErrorMessage( 0, "The vtk_type parameter for VTK Binary files should be either VTK_TYPE_ASCII or VTK_TYPE_BINARY.", "PrintMeshVTK", 1 );
  return;
}

/** 
  I blatantly copied this function from here: [draw.h](http://basilisk.fr/src/draw.h) */
static bool cfilter (Point point, scalar c, double cmin)
{
  double cmin1 = 4.*cmin;
  if (c[] <= cmin) {
    foreach_dimension()
      if (c[1] >= 1. - cmin1 || c[-1] >= 1. - cmin1)
	return true;
    return false;
  }
  if (c[] >= 1. - cmin) {
    foreach_dimension()
      if (c[1] <= cmin1 || c[-1] <= cmin1)
	return true;
    return false;
  }
  int n = 0;
  double min = HUGE, max = - HUGE;
  foreach_neighbor(1) {
    if (c[] > cmin && c[] < 1. - cmin && ++n >= (1 << dimension))
      return true;
    if (c[] > max) max = c[];
    if (c[] < min) min = c[];
  }
  return max - min > 0.5;
}

/** 
  Prints a VTK file with the current interface.
  NOTE: Currently, I am printing the interface in ASCII format.
  I will make the Binary version soon. */
void PrintInterfaceVTK(int n, double time)
{
  double fmin = 1e-3; // do not reconstruct fragments smaller than this

  // === Counting how many line segments we will draw and how many vertices
  int local_count_polys = 0;
  int local_count_vertices = 0;
  foreach(serial) {
    if( cfilter (point, f, fmin) ) {
      local_count_polys++;
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);
      #if dimension==2
        coord v[2];
        int m = facets (n, alpha, v);
      #else
        coord v[12];
        int m = facets (n, alpha, v, 1.1);
      #endif
      local_count_vertices += m;
    }
  }

  // === Alocating memory for the local vertices
  int *local_poly_count_vertices = (int *)malloc( local_count_polys*sizeof(int) );
  double *local_vertices_x = (double *)malloc( local_count_vertices*sizeof(double) );
  double *local_vertices_y = (double *)malloc( local_count_vertices*sizeof(double) );
  double *local_vertices_z = (double *)malloc( local_count_vertices*sizeof(double) );

  // === Calculating all polygons again (bad) and saving the vertices to the local arrays
  int index_polygon = 0, index_vertex = 0;
  foreach(serial) {
    if( cfilter (point, f, fmin) ) {
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);
      #if dimension==2
        coord v[2];
        int m = facets (n, alpha, v);
      #else
        coord v[12];
        int m = facets (n, alpha, v, 1.1);
      #endif
      local_poly_count_vertices[index_polygon++] = m;
      for( int i=0; i<m; i++ ) {
        local_vertices_x[index_vertex] = x + v[i].x*Delta;
        local_vertices_y[index_vertex] = y + v[i].y*Delta;
        #if dimension==2
          local_vertices_z[index_vertex] = 0.0;
        #else
          local_vertices_z[index_vertex] = z + v[i].z*Delta;
        #endif
        index_vertex++;
      }
    }
  }

  // === Counting how many polys and vertices we have in total between all processes
  int total_count_vertices = 0, total_count_polys = 0;
  #if _MPI
    MPI_Allreduce(&local_count_vertices, &total_count_vertices, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_count_polys, &total_count_polys, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #else
    total_count_vertices = local_count_vertices;
    total_count_polys = local_count_polys;
  #endif

  // === Sending all the vertices to the global vertex array in process id=0
  int *poly_count_vertices = NULL;
  double *vertices_x = NULL, *vertices_y = NULL, *vertices_z = NULL;
  #if _MPI
    if( pid()==0 ) {
      poly_count_vertices = (int *)malloc( total_count_polys*sizeof(int) );
      vertices_x = (double *)malloc( total_count_vertices*sizeof(double) );
      vertices_y = (double *)malloc( total_count_vertices*sizeof(double) );
      vertices_z = (double *)malloc( total_count_vertices*sizeof(double) );
    }

    MPI_Gather_Uneven(local_poly_count_vertices, local_count_polys, MPI_INT, poly_count_vertices, 0);
    MPI_Gather_Uneven(local_vertices_x, local_count_vertices, MPI_DOUBLE, vertices_x, 0);
    MPI_Gather_Uneven(local_vertices_y, local_count_vertices, MPI_DOUBLE, vertices_y, 0);
    MPI_Gather_Uneven(local_vertices_z, local_count_vertices, MPI_DOUBLE, vertices_z, 0);

    // === Releasing local memory
    free(local_poly_count_vertices);
    free(local_vertices_x);
    free(local_vertices_y);
    free(local_vertices_z);
  #else
    poly_count_vertices = local_poly_count_vertices;
    vertices_x = local_vertices_x;
    vertices_y = local_vertices_y;
    vertices_z = local_vertices_z;
  #endif

  // === From this point we will only do file-writing stuff. Only process ZERO will do it
  if( pid()!=0 )
    return;

  // === Opening the vtk file
  char nomeArq[900];
  sprintf(nomeArq, "%s/Interface-N%d.vtk", folder_name, n);
  FILE *arq = fopen(nomeArq, "wt");
  if( !arq ) {
    printf("\n\n PrintInterfaceVTK_2D: Problem opening file... \n\n");
    exit(0);
  }

  // === Writing the VTK header
  fprintf(arq, "# vtk DataFile Version 2.0\n");
  fprintf(arq, "INTERFACE. step %d time %lf\n", n, time);
  fprintf(arq, "ASCII\n");
  fprintf(arq, "DATASET POLYDATA\n");
  fprintf(arq, "POINTS %d float\n", total_count_vertices);

  
  // === Writing all the surface vertices
  for( index_vertex=0; index_vertex<total_count_vertices; index_vertex++ )
    fprintf(arq, "%e %e %e\n", vertices_x[index_vertex], vertices_y[index_vertex], vertices_z[index_vertex]);
  

  // === Writing the polygons conectivity
  int global_index_vertex = 0;
  #if dimension==2
    fprintf(arq, "LINES %d %d\n", total_count_polys, total_count_polys + total_count_vertices);
  #else
    fprintf(arq, "POLYGONS %d %d\n", total_count_polys, total_count_polys + total_count_vertices);
  #endif
  for( int index_polygon=0; index_polygon<total_count_polys; index_polygon++ ) {
    fprintf(arq, "%d", poly_count_vertices[index_polygon]);
    int index_vertex;
    for( index_vertex=0; index_vertex<poly_count_vertices[index_polygon]; index_vertex++ )
      fprintf(arq, " %d", global_index_vertex++);
    fprintf(arq, "\n");
  }


  // === Releasing memory
  free(poly_count_vertices);
  free(vertices_x);
  free(vertices_y);
  free(vertices_z);
  fclose(arq);
  return;
}

void swap_bs(double* xp, double* yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void InitializeInterfaceFromVTK(const char *file_name, double **out_centers_x, double **out_centers_y, int *out_num_centers)
{
  // === Opening the vtk file
  FILE *arq = fopen(file_name, "rt");
  if( !arq ) {
    printf("\n\n InitializeInterfaceFromVTK: Problem opening file... \n\n");
    exit(0);
  }

  // === Writing the VTK header
  double time;
  int total_count_vertices, n;
  fscanf(arq, "# vtk DataFile Version 2.0\n");
  fscanf(arq, "INTERFACE. step %d time %lf\n", &n, &time);
  fscanf(arq, "ASCII\n");
  fscanf(arq, "DATASET POLYDATA\n");
  fscanf(arq, "POINTS %d float\n", &total_count_vertices);
  
  // === Writing all the surface vertices
  double *all_vertices_x = (double *)malloc( total_count_vertices*sizeof(double) );
  double *all_vertices_y = (double *)malloc( total_count_vertices*sizeof(double) );
  double aux_vertex_z;
  for( int index_vertex=0; index_vertex<total_count_vertices; index_vertex++ )
    fscanf(arq, "%lf %lf %lf\n", &all_vertices_x[index_vertex], &all_vertices_y[index_vertex], &aux_vertex_z);
  


  // === Writing the polygons conectivity
  int total_count_polys, aux, index_center = 0;
  #if dimension==2
    fscanf(arq, "LINES %d %d\n", &total_count_polys, &aux);
  #else
    fscanf(arq, "POLYGONS %d %d\n", &total_count_polys, &aux);
  #endif

  // Allocating memory for the centers
  double *centers_x = (double *)malloc(total_count_polys*sizeof(double));
  double *centers_y = (double *)malloc(total_count_polys*sizeof(double));

  index_center = 0;
  for( int index_polygon=0; index_polygon<total_count_polys; index_polygon++ ) {
    int count_vertices_polygon, global_index_vertex[10];
    fscanf(arq, "%d", &count_vertices_polygon);
    int index_vertex;
    for( index_vertex=0; index_vertex<count_vertices_polygon; index_vertex++ )
      fscanf(arq, " %d", &global_index_vertex[index_vertex]);
    fscanf(arq, "\n");

    if( count_vertices_polygon!=2 ) {
      printf("Something wrong with number of vertices in a line: %d\n", count_vertices_polygon);
      exit(1);
    }

    centers_x[index_center] = 0.5*( all_vertices_x[global_index_vertex[0]] + all_vertices_x[global_index_vertex[1]] );
    centers_y[index_center] = 0.5*( all_vertices_y[global_index_vertex[0]] + all_vertices_y[global_index_vertex[1]] );

    if( centers_x[index_center]>0.005 ) {
      // printf("center: %lf %lf\n", centers_x[index_center], centers_y[index_center]);
      index_center++;
    }
  }

  int i, j;
  // Ordering by y-coordinate (just a bubble sort, dont care about efficiency here)
  n = index_center;
  bool swapped;
  for (i = 0; i < n - 1; i++) {
      swapped = false;
      for (j = 0; j < n - i - 1; j++) {
          if (centers_y[j] > centers_y[j + 1]) {
              swap_bs(&centers_y[j], &centers_y[j + 1]);
              swap_bs(&centers_x[j], &centers_x[j + 1]);
              swapped = true;
          }
      }

      // If no two elements were swapped by inner loop,
      // then break
      if (swapped == false)
          break;
  }

  free(all_vertices_x);
  free(all_vertices_y);

  *out_centers_x = centers_x;
  *out_centers_y = centers_y;
  *out_num_centers = index_center;
}

void calculate_energy_budget(double *kinetic, double *surface, double *elastic, double *dissipated, double Oh, double J, double regularization)
{
  double kinetic_temp = 0.0;
  double surface_temp = 0.0;
  double elastic_temp = 0.0;
  double dissipated_temp = 0.0;


  foreach(reduction(+:kinetic_temp) reduction(+:surface_temp) reduction(+:elastic_temp) reduction(+:dissipated_temp)) {

    // ==== Calculating kinetic energy
    kinetic_temp += f[]*(2.0*pi*y)*0.5*(sq(u.x[]) + sq(u.y[]))*sq(Delta);

    // ==== Calculating surface energy
    double fmin = 1e-03;
    if( cfilter (point, f, fmin) ) {
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);
      coord v[2];
      facets (n, alpha, v);
      
      double p1x = x + v[0].x*Delta;
      double p1y = y + v[0].y*Delta;
      double p2x = x + v[1].x*Delta;
      double p2y = y + v[1].y*Delta;

      double center_y = 0.5*(p1y + p2y);
      double center_x = 0.5*(p1x + p2x);

      if( center_x>1e-3 )
        surface_temp += 2.0*pi*center_y*sqrt( sq(p1x - p2x) + sq(p1y - p2y) );
    }

    // ==== Calculating elastic energy
    // HOW?

    // ==== Calculating viscously dissipated energy
    // Note: this is only the dissipated energy at the current timestep
    // Note: If you want the total dissipation over the simulation, you have to integrate this over time
    double dudx = (u.x[1, 0] - u.x[-1, 0])/(2.0*Delta);
    double dudy = (u.x[0, 1] - u.x[0, -1])/(2.0*Delta);
    double dvdx = (u.y[1, 0] - u.y[-1, 0])/(2.0*Delta);
    double dvdy = (u.y[0, 1] - u.y[0, -1])/(2.0*Delta);
    double axi_term = u.y[]/max(y, 1e-12);
    double norm_strain = sqrt( sq(dudx) + sq(dvdy) + sq(axi_term) + 0.5*sq(dudy + dvdx) );
    double bingham_viscosity = Oh + J/(2.0*norm_strain + 1e-10)*(1.0 - exp(-norm_strain/regularization));
    dissipated_temp += f[]*(2.0*pi*y)*sq(Delta)*bingham_viscosity*(2.0*dudx*dudx + 2.0*dvdy*dvdy + 2.0*axi_term + (dvdx  + dudy)*(dvdx  + dudy));
  }

  *kinetic = kinetic_temp;
  *surface = surface_temp;
  *elastic = elastic_temp;
  *dissipated += dt*dissipated_temp;

  return;
}

void **AllocateMatrix(int Linhas, int Colunas, size_t TamanhoElemento, void *ValorInicial)
{
    void **matriz;
    int i, j, byteAtualValor, byteAtualMatriz;
    unsigned char *charValorInicial, *charMatriz;

    charValorInicial = (unsigned char *)ValorInicial;

    matriz = (void **)malloc( Linhas*sizeof(void *) );
    for( i=0; i<Linhas; i++ )
        matriz[i] = malloc( Colunas*TamanhoElemento );

    //Inicializando todos os elementos com o valor inicial (soh funciona se for um tipo int)
    if( ValorInicial!=NULL ) {
        for(i=0; i<Linhas; i++) {
            charMatriz = (unsigned char *)matriz[i];
            byteAtualMatriz = 0;
            for(j=0; j<Colunas; j++) {
                for( byteAtualValor = 0; byteAtualValor<TamanhoElemento; byteAtualValor++ )
                    charMatriz[byteAtualMatriz++] = charValorInicial[byteAtualValor];
            }
        }
    }


    return matriz;
}

void DeallocateMatrix(void **Matriz, int Linhas, int Colunas)
{
    int i;

    for( i=0; i<Linhas; i++ )
        free(Matriz[i]);
    free(Matriz);
}

double uni_max_mesh_x = - 1.0, uni_max_mesh_y = -1.0;
void PrintMeshVTK_Uniform_Nico(int n, double time, double mesh_level, scalar *list_scalar_data, const char **list_scalar_names)
{
  double dx = 5.0/pow(2.0, mesh_level);
  foreach(reduction(max:uni_max_mesh_x) reduction(max:uni_max_mesh_y)) {
    double max_x = x + 0.5*Delta;
    double max_y = y + 0.5*Delta;
    if( max_x<=1.1 && max_y<=2.5 ) {
      uni_max_mesh_x = ((x+0.5*Delta)>uni_max_mesh_x) ? x + 0.5*Delta : uni_max_mesh_x;
      uni_max_mesh_y = ((y+0.5*Delta)>uni_max_mesh_y) ? y + 0.5*Delta : uni_max_mesh_y;
    }
  }
  int num_cells_x = round(uni_max_mesh_x/dx);
  int num_cells_y = round(uni_max_mesh_y/dx);

  double initial_value = 0.0;
  double **matrix_u = (double **)AllocateMatrix(num_cells_x, num_cells_y, sizeof(double), &initial_value);
  double **matrix_v = (double **)AllocateMatrix(num_cells_x, num_cells_y, sizeof(double), &initial_value);
  double **matrix_f = (double **)AllocateMatrix(num_cells_x, num_cells_y, sizeof(double), &initial_value);
  for(int i=0; i<num_cells_x; i++) {
    double x_interp = (i + 0.5)*dx;
    for(int j=0; j<num_cells_y; j++ ) {
      double y_interp = (j + 0.5)*dx;
      matrix_u[i][j] = interpolate(u.x, x_interp, y_interp);
      matrix_v[i][j] = interpolate(u.y, x_interp, y_interp);
      matrix_f[i][j] = interpolate(f, x_interp, y_interp);
    }
  }

  // From below this point we only continue in process 0
  if( pid() ) {
    DeallocateMatrix((void **)matrix_u, num_cells_x, num_cells_y);
    DeallocateMatrix((void **)matrix_v, num_cells_x, num_cells_y);
    DeallocateMatrix((void **)matrix_f, num_cells_x, num_cells_y);
    return;
  }

  // === Opening the vtk file
  char nomeArq[900];
  sprintf(nomeArq, "%s/UniformMesh-N%d.vtk", folder_name, n);
  FILE *arq = fopen(nomeArq, "wt");
  if( !arq ) {
    printf("\n\n PrintMeshVTK_Uniform_Nico: Problem opening file... \n\n");
    exit(0);
  }

  //Cabecalho
  int i, j;
  fprintf(arq, "# vtk DataFile Version 2.0\nStructured Grid Example\nASCII\n");
  fprintf(arq, "DATASET RECTILINEAR_GRID\n");
  fprintf(arq, "DIMENSIONS %d %d 1\n", num_cells_x+1, num_cells_y+1);
  fprintf(arq, "X_COORDINATES %d float\n", num_cells_x+1);
  for( i=0; i<=num_cells_x; i++ )
    fprintf(arq, "%lf ", 0.0 + i*dx);
  fprintf(arq, "\nY_COORDINATES %d float\n", num_cells_y + 1);
  for( i=0; i<=num_cells_y; i++ )
    fprintf(arq, "%lf ", 0.0 + i*dx);
  fprintf(arq, "\nZ_COORDINATES %d float\n", 1);
  fprintf(arq, "0.0\n\n");

  //Vetores velocidade (estou imprimindo valores de U e V nos vertices das celulas fazendo interpolacao)


  /* Velocity scalars */
  fprintf(arq, "CELL_DATA %d\n", (num_cells_x)*(num_cells_y));
  fprintf(arq, "VECTORS velocity float\n");
  for(j=0; j<num_cells_y; j++)
    for(i=0; i<num_cells_x; i++)
      fprintf(arq, "%e %e 0.0000\n", matrix_u[i][j], matrix_v[i][j]);

  fprintf(arq, "SCALARS fractions float 1\n");
  fprintf(arq, "LOOKUP_TABLE default\n");
  for(j=0; j<num_cells_y; j++)
    for(i=0; i<num_cells_x; i++)
      fprintf(arq, "%e\n", matrix_f[i][j]);
  fprintf(arq, "\n");

  fclose(arq);
  DeallocateMatrix((void **)matrix_u, num_cells_x, num_cells_y);
  DeallocateMatrix((void **)matrix_v, num_cells_x, num_cells_y);
  DeallocateMatrix((void **)matrix_f, num_cells_x, num_cells_y);
  return;
}



