/** 
## Short intro
This page is to explain the equidistant averaging function. Notes on convergence in time and space can be respectively found [here](equidistant_time.c) and [here](equidistant_space.c).
Note that the convergence in space is carried out only for a constant high resolution grid, since the errors of the low grid values will dominate the higher resolution ones. Having a combination of resolutions is left as an exercise. 

## The functions 


/** NOTE: Very much in progress.

These functions are geared towards equidistant diagnostics. The level for which diagnostics are done can be specified. Values will be linearly interpolated if necesarry. 

So far this function works with MPI and only works in 3D. 


Short description of what they do:

//void equi_diag(scalar s, int lvl, int diagi) {
    adds values to the equifield array
}

//void equi_output_ascii (scalar s, char * fname, int lvl, int diagi) {
    outputs equifield array to ascii file fname in the format: 
    place = x*n*n + y*n + z
    data is devided by diagi to average
}

//void equi_output_binary(scalar s, char * fname, int lvl, int diagi)i {
    outputs equifield to file binary fname in the format: 
    place = x*n*n + y*n + z
}

//void equi_import_binary (scalar s, char * fname, int lvl) { 
    imports data written by the equi_output_binary(...) function
    data is read into scalar s, no interpolation is used.
}


TODO implement 2D
TODO walking average since values in equifield could get big, or divide before summing?
TODO multiple fields at once? 
*/

double *equifield = NULL; 	// The diagnostics field

int equi_diag (scalar s, int lvl, int diagi){  
    int len = 1;		
    int n = 1<<lvl;
    if(!equifield) {
        equifield = calloc(len*n*n*n, sizeof(double)); // Init with zeros
    } 	

    double dDelta = 0.9999999*L0/n;  

    restriction({s}); 		// Make sure that coarse level values exist
    boundary({s}); 		// Update boundaries
    
    // Iterate over diag level of leaf cell
    foreach_level_or_leaf(lvl) {
  	if(level < lvl){
	    int pn = 1<<(lvl-level);
	    for(int i = 0; i < pn; i++) {
		for(int j = 0; j < pn; j++) {
		    for(int k = 0; k < pn; k++) {
			
			double ctrans = dDelta/2. - Delta/2.;	  // Translation to corner of cell
			double xx = x + ctrans + i*dDelta;   	  // Coordinates in equidistant grid
			double yy = y + ctrans + j*dDelta;
 			double zz = z + ctrans + k*dDelta;

		        int ii = round((xx - dDelta/2.)/dDelta);  // Iteration in equidistant grid
	         	int jj = round((yy - dDelta/2.)/dDelta);
            	  	int kk = round((zz - dDelta/2.)/dDelta);

	    		int place = ii*n*n + jj*n + kk*len;  	  // Location in equidistant grid 
			double temp = interpolate(s, xx, yy, zz); // Linear interpolation 

 			equifield[place] += temp;
		    }
		}
	    }

        } else {
 
            int ii = round((x - dDelta/2.)/dDelta); // x adaptive is x equidistant 
       	    int jj = round((y - dDelta/2.)/dDelta); // Derive place in equidistant grid
     	    int kk = round((z - dDelta/2.)/dDelta);

	    int place = ii*n*n + jj*n + kk*len;
	    equifield[place] += s[];
	}
    }
    return diagi + 1;
}


/** 
Output equidistant double field to double binary file 
*/
void equi_output_binary (scalar s, char * fname, int lvl, int diagi) {    
    int len = 1;
    int n = 1<<lvl;

    if(pid() == 0){
    FILE * fpm = fopen(fname, "w");     
// Make sure all thread values are summed up using MPI
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, equifield, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
@endif
    int sizef = sizeof(double);
    int len = 1;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
		int place = i*n*n + j*n + k*len;
		double temp = (double)(equifield[place]/((double)diagi)); // Divide by number of additions
	        fseek(fpm, place*sizef, SEEK_SET);
    		fwrite(&temp, sizef, 1, fpm); ;
	    }
	}    
    }
    fclose(fpm);	
    }
@if _MPI
    else // slave
    MPI_Reduce (equifield, NULL, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
@endif
}

/** 
Output equidistant double field to ascii file 
*/
void equi_output_ascii (scalar s, char * fname, int lvl, int diagi) {
    int len = 1;
    int n = 1<<lvl;

    if(pid() == 0){
    FILE * fpm = fopen(fname, "w"); 
// Make sure all thread values are summed up using MPI
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, equifield, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
@endif
    int len = 1;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
		int place = i*n*n + j*n + k*len;
		double temp = (equifield[place]/((double)diagi)); // Divide by number of additions
	        fprintf(fpm, "%g\t", temp);
	    }
	}    
    }
    fclose(fpm);	
    }
@if _MPI
    else // slave
    MPI_Reduce (equifield, NULL, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
@endif
}

/** 
Importing binary doubles from file written by equi_output_binary(...) 
It is a bit bodged together but for my application it seems fast enough.

Array equifield is filled with the imported data after which a scalar field is filled.
*/
void equi_import_binary (scalar s, char * fname, int lvl){  
    int len = 1;		
    int n = 1<<lvl;

    free(equifield);	// Make sure equifield is reset 
    equifield = NULL;	
    
    if(!equifield) {
        equifield = calloc(len*n*n*n, sizeof(double)); // Init with zeros
    } 	

    double dDelta = 0.9999999*L0/n;  

    if(pid() == 0){
    FILE * fpm = fopen(fname, "r");     
    int sizef = sizeof(double);
    // Fill equifield array
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
		int place = i*n*n + j*n + k*len;
		fseek(fpm, place*sizef, SEEK_SET);
    		fread(&equifield[place], sizef, 1, fpm); ;
	    }
	}    
    }
    fclose(fpm);	
    }

@if _MPI 
    MPI_Allreduce (MPI_IN_PLACE, equifield, len*n*n*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
@endif
	

    // Filling each cell with respective value
    foreach() {
  	if(level < lvl){
 	    double aveval = 0.;
	    int pn = 1<<(lvl-level);
	    for(int i = 0; i < pn; i++) {
		for(int j = 0; j < pn; j++) {
		    for(int k = 0; k < pn; k++) {
			
			double ctrans = dDelta/2. - Delta/2.;	// Translation to corner of cell
			double xx = x + ctrans + i*dDelta;   	// Coordinates in equidistant grid
			double yy = y + ctrans + j*dDelta;
 			double zz = z + ctrans + k*dDelta;

		        int ii = round((xx - dDelta/2.)/dDelta); // Iteration in equidistant grid
	         	int jj = round((yy - dDelta/2.)/dDelta);
            	  	int kk = round((zz - dDelta/2.)/dDelta);

	    		int place = ii*n*n + jj*n + kk*len;	// Location in equidistant grid 

 			aveval += equifield[place];
		    }
		}
	    }
	    s[] = 1.*aveval/cube(pn);

        } else {
 
            int ii = round((x-dDelta/2)/dDelta); 
       	    int jj = round((y-dDelta/2)/dDelta); 
     	    int kk = round((z-dDelta/2)/dDelta);

	    int place = ii*n*n + jj*n + kk*len;
	    s[] = 1.* equifield[place];
	}
    }
}





