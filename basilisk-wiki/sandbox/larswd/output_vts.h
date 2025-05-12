/*
This function outputs the entire computational domain from a multilayer solver to a vts file which then can be read using Paraview. 
Is compatible with multigrid and OpenMP (Thanks to Øystein Lande, whose code this function was developed from.) 
*/
    

void output_vts_ascii_all_layers(FILE* fp, scalar* list, int N)
{

    
    int nthreads_ = omp_get_num_threads();
    omp_set_num_threads(1); 
    
        // MULTIGRID
    
 fputs("<?xml version=\"1.0\"?>\n<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, N, 0, nl);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, N, 0, nl);
#endif

#if dimension == 2
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, nl);
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, nl);
#endif

    // Loop over velocity data and store kinematics in cell vector stucture
    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");
    for(int i = nl-1; i >= 0; i--){
    foreach() {
#if dimension == 1
        fprintf(fp, "%g %g 0.\n", u.x[0,0,i], w[0,0,i]);
#endif
#if dimension == 2
        fprintf(fp, "%g %g %g\n", u.x[0,0,i], u.y[0,0,i], w[0,0,i]);
#endif
    }
    }
    // Dummy text to get correct amount of layers
    foreach() {
        fprintf(fp, "0 0 0.\n");
    }
    fputs("\t\t\t\t </DataArray>\n", fp);


    // loop over all scalars in scalarlist and store values as cell data
    for (scalar s in list) {
        if (strcmp(s.name, "eta") == 0) {
            fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
            for(int i = nl-1; i >= 0; i--){
                foreach() {
                if (h[] > dry) {
                    fprintf(fp, "%g\n", s[0,0,i]);
                }
                else {
                    fprintf(fp, "nan\n");
                }
            }
        }
        fputs("\t\t\t\t </DataArray>\n", fp);
        }
        else {
            fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
            for(int i = nl-1; i >= 0; i--){
            foreach() {
                fprintf(fp, "%g\n", s[0,0,i]);
            }
            }
            foreach(){
                fprintf(fp, "0\n");
            }
            fputs("\t\t\t\t </DataArray>\n", fp);
        }   
    }   

    fputs("\t\t\t </CellData>\n", fp);
    // Coordinates 
    fputs("\t\t\t <Points>\n", fp);
#if dimension == 1
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
        double* zcorr = (double *)malloc((N+1)*sizeof(double));
        for (int j = 0; j <= N; j++){
            zcorr[j] = 0;
        }
        int k;
        for(int i = nl-1; i >= 0; i--){
            k = 0;
            foreach_vertex(serial) {
                fprintf(fp, "%12.4f %12.4f 0.\n", x,zcorr[k],0 );
                k++;
            }
            zcorr = zcorr - h[0,0,i];
        }
        foreach_vertex(serial){
                fprintf(fp, "%12.4f %12.4f 0.\n", x, zcorr);
        }
        free(zcorr);
#endif
#if dimension == 2
    fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
        
        double* zcorr = (double*)malloc((N+1)*(N+1)*sizeof(double));
        for (int j = 0; j < (N+1)*(N+1); j++){
            zcorr[j] = 0;
        }
        int j, k;
        for(int i = nl-1; i >= 0; i--){
            j = 0;
            foreach_vertex(serial) {
                fprintf(fp, "%12.4f %12.4f %12.4f\n", x,y, zcorr[j]);
                
                if (h[0,0,i] < 1e-3){
                    zcorr[j] = zcorr[j] - h[-1,-1,i];
                } else {
                    zcorr[j] = zcorr[j] - h[0,0,i];
                }
                j++;
            }
        }
        j = 0;
        foreach_vertex(serial){
                fprintf(fp, "%12.4f %12.4f %12.4f\n", x,y, zcorr[j]);
                j++;
        }
        free(zcorr);
#endif

    fputs("\t\t\t\t </DataArray>\n", fp);
    fputs("\t\t\t </Points>\n", fp);
    fputs("\t\t </Piece>\n", fp);

    // write time value
    fprintf(fp, "\t\t <FieldData> \n");
    fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray > \n");
    fprintf(fp, "\t\t </FieldData> \n");
    fputs("\t </StructuredGrid>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);

    omp_set_num_threads(nthreads_);

}