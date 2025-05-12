/*
Simple function for writing multilayer data to vts format which can be viewed in paraview.
Made by: Oystein Lande
*/


double dt_start = 0;
scalar id_boundary;
extern scalar w;
extern scalar phi;
extern int num_omp;

double nthreads_ = 1;
double gravity_ = 1.;
// Define alternative traversal order to match with the x ordering of vts (basilisk uses z ordering by default)

@def foreach_xorder()
OMP_PARALLEL() {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = { 0 };
    point.level = depth(); point.n = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
        for (_k = GHOSTS; _k < point.n + GHOSTS; _k++) {
#if dimension == 1
            point.i = _k;
#endif
#if dimension > 1
            point.j = _k;
            for (point.i = GHOSTS; point.i < point.n + GHOSTS; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_xorder()
#if dimension > 1
                }
#endif
}
}
@


@def foreach_face_xorder()
OMP_PARALLEL() {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = { 0 };
    point.level = depth(); point.n = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
        for (_k = GHOSTS; _k <= point.n + GHOSTS; _k++) {
#if dimension == 1
            point.i = _k;
#endif

#if dimension > 1
            point.j = _k;
            for (point.i = GHOSTS; point.i <= point.n + GHOSTS; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_face_xorder()
#if dimension > 1
                }
#endif
        }
}
@

@def foreach_vertex_xorder()
foreach_face_xorder() {
    x -= Delta / 2.;
#if dimension > 1  
    y -= Delta / 2.;
#endif
    @
    @define end_foreach_vertex_xorder() } end_foreach_face_xorder()





@def foreach_xorder()
OMP_PARALLEL() {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = { 0 };
    point.level = depth(); point.n = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
        for (_k = GHOSTS; _k < point.n + GHOSTS; _k++) {
#if dimension == 1
            point.i = _k;
#endif
#if dimension > 1
            point.j = _k;
            for (point.i = GHOSTS; point.i < point.n + GHOSTS; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_xorder()
#if dimension > 1
                }
#endif
}
}
@


@def foreach_face_xorder()
OMP_PARALLEL() {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = { 0 };
    point.level = depth(); point.n = 1 << point.level;
    int _k;
    OMP(omp for schedule(static))
        for (_k = GHOSTS; _k <= point.n + GHOSTS; _k++) {
#if dimension == 1
            point.i = _k;
#endif

#if dimension > 1
            point.j = _k;
            for (point.i = GHOSTS; point.i <= point.n + GHOSTS; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_face_xorder()
#if dimension > 1
                }
#endif
        }
}
@

@def foreach_vertex_xorder()
foreach_face_xorder() {
    x -= Delta / 2.;
#if dimension > 1  
    y -= Delta / 2.;
#endif
    @
    @define end_foreach_vertex_xorder() } end_foreach_face_xorder()



@def foreach_xorder_boundary()
OMP_PARALLEL() {
int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
Point point = { 0 };
point.level = depth(); point.n = 1 << point.level;
int _k;
OMP(omp for schedule(static))
    for (_k = 0; _k < point.n + 2*GHOSTS; _k++) {
#if dimension == 1
                point.i = _k;
#endif
#if dimension > 1
                point.j = _k;
                for (point.i = 0; point.i < point.n + 2*GHOSTS; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_xorder_boundary()
#if dimension > 1
                }
#endif
            }
    }
    @


@def foreach_face_xorder_boundary()
OMP_PARALLEL() {
int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
Point point = { 0 };
point.level = depth(); point.n = 1 << point.level;
int _k;
OMP(omp for schedule(static))
    for (_k = 0; _k <= point.n + 2*GHOSTS; _k++) {
#if dimension == 1
                point.i = _k;
#endif

#if dimension > 1
                point.j = _k;
                for (point.i = 0; point.i <= point.n + 2*GHOSTS+1; point.i++)
                {
#endif
                    POINT_VARIABLES
                        @
                        @def end_foreach_face_xorder_boundary()
#if dimension > 1
                }
#endif
            }
    }
    @

        @def foreach_vertex_xorder_boundary()
        foreach_face_xorder_boundary() {
        x -= Delta / 2.;
#if dimension > 1  
        y -= Delta / 2.;
#endif
        @
        @define end_foreach_vertex_xorder_boundary() } end_foreach_face_xorder_boundary()

    
event defaults(i = 0)
{
    assert(nl > 0);
    id_boundary = new scalar[nl];
    reset({id_boundary}, 0.);
}




void output_vts_bin_all_layers_multivar(scalar * list, FILE* fp)
{

    #if defined(_OPENMP)
        //int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
    #endif
        
        int count = 0;
        // MULTIGRID

        fputs("<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);

#if dimension == 1
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n", 0, N, 0, nl);
#endif

#if dimension == 2
    fprintf(fp, "\t <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, nl);
#endif

    
    // write time value
    fprintf(fp, "\t\t <FieldData> \n");
    fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", "TimeValue", t + dt_start, t + dt_start);
    fprintf(fp, "\t\t\t %.3f \n", t + dt_start);
    fprintf(fp, "\t\t\t </DataArray > \n");
    fprintf(fp, "\t\t </FieldData> \n");
    
#if dimension == 1
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d  0 0\">\n", 0, N, 0, nl);
#endif

#if dimension == 2
    fprintf(fp, "\t\t <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, N, 0, N, 0, nl);
#endif

    // Loop over velocity data and store kinematics in cell vector stucture
    fputs("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
    for (scalar s in list) {
        fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"1\" Name=\"%s\" format=\"appended\" offset=\"%d\">\n", s.name, count);
        fputs("\t\t\t\t </DataArray>\n", fp);
#if dimension == 1
        count += ((N * nl) + 1) * 8;
#endif
#if dimension == 2
        count += ((N * N * nl) + 1) * 8;
#endif
    }
    // append velocity
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"appended\" offset=\"%d\">\n", count);
    fputs("\t\t\t\t </DataArray>\n", fp);
    fputs("\t\t\t </CellData>\n", fp);
    

    // Coordinates 
    fputs("\t\t\t <Points>\n", fp);
#if dimension == 1
    count += ((N * nl * 3) + 1) * 8;
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\">\n", count);
#endif
#if dimension == 2
    count += ((N * N * nl * 3) + 1) * 8;
    fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\">\n", count);
#endif

    fputs("\t\t\t\t </DataArray>\n", fp);
    fputs("\t\t\t </Points>\n", fp);
    fputs("\t\t </Piece>\n", fp);

    fputs("\t </StructuredGrid>\n", fp);
    fputs("\t <AppendedData encoding=\"raw\">\n", fp);
    fputs("_", fp);
    
    unsigned long long block_len = 8;
    
    // Append scalars
    for (scalar s in list) {
#if dimension == 1
    block_len = 8 * N * nl;
#endif
#if dimension == 2
    block_len = 8 * N * N * nl;
#endif
        fwrite(&block_len, sizeof(unsigned long long), 1, fp);
        foreach_layer() {
            foreach_xorder() {
                //double tt = phi[]+gravity_*eta[]; 
                //fwrite(&tt, sizeof(double), 1, fp);
                fwrite(&val(s), sizeof(double), 1, fp);
            }
        }
    }
    // append velocity  
#if dimension == 1
    block_len = 8 * 3 * N * nl;
    double dummy = 0.;
#endif
#if dimension == 2
    block_len = 8 * 3 * N * N * nl;
#endif
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);
    foreach_layer() {
        foreach_xorder() {
#if dimension == 1
            fwrite(&u.x[], sizeof(double), 1, fp);
            fwrite(&w[], sizeof(double), 1, fp);
            fwrite(&dummy, sizeof(double), 1, fp);
#endif
#if dimension == 2
            fwrite(&u.x[], sizeof(double), 1, fp);
            fwrite(&u.y[], sizeof(double), 1, fp);
            fwrite(&w[], sizeof(double), 1, fp);
#endif
        }

    }
    // append points
    // Coordinates 
#if dimension == 1
    block_len = 8 * 3 * (N + 1) * (nl + 1);
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);
    // first do the bottom coordinates stored in zb
    scalar hh[];
    foreach_vertex_xorder() {
        hh[] = (zb[] + zb[-1]) / 2.;
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&hh[], sizeof(double), 1, fp);
        fwrite(&dummy, sizeof(double), 1, fp);
    }

    foreach_layer() {
        //for (scalar h_layer in hl){
        foreach_vertex_xorder() {
            hh[] += (h[] + h[-1]) / 2.;
            fwrite(&x, sizeof(double), 1, fp);
            fwrite(&hh[], sizeof(double), 1, fp);
            fwrite(&dummy, sizeof(double), 1, fp);
        }
    }
#endif
#if dimension == 2
    block_len = 8 * 3 * (N + 1) * (N + 1) * (nl + 1);
    fwrite(&block_len, sizeof(unsigned long long), 1, fp);

    // first do the bottom coordinates stored in zb
    scalar hh[];
    foreach_vertex_xorder() {
        hh[] = (zb[] + zb[-1] + zb[0, -1] + zb[-1, -1]) / 4.;
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&y, sizeof(double), 1, fp);
        fwrite(&hh[], sizeof(double), 1, fp);
    }

    foreach_layer() {
        //for (scalar h_layer in hl){
        foreach_vertex_xorder() {
            hh[] += (h[] + h[-1] + h[0, -1] + h[-1, -1]) / 4.;
            fwrite(&x, sizeof(double), 1, fp);
            fwrite(&y, sizeof(double), 1, fp);
            fwrite(&hh[], sizeof(double), 1, fp);
        }
    }
#endif


    fputs("\t\n", fp);
    fputs("\t </AppendedData>\n", fp);
    fputs("</VTKFile>\n", fp);
    fflush(fp);

    #if defined(_OPENMP)
        omp_set_num_threads(num_omp);
    #endif
}





