/**
## Write the left BC from a precurser run
If we need the inflow BC from a precurser run, we write the left BC of a variable s to the disk using MPI and read the data to new simulation. By filling in left BC of a new variable in new simulation, we can prescribe: b[left] = s[-1,0]; This function is adapted from [Antoon's toolbox](http://www.basilisk.fr/sandbox/Antoonvh/muhhtob.h)
*/

int write_inflow(scalar s, char *fname, int dlevel){
#if !_MPI
        FILE *fpm;
        fpm = fopen(fname, "wb");
#elif _MPI
        MPI_File fpm;
        MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fpm);
#endif
        int siz = 8;
        int g = 0;
        int o = -1 - BGHOSTS;
        int CR = (1 << dlevel);
        foreach_boundary(left, reduction(+:g))
        {
                // first z direction, then y direction
                g = point.k + o + CR * (point.j + o);
#if !_MPI
                fseek(fpm, g * siz, SEEK_SET);
                fwrite(&s[-1, 0], siz, 1, fpm);
#elif _MPI
                MPI_File_seek(fpm, g * siz, MPI_SEEK_SET);
                MPI_Status status;
                MPI_File_write(fpm, &s[-1, 0], 1, MPI_DOUBLE, &status);
#endif
        }
#if _MPI
        MPI_File_close(&fpm);
#endif
        return 0;
}

/**
## Read the left BC from a precurser run
*/

int read_inflow(scalar s, char *fname, int dlevel)
{
#if !_MPI
        FILE *fpm;
        if (!(fpm = fopen(fname, "rb")))
        {
                fprintf(stderr, "Cound not open file with name '%s'\n", fname);
                return 1;
        }
#elif _MPI
        MPI_File fpm;
        MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpm);
#endif
        int siz = 8;
        int g = 0;
        int o = -1 - BGHOSTS;
        int CR = (1 << dlevel);
        foreach_boundary(left, reduction(+:g))
        {
                // first z direction, then y direction
                g = point.k + o + CR * (point.j + o);
#if !_MPI
                fseek(fpm, g * siz, SEEK_SET);
                fread(&s[-1, 0], siz, 1, fpm);
#elif _MPI
                MPI_File_seek(fpm, g * siz, MPI_SEEK_SET);
                MPI_Status status;
                MPI_File_read(fpm, &s[-1, 0], 1, MPI_DOUBLE, &status);
#endif
        }
#if _MPI
        MPI_File_close(&fpm);
#endif
        return 0;
}

/**
Since the idealized percurser run is normally prescribed with periodic boundary condition and the ghost points are not available, we take the first slice from left hand side. This function can take the slice from tree grid with the price of refine left boundary condition. 
*/
void write_inflow_pbc(scalar s, char *fname, int dlevel)
{
#if !_MPI
    FILE *fpm;
    fpm = fopen(fname, "wb");
#elif _MPI
    MPI_File fpm;
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fpm);
#endif
    int siz = 8;
    int g = 0;
    int o = -1 - BGHOSTS;
    int CR = (1 << dlevel);
    // here i cheat by refine the grid near the left boundary
    refine(x <= Delta && level == (dlevel-1));
    foreach(reduction(+:g))
    {
        if ((point.i+o) == 0){
        // first z direction, then y direction
        g = point.k + o + CR * (point.j + o);
#if !_MPI
        fseek(fpm, g * siz, SEEK_SET);
        fwrite(&s[], siz, 1, fpm);
#elif _MPI
        MPI_File_seek(fpm, g * siz, MPI_SEEK_SET);
        MPI_Status status;
        MPI_File_write(fpm, &s[], 1, MPI_DOUBLE, &status);
#endif
        }
    }
#if _MPI
    MPI_File_close(&fpm);
#endif
}

