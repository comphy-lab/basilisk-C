/**
These altarnatives to `dump()` are based on the existing `dump()`
functionality and the existing helper functions are re-used.
 */
#include "output.h"

/**
`dump2("name")` lets each thread write its onw file `name-`$x$. with
$x$ the `pid()` with quite possibly some leading zeros (for upto
100000 cores). Such that,

~~~bash
cat name-* > name
~~~

can be used to concatenate the files for e.g. `restore()` and `bview`. 
 */


trace
void dump2 (struct Dump p)
{
  FILE * fp;
  char * file = p.file;
  if (file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }
  char * name = NULL;
  if (file) {
    char pd[97];
    sprintf(pd, "-%05d", pid());
    name = (char *) malloc (strlen(file) + 99);
    strcpy (name, file);
    strcat (name, pd);
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  assert (fp);
  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size[];
  scalar * list = list_concat ({size}, dlist); free (dlist);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
			       dump_version };
#if MULTIGRID_MPI
  for (int i = 0; i < dimension; i++)
    (&header.n.x)[i] = mpi_dims[i];
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  if (pid() == 0)
    dump_header (fp, &header, list);

  subtree_size (size, false);
  foreach_cell() {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (is_local(cell)) {
      if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
	perror ("dump(): error while writing flags");
	exit (1);
      }
      for (scalar s in list)
	if (fwrite (&s[], sizeof(double), 1, fp) < 1) {
	  perror ("dump(): error while writing scalars");
	  exit (1);
	}
    }
      if (is_leaf(cell))
	continue;
  }
  free (list);
  if (file) {
    fclose (fp);
    free (name);
  }
}
/**
   `dump3()` shoud work well 
 */

trace
void dump3 (struct Dump p)
{
  FILE * fh = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;
#if _MPI
  if (fh != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }
#endif
  char * name = NULL;
  if (file) {
    name = (char *) malloc (strlen(file) + 2);
    strcpy (name, file);
    if (p.buffered)
      strcat (name, "~");
    if ((fh = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);    
    }
  }
  if (fh == NULL) {
    perror (name);
    exit (1);    
  }
  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size[];
  scalar * list = list_concat ({size}, dlist); free (dlist);
  struct DumpHeader header = {t, list_len(list), iter, depth(), npe(),
			      dump_version};
 
#if MULTIGRID_MPI
  for (int i = 0; i < dimension; i++)
    (&header.n.x)[i] = mpi_dims[i];
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  if (pid() == 0)
    dump_header (fh, &header, list);

  /**
For parallel writing, each thread starts with an offset.
  */
#if _MPI
  long offset = 0;
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  for (scalar s in list)
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(s.name);
  scalar index = {-1};
  index = new scalar;
  z_indexing (index, false);
  foreach_cell() {
    if (is_local(cell)) {
      if (offset == (long)0)
	offset = sizeofheader + index[]*cell_size;
      break;
    }
  }
  delete ({index});
  fseek (fh, offset, SEEK_SET);
#endif
  /**
     To help encode the tree-grid structure, the size of each subtree
     is computed and dumped for each cell.
  */
  subtree_size (size, false);
  foreach_cell() {
    if (is_local(cell)) {
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      for (scalar s in list)
	fwrite (&s[], 1, sizeof(double), fh);
    }
    if (is_leaf(cell))
      continue;
  }
  free (list);
  if (file) {
    fclose (fh);
    if (p.buffered && pid() == 0)
      rename (name, file);
    free (name);
  }
}
/**
`dump4()` employs functions for output defined by the MPI library.

So-called 'hints' may be passed by modifying `MPI_INFO` for optimal
parallel output performance.
 */

#if _MPI
MPI_Info MPI_INFO;

static void MPI_dump_header (MPI_File fp, struct DumpHeader * header, scalar * list) {
  MPI_Status status;
  MPI_File_write (fp, header, sizeof(struct DumpHeader), MPI_BYTE , &status);
  for (scalar s in list) {
    unsigned len = strlen (s.name);
    MPI_File_write (fp, &len, 1, MPI_UNSIGNED, &status);
    MPI_File_write (fp, s.name, len, MPI_CHAR, &status);
  }
  double o[4] = {X0, Y0, Z0, L0};
  MPI_File_write (fp, o, 4, MPI_DOUBLE, &status);
}

trace
void MPI_dump (struct Dump p) {
  if (p.fp) {
    fprintf (ferr, "You should provide a file name instead of a pointer\n");
    exit (1);
  }
  MPI_File  fh;
  char * file = p.file;
  if (file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI_dump\n");
    exit(1);
  }
  char name[strlen(file) + 2];
  strcpy (name, file);
  if (MPI_INFO == NULL) 
    MPI_File_open (MPI_COMM_WORLD, name,
		   MPI_MODE_WRONLY | MPI_MODE_CREATE,
		   MPI_INFO_NULL, &fh);
  else
    MPI_File_open (MPI_COMM_WORLD, name,
		   MPI_MODE_WRONLY | MPI_MODE_CREATE,
		   MPI_INFO, &fh);
  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size[];
  scalar * list = list_concat ({size}, dlist);
  free (dlist);
  struct DumpHeader header = {t, list_len(list), iter, depth(),
			      npe(), dump_version };
#if MULTIGRID_MPI
  for (int i = 0; i < dimension; i++)
    (&header.n.x)[i] = mpi_dims[i];
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  if (pid() == 0)
    MPI_dump_header (fh, &header, list);
  scalar index = {-1};
  index = new scalar;
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double); 
  for (scalar s in list)
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(s.name);
  subtree_size (size, false);
  long offset = 0;
  foreach_cell() { //We find the index of the first local cell
    if (is_local(cell)) {
      if (offset == (long)0)
	offset = sizeofheader + index[]*cell_size;
      break; 
    }
  }
  MPI_File_seek (fh, offset , MPI_SEEK_SET);
  MPI_Status status;
  foreach_cell() {
    if (is_local(cell)) {
      unsigned flags = is_leaf(cell) ? leaf : 0;
      MPI_File_write(fh, &flags,1,MPI_UNSIGNED, &status);
      for (scalar s in list)
	MPI_File_write(fh, &s[],1,MPI_DOUBLE, &status);
    }
    if (is_leaf(cell))
      continue;
  }
  delete ({index});
  free (list);
  MPI_File_close (&fh);
}
#endif
