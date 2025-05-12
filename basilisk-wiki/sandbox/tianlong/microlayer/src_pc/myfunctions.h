#ifndef _MY_FUNCTIONS_H
#define _MY_FUNCTIONS_H

void outputFacetsParallel(char * filename)
{
  char names[80];
  sprintf(names,"interfaces%d", pid());
  FILE * fp = fopen (names,"w");
  output_facets(f,fp);
  fclose(fp);
#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  char command[80];
  sprintf(command, "cat interfaces* > %s", filename);
  system(command);
#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  sprintf(command, "rm interfaces*");
  system(command);
}
#if USE_MY_SOLID
void outputWallTemperature(char * filename, scalar Tout)
{
  FILE *out = fopen(filename, "w");
  extern int Ncell;
  extern double Lsize;
  double delta = Lsize / (double)(Ncell);
  double ystart = 0.5 * delta;
  double yend = ystart + Lsize - SOLID_LEN_y - delta;
  int len = 0;
  for (double pp = ystart; pp <= yend; pp += delta)
  {
    ++len;
  }

  double *TT = (double *)calloc(len, sizeof(double));
  double *pp = (double *)calloc(len, sizeof(double));

  for (int ii = 0; ii < len; ++ii)
  {
    double pos = ystart + ii * delta;
    double Tmp = interpolate(Tout, 0, pos);
    Tmp = Tmp == nodata ? 0.0 : Tmp;
    pp[ii] = pos;
    TT[ii] = Tmp;
  }

#if _MPI
  MPI_Allreduce(MPI_IN_PLACE, TT, len, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  for (int ii = 0; ii < len; ++ii)
  {
        if (pid() == 0)
        {
          fprintf(out, "%g %g %g %g\n", pp[ii], TT[ii], delta_q(pp[ii], t), t);
          fflush(out);
        }
  }

#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  free(TT);
  free(pp);
  fclose(out);
}
void outputWallCentral(char * filename, scalar Tout)
{
  FILE *out = fopen(filename, "w");
  extern int Ncell;
  extern double Lsize;
  double delta = Lsize / (double)(Ncell);
  double xstart = 0.5 * delta;
  double xend = xstart + Lsize - SOLID_LEN_x - delta;
  int len = 0;
  for (double pp = xstart; pp <= xend; pp += delta)
  {
    ++len;
  }

  double *TT = (double *)calloc(len, sizeof(double));
  double *pp = (double *)calloc(len, sizeof(double));

  for (int ii = 0; ii < len; ++ii)
  {
    double pos = xstart + ii * delta;
    double Tmp = interpolate(Tout, pos, 0.5 * delta);
    Tmp = Tmp == nodata ? 0.0 : Tmp;
    pp[ii] = pos;
    TT[ii] = Tmp;
  }

#if _MPI
  MPI_Allreduce(MPI_IN_PLACE, TT, len, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  for (int ii = 0; ii < len; ++ii)
  {
        if (pid() == 0)
        {
          fprintf(out, "%g %g %g %g\n", pp[ii], TT[ii], delta_q(pp[ii], t), t);
          fflush(out);
        }
  }

#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  free(TT);
  free(pp);
  fclose(out);
}
#endif
void setJumpVarAdaptation(scalar var, int num)
{
  if(num == 0) return;
  scalar new[];
  foreach()
  {
    new[] = var[];
  }

  for (int ii = 0; ii < num * 2; ++ii)
  {
    foreach ()
    {
      if (var[] == 0. || var[] == 1.)
      {
        bool is_neigh = false;
        foreach_neighbor(1)
        {
          bool is_finest = point.level == grid->maxdepth;
          const double ferror = 1.0e-10;
          if (fabs(var[]) > ferror && fabs(var[]) < (1. - ferror))
          {
            is_neigh = true;
          }
          // if (fabs(var[]) > 0. && fabs(var[]) < 1. && is_finest)
          // {
          //   is_neigh = true;
          // }
        }
        if (is_neigh)
        {
          double rnumb = 0.7 * pow(-1, ii);
          new[] = rnumb;
        }
      }
    }
    foreach ()
    {
      var[] = new[];
    }
  }
}

//To avoid strange problems after reading the initial temperature profile
trace
bool restoreModified (struct Dump p)
{
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    return false;
  assert (fp);

  struct DumpHeader header;  
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }

#if TREE
  init_grid (1);
  foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }
  tree->dirty = true;
#else // multigrid
#if MULTIGRID_MPI
  if (header.npe != npe()) {
    fprintf (ferr,
	     "restore(): error: the number of processes don't match:"
	     " %d != %d\n",
	     header.npe, npe());
    exit (1);
  }
  dimensions (header.n.x, header.n.y, header.n.z);
  double n = header.n.x;
  int depth = header.depth;
  while (n > 1)
    depth++, n /= 2;
  init_grid (1 << depth);
#else // !MULTIGRID_MPI
  init_grid (1 << header.depth);
#endif
#endif // multigrid

  bool restore_all = (p.list == all);
  scalar * list = dump_list (p.list ? p.list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (list)) {
      fprintf (ferr,
	       "restore(): error: the list lengths don't match: "
	       "%ld (file) != %d (code)\n",
	       header.len - 1, list_len (list));
      exit (1);
    }
  }
  else { // header.version != 161020
    if (header.version != dump_version) {
      fprintf (ferr,
	       "restore(): error: file version mismatch: "
	       "%d (file) != %d (code)\n",
	       header.version, dump_version);
      exit (1);
    }
    
    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
	fprintf (ferr, "restore(): error: expecting len\n");
	exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
	fprintf (ferr, "restore(): error: expecting s.name\n");
	exit (1);
      }
      name[len] = '\0';

      if (i > 0) { // skip subtree size
	bool found = false;
	for (scalar s in list)
	  if (!strcmp (s.name, name)) {
	    input = list_append (input, s);
	    found = true; break;
	  }
	if (!found) {
	  if (restore_all) {
	    scalar s = new scalar;
	    free (s.name);
	    s.name = strdup (name);
	    input = list_append (input, s);
	  }
	  else
	    input = list_append (input, (scalar){INT_MAX});
	}
      }
    }
    free (list);
    list = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin (o[0], o[1], o[2]);
    size (o[3]);
  }

#if MULTIGRID_MPI
  long cell_size = sizeof(unsigned) + header.len*sizeof(double);
  long offset = pid()*((1 << dimension*(header.depth + 1)) - 1)/
    ((1 << dimension) - 1)*cell_size;
  if (fseek (fp, offset, SEEK_CUR) < 0) {
    perror ("restore(): error while seeking");
    exit (1);
  }
#endif // MULTIGRID_MPI
  
  scalar * listm = is_constant(cm) ? NULL : (scalar *){fm};
#if TREE && _MPI
  restore_mpi (fp, list);
#else
  foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    // skip subtree size
    fseek (fp, sizeof(double), SEEK_CUR);
    for (scalar s in list) {
      double val;
      if (fread (&val, sizeof(double), 1, fp) != 1) {
	fprintf (ferr, "restore(): error: expecting a scalar\n");
	exit (1);
      }
      if (s.i != INT_MAX)
	s[] = val;
    }
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, listm, 0, NULL);
    if (is_leaf(cell))
      continue;
  }
  for (scalar s in all)
    s.dirty = true;
#endif
  
  scalar * other = NULL;
  for (scalar s in all)
    if (!list_lookup (list, s) && !list_lookup (listm, s))
      other = list_append (other, s);
  reset (other, 0.);
  free (other);
  
  free (list);
  if (file)
    fclose (fp);

  // // the events are advanced to catch up with the time  
  // while (iter < header.i && events (false))
  //   iter = inext;
  // events (false);
  // while (t < header.t && events (false))
  //   t = tnext;
  // t = header.t;
  // events (false);
  
  return true;
}
#endif //_MY_FUNCTIONS_H