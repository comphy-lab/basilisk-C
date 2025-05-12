/**
Sometimes it may be usefull to extend the domain to get increased statistics. 
This functions copies a grid and solution $2^{dimension-1}$ times and extends it over the bottom boundary.

See the [example](enlarge_example.c) for reference.
*/

void enlarge_hor(scalar * list){
  double * data;
  unsigned * flags;
  long unsigned int n = 0, j = 0, jj = 0;
  int len = list_len(list);
  foreach_cell(reduction(+:n))
    n++;
  data = malloc (n*len*sizeof(double));
  flags = malloc (n*sizeof(unsigned));
  foreach_cell(){ //"Dump" to arrays
    flags[jj++] = is_leaf(cell) ? leaf : 0;
    for (scalar s in list)
      data[j++] = s[];
    if (is_leaf(cell)) // No halo data
      continue;
  }
  unrefine(level >= 1);
  foreach_cell(){
    if (point.level >= 1 && y < Y0 + L0/2.){
      if (point.level == 1)
	j = jj = 0;
      for (scalar s in list)
	s[] = data[j++];
      if (!(flags[jj++] & leaf) && is_leaf(cell))
	refine_cell(point, list, 0, NULL);
      if (is_leaf(cell))// Still No halos
	continue;
    }
  }
  free(data); free(flags);
  L0 *= 2.; X0 *= 2.; Z0 *= 2.;
}
/**
This very similar function extends the domain in all `dimension` directions.
See this [example](turbulence2D.c)
*/
void enlarge_all(scalar * list){
  double * data;
  unsigned * flags;
  long unsigned int n = 0, j = 0, jj = 0;
  int len = list_len(list);
  foreach_cell(reduction(+:n))
    n++;
  data = malloc (n*len*sizeof(double));
  flags = malloc (n*sizeof(unsigned));
  foreach_cell(){ //"Dump" to arrays
    flags[jj++] = is_leaf(cell) ? leaf : 0;
    for (scalar s in list)
      data[j++] = s[];
    if (is_leaf(cell)) // No halo data
      continue;
  }
  unrefine(level >= 1);
  foreach_cell(){
    if (point.level >= 1){ // Do not start at root cell
      if (point.level == 1)
	j = jj = 0;
      for (scalar s in list)
	s[] = data[j++];
      if (!(flags[jj++] & leaf) && is_leaf(cell))
	refine_cell(point, list, 0, NULL);
      if (is_leaf(cell))// Still No halos
	continue;
    }
  }
  free(data); free(flags);
  foreach_coarse_level(0){ //Also restrict the root cell
    for (scalar s in list){
      double sum = 0;
      foreach_child()
	sum += s[];
      s[] = sum /(double)(1 << dimension);
    }
  }
  L0 *= 2.; X0 *= 2.; Y0 *= 2.; Z0 *= 2.;
}
/**
Child points   ->  rotated
1st = -1 -1 -1 -> -1 -1 -1 = 1st  
2nd = -1 -1  1 ->  1 -1 -1 = 5th  
3rd = -1  1 -1 -> -1 -1  1 = 2nd   
4th = -1  1  1 ->  1 -1  1 = 6th  
5th =  1 -1 -1 -> -1  1 -1 = 3rd  
6th =  1 -1  1 ->  1  1 -1 = 7th  
7th =  1  1 -1 -> -1  1  1 = 4th  
8th =  1  1  1 ->  1  1  1 = 8th  
 */
int rotate_all[8] = {0, 4, 1, 5, 2, 6, 3, 7}, ind[8][3], * child_nr;

void write_data (Point point, double * data, unsigned * flags, scalar * list,
		 long unsigned int j[1], long unsigned int jj[1]) {
  flags[jj[0]++] = is_leaf(cell) ? leaf : 0;
  for (scalar s in list)
    data[j[0]++] = s[];
}

Point goto_sibling (Point point) {
  int old = rotate_all[child_nr[point.level]- 1];
  int new = rotate_all[child_nr[point.level]];
  point.i += ind[new][0] - ind[old][0];
  point.j += ind[new][1] - ind[old][1];
  point.k += ind[new][2] - ind[old][2];
  return point;
}

void swap_direction (scalar * list) {
  double * data;
  unsigned * flags;
  long unsigned int n = 0, j = 0, jj = 0;
  foreach_cell(reduction(+:n))
    n++;
  child_nr = (int*)      calloc (depth() + 1, sizeof(int));
  data     = (double*)   malloc (n*list_len(list)*sizeof(double));
  flags    = (unsigned*) malloc (n*sizeof(unsigned));
  for (int i = 0; i < 2 ; i++)
    for (int j = 0; j < 2 ; j++)
      for (int k = 0; k < 2 ; k++)
	ind[jj][0] = i, ind[jj][1] = j, ind[jj++][2] = k;
  jj = 0;
  /**
We will traverse the tree and keep track of the number of iterated
cells for the current traversal at each level (`child_nr[]`).
  */
  Point point; point.i = point.j = point.k = GHOSTS, point.level = 0;
  write_data (point, data, flags, list, &j, &jj);
  child_nr[point.level] = 8; //We are done at level = 0.
  
  while (1) {
    if (!is_leaf(cell)) {//non-leaf: Go to child
      point.level++;
      point.i +=  point.i - GHOSTS + ind[rotate_all[0]][0];
      point.j +=  point.j - GHOSTS + ind[rotate_all[0]][1];
      point.k +=  point.k - GHOSTS + ind[rotate_all[0]][2];
    } else { //Leaf: Go to a sibling ...
      if (child_nr[point.level] < 8) 
	point = goto_sibling (point);
      else { //... or a parent
	while (child_nr[point.level] == 8 && point.level > 0) {
	  child_nr[point.level] = 0;
	  point.i = (point.i + GHOSTS)/2;
	  point.j = (point.j + GHOSTS)/2;
	  point.k = (point.k + GHOSTS)/2;
	  point.level-- ;
	}
	if (point.level > 0)	// Goto the next sibling at this coarser level
	  point = goto_sibling (point);
	else                    // Back the root level...
 	  break;
      }
    }
    child_nr[point.level]++; //Mark this child as iterated
    write_data (point, data, flags, list, &j, &jj);
  }
  
  //"Restore" from the arrays
  unrefine(level > 0);
  foreach_cell() {
    if (point.level == 0)
      j = jj = 0;
    for (scalar s in list)
      s[] = data[j++];
    if (!(flags[jj++] & leaf) && is_leaf(cell))
      refine_cell(point, list, 0, NULL);
    if (is_leaf(cell)) //No halos
      continue;
  }
  free (flags);
  free (data);
  free (child_nr);
}

	   


      
