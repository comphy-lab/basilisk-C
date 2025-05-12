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
int swap_array[8] = {0, 4, 1, 5, 2, 6, 3, 7};
int ind[8][3];

void swap_direction (scalar * list) {
  double * data;
  unsigned * flags;
  long unsigned int n = 0, j = 0, jj = 0;
  foreach_cell(reduction(+:n))
    n++;
  data  = malloc (n*list_len(list)*sizeof(double));
  flags = malloc (n*sizeof(unsigned));
  for (int i = 0; i < 2 ; i++)
    for (int j = 0; j < 2 ; j++)
      for (int k = 0; k < 2 ; k++) {
	ind[jj][0] = i, ind[jj][1] = j, ind[jj][2] = k;
	jj++;
      }
  jj = 0;
  /**
We will traverse the tree and keep track of the number of iterated
cells at each level. Per traversal, the maximum is eight.  
   */
  int child_nr[depth() + 1];
  for (int d = 0; d <= depth(); d++)
    child_nr[d] = 0;
  
  /**
     We write our own cell iterator
  */
  Point point; point.i = point.j = point.k = GHOSTS, point.level = 0;
  
  /**
     The root is not rotated...
    */
  child_nr[point.level] = 8; //We are done at level = 0.
  write_data (point, list, data, flags, &j, &jj);
  bool next_cell = true;
  while (next_cell) {
    if (!is_leaf(cell)) { //Go to child
      point.level++;
      point.i =  2*point.i - GHOSTS + ind[swap_array[0]][0];
      point.j =  2*point.j - GHOSTS + ind[swap_array[0]][1];
      point.k =  2*point.k - GHOSTS + ind[swap_array[0]][2];
      child_nr[point.level]++; //Mark as being iterated
    } else { //a leaf: We go to a sibling or a parent...
      if (child_nr[point.level] < 8) { //Goto sibling
	int old = swap_array[child_nr[point.level] - 1];
	int new = swap_array[child_nr[point.level]];
	point.i += ind[new][0] - ind[old][0];
	point.j += ind[new][1] - ind[old][1];
	point.k += ind[new][2] - ind[old][2];
	child_nr[point.level]++;
      } else {
	while (child_nr[point.level] == 8 && point.level > 0) { //reset and move down
	  child_nr[point.level] = 0;
	  point.i = (point.i + GHOSTS)/2;
	  point.j = (point.j + GHOSTS)/2;
	  point.k = (point.k + GHOSTS)/2;
	  point.level-- ;
	}
	if (point.level > 0) {	//Goto the next sibling at this coarser level
	  int old = swap_array[child_nr[point.level] - 1];
	  int new = swap_array[child_nr[point.level]];
	  point.i += ind[new][0] - ind[old][0];
	  point.j += ind[new][1] - ind[old][1];
	  point.k += ind[new][2] - ind[old][2];
	  child_nr[point.level]++;
	} else
	  next_cell = false;
      }
    }
    if (next_cell) { //Write the data
      flags[jj++] = is_leaf(cell) ? leaf : 0;
      for (scalar s in list)
	data[j++] = s[];
    }
  }

  //"Restore"
  
  unrefine(level >= 1);
  foreach_cell(){
    if (point.level == 0)
      j = jj = 0;
    for (scalar s in list)
      s[] = data[j++];
    if (!(flags[jj++] & leaf) && is_leaf(cell))
      refine_cell(point, list, 0, NULL);
    if (is_leaf(cell))// Still No halos
      continue;
  }
}

	   


      
