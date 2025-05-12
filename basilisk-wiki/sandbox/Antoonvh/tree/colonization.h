/**
# Algoritmic botany

Space colonization algortihm of Runions et al. (2007). 

## The relevant variables (and types)

We want to have a list of branch sections that define a tree (or a
root?).
 */
#include "treegen.h"
/**
These are based on a set of tree nodes, consisting of the centered
location, their index, the index of the parent, and possible child
node indices.
 */
typedef struct{
  double x;
  double y;
  double z;
  int j;    //index
  int p;    //Parent index
  int nrc;  //nr of children (0, 1 or 2)
  int c[2]; //child indices
} Tnode;
/**
## User interface

These `tnodes` are computed from a set of attraction-point
coordinates, following the space colonisation alghorithm of Runions et
al. (2017). The function returns the number of computed tree nodes.
*/

struct Cnodes {     //Input structure for node computation
  coord * atr;      //Attraction point coordinates
  int na;           //Number of attraction points
  Tnode ** tnodes;  //The tree nodes (input and output)
  int tni;          //Number of initial tree nodes
  double di;        //radius of influence
  double D;         //displacement length 
  double dk;        //kill distance
  coord g;          //Biassing vector 
  int hnr;          //Halting number of tree nodes
 
  void (* myfun) (Tnode * tnodes, int tn, coord * atr, int na); //Analize growth
}; 
 
int coloninze (struct Cnodes C); //prototype

/**
## Helper functions

The algorithm employs the following helper functions:

* Couple influencing atr. points to nearest node, mark influenced
  nodes and mark atr. points to be removed. Returns the number of
  influenced nodes.
* Compute the vectors for each influencing attraction point.
* Add new nodes to the tree based on a list of vectors and the
  displacement distance.
* Kill nearby attraction points.

Starting with an $\mathcal{O}\left(\mathtt{na}\times\mathtt{tn}\right)$
approach:

*/
trace
int get_ind_atr (Tnode * nodes, int tn,          //tree nodes
		 coord * atr, int na, double di, //atr. points and d_i
		 int * inda,                     //nearest node foreach atr. p.
		 int * indn,                     //nr atr p. per node (if zero array input)
		 int * toremove, int * rp, double dk) {//remove indices and d_k  
  int ni = 0;    //Counter for nodes being influenced (i.e. new nodes)
  rp[0] = 0;     //counter for atr. points to remove 
  di = sq(di);   //squared distances are used
  dk = sq(dk);
  for (int a = 0;  a < na; a++) {    //for all atr. points
    double dm = HUGE;
    int nn = -1;
    for (int n = 0; n < tn; n++) {   //For all tree nodes
      double d = 0;
      foreach_dimension() 
	d += sq(nodes[n].x - atr[a].x);
      if (d < dm && nodes[n].nrc != 2) {
	dm = d;
	nn = n;
      }
    }
    if (dm < dk) {
      toremove[rp[0]] = a;
      rp[0]++;
    }
    if (dm < di && nn != -1) {
      indn[nn]++;   
      inda[a] = nn;
      if (indn[nn] == 1)
	ni++;
    } else
      inda[a] = -1; //Not influencing, maybe killed...
  }
  return ni;
}

trace
int get_vectors (Tnode * nodes, int tn, coord * atr, int * inda, coord * ni, int na) { 
  int nr = 0;
  for (int a = 0; a < na; a++)  //Foreach atr. points
    if (inda[a] != -1){         //that is influencing node inda[a]
      foreach_dimension()
	ni[a].x = atr[a].x - nodes[inda[a]].x;
      normalize (&ni[a]);
    }
  return nr;
}

trace
int add_nodes (Tnode * nodes,int * indn, int tn, coord * ni, int * inda, int na, double D) {
  coord * V = calloc (tn, sizeof(coord));
  for (int a = 0; a < na; a++) {  //Foreach atr. points
    if (inda[a] != -1) {          //that is influencing node inda[a]
      foreach_dimension() {
	V[inda[a]].x += ni[a].x;
      }
    }
  }
  int nn = 0; //new node counter
  for (int n = 0; n < tn; n++) { //Foreach node
    if (indn[n] > 0) {           //That is being influenced by indn[n] atr. points
      normalize (&V[n]);
      foreach_dimension()
	nodes[tn + nn].x = nodes[n].x + D*V[n].x;
      nodes[tn + nn].p = n;
      nodes[n].c[nodes[n].nrc++] = tn + nn;
      nodes[tn + nn].nrc = 0;
      nn++;
    }
  }
  free (V);
  return nn;
}

trace
int kill_atr (coord * atr, int na, int * toremove, int rp) {
  int j = 0;
  int m = 0;
  while (j < na - rp) {
    while (m < rp ? j + m == toremove[m] : 0)
      m++;
    while (m < rp ? j < na - rp && j + m != toremove[m] : j < na - rp) {
      atr[j]   =  atr[j + m];
      j++;
    }
  }
  return rp;
}

/**
## The implementation of the algorithm

The algorithm is readily implemented with the help of the helper
functions.
 */
int colonize (struct Cnodes C) {
  if (!C.hnr)
    C.hnr = INT_MAX;
  int tn = C.tni, na = C.na;
  while (tn < C.hnr && na > 0) {      
    int inda[C.na], toremove[C.na];
    int * indn = calloc(tn, sizeof(int));
    int rp = 0;
    int nn = get_ind_atr (*C.tnodes, tn,
			  C.atr, na, C.di,
			  inda,
			  indn,
			  toremove, &rp, C.dk);
    if (nn) {
      coord * ni = calloc (na, sizeof(coord));
      get_vectors (*C.tnodes, tn,  C.atr, inda, ni, na);
      *C.tnodes = realloc (*C.tnodes, (tn + nn)*sizeof(Tnode));
      add_nodes (*C.tnodes, indn, tn, ni, inda, na, C.D);
      free (ni);
      free (indn);
      tn += nn;
    } else  //done
      break;
    if (rp)
      na -= kill_atr (C.atr, na, toremove, rp);
    if (C.myfun)
      C.myfun (*C.tnodes, tn, C.atr, na);
  }	
  return tn; 
}

/**
## Utilities

Further utilities include:

* A function that computes the path length to the farrest leaf (growing age).
* Create a `Branch` list from nodes. 

A recursive approach is used to get the distance to leaves.

*/

int path_length (Tnode * nodes, int n, int * path) {
  if (nodes[n].nrc == 0) // a leaf
    return path[n] = 0;
  if (nodes[n].nrc == 1) {//Follow children
    int ind = nodes[n].c[0];
    if (path[ind] != -1)
      return path[n] = path[ind] + 1;
    else
      return path[n] = path_length(nodes, ind, path) + 1;  
  }
  if (nodes[n].nrc == 2) {// Follow longest branch
    int ind1 = nodes[n].c[0];
    int ind2 = nodes[n].c[1];
    return path[n] = max(path_length(nodes, ind1, path), path_length(nodes, ind2, path)) + 1;
  }
  fprintf (stderr, "#This should not happen...\n");
  return 0;
}
trace
void path_lengths (Tnode * nodes, int * path, int tn) {
  for (int n = 0; n < tn; n++) 
    path[n] = -1;
  for (int n = 0; n < tn; n++) //There could be many roots.. 
    path[0] = path_length (nodes, 0, path);
}

void putdata (Tnode * nodes,int n, Branch * l,int b, int j, int * d, double Rstem) {
  foreach_dimension() {
    l[b].start.x = nodes[n].x;
    l[b].end.x = nodes[nodes[n].c[j]].x;
  }
  // For the radius is a fraction of the stem size
  l[b].R = Rstem*max(sqrt((double)d[nodes[n].c[j]]/d[0]), 1./4.);
}
trace
int nodestotree (Tnode * nodes, int tn, Branch ** l, double rstem) {
  int s = 0;
  for (int n = 0; n < tn; n++)
    s += nodes[n].nrc;
  *l = realloc (*l, s*sizeof(Branch));
  int d[tn];
  path_lengths (nodes, d, tn);
  int b = 0;
  for (int n = 0; n < tn; n++) {
    for (int j = 0; j < nodes[n].nrc; j++) {
      putdata (nodes, n, *l, b, j, d, rstem);
      b++;
    }
  }
  return s;
}


  
 

/**
## Reference

Runions, Adam, Brendan Lane, and Przemyslaw Prusinkiewicz. *Modeling
Trees with a Space Colonization Algorithm.* NPH 7 (2007): 63-70.
 */
