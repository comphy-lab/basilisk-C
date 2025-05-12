/**
# Generate a fractal tree

   Stem parameters:
 */
coord stem_pos; 
double stem_rad = 1.; //Base size
/**
   Length-over-radius ratio, and number of fractal levels. The stem is at
level = 0.
 */
double LoR = 20., LoR_spread = 5.; //LoR_i = LoR + noise()*LoR_spread
int levels = 1;
/**
Branch parameters: Number of branches will be rounded to the nearest
integer after the spreading randomization. The `new_angle` parameter
is the angle relative to the parent branch (100% randomized rotation).
 */
double fractal_factor = 0.7, fractal_factor_spread = 0.1; 
double split_branches = 2.1, split_branches_spread = 0.5; //Starting at the end of parent  
double side_branches  = 4.1, side_branches_spread = 0.5;  //Starting along the parent
double side_start = 0.65, side_start_spread = .3;         //Relative along-parent position 

double new_angle = 40.*pi/180., new_angle_spread = 10.*pi/180.; 

/**
That concludes the inputs.
*/

typedef struct {
  coord start, end;
  double R;
  int parent;
} Branch;

unsigned int nb = 0; //Number of segments

coord cross (coord a, coord b) {
  coord new;
  new.x = a.y*b.z - a.z*b.y;
  new.y = a.z*b.x - a.x*b.z;
  new.z = a.x*b.y - a.y*b.x;
  return new;
}

void normalize3 (coord * n)
{
  double norm = 0.;
  foreach_dimension(3)
    norm += sq(n->x);
  norm = sqrt(norm);
  foreach_dimension(3)
    n->x /= norm;
}

unsigned int add_bones (Branch * branches, unsigned int p, unsigned int j) {
  Branch par = branches[p];
  int indp = j;
  int sidenr  = (int)round(side_branches + side_branches_spread*noise()) + 0.5;
  double rndm_ngl = pi*noise();
  
  coord  randomvec = {noise(), noise(), noise()};
  for (int m = 0; m < sidenr; m++) {
    j++;
    branches[j].parent = indp;
    double sd_strt = side_start + noise()*side_start_spread;
    foreach_dimension(3) 
      branches[j].start.x = par.start.x + sd_strt*(par.end.x-par.start.x);
    double frcl_fctr = fractal_factor + fractal_factor_spread*noise();
    branches[j].R = par.R*frcl_fctr;
    double length = (LoR + noise()*LoR_spread)*branches[j].R;
    double ngl = new_angle + new_angle_spread*noise();
    rndm_ngl += 2.*pi/((double)sidenr) + new_angle_spread*noise();
    coord vec;
    foreach_dimension(3) {
      vec.x = (par.end.x - par.start.x);
    }
/**
We need to rotate this vector with angle `ngl` and `rndm_ngl`
(pitch and yaw)

This was helpfull:

tparker, *Rotate a vector by a randomly oriented angle*, URL (version:
2017-09-30): [https://scicomp.stackexchange.com/q/27969]()
*/
    normalize3 (&vec);
    coord vecx = cross (vec, randomvec);
    normalize3 (&vecx);
    coord vecy = cross (vec, vecx);
    coord vecn;
    foreach_dimension(3)
      vecn.x = (sin(ngl)*cos(rndm_ngl)*vecx.x +
		sin(ngl)*sin(rndm_ngl)*vecy.x +
		cos(ngl)*vec.x);
    foreach_dimension(3)
      branches[j].end.x = branches[j].start.x + vecn.x*length;
    foreach_dimension(3)
      branches[j].start.x += vecn.x*par.R; //away from parent center line
  }
    int splitnr = (int)round(split_branches + split_branches_spread*noise()) + 0.5;
    rndm_ngl = pi*noise();
    for (int m = 0; m < splitnr; m++) {
      j++;
    branches[j].parent = indp;
    branches[j].start = par.end;
    double frcl_fctr = fractal_factor + fractal_factor_spread*noise();
    branches[j].R = par.R*frcl_fctr;
    double length = (LoR + noise()*LoR_spread)*branches[j].R;
    double ngl = new_angle + new_angle_spread*noise();
    rndm_ngl += 2*pi/((double)splitnr) + new_angle_spread*noise();
    coord vec;
    foreach_dimension(3)
      vec.x = (par.end.x - par.start.x);
    normalize3 (&vec);
    coord vecx = cross (vec, randomvec);
    normalize3 (&vecx);
    coord vecy = cross (vec, vecx);
    coord vecn;
    foreach_dimension(3)
      vecn.x = (sin(ngl)*cos(rndm_ngl)*vecx.x +
		cos(ngl)*sin(rndm_ngl)*vecy.x +
		cos(ngl)*vec.x);
    foreach_dimension(3)
      branches[j].end.x = branches[j].start.x + vecn.x*length;
    foreach_dimension(3)
      branches[j].start.x += vecn.x*par.R; //away from parent center line
  }
  return splitnr + sidenr;
}

Branch * tree_skeleton () {
  Branch * branches;
  //array or branches
  int maxb = (round(split_branches + split_branches_spread) +
	      round(side_branches  + side_branches_spread));
  unsigned int maxnl = (unsigned int)(pow(maxb, levels + 1) + 3);
  branches  = (Branch*) malloc(sizeof(Branch)*maxnl);
  // init stem
  int j = 0;
  branches[j].start = stem_pos;
  branches[j].R = stem_rad;
  coord ey = {0, 1, 0};
  double length = (LoR + noise()*LoR_spread)*branches[j].R;
  foreach_dimension(3)
    branches[j].end.x = stem_pos.x + length*ey.x;
  branches[j].parent = 0;
  nb ++;      
  int pc = 0; 
  
  for (int l = 0 ; l < levels; l++) {
    // At this level, all branches are parents
    while (pc < nb) {
      int n = add_bones (branches, pc, j); //parent pc to child j 
      pc++;
      j += n;
    }
    nb = j + 1;
  }
  return branches;	      
}

/**
## Skeleton to volume fractions

 */
#include "PointTriangle.h" //PointSegmentDistance exists!
#include "fractions.h"
struct tree_int {
  Branch * branches;
  scalar c;
  face vector fs;
  scalar J;
  double smooth;
  scalar * alist;  //For adaptation
  double * crit;   //...
  int minlevel;
  int maxlevel;
  scalar * ulist;
  int stop;
  double stopc;
};

#define SMOOTHR ((d - branches[j].R) > 0 ?exp (-sq((d - branches[j].R)/(branches[j].R*ti.smooth))) : 1)
//#define SMOOTHR ((d - branches[j].R) > 0 ? ((d - branches[j].R)/branches[j].R < pi/2 ? (cos (d - branches[j].R)/branches[j].R) : 0) : 1)
trace
void tree_interface (struct tree_int ti) {
  Branch * branches = ti.branches;
  vertex scalar phi[], phiR[];
  scalar J = automatic (ti.J);
  
  if (!ti.smooth) {
    foreach_vertex() {
      phi[] = HUGE;
      for (int j = 0; j < nb; j++) {
	coord pos = (coord){x, y, z};
	coord p0 = branches[j].start;
	coord p1 = branches[j].end;
	coord Sc;
	double  sP;
	double d = sqrt(PointSegmentDistance (&pos, &p0, &p1, &Sc, &sP))
	  - branches[j].R;
	if (d < phi[]) {
	  phi[] = d;
	  if (ti.J.i)
	    J[] = j;
	}
      }
    }
  } else {
    foreach_vertex() {
      phi[] = 0;
      double mind = HUGE, phiR = 0;
      for (int j = 0; j < nb; j++) {
	coord pos = (coord){x, y, z};
	coord p0 = branches[j].start;
	coord p1 = branches[j].end;
	coord Sc;
	double sP;
	double d = sqrt(PointSegmentDistance (&pos, &p0, &p1, &Sc, &sP));
	
	phi[] += d > 0 ? sq(1/d) : HUGE ;
	phiR  += SMOOTHR*branches[j].R;
	if (ti.J.i)
	  if (d < mind) {
	    mind = d;
	    J[] = j;
	  }
      }
      phi[] = sqrt(1./phi[]) - phiR;
    }
  }
  boundary ({phi});
  fractions (phi, ti.c, ti.fs);
}


/**
# Skeleton to fractions with refinement 

 */

void refine_vertex (Point point, scalar s) {
  for (int i = 0; i < 2; i++) 
    for (int j = 0; j < 1 + (dimension > 1); j++) 
      for (int k = 0; k < 1 + (dimension > 2); k++) { 
	fine (s, 2*i, 2*j, 2*k) = s[i, j, k];
	fine (s, 2*i + 1, 2*j    , 2*k)     = nodata; //1
#if (dimension > 1)
	fine (s, 2*i    , 2*j + 1, 2*k)     = nodata; //1
	fine (s, 2*i + 1, 2*j + 1, 2*k)     = nodata; //2
#if (dimension > 2)
	fine (s, 2*i    , 2*j    , 2*k + 1) = nodata; //1
	fine (s, 2*i + 1, 2*j    , 2*k + 1) = nodata; //2
	fine (s, 2*i    , 2*j + 1, 2*k + 1) = nodata; //2
	fine (s, 2*i +1 , 2*j + 1, 2*k + 1) = nodata; //3
#endif 
#endif
      }
}

void tree_interface_adapt (struct tree_int ti) {
  Branch * branches = ti.branches;
  vertex scalar phi[];
  phi.refine = refine_vertex;
  foreach_vertex()
    phi[] = nodata;

  scalar * update = NULL;
  if (!ti.ulist)
    update = list_copy (all);
  else {
    update = list_copy (ti.ulist);
    update = list_add(update, phi);
  }
  
  int it = 0, computed;
  do {
    computed = 0;
    if (ti.stopc)
      ti.stop = grid->tn/(int)ti.stopc;
    if (!ti.smooth) {
      foreach_vertex(reduction (+:computed)) {
	if (phi[] == nodata) {
	  computed++;
	  for (int j = 0; j < nb; j++) {
	    coord pos = (coord){x, y, z};
	    coord p0 = branches[j].start;
	    coord p1 = branches[j].end;
	    coord Sc;
	    double  sP;
	    double d = sqrt(PointSegmentDistance (&pos, &p0, &p1, &Sc, &sP))
	      - branches[j].R;
	    if (d < phi[]) 
	      phi[] = d;
	  }
	}
      }
    } else {
      foreach_vertex(reduction (+:computed)) {
	if (phi[] == nodata) {
	  computed++;
	  phi[] = 0;
	  double phiR = 0;
	  for (int j = 0; j < nb; j++) {
	    coord pos = (coord){x, y, z};
	    coord p0 = branches[j].start;
	    coord p1 = branches[j].end;
	    coord Sc;
	    double sP;
	    double d = sqrt(PointSegmentDistance (&pos, &p0, &p1, &Sc, &sP));
	    phi[] += d > 0 ? sq(1/d) : HUGE;
	    phiR  += SMOOTHR*branches[j].R;
	  }
	  phi[] = sqrt(1./phi[]) - phiR;
	}
      }
    }
    boundary ({phi});
    fractions (phi, ti.c, ti.fs);
    scalar c = ti.c;
    face vector fss = ti.fs;
    boundary ({c, fss});
#if EMBED
    if (fractions_cleanup (ti.c, ti.fs))
      boundary ({c, fss});
#endif
    
    if (pid() == 0)
      printf ("# Tree Reconstruction Adapt Iteration %d: %ld cells\n"
	      "#      and %d newly computed vertices\n", ++it, grid->tn, computed);
  } while (adapt_wavelet (ti.alist, ti.crit, ti.maxlevel, ti.minlevel, update).nf > ti.stop);
  scalar c = ti.c;
  face vector fss = ti.fs;
  boundary ({c, fss});
  // Last time for R
  if (ti.J.i) {
    tree_interface (ti);
    boundary ({c, fss});
#if EMBED
    if (fractions_cleanup (ti.c, ti.fs))
      boundary ({c, fss});
#endif
  }
}
