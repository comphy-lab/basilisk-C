/**
# Attempt to Create simply connected regions for domain decomposition

![Domain decomposition for various theads, and varying amount of wiggle room](simply-connected/dd.mp4)

We may tune the load MPI balancing a bit, hoping for simply connected
regions. For the routine to be easy to evaluate and "MPI compatible",
it is not perfect. It only reduces the chance of broken domains. The
concept is presented here, using a serial routine.

## Identifying the locations where a domain maybe "cut"

When a cell is the top-right(-front) child from a lineage of top-right
parents (if you know what I mean), it is likely that the next cell in
the N-order is "far". As such, we compute the "lineage
depth". Further, we check if the next cell is vertically down (`vert`
= true), or horizontally left (`vert` = false). As the vertical jumps
are bigger.
*/

int lineage_depth (Point point, bool * vert) {
  int tzx = __builtin_ctz(point.i - BGHOSTS);
  int tzy = __builtin_ctz(point.j - BGHOSTS);
  if (tzy > tzx)
    *vert = true;
  else
    *vert = false;
  return min(tzx, tzy);
}
/**
The lineage depth is part of domain-splitting likelyness
grey-scale. The leafs' level and the area (volume) spanned by the
mpi domain should be taken into account. We may relate level and area in a simple manner.
 */
double area_equivalent_level (double area) {
  int divisions = 1 << dimension;//2, 4, or 8
  return -log(area/sq(L0))/log(divisions);
}
/**
We can now check if a cell is a potential "domain breaker"

![Domain breaker points on a $64^2$ mesh with 17 cores](simply-connected/db.png)
 */
bool domain_breaker (Point point, double area) {
  assert (area > 0);
  bool vert;
  int ld = lineage_depth (point, &vert);
  if (ld) {
    if ((level - ld) <= area_equivalent_level(area) + (vert == true) + 0) // 0 can be tuned at the expense 
      return true;
  }
  return false;
}
/**
## Deciding to unblance, hoping it aid to simply connect domains.

We allow the domain to be unblanced by some fraction, this
`wiggle_room`, is needed to make domains start and end at the
potential domain_breaker locations. It could be that the wiggle room
is too little to do anything about it. It could also be that the unbalance does not fix broken domain.
 */
double wiggle_room = 0.1; //Load in a domain can deviate up to 20% from the balanced distribution. 

/**
The location of the domain cut is marked with the cumulative load in its last cell
 */
double load_for_pid(int pid, int n, scalar load_cumu, double max_load, double area) {
  double tg = (pid + 1)*max_load/n;
  double min_diff = HUGE;
  double nl_breaker = HUGE;
  double load_balanced = HUGE;
  double min_diff_bal = HUGE;
  foreach() {
    if (domain_breaker(point, area)) {
      if (min_diff > fabs(load_cumu[] - tg)) {
	min_diff = fabs(load_cumu[] - tg);
	nl_breaker = load_cumu[];
      }
    }
    if (min_diff_bal > fabs(load_cumu[] - tg)) {
      min_diff_bal = fabs(load_cumu[] - tg);
      load_balanced = load_cumu[];
    }
  }
  if (min_diff < wiggle_room*max_load/n) {
    return nl_breaker;
  }
  return load_balanced;
}
/**
## Setup 
 */
#include "view.h"
#include "utils.h"
#include "tag.h"
int main() {
  X0 = Y0 = -L0/2.;
  init_grid (N);
  // Scalar field for the cummulative load, and the domain decomposition
  scalar zi[], dd[];
  double zv = 0;
  int n = 17;
  foreach() {
    zi[] = zv;
    zv++;
    // Check domain breaker points
    dd[] = domain_breaker(point, sq(L0)/n);
  }
  output_ppm (dd, file = "db.png", n = 512);
  // Loop over thread count and wiggle room range
  for (int n = 7; n < 16; n += 2) {
    for (wiggle_room = 0.1; wiggle_room < 0.35; wiggle_room += 0.01) {
      for (int i = n - 1; i >= 0; i--) {
	double ldp = load_for_pid(i, n, zi, zv, L0/sq(L0)/n);
	foreach() {
	  if (zi[] <= ldp)
	    dd[] = i;
	}
      }

      squares ("dd", min = 0, max = n -1);
      char str[99];
      sprintf (str, "npe() = %d", n);
      draw_string (str, size = 30);
      sprintf (str, "wiggle room = %g %%", wiggle_room*100);
      draw_string (str, pos = 1, size = 25);

      // Count broken domains.
      int tot = 0;
      scalar t[];
      for (int i = 0; i < n; i++) {
	foreach() 
	  t[] = (fabs(dd[] - i) < 1e-3) ? 1: 0;
	if (tag(t) > 1) 
	  tot++;
      }
      sprintf (str, "%d broken domain(s)", tot);
      draw_string (str, size = 30, pos = 3);
      
      save ("dd.mp4", opt = "-r 2");
    }
  }
}
