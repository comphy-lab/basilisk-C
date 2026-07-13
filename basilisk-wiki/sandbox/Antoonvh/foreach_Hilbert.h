/**
# Visit cells along a Hilbert-curve
*/
#define STACKSIZEH 20

#define _pushH(level, ig, jg, figure, next_child) {	  \
    _s++;						  \
    stack[_s].l = level;				  \
    stack[_s].i = ig; stack[_s].j = jg;			  \
    stack[_s].fig = figure;				  \
    stack[_s].nc = next_child;				  \
  }
#define _popH() {					  \
    point.level = stack[_s].l;				  \
    point.i = stack[_s].i; point.j = stack[_s].j;	  \
    figure = stack[_s].fig;				  \
    next_child = stack[_s].nc;				  \
    _s--;						  \
  }
/**
## Hilbert-indexing

We use binary codes to describe the x and y offset of the cildren for the four patterns.
 */
// 7 = 99;  // 0b01100011; // xxxxyyyy
// n = 54;  // 0b00110110;
// c = 156; // 0b10011100;
// u = 201; // 0b11001001;
unsigned char H_figs[4] = {99, 54, 156, 201}; // for fig index , 0, 1, 2, 3

#define IDX (2*point.i - GHOSTS + ((H_figs[figure] >> (7 - next_child)) & 1))
#define IDY (2*point.j - GHOSTS + ((H_figs[figure] >> (3 - next_child)) & 1))

char child_figs_n[4][4] =  {{1, 0, 0, 3}, // seq of figures of children in a fig 7 sequence
			    {0, 1, 1, 2}, // seq of figures of children in a fig n sequence
			    {3, 2, 2, 1}, // etc.
			    {2, 3, 3, 0}};  
char root_fig = 1; // 1 = n
/**
Visit cells on first encounter
 */
macro2 foreach_cell_root_Hilbert(Point root)   {
  Point point = {0};
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  struct {int i, j;
    unsigned char l, fig, nc;} stack[STACKSIZEH];
  int _s = -1;
  _pushH (0, root.i, root.j, root_fig, 0);
  while (_s >= 0) {
   unsigned char figure, next_child;
    _popH();
    if (!allocated (0))
      continue;
    if (next_child == 0)  //first visit
      {...}
    if (point.level < grid->depth) { // Maybe visit children...
      if (next_child < 3) // store state of parent
	_pushH (point.level, point.i, point.j, figure, next_child + 1);
      _pushH (point.level + 1, IDX, IDY, child_figs_n[figure][next_child], 0);
    }
  }
}


macro2 foreach_cell() {
  {
#if dimension == 1
    Point root = {GHOSTS,0};
#elif dimension == 2
    Point root = {GHOSTS,GHOSTS,0};
#else // dimension == 3
    Point root = {GHOSTS,GHOSTS,GHOSTS,0};
#endif
    foreach_cell_root_Hilbert (root)
      {...}
  }
}

macro2 foreach_cell_all() {
  {
  Point root = {0};
  for (root.i = GHOSTS*Period.x; root.i <= GHOSTS*(2 - Period.x); root.i++)
#if dimension >= 2
    for (root.j = GHOSTS*Period.y; root.j <= GHOSTS*(2 - Period.y); root.j++)
#endif
#if dimension >= 3
      for (root.k = GHOSTS*Period.z; root.k <= GHOSTS*(2 - Period.z); root.k++)
#endif
	foreach_cell_root_Hilbert (root)
	  {...}
  }
}
/**
Visit cells on last encounter
 */
struct {int i, j;
      unsigned char l, fig, nc;} stack[STACKSIZEH];
macro2 foreach_cell_post_root(bool condition,Point root)
{
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};
    //struct {int i, j;
    // unsigned char l, fig, nc;} stack[STACKSIZEH];
    int _s = -1;
    _pushH (0, root.i, root.j, root_fig, 0);
    while (_s >= 0) {
      unsigned char figure, next_child;
      _popH();
      if (next_child > 3) {
	{...}
	continue;
      }
      if (!allocated (0))  // We have moved too far down the tree
	continue; // Should we tell the stack?
      _pushH (point.level, point.i, point.j, figure, next_child + 1);
      if (condition && point.level < grid->depth)
	_pushH (point.level + 1, IDX, IDY, child_figs_n[figure][next_child], 0);
    }
  }
}

macro2 foreach_cell_post(bool condition)
  {
    {
#if dimension == 1
    Point root = {GHOSTS,0};
#elif dimension == 2
    Point root = {GHOSTS,GHOSTS,0};
#else // dimension == 3
    Point root = {GHOSTS,GHOSTS,GHOSTS,0};
#endif
    foreach_cell_post_root(condition, root)
      {...}
    }
  }

macro2 foreach_cell_post_all(bool condition)
{
  {
  Point root = {0};
  for (root.i = 0; root.i <= 2*GHOSTS; root.i++)
#if dimension >= 2
    for (root.j = 0; root.j <= 2*GHOSTS; root.j++)
#endif
#if dimension >= 3
      for (root.k = 0; root.k <= 2*GHOSTS; root.k++)
#endif
	foreach_cell_post_root (condition, root)
	  {...}
  }
}

macro2 foreach_leaf()
{
  foreach_cell()
    if (is_leaf (cell)) {
      if (is_active(cell) && is_local(cell)) 
	{...}
      continue;
    }
}

macro2 foreach_cell_restore (ivec d = Dimensions, int rootlevel = 0)
{
  if (!d.x) d.x = 1;
  if (!d.y) d.y = 1;
  for (int ox = 0; ox < d.x; ox++)
    for (int oy = 0; oy < d.y; oy++) {
      Point root = {GHOSTS + ox, GHOSTS + oy, rootlevel};
      foreach_cell_root_Hilbert (root)
	{...}
    }
}

