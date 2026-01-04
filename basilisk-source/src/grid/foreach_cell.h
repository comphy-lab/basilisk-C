/* ===============================================================
 *                    Tree traversal
 * recursive() below is for reference only. The macro
 * foreach_cell() is a stack-based implementation of
 * recursive(). It is about 12% slower than the recursive
 * version and 60% slower than simple array traversal.
 *
 * This article was useful:
 * http://www.codeproject.com/Articles/418776/How-to-replace-recursive-functions-using-stack-and
 *
 * =============================================================== */

#define _BOTTOM (2*point.j - GHOSTS)
#define _TOP    (_BOTTOM + 1)
#define _LEFT   (2*point.i - GHOSTS)
#define _RIGHT  (_LEFT + 1)
#define _BACK   (2*point.k - GHOSTS)
#define _FRONT  (_BACK + 1)

#if 0
void recursive (Point point)
{
  if (point.level == grid->depth) {
    /* do something */
  }
  else {
    Point p1 = point; p1.level = point.level + 1;
    p1.i = _LEFT;  p1.j = _TOP;    recursive (p1);
    p1.i = _RIGHT; p1.j = _TOP;    recursive (p1);
    p1.i = _LEFT;  p1.j = _BOTTOM; recursive (p1);
    p1.i = _RIGHT; p1.j = _BOTTOM; recursive (p1);
  }
}
#endif

#define STACKSIZE 20
#if dimension == 1
#define _push(b,c,d,e,f)			  \
  { _s++; stack[_s].l = b;			  \
    stack[_s].i = c;				  \
    stack[_s].stage = f; }
#define _pop()					  \
  { point.level = stack[_s].l;			  \
    point.i = stack[_s].i;			  \
    stage = stack[_s].stage; _s--; }
#elif dimension == 2
#define _push(b,c,d,e,f)			  \
  { _s++; stack[_s].l = b;			  \
    stack[_s].i = c; stack[_s].j = d;		  \
    stack[_s].stage = f; }
#define _pop()					  \
  { point.level = stack[_s].l;			  \
    point.i = stack[_s].i; point.j = stack[_s].j; \
    stage = stack[_s].stage; _s--; }
#else // dimension == 3
#define _push(b,c,d,e,f)						\
  { _s++; stack[_s].l = b;						\
    stack[_s].i = c; stack[_s].j = d; stack[_s].k = e;			\
    stack[_s].stage = f; }
#define _pop()								\
  { point.level = stack[_s].l;						\
    point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; \
    stage = stack[_s].stage; _s--; }
#endif

macro2 foreach_cell_root (Point root)
{
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};
#if dimension == 1
    struct { int l, i, stage; } stack[STACKSIZE];
#elif dimension == 2
    struct { int l, i, j, stage; } stack[STACKSIZE];
#else // dimension == 3
    int kg = 0; NOT_UNUSED(kg);
    struct { int l, i, j, k, stage; } stack[STACKSIZE];
#endif
    int _s = -1;
    _push (root.level, root.i, root.j, root.k, 0);
    while (_s >= 0) {
      int stage;
      _pop();
      if (!allocated (0,0,0))
	continue;
      switch (stage) {
      case 0: {

	{...}

	if (point.level < grid->depth) {
	  _push (point.level, point.i, point.j, point.k, 1);
	  _push (point.level + 1, _LEFT, _BOTTOM, _BACK, 0);
	}
	break;
      }
#if dimension == 1
      case 1: _push (point.level + 1, _RIGHT, _BOTTOM, _BACK, 0); break;
#elif dimension == 2
      case 1: _push (point.level, point.i, point.j, point.k, 2);
	_push (point.level + 1, _LEFT,  _TOP, _BACK, 0); break;
      case 2: _push (point.level, point.i, point.j, point.k, 3);
	_push (point.level + 1, _RIGHT, _BOTTOM, _BACK, 0); break;
      case 3: _push (point.level + 1, _RIGHT, _TOP, _BACK, 0); break;
#elif dimension == 3
      case 1: _push (point.level, point.i, point.j, point.k, 2);
	_push (point.level + 1, _LEFT, _BOTTOM, _FRONT, 0); break;
      case 2: _push (point.level, point.i, point.j, point.k, 3);
	_push (point.level + 1, _LEFT, _TOP, _BACK, 0); break;
      case 3: _push (point.level, point.i, point.j, point.k, 4);
	_push (point.level + 1, _LEFT, _TOP, _FRONT, 0); break;
      case 4: _push (point.level, point.i, point.j, point.k, 5);
	_push (point.level + 1, _RIGHT, _BOTTOM, _BACK, 0); break;
      case 5: _push (point.level, point.i, point.j, point.k, 6);
	_push (point.level + 1, _RIGHT, _BOTTOM, _FRONT, 0); break;
      case 6: _push (point.level, point.i, point.j, point.k, 7);
	_push (point.level + 1, _RIGHT, _TOP, _BACK, 0); break;
      case 7: _push (point.level + 1, _RIGHT, _TOP, _FRONT, 0); break;
#endif
      }
    }
  }
}

macro2 foreach_cell()
{
  {
#if dimension == 1
    Point root = {GHOSTS,0};
#elif dimension == 2
    Point root = {GHOSTS,GHOSTS,0};
#else // dimension == 3
    Point root = {GHOSTS,GHOSTS,GHOSTS,0};
#endif
    foreach_cell_root (root)
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
	  foreach_cell_root (root)
	    {...}
  }
}

macro2 foreach_cell_post_root (bool condition, Point root)
{
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point = {0};
#if dimension == 1
    struct { int l, i, stage; } stack[STACKSIZE];
#elif dimension == 2
    struct { int l, i, j, stage; } stack[STACKSIZE];
#else // dimension == 3
    int kg = 0; NOT_UNUSED(kg);
    struct { int l, i, j, k, stage; } stack[STACKSIZE];
#endif
    int _s = -1;
    _push (0, root.i, root.j, root.k, 0); /* the root cell */
    while (_s >= 0) {
      int stage;
      _pop();
      if (!allocated (0,0,0))
	continue;
      switch (stage) {
      case 0: {
	if (point.level == grid->depth) {
	  _push (point.level, point.i, point.j, point.k, 8);
	}
	else {
	  _push (point.level, point.i, point.j, point.k, 1);
	  if (condition)
	    _push (point.level + 1, _LEFT, _BOTTOM, _BACK, 0);
	}
	break;
      }
#if dimension == 1
      case 1:
	_push (point.level, point.i, point.j, point.k, 2);
	if (condition)
	  _push (point.level + 1, _RIGHT, _BOTTOM, _BACK, 0);
	break;
#elif dimension == 2
      case 1:
	_push (point.level, point.i, point.j, point.k, 2);
	if (condition)
	  _push (point.level + 1, _LEFT, _TOP, _BACK, 0);
	break;
      case 2:
	_push (point.level, point.i, point.j, point.k, 3);
	if (condition)
	  _push (point.level + 1, _RIGHT, _BOTTOM, _BACK, 0);
	break;
      case 3:
	_push (point.level, point.i, point.j, point.k, 4);
	if (condition)
	  _push (point.level + 1, _RIGHT, _TOP, _BACK, 0);
	break;
#elif dimension == 3
      case 1:
	_push (point.level, point.i, point.j, point.k, 2);
	if (condition)
	  _push (point.level + 1, _LEFT, _BOTTOM, _FRONT, 0);
	break;
      case 2:
	_push (point.level, point.i, point.j, point.k, 3);
	if (condition)
	  _push (point.level + 1, _LEFT, _TOP, _BACK, 0);
	break;
      case 3:
	_push (point.level, point.i, point.j, point.k, 4);
	if (condition)
	  _push (point.level + 1, _LEFT, _TOP, _FRONT, 0);
	break;
      case 4:
	_push (point.level, point.i, point.j, point.k, 5);
	if (condition)
	  _push (point.level + 1, _RIGHT, _BOTTOM, _BACK, 0);
	break;
      case 5:
	_push (point.level, point.i, point.j, point.k, 6);
	if (condition)
	  _push (point.level + 1, _RIGHT, _BOTTOM, _FRONT, 0);
	break;
      case 6:
	_push (point.level, point.i, point.j, point.k, 7);
	if (condition)
	  _push (point.level + 1, _RIGHT, _TOP, _BACK, 0);
	break;
      case 7:
	_push (point.level, point.i, point.j, point.k, 8);
	if (condition)
	  _push (point.level + 1, _RIGHT, _TOP, _FRONT, 0);
	break;	
#endif
      default:
	{...}

      }
    }
  }
}

macro2 foreach_cell_post (bool condition)
{
  {
#if dimension == 1
    Point root = {GHOSTS,0};
#elif dimension == 2
    Point root = {GHOSTS,GHOSTS,0};
#else // dimension == 3
    Point root = {GHOSTS,GHOSTS,GHOSTS,0};
#endif
    foreach_cell_post_root (condition, root)
      {...}
  }
}

macro2 foreach_cell_post_all (bool condition)
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

#if dimension == 1
macro2 foreach_cell_restore (ivec d = Dimensions, int rootlevel = 0)
{
  for (int ox = 0; ox < d.x; ox++) {
    Point root = {GHOSTS + ox, rootlevel};
    foreach_cell_root (root)
      {...}
  }
}
#elif dimension == 2
macro2 foreach_cell_restore (ivec d = Dimensions, int rootlevel = 0)
{
  for (int ox = 0; ox < d.x; ox++)
    for (int oy = 0; oy < d.y; oy++) {
      Point root = {GHOSTS + ox, GHOSTS + oy, rootlevel};
      foreach_cell_root (root)
	{...}
    }
}
#elif dimension == 3
macro2 foreach_cell_restore (ivec d = Dimensions, int rootlevel = 0)
{
  for (int ox = 0; ox < d.x; ox++)
    for (int oy = 0; oy < d.y; oy++)
      for (int oz = 0; oz < d.z; oz++) {
	Point root = {GHOSTS + ox, GHOSTS + oy, GHOSTS + oz, rootlevel};
	foreach_cell_root (root)
	  {...}
      }
}
#endif // dimension == 3
