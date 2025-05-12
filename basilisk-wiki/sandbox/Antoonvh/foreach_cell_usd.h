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
/**
We switch the definitions for the _TOP and _BOTTOM children:

Also every function and definition gets a new name:
 */
#include "grid/foreach_cell.h"

#define _TOP2    (2*point.j-GHOSTS)
#define _BOTTOM2 (_TOP2 + 1)


@def foreach_cell_root2(root)
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point;
#if dimension == 1
    struct { int l, i, stage; } stack[STACKSIZE];
#elif dimension == 2
    struct { int l, i, j, stage; } stack[STACKSIZE];
#else // dimension == 3
    int kg = 0; NOT_UNUSED(kg);
    struct { int l, i, j, k, stage; } stack[STACKSIZE];
#endif
    int _s = -1;
    _push (0, root.i, root.j, root.k, 0);
    while (_s >= 0) {
      int stage;
      _pop();
      if (!allocated (0))
	continue;
      switch (stage) {
      case 0: {
	POINT_VARIABLES;
	/* do something */
@
@def end_foreach_cell_root2()
        if (point.level < grid->depth) {
	  _push (point.level, point.i, point.j, point.k, 1);
          _push (point.level + 1, _LEFT, _BOTTOM2, _BACK, 0);
        }
        break;
      }
#if dimension == 1
      case 1: _push (point.level + 1, _RIGHT, _BOTTOM2, _BACK, 0); break;
#elif dimension == 2
      case 1: _push (point.level, point.i, point.j, point.k, 2);
	      _push (point.level + 1, _LEFT,  _TOP2, _BACK, 0); break;
      case 2: _push (point.level, point.i, point.j, point.k, 3);
	      _push (point.level + 1, _RIGHT, _BOTTOM2, _BACK, 0); break;
      case 3: _push (point.level + 1, _RIGHT, _TOP2, _BACK, 0); break;
#elif dimension == 3
      case 1: _push (point.level, point.i, point.j, point.k, 2);
	      _push (point.level + 1, _LEFT, _BOTTOM2, _FRONT, 0); break;
      case 2: _push (point.level, point.i, point.j, point.k, 3);
	      _push (point.level + 1, _LEFT, _TOP2, _BACK, 0); break;
      case 3: _push (point.level, point.i, point.j, point.k, 4);
	      _push (point.level + 1, _LEFT, _TOP2, _FRONT, 0); break;
      case 4: _push (point.level, point.i, point.j, point.k, 5);
	      _push (point.level + 1, _RIGHT, _BOTTOM2, _BACK, 0); break;
      case 5: _push (point.level, point.i, point.j, point.k, 6);
	      _push (point.level + 1, _RIGHT, _BOTTOM2, _FRONT, 0); break;
      case 6: _push (point.level, point.i, point.j, point.k, 7);
	      _push (point.level + 1, _RIGHT, _TOP2, _BACK, 0); break;
      case 7: _push (point.level + 1, _RIGHT, _TOP2, _FRONT, 0); break;
#endif
      }
    }
  }
@

@def foreach_cell2() {
#if dimension == 1
  Point root = {GHOSTS,0};
#elif dimension == 2
  Point root = {GHOSTS,GHOSTS,0};
#else // dimension == 3
  Point root = {GHOSTS,GHOSTS,GHOSTS,0};
#endif
  foreach_cell_root2 (root)
@
@define end_foreach_cell2() end_foreach_cell_root2() }

@def foreach_cell_all2() {
  Point root = { .level = 0 };
  for (root.i = GHOSTS*Period.x; root.i <= GHOSTS*(2 - Period.x); root.i++)
#if dimension >= 2
    for (root.j = GHOSTS*Period.y; root.j <= GHOSTS*(2 - Period.y); root.j++)
#endif
#if dimension >= 3
      for (root.k = GHOSTS*Period.z; root.k <= GHOSTS*(2 - Period.z); root.k++)
#endif
	foreach_cell_root2 (root)
@
@define end_foreach_cell_all2() end_foreach_cell_root2() }

@def foreach_cell_post_root2(condition, root)
  {
    int ig = 0, jg = 0;	NOT_UNUSED(ig); NOT_UNUSED(jg);
    Point point;
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
      if (!allocated (0))
	continue;
      switch (stage) {
      case 0: {
        POINT_VARIABLES;
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
      default: {
        POINT_VARIABLES;
	/* do something */
@
@def end_foreach_cell_post_root2()
      }
      }
    }
  }
@

@def foreach_cell_post2(condition)
  {
#if dimension == 1
    Point root = {GHOSTS,0};
#elif dimension == 2
    Point root = {GHOSTS,GHOSTS,0};
#else // dimension == 3
    Point root = {GHOSTS,GHOSTS,GHOSTS,0};
#endif
    foreach_cell_post_root2(condition, root)
@
@define end_foreach_cell_post2() end_foreach_cell_post_root2() }

@def foreach_cell_post_all2(condition) {
  Point root = { .level = 0 };
  for (root.i = 0; root.i <= 2*GHOSTS; root.i++)
#if dimension >= 2
    for (root.j = 0; root.j <= 2*GHOSTS; root.j++)
#endif
#if dimension >= 3
      for (root.k = 0; root.k <= 2*GHOSTS; root.k++)
#endif
	foreach_cell_post_root2 (condition, root)
@
@define end_foreach_cell_post_all2() end_foreach_cell_post_root2() }

@def foreach_leaf2() foreach_cell2()
  if (is_leaf (cell)) {
    if (is_active(cell) && is_local(cell)) {
@
@define end_foreach_leaf2() } continue; } end_foreach_cell2()
