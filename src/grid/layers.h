/**
# Multiple "layers"

This adds support for multiple "layers" i.e. a constant number `nl` of
contiguous scalar fields. */

int nl = 1;

/**
We redefine two types of block traversals. An "outer" traversal,
usable outside foreach() loops, which uses a global `_layer` index... */

int _layer = 0;

#undef foreach_block
@define foreach_block() for (_layer = 0; _layer < nl; _layer++)
@define end_foreach_block() _layer = 0;
  
/**
... and an "inner" traversal, usable within foreach() loops, which
uses the `point.l` index as local layer index. */
  
@define foreach_block_inner() for (point.l = 0; point.l < nl; point.l++)
@define end_foreach_block_inner() point.l = 0;
  
/**
We also redefine the "per field" (inner) traversal. */

#undef foreach_blockf
@def foreach_blockf(_f)
  for (point.l = 0; point.l < _attribute[_f.i].block; point.l++)
@
@define end_foreach_blockf() point.l = 0;
  
/**
The two indices are combined to access field values. In practice only
one of the indices is used. */
  
@undef _index
@def _index(a,m)
  (a.i + (_layer + point.l + m < _attribute[a.i].block ?
	  _layer + point.l + m : 0))
@

@undef val
@define val(a,k,p,m) data(k,p,m)[_index(a,m)]

/**
foreach_layer() is just an alias for foreach_block(). */
  
#define foreach_layer() foreach_block()
