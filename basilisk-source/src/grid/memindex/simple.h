/**
# Simple multi-dimensional arrays

This is the simplest implementation of a full multidimensional array,
without any optimisation for sparsity.

The interface is the same as that in [range.h]() which implements
"range sparse arrays".

The `Memindex` structure defines multi-dimensional arrays. The
interface is that used by [/src/grid/tree.h](). */

struct _Memindex {
  int len;
#if dimension == 1
  char ** b;
#elif dimension == 2  
  char *** b;
#else // dimension == 3
  char **** b;
#endif
};

#define Memindex struct _Memindex *

/**
The `mem_data()` macros return the data stored at a specific
(multidimensional) index. It assumes that the index is allocated. This
can be checked with `mem_allocated()`. */

#if dimension == 1
inline static
const bool mem_allocated (const Memindex m, const int i) {
  return (i >= 0 && i < m->len && m->b[i]);
}
#define mem_data(m,i) ((m)->b[i])
#elif dimension == 2
inline static
const bool mem_allocated (const Memindex m, const int i, const int j) {
  return (i >= 0 && i < m->len && m->b[i] &&
	  j >= 0 && j < m->len && m->b[i][j]);
}
#define mem_data(m,i,j) ((m)->b[i][j])
#else // dimension == 3
inline static
const bool mem_allocated (const Memindex m, const int i, const int j, const int k) {
  return (i >= 0 && i < m->len && m->b[i] &&
	  j >= 0 && j < m->len && m->b[i][j] &&
	  k >= 0 && k < m->len && m->b[i][j][k]);
}
#define mem_data(m,i,j,k) ((m)->b[i][j][k])
#endif // dimension == 3

/**
The `mem_new()` function returns a new (empty) `Memindex`. */

Memindex mem_new (int len)
{
  Memindex m = calloc (1, sizeof (struct _Memindex));
  return m;
}

/**
The `mem_destroy()` function frees all the memory allocated by a given
`Memindex`. */

void mem_destroy (Memindex m, int len)
{
#if dimension > 1
  for (int i = 0; i < len; i++)
    if (m->b[i]) {
  #if dimension > 2
      for (int j = 0; j < len; j++)
	if (m->b[i][j])
	  free (m->b[i][j]);
  #endif // dimension > 2
      free (m->b[i]);
    }
#endif // dimension > 1
  if (m->b)
    free (m->b);
  free (m);
}

/**
The `mem_assign()` function assigns a (pointer) value to a given index. */

void mem_assign (Memindex m, int i, int len, void * b)
{
  if (!m->b) {
    m->b = calloc (len, sizeof(char *));
    m->len = len;
  }
  mem_data(m,i) = b;
}
#elif dimension == 2
void mem_assign (Memindex m, int i, int j, int len, void * b)
{
  if (!m->b) {
    m->b = malloc (len*sizeof(char *));
    m->len = len;
  }
  if (!m->b[i])
    m->b[i] = calloc (len, sizeof(char *));
  mem_data(m,i,j) = b;
}
#else // dimension == 3
void mem_assign (Memindex m, int i, int j, int k, int len, void * b)
{
  if (!m->b) {
    m->b = malloc (len*sizeof(char *));
    m->len = len;
  }
  if (!m->b[i])
    m->b[i] = calloc (len, sizeof(char *));
  if (!m->b[i][j])
    m->b[i][j] = calloc (len, sizeof(char *));
  mem_data(m,i,j,k) = b;
}
#endif // dimension == 3

/**
The `mem_free()` function frees a given index. */

#if dimension == 1
void mem_free (Memindex m, int i, int len, void * b)
{
  mem_data(m,i) = NULL;
}
#elif dimension == 2
void mem_free (Memindex m, int i, int j, int len, void * b)
{
  mem_data(m,i,j) = NULL;
}
#else // dimension == 3
void mem_free (Memindex m, int i, int j, int k, int len, void * b)
{
  mem_data(m,i,j,k) = NULL;
}
#endif // dimension == 3

/**
The `foreach_mem()` macro traverses every `_i` allocated elements of
array `_m` taking into account a periodicity of `_len` (and ghost
cells). */

macro foreach_mem (Memindex index, int len, int _i) {
  Memindex _m = index;
  int _len = len;
  Point point = {0};
  for (point.i = max(Period.x*GHOSTS, 0);
       point.i < min(_len - Period.x*GHOSTS, _len);
       point.i += _i)
    if (_m->b[point.i])
#if dimension > 1
      for (point.j = max(Period.y*GHOSTS, 0);
	   point.j < min(_len - Period.y*GHOSTS, _len);
	   point.j += _i)
	if (_m->b[point.i][point.j])
#if dimension > 2
	  for (point.k = max(Period.z*GHOSTS, 0);
	       point.k < min(_len - Period.z*GHOSTS, _len);
	       point.k += _i)
	    if (_m->b[point.i][point.j][point.k])
#endif // dimension > 2
#endif // dimension > 1
	      {...}
}
