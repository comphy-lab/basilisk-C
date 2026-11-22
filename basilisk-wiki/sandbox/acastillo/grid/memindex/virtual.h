/**
# Multi-dimensional arrays using virtual memory

This implementation uses the Unix/Linux virtual memory subsystem to
more efficiently store sparse multidimensional arrays.

The `Memindex` structure defines multi-dimensional arrays. The
interface is that used by [/src/grid/tree.h](). */

#include <sys/mman.h>

struct _Memindex {
  ivec len;
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
  return (i >= 0 && i < m->len.x && m->b[i]);
}
#define mem_data(m,i) ((m)->b[i])
#elif dimension == 2
inline static
const bool mem_allocated (const Memindex m, const int i, const int j) {
  return (i >= 0 && i < m->len.x && m->b[i] &&
	  j >= 0 && j < m->len.y && m->b[i][j]);
}
#define mem_data(m,i,j) ((m)->b[i][j])
#else // dimension == 3
inline static
const bool mem_allocated (const Memindex m, const int i, const int j, const int k) {
  return (i >= 0 && i < m->len.x && m->b[i] &&
	  j >= 0 && j < m->len.y && m->b[i][j] &&
	  k >= 0 && k < m->len.z && m->b[i][j][k]);
}
#define mem_data(m,i,j,k) ((m)->b[i][j][k])
#endif // dimension == 3

/**
The `mem_new()` function returns a new (empty) `Memindex`. */

Memindex mem_new (ivec len)
{
  Memindex m = calloc (1, sizeof (struct _Memindex));
  m->len = len;
  return m;
}

/**
The `mem_destroy()` function frees all the memory allocated by a given
`Memindex`. */

void mem_destroy (Memindex m, ivec len)
{
  if (m->b) {
#if dimension > 1
    for (int i = 0; i < len.x; i++)
      if (m->b[i]) {
  #if dimension > 2
        for (int j = 0; j < len.y; j++)
	  if (m->b[i][j])
	    munmap (m->b[i][j], m->len.z*sizeof (char *));
  #endif // dimension > 2
        munmap (m->b[i], m->len.y*sizeof (char *));
      }
#endif // dimension > 1
    munmap (m->b, m->len.x*sizeof (char *));
  }
  free (m);
}

/**
The `mem_assign()` function assigns a (pointer) value to a given index. */

#define mymap(len) mmap (NULL, len*sizeof (char *), PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0)

#if dimension == 1
void mem_assign (Memindex m, int i, ivec len, void * b)
{
  if (!m->b) {
    m->b = mymap (len.x);
    m->len = len;
  }
  mem_data(m,i) = b;
}
#elif dimension == 2
void mem_assign (Memindex m, int i, int j, ivec len, void * b)
{
  assert (b != NULL);
  if (!m->b) {
    m->b = mymap (len.x);
    m->len = len;
  }
  if (!m->b[i])
    m->b[i] = mymap (len.y);
  mem_data(m,i,j) = b;
}
#else // dimension == 3
void mem_assign (Memindex m, int i, int j, int k, ivec len, void * b)
{
  if (!m->b) {
    m->b = mymap (len.x);
    m->len = len;
  }
  if (!m->b[i])
    m->b[i] = mymap (len.y);
  if (!m->b[i][j])
    m->b[i][j] = mymap (len.z);
  mem_data(m,i,j,k) = b;
}
#endif // dimension == 3

/**
The `mem_free()` function frees a given index. */

static bool mem_free1d (char ** const b, const int i, int len)
{
 // Assumes 4k pages = 512*sizeof(char *) = 512*8 = 4096
  char ** const start = b + (i/512)*512;
  char ** const end = b + len < start + 512 ? b + len : start + 512;
  for (char ** j = start; j < end; j++)
    if (*j)
      return false;
  if (madvise (start, 4096, MADV_DONTNEED))
    perror ("virtual.h: madvise: ");
  return true;
}

#if dimension == 1
void mem_free (Memindex m, int i, ivec len, void * b)
{
  if (m->b[i]) {
    m->b[i] = NULL;
    mem_free1d (m->b, i, len.x);
  }
}
#elif dimension == 2
void mem_free (Memindex m, int i, int j, ivec len, void * b)
{
  if (m->b[i][j]) {
    m->b[i][j] = NULL;
    if (mem_free1d (m->b[i], j, len.y))
      mem_free1d ((char **)m->b, i, len.x);
  }
}
#else // dimension == 3
void mem_free (Memindex m, int i, int j, int k, ivec len, void * b)
{
  if (m->b[i][j][k]) {
    m->b[i][j][k] = NULL;
    if (mem_free1d (m->b[i][j], k, len.z))
      if (mem_free1d ((char **)m->b[i], j, len.y))
	mem_free1d ((char **)m->b, i, len.x);
  }
}
#endif // dimension == 3

/**
The `foreach_mem()` macro traverses every `_i` allocated elements of
array `_m` taking into account a periodicity of `_len` (and ghost
cells). */

macro foreach_mem (Memindex index, ivec len, int _i) {
  Memindex _m = index;
  ivec _len = len;
  Point point = {0};
  if (_m->b)
    for (point.i = max(Period.x*GHOSTS, 0);
         point.i < min(_len.x - Period.x*GHOSTS, _len.x);
         point.i += _i)
      if (_m->b[point.i])
#if dimension > 1
        for (point.j = max(Period.y*GHOSTS, 0);
	     point.j < min(_len.y - Period.y*GHOSTS, _len.y);
	     point.j += _i)
	  if (_m->b[point.i][point.j])
#if dimension > 2
	    for (point.k = max(Period.z*GHOSTS, 0);
	         point.k < min(_len.z - Period.z*GHOSTS, _len.z);
	         point.k += _i)
	      if (_m->b[point.i][point.j][point.k])
#endif // dimension > 2
#endif // dimension > 1
	        {...}
}
