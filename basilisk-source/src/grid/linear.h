/* Linear quadtree implementation based on 
 * 
 * K. Aizawa, K. Motomura, S. Kimura, R. Kadowaki and J. Fan
 * "Constant time neighbor finding in quadtrees"
 * ISCCSP 2008, Malta, 12-14 March 2008.
 * http://www.lcad.icmc.usp.br/~jbatista/procimg/quadtree_neighbours.pdf
 *
 * This uses a Z-grid ordering
 */

#define GRIDNAME "Linear quadtree"

#include <stdio.h>
#include <assert.h>

#define _I (quad_x(point.i))
#define _J (quad_y(point.i))

typedef struct {
  Data * d;
  int i, level, depth;
} Point;

#define data(k,l) point.d[quad_neighbor(point.i, k, l)]

/* the maximum level */
static int quad_r;

static int size (int l)
{
  return ((1 << (2*(l + 1))) - 1)/3;
}

int level (int p)
{
  p = 3*p + 1;
  int l = 0;
  while (p > 0) { p /= 4; l++; }
  return l - 1;
}

int code (int p, int l)
{
  return (p - size (l - 1)) << 2*(quad_r - l);
}

int index (int code, int l)
{
  return size(l - 1) + (code >> 2*(quad_r - l));
}

int quad_x (int p)
{
  int n = code (p, level(p)), a = 0, m = 1, b = 1;
  for (int i = 0; i < 2*quad_r - 1; i += 2, m *= 4, b *= 2)
    a += ((m & n) != 0)*b;
  return a;
}

int quad_y (int p)
{
  int n = code (p, level(p)), a = 0, m = 2, b = 1;
  for (int i = 1; i < 2*quad_r; i += 2, m *= 4, b *= 2)
    a += ((m & n) != 0)*b;
  return a;
}

int repeat (int a)
{
  int s = 0;
  for (int i = 0; i < quad_r; i++, a *= 4)
    s += a;
  return s;
}

static int quad_left, quad_right = 1, quad_top = 2, quad_bottom;

#define quad(n, d) ((((n|quad_bottom)+(d&quad_left))&quad_left)|(((n|quad_left)+(d&quad_bottom))&quad_bottom))

static int quad_id[3][3];

int quad_neighbor (int p, int i, int j)
{
  int d = quad_id[i+1][j+1];
  if (d == 0) return p;
  int l = level (p);
  int n = code (p, l);
  d <<= (2*(quad_r - l));
  return index (quad(n, d), l);
}

int quad_neighbor_finest (int p, int i, int j)
{
  int d = quad_id[i+1][j+1];
  if (d == 0) return p;
  int s = size (quad_r - 1);
  int n = p - s;
  return s + quad(n,d);
}

void * quadtree (int r, size_t s)
{
  quad_r = r;
  void * q = malloc (s*size (r));
  quad_left = repeat (1);
  quad_bottom = repeat (2);
  quad_id[0][2] = quad_top|quad_left;    
  quad_id[1][2] = quad_top;    
  quad_id[2][2] = quad_top|quad_right;
  quad_id[0][1] = quad_left;        quad_id[1][1] = 0;      quad_id[2][1] = quad_right;
  quad_id[0][0] = quad_bottom|quad_left; quad_id[1][0] = quad_bottom; 
  quad_id[2][0] = quad_bottom|quad_right;
  return q;
}

void * init_grid (int n)
{
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (stderr, "quadtree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  void * q = quadtree (depth, sizeof (Data));
  return q;
}

void free_grid (void * m)
{
  free (m);
}

#define STACKSIZE 20
#define _push(c,d)							\
  { _s++; stack[_s].i = c; stack[_s].stage = d; }
#define _pop(c,d)							\
  { c = stack[_s].i; d = stack[_s].stage; _s--; }

#define foreach(grid) \
{									  \
  struct { int i, stage; } stack[STACKSIZE]; int _s = -1; /* the stack */ \
  Point point;								\
  point.d = grid;							\
  point.depth = quad_r;							\
  _push (0, 0); /* the root cell */					\
  while (_s >= 0) {							\
    int stage;								\
    _pop (point.i, stage);						\
    if (!stage) {							\
      point.level = level (point.i);					\
      if (point.level < point.depth) {					\
        _push (point.i, 1);						\
	_push (4*point.i + 1, 0);					\
      }									\
      else {								\
        VARIABLES;							\
	/* do something */

#define end_foreach()				                        \
      }									\
    }								        \
    else {								\
      if (stage < 3)							\
	_push (point.i, stage + 1);					\
      _push (4*point.i + stage + 1, 0);					\
    }									\
  }									\
}

#if 0
int main (int argc, char * argv[])
{
  void * grid = init_grid (atoi (argv[1]));
  foreach (grid) {
    printf ("%d %d %d\n", I, J, point.level);
  } end_foreach();
  free_grid (grid);
}
#endif

#if 0
int main (int argc, char * argv[])
{
  int r = atoi(argv[1]);
  void * q = quadtree (r, sizeof (double));
  for (int p = size (r - 1); p < size (r); p++) {
    int q = neighbor (p, quad_bottom | quad_left, r);
    printf ("%d %d %d ( %d , %d ) | ( %d , %d ) ",
	    p, index (code (p, r), level (p), r), level (p),
	    x(p, r), y(p, r),
	    x(q, r), y(q, r));
    printf ("\n");
  }
  return 0;
}
#endif
