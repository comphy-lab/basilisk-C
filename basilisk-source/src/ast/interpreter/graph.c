typedef struct _Node Node;

typedef struct {
  Node * n;
  float c;
} Edge;

struct _Node {
  Key * k;
  Edge * parent, * child;
};

KHASH_MAP_INIT_INT(GRAPH, Node *)
  
typedef struct {
  khash_t(GRAPH) * hash;
} Graph;

static inline
bool identical_nodes (const Node * a, const Node * b)
{
  return a == b;
}

static
Node * get_node (Graph * g, Key * key)
{
  khiter_t k = kh_get (GRAPH, g->hash, key->j);
  if (k != kh_end (g->hash))
    return kh_value (g->hash, k);
  Node * n = calloc (1, sizeof (Node));
  n->k = key;
  int ret;
  k = kh_put (GRAPH, g->hash, key->j, &ret);
  assert (ret > 0);
  kh_value (g->hash, k) = n;
  return n;
}

static
int nedge (const Edge * list)
{
  int n = 0;
  if (list) for (const Edge * e = list; e->n; e++, n++);
  return n;
}

static
float remove_edge (Edge * list, Node * b)
{
  if (list)
    for (Edge * e = list; e->n; e++)
      if (identical_nodes (e->n, b)) {
	float c = e->c;
	for (; e->n; e++)
	  *e = *(e + 1);
	return c;
      }
  return 0.;
}

static
Edge * add_mono_edge (Edge * list, Node * b, float c)
{
  if (list) {
    for (Edge * e = list; e->n; e++)
      if (identical_nodes (e->n, b)) {
	e->c += c;
	if (!e->c)
	  for (; e->n; e++)
	    *e = *(e + 1);
	return list;
      }
  }
  int n = nedge (list);
  list = realloc (list, (n + 2)*sizeof(Edge));
  list[n + 1].n = NULL;
  list[n].n = b;
  list[n].c = c;
  return list;
}

static
void add_edge (Node * a, Node * b, float c)
{
  assert (fabs (c) < 100);
  a->child = add_mono_edge (a->child, b, c);
  b->parent = add_mono_edge (b->parent, a, c);
}

static Key * matrix_key (const System * s, int i, int j)
{
  Dimension ** r = s->r + i;
  if (!r[0]->c)
    return NULL;
  foreach_key (r[0], c)
    if (c->j == j)
      return c;
  return NULL;
}

static
int leftmost (const Dimension * d)
{
  int left = INT_MAX;
  for (Key ** c = d->c; c[0]; c++)
    if (c[0]->j < left)
      left = c[0]->j;
  return left;
}

Graph * system_to_graph (const System * s)
{
  Graph * g = malloc (sizeof (Graph));
  g->hash = kh_init (GRAPH);
  int n = 0;
  foreach_constraint (s, i) {
    if (i->c) {
      int left = leftmost (i);
      Key * kleft = matrix_key (s, n, left);
      float coef = - matrix (s, n, left);
      Node * node = get_node (g, kleft);
      foreach_key (i, c)
	if (c->j != left)
	  add_edge (get_node (g, c), node, matrix (s, n, c->j)/coef);
    }
    n++;
  }
  return g;
}

void check_node (const Node * n)
{
  if (n->child)
    for (const Edge * e = n->child; e->n; e++) {
      assert (e->n->parent);
      int np = 0;
      for (const Edge * f = e->n->parent; f->n; f++)
	if (identical_nodes (f->n, n))
	  np++;
      assert (np == 1);
    }
  if (n->parent)
    for (const Edge * e = n->parent; e->n; e++) {
      assert (e->n->child);
      int np = 0;
      for (const Edge * f = e->n->child; f->n; f++)
	if (identical_nodes (f->n, n))
	  np++;
      assert (np == 1);
    }
}

static
void remove_node_internal (Graph * g, Node * node)
{
  khiter_t k = kh_get (GRAPH, g->hash, node->k->j);
  kh_del (GRAPH, g->hash, k);
  free (node->child);
  free (node->parent);
  free (node);  
}

static
void edge_collapse_parent (Graph * g, Node * node, Node * parent)
{
  assert (remove_edge (parent->child, node));
  if (node->child)
    for (Edge * e = node->child; e->n; e++) {
      assert (remove_edge (e->n->parent, node));
      add_edge (parent, e->n, e->c*node->parent->c);
    }
  check_node (parent);
  remove_node_internal (g, node);
}

static
void edge_collapse_child (Graph * g, Node * node, Node * child)
{
  float c = remove_edge (child->parent, node);
  assert (c);
  if (node->parent)
    for (Edge * e = node->parent; e->n; e++) {
      assert (remove_edge (e->n->child, node));
      add_edge (e->n, child, e->c*c);
    }
  if (node->child)
    for (Edge * e = node->child; e->n; e++)
      if (!identical_nodes (e->n, child)) {
	assert (remove_edge (e->n->parent, node));
	add_edge (child, e->n, e->c*c);
      }
  check_node (child);
  remove_node_internal (g, node);
}

static
void remove_node (Graph * g, Node * node)
{
  if (node->child)
    for (Edge * e = node->child; e->n; e++)
      assert (remove_edge (e->n->parent, node));
  if (node->parent)
    for (Edge * e = node->parent; e->n; e++)
      assert (remove_edge (e->n->child, node));
  remove_node_internal (g, node);
}

int graph_simplify (Graph * g)
{
  int ns = 0;
  for (khiter_t k = kh_begin (g->hash); k != kh_end (g->hash); ++k)
    if (kh_exist (g->hash, k)) {
      Node * node = kh_value (g->hash, k);
      check_node (node);
      if (!atof (ast_right_terminal (node->k->parent)->start)) {
#if 1	
	if (nedge (node->child) == 0) {
	  remove_node (g, node);
	  ns++;
	}
#endif
#if 1
	else if (nedge (node->parent) == 1) {
	  edge_collapse_parent (g, node, node->parent->n);
	  ns++;
	}
#endif
      }
#if 1      
      else if (nedge (node->parent) == 1 &&
	       !atof (ast_right_terminal (node->parent->n->k->parent)->start)) {
	edge_collapse_child (g, node->parent->n, node);
	ns++;	  	
      }
#endif
    }
  return ns;
}

static
void print_node_label (const Key * k, FILE * fp)
{
  fprintf (fp, "  %d [label=\"", k->j);
  print_key_label (k, '\n', fp, LINENO);
#if 0  
  if (should_be_dimensionless (k))
    fputs ("\" style=filled fillcolor=\"aquamarine", fp);
#endif
  fprintf (fp, " %d\"]\n", k->j);
}

void graph_dot (const Graph * g, FILE * fp)
{
  fputs ("digraph mygraph {\n", fp);
  int key;
  Node * node;
  kh_foreach (g->hash, key, node, {
      if (node->child || node->parent) // do not print "orphaned" nodes
	print_node_label (node->k, fp);
      if (node->child)
	for (Edge * e = node->child; e->n; e++)
	  fprintf (fp, "  %d -> %d [label=\"%g\"]\n",
		   key, e->n->k->j, e->c);
#if 0  
      if (node->parent)
	for (Edge * e = node->parent; e->n; e++)
	  fprintf (fp, "  %d -> %d [label=\"p%g\"]\n",
		   e->n->k->index, key, e->c);
#endif
    });
  fputs ("}\n", fp);
}
