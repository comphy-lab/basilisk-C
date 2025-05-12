/**
# Dimensional analysis 

This file contains the code implementing dimensional analysis in
Basilisk.

We will assume in what follows that the reader is familiar with the
basic concepts of dimensional analysis and of the conventions used to
represent dimensions within Basilisk. These are introduced in the
[Dimension tutorial](/Tutorial.dimensions).

## Introduction

The primary goal of the following code is to check that all
(arithmetic) operations done by a program written in Basilisk C are
dimensionally consistent.

Interestingly, while dimensional consistency is clearly a cornerstone
of physics, there does not seem to be a widely used language which
ensure this consistency at the programming level. The code below is an
attempt to fill this gap, at least in the context of code written with
Basilisk.

A reasonably simple way to impose this consistency would be to
associate different types to arithmetic variables, depending on their
dimensions. For example instead of writing

~~~c
double L = 10, T = 4, H = 30, G = 9.81, U = L/T + sqrt(G*H);
~~~

one would write

~~~c
Length L = 10, H = 30;
Time T = 4;
Acceleration G = 9.81;
Speed U = L/T + sqrt(G*H);
~~~

the compiler would then be able to check that operations (such as the
addition above) only happen between compatible types (i.e. Speed in
the example).

An important drawback of this approach (which has been used in the
past) is that it obviously requires major code rewriting, as well as
important constraints on the way code is written. In particular, all
dimensions must be known and correctly imposed beforehand: writing
code which is "dimension-independent" is not possible anymore.

It also seems obvious that associating dimensions with variables (via
their types or otherwise) is not the most general approach. Dimensions
(and units) should be associated with numerical values (i.e. with the
content of variables). One could for example design CPU chips which
would carry out operations on floating point numbers with associated
dimensions. This chip could then return errors (i.e. a "Floating point
exception") in the case of operations done on values with incompatible
dimensions. Building a new special-purpose chip is clearly not a very
practical solution however. Furthermore this approach would still
require specifying the dimensions of every numerical value used as
input to the program, which could also be quite cumbersome.

The approach taken here, called "Dimensional Inference", also relies
on associating dimensions with values (and not variables or types) and
has been known since at least the work of [Wand and O'Keefe,
1991](#wand1991). The principle is simple:

1) collect all "dimensional constraints" i.e. relationships necessary
to impose consistency when performing a given arithmetic operation. By
definition, these constraints will always involve the dimensions of
the input constants.

2) this collection of constraints linearly relates the dimensions of
all the numerical constants used in the program i.e. they form a
linear system of equations.

3) try to solve this system. There are three main cases: a) if the
system does not have a solution, then there is a dimensional
inconsistency somewhere. b) if the system has more than one solution
it requires additional constraints (i.e. directly setting the
dimensions of more input constants). c) if the system has a solution
then the dimensions of all input constants can be uniquely determined.

A major advantage of this approach is that the dimensions of most
constants can usually be "inferred" from the explicit specifications
of only a few input constants. Minimal code modification is thus
necessary to ensure full dimensional consistency.

## References

~~~bib
@inproceedings{wand1991,
  title={Automatic Dimensional Inference.},
  author={Wand, Mitchell and O'Keefe, Patrick},
  booktitle={Computational Logic-Essays in Honor of Alan Robinson},
  pages={479--483},
  year={1991},
  organization={Citeseer}
}

@article{sandberg2003,
  title={Automatic dimensional consistency checking for 
         simulation specifications},
  author={Sandberg, Mikael and Persson, Daniel and Lisper, Bj{\"o}rn},
  journal={Group},
  volume={1},
  number={2},
  pages={M1},
  year={2003}
}

@inproceedings{farrimond2007,
  title={Dimensional Inference Using Symbol Lives},
  author={Farrimond, Brian and Collins, John},
  booktitle={SETP},
  pages={152--159},
  year={2007}
}
~~~

## Implementation

The primary reason for associating dimensions with types is that this
allows for "static analysis" of dimensional consistency i.e. the
compiler can check the consistency directly. When dimensions are
associated with values, the only way to ensure consistency is to "run"
the code ("dynamic analysis"), at least partially, and check that the
operations performed involve values which are dimensionally
consistent.

To do so, the implementation below uses the
[interpreter](interpreter.c) to "run" the code being analysed. Of
course running the full code within an interpreter just to check for
consistency would be too expensive, so a simplified version of the
code is run. The major simplifications are: the mesh contains a single
grid point and the number of times loops are iterated is limited. This
includes of course the number of timesteps which are performed. This
"run" is performed at compilation time and usually completes in a
reasonable amount of time, compared to the total compilation time.

During this runtime phase, all arithmetic operations involving types
to which dimensions can be associated (essentially floating point
numbers, see below for details) are intercepted and the corresponding
dimensional constraints are added to the linear system. Some obvious
dimensional errors can also be detected directly within this first
pass of analysis.

A potential problem with the simplifications above is that they may
miss some branches which will be taken when the full code is run. For
example, an upwind scheme may take a branch which will depend on the
sign of the value of a field. If a single grid point is used, this
field will have a single value and only one of the branches will be
taken. If this happens, dimensional analysis will not be
exhaustive. To avoid this, the interpreter implements a sophisticated
scheme where values can be specified as "unset" i.e. having an
undefined value. If a condition depends on such an "unset" value, the
interpreter makes sure to take both branches and to properly propagate
"unset" values through their chains of dependencies.

Note that this approach imposes some constraints on the dimensions of
values defined within conditional branches.

More details on implementation are given in the documentation
associated with each section of the code below.

## Dimension tests

Test cases for dimensional analysis are in
`/src/ast/interpreter/dimension-tests/test*.c` they are run when doing

~~~bash
cd src/ast/interpreter/
make check
~~~

They also include some documentation of various "corner cases".

# Data structures */

#if 0
# define DEBUG(...) (__VA_ARGS__)
#else
# define DEBUG(...)
#endif

typedef struct _List List;
typedef struct _Key Key;
typedef struct _Dimension Dimension;
typedef struct _System System;

struct _List {
  Dimension * d;
  List * next;
};

struct _Key {
  Ast * parent;
  char * field, * label;
  int j, refs, used;
  List * dimensions;
};

static
void key_add_dimension (Key * k, Dimension * d, Stack * stack)
{
  k->refs++;
  if (!k->dimensions) {
    k->dimensions = static_calloc (1, sizeof (List), stack);
    k->dimensions->d = d;
    return;
  }
  List * l = k->dimensions;
  while (l->next) l = l->next;
  l->next = static_calloc (1, sizeof (List), stack);
  l->next->d = d;
}

static
bool key_remove_dimension (Key * k, const Dimension * d)
{
  for (List * l = k->dimensions, * prev = NULL; l; prev = l, l = l->next)
    if (l->d == d) {
      if (!prev)
	k->dimensions = l->next;
      else
	prev->next = l->next;
      k->refs--;
      return true;
    }
  return false;
}

/**
The dimension is
$$
\mathbf{a} + \sum b_i [c_i]
$$
*/

struct _Dimension {
  double * a, * b;
  Key ** c;
  const Ast * origin;
  int row;
  char * warn;
  System * s;
};

/**
The system of dimensional constraints. */

struct _System {
  Dimension ** r;
  Allocator * alloc;
  Key ** index;
  int m, n;
  FILE * output;
  int redundant, finite, lineno, warn;
  bool dimensionless;
};

/**
A special value for [*]. */

static
Dimension * const dimension_any = (Dimension *) 128;

/**
# Outputs */

const char * ast_file_crop (const char * file)
{
  int len = strlen (BASILISK);
  if (strlen (file) > len && !strncmp (file, BASILISK, len))
    return file + len - 4;
  return file;
}

static inline
bool is_identifier (char c)
{
  return c == '_' || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

/**
The function below just replaces 'val(s,0,1,0)' with 's[0,1]' and
'_attribute[s.i].' with 's.'. */

static
char * simplified_expression (const char * string)
{
  char * smp = malloc (strlen (string) + 1), * fp = smp;
  const char * s = string;
  while (*s != '\0') {
    if (!strncmp (s, "val(", 4) && (s == string || !is_identifier (*(s - 1)))) {
      s += 4;
      while (*s != ',') *fp++ = *s++;
      s++;
      *fp++ = '[';
      int para = 1;
      const char * end = s;
      while (para) {
	if (*end == '(') para++;
	if (*end == ')') para--;
	end++;
      }
      const char * start = end - 2;
      while (*start == '0' || *start == ',') start--;
      while (s <= start) *fp++ = *s++;
      *fp++ = ']';
      s = end - 1;
    }
    else if (!strncmp (s, "_attribute[", 11) && (s == string || !is_identifier (*(s - 1)))) {
      s += 11;
      while (strncmp (s, ".i]", 3)) *fp++ = *s++;
      s += 2;
    }
    else
      *fp++ = *s;
    s++;
  }
  *fp = *s;
  return smp;
}

static
void print_simplified_expression (const char * string, FILE * fp)
{
  char * s = simplified_expression (string);
  fputs (s, fp);
  free (s);
}

#define EXPR_MAGIC "10293847566574839201."

enum {
  LINENO   = 1,
  CARRIAGE = 1 << 1,
  INDEX    = 1 << 2,
  NORIGIN  = 1 << 3,
  REFS     = 1 << 4
};

static
Ast * get_explicit_dimension (const Ast * n)
{
  Ast * array = ast_ancestor (ast_parent (n, sym_primary_expression), 2);
  if (ast_schema (array, sym_array_access,
		  0, sym_postfix_expression,
		  0, sym_primary_expression) &&
      array->child[3])
    return array;
  return NULL;
}

static
Ast * parent_expression (const Key * c)
{
  Ast * n = c->parent;
  Ast * parent = ast_parent (n, sym_assignment_expression);
  if (ast_evaluate_constant_expression (parent) < DBL_MAX) {
    if (ast_schema (ast_ancestor (parent, 2), sym_init_declarator,
		    2, sym_initializer))
      n = ast_ancestor (parent, 2);
    else if (ast_schema (parent->parent, sym_assignment_expression,
			 1, sym_assignment_operator,
			 0, token_symbol ('=')))
      n = parent->parent;
  }
  Ast * array = get_explicit_dimension (n);
  if (array)
    return array;
  return n;
}

static
char * key_label (const Key * c, char sep, int flags)
{
  Ast * n = parent_expression (c);
  AstTerminal * t = ast_left_terminal (n);
  if (t->start && !strcmp (t->start, EXPR_MAGIC))
    n = ast_parent (n, sym_additive_expression)->child[0];
  const char * file = ast_file_crop (t->file);
  int len = strlen (file) + (c->field ? strlen(c->field): 0) + 30;
  char * label = malloc (len);
  snprintf (label, len, "%s:%d:%c'%s%s",
	    ast_file_crop (t->file), (flags & LINENO) ? t->line : 0, sep,
	    c->field && c->field[0] != '\0' ? c->field : "",
	    c->field && c->field[0] != '\0' ? "[] = " : "");
  char * s = ast_str_append (n, NULL);
  if (s) {
    char * s1 = simplified_expression (ast_crop_before (s));
    str_append (label, s1, "\'");
    free (s), free (s1);
  }
  else
    str_append (label, "0'");
  return label;
}

static
void print_key_label (const Key * c, char sep, FILE * fp, int flags)
{
  char * label = key_label (c, ' ', true);
  fputs (label, fp);
  free (label);
  if (flags & INDEX)
    fprintf (fp, " %d", c->j);
  if (flags & REFS)
    fprintf (fp, " %d", c->refs);
#if 0
  fprintf (fp, " %d", c->flags);
#endif
#if 0
  fprintf (fp, " %d", c->refs);
#endif
}

#define foreach_key(a, b) if (a->c) for (Key ** _p = a->c, * b = *_p; b; b = *++_p)

#define END 1e10

static
void dimension_print (const Dimension * d, FILE * fp, double coef, int flags)
{
  if (d == dimension_any) {
    fputs ("[*]", fp);
    return;
  }
  if (!d || (!d->a && !d->b)) {
    fputs ("[0]", fp);
    return;
  }
  if (d->a) {
    fputc ('[', fp);
    for (double * i = d->a; *i < END; i++) {
      if (*i)
	fprintf (fp, "%g", (*i)*coef);
      else
	fputs ("0", fp);
      if (*(i + 1) < END)
	fputc (',', fp);
    }
    fputc (']', fp);
  }
  if (d->b) {
    double * i = d->b;
    foreach_key (d, c) {
      //      assert ((*i) != 0.);
      fprintf (fp, " %s ", (*i)*coef > 0. ? "+" : "-");
      if (fabs ((*i)*coef) != 1.)
	fprintf (fp, "%g*[", fabs ((*i)*coef));
      else
	fputc ('[', fp);
      print_key_label (c, ' ', fp, flags);
      fputc (']', fp);
      i++;
    }
  }
}

/**
# Various utility functions */

static inline
bool can_have_dimension (const Value * v) {
  return v && !v->pointer && (v->type->sym == sym_LONG ||
			      v->type->sym == sym_DOUBLE ||
			      v->type->sym == sym_FLOAT);
}

static inline
bool can_assign_dimension (const Value * v) {
  return v && !v->pointer && (v->type->sym == sym_DOUBLE ||
			      v->type->sym == sym_FLOAT);
}

static
int dimension_type_size (const Ast * type)
{
  int size = ast_base_type_size (type);
  switch (type->sym) {
  case sym_LONG: case sym_FLOAT: case sym_DOUBLE:
    // case sym_COMPLEX: // fixme: need to be added below and tested
    size += sizeof (Ast *);
    break;
  }
  return size;
}

static inline
int unknowns (const Dimension * d)
{
  int n = 0;
  if (d && d != dimension_any)
    foreach_key (d, c)
      n++;
  return n;
}

static
Ast * binary_expression_parent (Ast * n)
{
  n = n->parent;
  while (n->child &&
	 (!n->child[1] ||
	  n->sym == sym_primary_expression ||
	  n->child[0]->sym == sym_unary_operator))
    n = n->parent;
  for (int * op = (int[]){
      sym_multiplicative_expression,
      sym_additive_expression,
      sym_shift_expression,
      sym_relational_expression,
      sym_equality_expression,
      -1 }; *op > 0; op++) {
    if (ast_schema (n, *op, 0, *op))
      return n;
    if (ast_schema (n->parent, *op, 0, *op))
      return n->parent;
  }
  if (ast_schema (n, sym_assignment_expression,
		  1, sym_assignment_operator))
    return n;
  return NULL;
}

static
int unique_unknowns (const Dimension * a, const Dimension *b)
{
  if (unknowns (b) > unknowns (a)) {
    const Dimension * tmp = a; a = b; b = tmp;
  }
  int nt = 0;
  if (a && a != dimension_any && a->c) {
    nt = unknowns (a);
    if (b && b != dimension_any && b->c)
      foreach_key (b, c) {
	Key ** found;
	for (found = a->c; found[0] && found[0] != c; found++);
	if (!found[0])
	  nt++;
      }
  }
  return nt;
}

static inline
int dimensions (const Dimension * d)
{
  int n = 0;
  if (d && d != dimension_any && d->a)
    for (double * a = d->a; *a != END; a++, n++);
  return n;
}

#define value_dimension(v) (*((Dimension **)(((char *) (v)->data.p) + ast_base_type_size ((v)->type))))

static void constraint_print (Dimension * d, FILE * fp, int flags);

static
Dimension * dimension_zero (Allocator * alloc, Ast * origin)
{
  Dimension * d = allocate (alloc, sizeof (Dimension));
  d->a = NULL;
  d->b = NULL;
  d->c = NULL;
  d->origin = origin;
  d->warn = NULL;
  return d;
}

#define ROUNDOFF 1e-12

static
Dimension * dimension_simplify (Dimension * d)
{
  int n = dimensions (d);
  double end = END;
  for (int i = n - 1; i >= 0; i--)
    if (fabs (d->a[i]) <= ROUNDOFF)
      d->a[i] = end;
    else
      end = 0.;
  if (d->a && d->a[0] == END)
    d->a = NULL;

  n = unknowns (d);
  for (int i = 0; i < n; i++)
    if (fabs (d->b[i]) <= ROUNDOFF) {
      for (int j = i; j < n - 1; j++)
	d->b[j] = d->b[j + 1];
      for (int j = i; j < n; j++)
	d->c[j] = d->c[j + 1];
      n--, i--;
    }
  if (n == 0)
    d->b = NULL, d->c = NULL;
  
  return d;
}

static
void add_constraint (System * s, Dimension * constraint, Stack * stack);

static
Dimension * new_dimension (Ast * n, Stack * stack)
{
  Allocator * alloc = stack_static_alloc (stack);
  Dimension * dimension = allocate (alloc, sizeof (Dimension));
  dimension->a = NULL;
  dimension->b = allocate (alloc, sizeof (double));
  dimension->b[0] = 1.;
  dimension->c = allocate (alloc, 2*sizeof (Key *));
  dimension->c[0] = allocate (alloc, sizeof (Key));
  Ast * prev = n;
  for (Ast ** c = n->parent->child; *c; c++)
    if (c[0]->sym == n->sym) {
      assert (n == prev);
      n = c[0];
      break;
    }
  *dimension->c[0] = (Key) { .parent = n };
  dimension->c[1] = NULL;
  dimension->origin = n;
  dimension->warn = NULL;
  return dimension;
}

static
bool in_c_file (const Ast * n)
{
  AstTerminal * t = ast_left_terminal (n);
  if (!t->file)
    return false;
  int len = strlen (t->file);
  return !strcmp (t->file + len - 2, ".c");
}

static
Dimension * get_dimension (Value * v, Stack * stack)
{
  if (!can_have_dimension (v))
    return NULL;
  
  if (value_dimension (v))
    return value_dimension (v);

  if (is_constant_expression (v)) {
    Ast * n = (Ast *) v;
    
    /**
    Character constants have no dimension. */
    
    if (n->sym == sym_I_CONSTANT && ast_terminal (n)->start[0] == '\'')
      return (value_dimension (v) = NULL);

    /**
    Read an explicitly-set dimension. */
    
    Ast * array = get_explicit_dimension (n);
    if (array) {
      Ast * expr = array->child[2];
      if (expr->sym != sym_expression)
	return (value_dimension (v) = dimension_any);
      Dimension * dimension = new_dimension (array, stack);
      Allocator * alloc = stack_static_alloc (stack);
      Dimension * rhs = allocate (alloc, sizeof (Dimension));
      rhs->b = NULL;
      int m = 0;
      foreach_item (expr, 2, d)
	m++;
      rhs->a = allocate (alloc, (m + 1)*sizeof (double));
      m = 0;
      foreach_item_r (expr, sym_assignment_expression, d) {
	Value * v = run (d, stack);
	bool error;
	rhs->a[m++] = value_double (v, &error, stack);
      }
      rhs->a[m] = END;
      rhs = dimension_simplify (rhs);
      if (rhs->a) {
	dimension->a = rhs->a;
	dimension->b[0] = - 1.;
      }
      add_constraint (interpreter_get_data (stack), dimension, stack);
      dimension->c[0]->used = 0;
      rhs->a = NULL;
      rhs->b = allocate (alloc, sizeof (double));
      rhs->b[0] = 1.;
      rhs->c = dimension->c;
      return (value_dimension (v) = rhs);
    }

    /**
    Non-zero constants in multiplicative expressions are dimensionless. */

    Ast * parent = binary_expression_parent (n);
    if (parent && (parent->sym == sym_multiplicative_expression ||
		   (parent->sym == sym_assignment_expression &&
		    (parent->child[1]->child[0]->sym == sym_MUL_ASSIGN ||
		     parent->child[1]->child[0]->sym == sym_DIV_ASSIGN))) &&
	(!ast_terminal(n) || atof (ast_terminal(n)->start)))
      return (value_dimension (v) = dimension_zero (stack_static_alloc (stack), (Ast *) v));
    
    /**
    In general, return a new unknown dimension for this constant. */

    return (value_dimension (v) = new_dimension (n, stack));
  }
  
  /**
  Array members without a dimension have the same dimension as array
  member zero. */
  
  if (ast_schema ((Ast *) v, sym_array_access)) {
    Ast * n = (Ast *) v;
    assert (n->child[2]);
    Value * a = run (n->child[0], stack);
    assert (a);
    Value * value = array_member_value (n, a, 0, false, stack);
    return (value_dimension (v) = value_dimension (value));
  }

  /**
  Everything else is dimensionless. */

  Dimension * d = dimension_zero (stack_static_alloc (stack), (Ast *) v);
  if ((value_flags (v) & unset) && v->type->sym != sym_LONG &&
      ((System *)interpreter_get_data (stack))->output) {
    AstTerminal * t = ast_left_terminal ((Ast *) v);
    char * s = ast_str_append ((Ast *) v, NULL), warn[200];
    s = simplified_expression (ast_crop_before (s));
    Ast * n = (Ast *) v, * call = ast_function_call_identifier ((Ast *) v);
    if (call && !strcmp (ast_terminal (call)->start, "val")) {
      call = NN(n, sym_function_call,
		NN(n, sym_postfix_expression,
		   NN(n, sym_primary_expression,
		      NA(n, sym_IDENTIFIER, "_field_name"))),
		NCA(n, "("),
		NN(n, sym_argument_expression_list,
		   ast_copy (ast_find ((Ast *) v, sym_argument_expression_list_item))),
		NCA(n, ")"));
      Value * vname = run (call, stack);
      ast_destroy (call);
      char * s1 = NULL;
      str_append (s1, value_data (vname, void *), strchr (s, '['), NULL);
      free (s);
      s = s1;
    }
    snprintf (warn, 199,
	      "%s:%d: warning: '%s' is unset: assuming it has dimension [0]\n",
	      ast_file_crop (t->file), t->line, s);
    free (s);
    int len = strlen (warn) + 1;
    d->warn = allocate (stack_static_alloc (stack), len*sizeof (char));
    strcpy (d->warn, warn);
  }
  return (value_dimension (v) = d);
}

static
void value_print_dimension (const Value * v, bool zero, int flags, FILE * fp)
{
  if (can_have_dimension (v)) {
    Dimension * dimension = value_dimension (v);
    if (dimension) {
      fputc (' ', fp);
      dimension_print (dimension, fp, 1., flags);
    }
    else if (zero)
      fputs (" [0]", fp);
  }
}

static
void value_print_dimension_hook (const Value * v, FILE * fp)
{
  value_print_dimension (v, false, LINENO, fp);
}

static
Dimension * dimensions_multiply (const Ast * origin, Allocator * alloc,
				 const Dimension * da, const Dimension * db,
				 double v)
{

  /**
  Multiplicative expressions between quantities with any dimension has
  any dimension. */
  
  if (da == dimension_any || db == dimension_any)
    return dimension_any;
  
  Dimension * dimension = allocate (alloc, sizeof (Dimension));
  dimension->origin = origin;
  
  int na = dimensions (da), nb = dimensions (db), n = na > nb ? na : nb;
  dimension->a = n ? allocate (alloc, sizeof (double)*(n + 1)) : NULL;
  if (n > 0)
    dimension->a[n] = END;
  for (int i = 0; i < n; i++)
    dimension->a[i] = 0.;
  for (int i = 0; i < na; i++)
    dimension->a[i] += da->a[i];
  for (int i = 0; i < nb; i++)
    dimension->a[i] += v*db->a[i];

  n = unique_unknowns (da, db);
  if (!n) {
    dimension->b = NULL;
    dimension->c = NULL;
    return dimension_simplify (dimension);
  }
  
  dimension->b = allocate (alloc, sizeof (double)*n);
  dimension->c = allocate (alloc, sizeof (Key *)*(n + 1));
  dimension->c[n] = NULL;

  double va = 1., vb = v;
  na = unknowns (da), nb = unknowns (db);
  if (nb > na) {
    const Dimension * d = da; da = db; db = d;
    va = v, vb = 1.;
    int n = na; na = nb; nb = n;
  }
  n = 0;
  for (int i = 0; i < na; i++)
    dimension->b[n] = va*da->b[i], dimension->c[n++] = da->c[i];
  for (int i = 0; i < nb; i++) {
    int j;
    for (j = 0; j < n && dimension->c[j] != db->c[i]; j++);
    if (j < n)
      dimension->b[j] += vb*db->b[i];
    else
      dimension->b[n] = vb*db->b[i], dimension->c[n++] = db->c[i];
  }
  return dimension_simplify (dimension);
}

static
void * not_homogeneous (Ast * n, Value * va, Value * vb, Stack * stack)
{
  AstTerminal * t = ast_left_terminal (n);
  System * s = interpreter_get_data (stack);
  const char * msg = s->warn ? "warning" : "error";
  fprintf (stderr, "%s:%d: %s: '", t->file, t->line, msg);
  print_simplified_expression (ast_crop_before (ast_str_append (n, NULL)), stderr);
  fputs ("' is not an homogeneous expression\n", stderr);
  t = ast_left_terminal ((Ast *)va);
  fprintf (stderr, "%s:%d: %s: '", t->file, t->line, msg);
  print_simplified_expression (ast_crop_before (ast_str_append ((Ast *)va, NULL)), stderr);
  fputs ("' has dimensions", stderr);
  value_print_dimension (va, true, LINENO, stderr);
  fputs ("\n", stderr);
  t = ast_left_terminal ((Ast *)vb);
  fprintf (stderr, "%s:%d: %s: '", t->file, t->line, msg);
  print_simplified_expression (ast_crop_before (ast_str_append ((Ast *)vb, NULL)), stderr);
  fputs ("' has dimensions", stderr);
  value_print_dimension (vb, true, LINENO, stderr);
  fputs ("\n", stderr);
  return NULL;
}

/**
# System of dimensional constraints */

static
System * system_new()
{
  System * system = calloc (1, sizeof (System));
  system->r = calloc (1, sizeof (Dimension *));
  system->dimensionless = true;
  system->lineno = true;
  return system;
}

static
void system_destroy (System * system)
{
  free (system->r);
  free (system->index);
  free (system);
}

static void constraint_print (Dimension * d, FILE * fp, int flags)
{
  if (d != dimension_any) {
    if (d->origin && !(flags & NORIGIN)) {
      AstTerminal * t = ast_left_terminal (d->origin);
      assert (t->file);
      char * s = ast_str_append (d->origin, NULL);
      fprintf (fp, "%s:%d: '", t->file ? ast_file_crop (t->file) : "null", t->line);
      print_simplified_expression (ast_crop_before (s), fp);
      fprintf (fp, "'%s", (flags & CARRIAGE) ? "\n\t└─ " : ": ");
      free (s);
    }
    if (unknowns (d) == 1) {
      fputc ('[', fp);
      print_key_label (*d->c, ' ', fp, flags | LINENO);
      fputs ("] = ", fp);
      double * b = d->b;
      assert (b[0]);
      double coef = - 1./b[0];;
      d->b = NULL;
      dimension_print (d, fp, coef, flags);
      d->b = b;
      fputc ('\n', fp);
    }
    else {
      dimension_print (d, fp, 1., flags);
      fputs (" = [0]\n", fp);
    }
  }
}

#define foreach_constraint(s, i) if (s->r) for (Dimension ** _p = s->r, * i = *_p; i; i = *++_p)

void system_print (System * s, FILE * fp)
{
  foreach_constraint (s, r)
    constraint_print (r, fp, LINENO | INDEX | REFS);
}

static
bool system_append (System * s, Dimension * constraint)
{
  int n = 0;
  foreach_constraint(s, r) {
    if (r == constraint)
      return false;
    n++;
  }
  s->r = realloc (s->r, (n + 2)*sizeof (Dimension *));
  s->r[n] = constraint;
  constraint->s = s;
  s->r[n + 1] = NULL;
  if (constraint->a)
    s->dimensionless = false;
  return true;
}

static
void add_constraint (System * s, Dimension * constraint, Stack * stack)
{
  if (!constraint || (!constraint->a && !constraint->b))
    return;
  if (stack_verbosity (stack) > 2)
    constraint_print (constraint, stderr, LINENO);
  if (system_append (s, constraint)) {
    foreach_key (constraint, c) {
      key_add_dimension (c, constraint, stack);
      c->used = 1;
    }
  }
}

static Dimension * row_replacement (Dimension * d, Dimension * with, Allocator * alloc)
{
  if (!d->c)
    return NULL;
  Key * p = with->c[0], ** c;
  double * j = d->b;
  for (c = d->c; c[0] && c[0] != p; c++, j++);
  if (!c[0])
    return NULL; // p not in d
  
  assert (with->b[0]);
  return dimensions_multiply (d->origin, alloc, d, with, - j[0]/with->b[0]);
}

static inline
double column (const Dimension * r, int j)
{
  if (!r->c)
    return 0.;
  double * b = r->b;
  foreach_key (r, c) {
    if (c->j == j)
      return *b;
    b++;
  }
  return 0.;
}

static inline
double matrix (const System * s, int i, int j)
{
  return column (s->r[i], j);
}

static
int compare_unknowns (const void * pa, const void * pb)
{
  Dimension * const * a = pa, * const * b = pb;
  int na = unknowns (a[0]), nb = unknowns (b[0]);
  if (na > nb)
    return 1;
  if (na < nb)
    return -1;
  return pa > pb ? 1 : -1;
}

static
void system_index (System * s)
{
  s->m = 0;
  foreach_constraint (s, r) {
    s->m++;
    foreach_key (r, c)
      c->j = -1;
  }
    
  /**
  Sort equations by increasing order of number of unknowns so that
  pivoting always uses the simplest pivot. This changes the pivoted
  system when it is row-deficient (i.e. unconstrained) but should not
  change the result otherwise (to within round-off errors). */
  
  qsort (s->r, s->m, sizeof (Dimension *), compare_unknowns);
  
  s->m = 0, s->n = 0;
  int row = 0;
  foreach_constraint (s, r) {
    r->row = row++;
    if (r->c || r->a) {
      s->m++;
      foreach_key (r, c)
	if (c->j == -1) {
	  c->j = s->n++;
	  s->index = realloc (s->index, s->n*sizeof (Key *));
	  s->index[s->n - 1] = c;
	}
    }
  }
}

#define WRITE_GRAPH 0

#if WRITE_GRAPH
# include "graph.c"

void system_write_graph (System * s, const char * name)
{
  FILE * fp = fopen (name, "w");
  if (!fp) {
    perror (name);
    return;
  }
  system_index (s);
  Graph * g = system_to_graph (s);
  // while (graph_simplify (g));
  //  graph_simplify (g);
  graph_dot (g, fp);
  fclose (fp);
  /**
  gnuplot> plot '/tmp/out' matrix with image 
  */
}

void system_write_matrix (System * s, const char * name)
{
  FILE * fp = fopen (name, "w");
  if (!fp) {
    perror (name);
    return;
  }
  int i = 0;
  foreach_constraint (s, r) {
    double * b = r->b;
    foreach_key (r, c) {
      fprintf (fp, "%d %d %g ", c->j, i, *b);
      print_key_label (c, ' ', fp, LINENO | INDEX);
      fputc ('\n', fp);
      b++;
    }
    i++;
  }
  fclose (fp);
}
#endif

static
void set_system (Dimension * j, System * sub)
{
  if (system_append (sub, j)) {
    foreach_key (j, c)
      for (List * l = c->dimensions; l; l = l->next)
	if (!l->d->s)
	  set_system (l->d, sub);
  }
}

static
bool remove_dimensionless_subsystems (System * s)
{
  bool found = false;
  foreach_constraint (s, j)
    j->s = NULL;
  foreach_constraint (s, j)
    if (!j->s && j->c && j->c[1]) {
      System * su = system_new();
      su->alloc = s->alloc;
      set_system (j, su);
      if (su->dimensionless) {
	found = true;
	foreach_constraint (su, j)
	  j->a = j->b = NULL, j->c = NULL;
      }
      system_destroy (su);
    }
  return found;
}

static inline
void set_row (System * s, int row, Dimension * d)
{
  s->r[row] = d; d->row = row;
}

static bool system_pivot (System * s, Stack * stack)
{
  system_index (s);

  if (!s->m)
    return true;
  
  /** 
  Gaussian elimination, algorithm directly adapted from
  [wikipedia](https://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode). */ 

  int h = 0; /* Initialization of the pivot row */
  int k = 0; /* Initialization of the pivot column */
  while (h < s->m && k < s->n) {
    
    /**
    Find the k-th pivot: we just take the first non-zero value. This
    relies on the equations already being sorted in order of
    increasing complexity (see system_index() above). */
    
    int i_piv = 0;
    double piv = 0.;
    for (int i = h; i < s->m && !piv; i++) {
      double f = matrix (s, i, k);
      if (f)
	piv = f, i_piv = i;
    }
    if (!piv) {
      /* No pivot in this column */
      DEBUG (fprintf (stderr, "no pivot %d: ", k),
	     print_key_label (s->index[k], ' ', stderr, LINENO | INDEX),
	     fputc ('\n', stderr));
      k++; // pass to next column
    }
    else {
      /* swap rows(h, i_min) */
      Dimension * r = s->r[h];
      set_row (s, h, s->r[i_piv]);
      set_row (s, i_piv, r);
      /* Do for all rows below pivot: */
      int nd = s->index[k]->refs;
      if (nd > 1) {
	Dimension * list[nd], ** i = list;
	for (List * _l = s->index[k]->dimensions; _l; _l = _l->next)
	  *i++ = _l->d;	
	for (int j = 0; j < nd; j++)
	  if (list[j]->row > h) {
	    Dimension * r = list[j];
	    double f = column (r, k)/piv;
	    /* Substract the scaled pivot row from the current row */
	    DEBUG (fprintf (stderr, "pivoting %d ", h), constraint_print (s->r[h], stderr, LINENO | INDEX | NORIGIN),
		   fprintf (stderr, "         in "), constraint_print (r, stderr, LINENO | INDEX | NORIGIN));
	    foreach_key(r, c)
	      assert (key_remove_dimension (c, r));
	    Dimension * d = dimensions_multiply (r->origin, s->alloc, r, s->r[h], - f);
	    foreach_key(d, c)
	      key_add_dimension (c, d, stack);
	    DEBUG (fprintf (stderr, "      gives "), constraint_print (d, stderr, LINENO | INDEX | NORIGIN));
	    /* If l.h.s. is zero and r.h.s. is not zero the system does not have a solution */
	    if (!d->c && d->a) {
	      AstTerminal * t = ast_left_terminal (s->r[h]->origin);
	      fprintf (stderr, "%s:%d: %s: the dimensional constraints below are not compatible\n",
		       t->file, t->line, s->warn ? "warning" : "error");
	      constraint_print (s->r[h], stderr, CARRIAGE | LINENO);
	      constraint_print (r, stderr, CARRIAGE | LINENO);
	      return false;
	    }
	    set_row (s, r->row, d);
	  }
      }
      /* Increase pivot row and column */
      h++, k++;
    }
  }

  /**
  Eliminate non-independent constraints */

  s->m = 0;
  for (Dimension ** r = s->r; *r; r++) {
    // constraint_print (r[0], stderr, LINENO | INDEX);
    if (r[0]->c)
      s->m++;
    else { // zero unknowns
      assert (!r[0]->a); // r.h.s. is zero
      r[0] = NULL;
      r++;
      for (; r[0]; r++)
	assert (!r[0]->c && !r[0]->a);
      break;
    }
  }

  return true;
}

static
Key ** system_unconstrained (const System * s)
{
  if (s->n == s->m)
    return NULL;
  Key ** unconstrained = NULL;
  int nu = 0;
  Dimension ** i;
  for (i = s->r; *i; i++);
  for (i--; i >= s->r; i--)
    if (unknowns (i[0]) > 1)
      foreach_key (i[0], c) {
	Key ** found;
	for (found = unconstrained; found && found[0] && found[0]->parent != c->parent; found++);
	if (!found || !found[0]) {
	  nu++;
	  unconstrained = realloc (unconstrained, (nu + 1)*sizeof (Key *));
	  unconstrained[nu - 1] = c;
	  unconstrained[nu] = NULL;
	}
      }
  return unconstrained;
}

static bool system_solve (System * s, Stack * stack)
{
#if WRITE_GRAPH  
  system_write_graph (s, "dimensions0.dot");
  { FILE * fp = fopen ("system0", "w"); system_print (s, fp), fclose (fp); }
#endif
  
  if (!system_pivot (s, stack))
    return false;
  
  DEBUG (fprintf (stderr, "@@ %d unknowns, %d constraints\n", s->n, s->m));
  if (s->m) {

#if WRITE_GRAPH
    { FILE * fp = fopen ("system1", "w"); system_print (s, fp), fclose (fp); }
#endif
    
    /**
    Backward substitution */

    for (Dimension ** i = s->r + s->m - 1; i >= s->r; i--)
      if (unknowns (i[0]) == 1 && i[0]->c[0]->refs > 1)
	for (Dimension ** j = s->r; j < i; j++) {
	  Dimension * r = row_replacement (*j, *i, s->alloc);
	  if (r) {
	    assert (r->c || !r->a);
	    DEBUG (fprintf (stderr, "replacing "), constraint_print (*i, stderr, LINENO | INDEX | NORIGIN),
		   fprintf (stderr, "in        "), constraint_print (*j, stderr, LINENO | INDEX | NORIGIN),
		   fprintf (stderr, "gives     "), constraint_print (r, stderr, LINENO | INDEX | NORIGIN));
	    *j = r;
	  }
	}
    
#if WRITE_GRAPH    
    system_write_graph (s, "dimensions2.dot");
    { FILE * fp = fopen ("system2", "w"); system_print (s, fp), fclose (fp); }
#endif

    if (remove_dimensionless_subsystems (s))
      system_index (s);

#if WRITE_GRAPH
    system_write_graph (s, "dimensions3.dot");
    { FILE * fp = fopen ("system3", "w"); system_print (s, fp), fclose (fp); }
#endif
    
  }

  if (s->output && s->lineno)
    fprintf (s->output, "%d constraints, %d unknowns\n", s->m, s->n);

  return true;
}

/**
# Dimension hooks to the interpreter */

static
Value * dimension_run (Ast * n, Stack * stack)
{
  if (!n || n == ast_placeholder)
    return NULL;

  switch (n->sym) {

  case sym_primary_expression: {
    Value * v = ast_run_node (n, stack);
    if (!v)
      return NULL;
    if (n->child[0]->sym == sym_constant &&
	(n->child[0]->child[0]->sym == sym_I_CONSTANT ||
	 n->child[0]->child[0]->sym == sym_F_CONSTANT))
      get_dimension (v, stack);
    else if (n->child[0]->sym == sym_IDENTIFIER &&
	     !strcmp (ast_terminal (n->child[0])->start, "_val_higher_dimension"))
      value_dimension (v) = dimension_any;
    return v;
  }

  /**
  All initializers of `double, float *` or `double, float []` arrays
  have the same dimensions. */

  case sym_init_declarator: {
    Ast * list, * type;
    if (!(list = ast_find (n, sym_initializer_list)) ||
	!(type = ast_find (ast_parent (n, sym_declaration), sym_declaration_specifiers,
			   0, sym_type_specifier,
			   0, sym_types)) ||
	(type->child[0]->sym != sym_DOUBLE && type->child[0]->sym != sym_FLOAT))
      return ast_run_node (n, stack);
    Ast * array = NULL;
    if (ast_schema (n, sym_init_declarator,
		    0, sym_declarator,
		    0, sym_direct_declarator,
		    1, token_symbol('['))) { // a `[]` array
      if (!(array = ast_schema (n, sym_init_declarator,
				0, sym_declarator,
				0, sym_direct_declarator,
				0, sym_direct_declarator,
				0, sym_generic_identifier,
				0, sym_IDENTIFIER)))
	return ast_run_node (n, stack);
    }
    else if (ast_schema (n, sym_init_declarator,
			 0, sym_declarator,
			 0, sym_pointer)) { // a `*` array
      if (!(array = ast_schema (n, sym_init_declarator,
				0, sym_declarator,
				1, sym_direct_declarator,
				0, sym_generic_identifier,
				0, sym_IDENTIFIER)))
	return ast_run_node (n, stack);
    }
    else
      return ast_run_node (n, stack);
    Value * value = ast_run_node (n, stack);
    int len = 0;
    foreach_item (list, 2, item) len++;
    if (len > 1) {
      char slen[20];
      snprintf (slen, 19, "%d", len);
      char * func = strdup ("_set_element_dimensions");
      if (type->child[0]->sym == sym_FLOAT)
	str_append (func, "_float");
      Ast * call =
	NN(n, sym_function_call,
	   NN(n, sym_postfix_expression,
	      NN(n, sym_primary_expression,
		 NA(n, sym_IDENTIFIER, func))),
	   NCA(n, "("),
	   NN(n, sym_argument_expression_list,
	      NN(n, sym_argument_expression_list,
		 NN(n, sym_argument_expression_list_item,
		    ast_attach (ast_new_unary_expression (n),
				ast_new_identifier (n, ast_terminal(array)->start)))),
	      NCA(n, ","),
	      NN(n, sym_argument_expression_list_item,
		 ast_new_constant (n, sym_I_CONSTANT, slen))),
	   NCA(n, ")"));
      free (func);
      // fprintf (stderr, "%d ", len); ast_print_tree (array, stderr, 0, 0, -1);
      ast_run_node (call, stack);
      ast_destroy (call);
    }
    return value;
  }
    
  default:
    return ast_run_node (n, stack);
    
  }
  return NULL;
}

static
Value * dimension_assign (Ast * n, Value * dst, Value * src, Stack * stack)
{
  if (can_assign_dimension (dst)) {
    Dimension * d = get_dimension (src, stack);
    if (!d)
      
      /**
      If the dimension is undefined we are (most probably) casting
      from values which cannot have dimensions i.e. 'int's etc. so we
      set the dimension to zero. */
      
      d = dimension_zero (stack_static_alloc (stack), (Ast *) src);

    value_dimension (dst) = d;
  }
  return dst;
}

static void warn_unset_zero (Dimension * d, Stack * stack)
{
  if (d && d->warn) {
    System * system = interpreter_get_data (stack);
    if (system->output)
      fputs (d->warn, system->output);
    d->warn = NULL;
  }
}

static
Dimension * homogeneous_dimensions (Ast * n, Value * va, Value * vb, Stack * stack)
{
  Dimension * da = get_dimension (va, stack), * db = get_dimension (vb, stack);
  if (da == dimension_any)
    return db;
  if (db == dimension_any)
    return da;

  warn_unset_zero (da, stack);
  warn_unset_zero (db, stack);
  
  /**
  We first consider the case where all dimensions are known and check
  that they are identical. */

  int nt = unique_unknowns (da, db);
  if (!nt) {
    if (dimensions (da) != dimensions (db))
      return not_homogeneous (n, va, vb, stack);
    else if (da && da->a)
      for (double * i = da->a, * j = db->a; *i < END; i++, j++)
	if (*i != *j)
	  return not_homogeneous (n, va, vb, stack);
    return da ? da : db;
  }

  /**
  One or both dimensions are unknown, we add the corresponding
  constraint to the system (if it is compatible). */
  
  Dimension * constraint = dimensions_multiply (n, stack_static_alloc (stack), da, db, -1.);
  if (!constraint->c && constraint->a)
    return not_homogeneous (n, va, vb, stack);
  add_constraint (interpreter_get_data (stack), constraint, stack);
  
  /**
  We return either the non-constant dimension or the dimension with
  the smallest number of unknowns. */
  
  Dimension * d;
  if (is_constant_expression (va))
    d = db;
  else if (is_constant_expression (vb))
    d = da;  
  else if (unknowns (da) < unknowns (db))
    d = da;
  else
    d = db;
  
  /**
  If d is NULL, we explicitly set the dimension to zero. */
  
  if (!d)
    d = dimension_zero (stack_static_alloc (stack), n);

  return d;
}

static inline
bool is_finite (const Value * v)
{
  if (!is_constant_expression (v))
    return true;
  double val =fabs (v->type->sym == sym_LONG ? value_data (v, long) :
		    v->type->sym == sym_FLOAT ? value_data (v, float) :
		    v->type->sym == sym_DOUBLE ? value_data (v, double) : 0.);
  if (val > 1e-30 && val < 1e30 && val != 1234567890)
    return true;
  return false;
}

static void dimension_error (Stack * stack)
{
  ((StackData *)stack_get_data (stack))->maxcalls = -1;
}

static
Value * dimension_binary_operation (Ast * n, Stack * stack, Value * a, Value * b, Value * value)
{
  if ((!can_assign_dimension (a) && !can_assign_dimension (b)) ||
      (is_constant_expression (a) && is_constant_expression (b) &&
       !get_explicit_dimension ((Ast *)a) && !get_explicit_dimension ((Ast *)b)))
    return value;

  int op = n->child[1]->sym;
  int mul = ((op == token_symbol('*') ||
	      ast_schema (n->child[1], sym_assignment_operator, 0, sym_MUL_ASSIGN)) ? 1 :
	     (op == token_symbol('/') ||
	      ast_schema (n->child[1], sym_assignment_operator, 0, sym_DIV_ASSIGN)) ? -1 :
	     0);
  if (op == token_symbol('+') || op == token_symbol('-') ||
      ast_schema (n->child[1], sym_assignment_operator, 0, sym_ADD_ASSIGN) ||
      ast_schema (n->child[1], sym_assignment_operator, 0, sym_SUB_ASSIGN)) {
    if (can_assign_dimension (value) &&
	!(value_dimension (value) = homogeneous_dimensions (n, a, b, stack)))
      dimension_error (stack);
  }
  else if (mul) {
    if (can_assign_dimension (value))
      value_dimension (value) = dimensions_multiply (n, stack_static_alloc (stack),
						     get_dimension (a, stack), get_dimension (b, stack), mul);
  }
  else if (can_have_dimension (a) && can_have_dimension (b) &&
	   !homogeneous_dimensions (n, a, b, stack))
    dimension_error (stack);
  return value;
}

static
Dimension * dimension_pow (Ast * origin, Stack * stack, Dimension * d, double e)
{
  if (d == dimension_any)
    return d;
  Dimension * r = dimension_zero (stack_static_alloc (stack), origin);
  if (e) {
    Allocator * alloc = stack_static_alloc (stack); 
    if (d->a) {
      int na = 0;
      for (double * a = d->a; *a != END; a++, na++);
      r->a = allocate (alloc, sizeof (double)*(na + 1));
      memcpy (r->a, d->a, sizeof (double)*(na + 1));
      for (double * a = r->a; *a != END; a++)
	*a *= e;
    }
    if (d->c) {
      int nc = 0;
      foreach_key (d, c) nc++;
      r->c = allocate (alloc, sizeof (Key *)*(nc + 1));
      memcpy (r->c, d->c, sizeof (Key *)*(nc + 1));
      r->b = allocate (alloc, sizeof (double)*(nc + 1));
      memcpy (r->b, d->b, sizeof (double)*nc);
      double * b = r->b;
      foreach_key (r, c)
	*b++ *= e;
    }
  }
  return r;
}

static
void print_params (Ast * call, Value ** params, const char * msg)
{
  Ast * list = call->child[2];
  int i = 0;
  foreach_item_r (list, sym_argument_expression_list_item, arg) {
    AstTerminal * t = ast_left_terminal (arg);
    fprintf (stderr, "%s:%d: %s: '%s' has dimensions",
	     t->file, t->line, msg, ast_crop_before (ast_str_append (arg, NULL)));
    value_print_dimension (params[i++], true, LINENO, stderr);
    fputc ('\n', stderr);
  }
}

static
Value * dimension_internal_functions (Ast * call, Ast * identifier, Value ** params, Stack * stack, Value * value)
{
  const char * name = ast_terminal (identifier)->start;
  if (!strcmp (name, "sqrt"))
    value_dimension (value) = dimension_pow (call, stack, get_dimension (params[0], stack), 0.5);
  else if (!strcmp (name, "atan2")) {
    Dimension * constraint = dimensions_multiply (call, stack_static_alloc (stack),
						  get_dimension (params[0], stack),
						  get_dimension (params[1], stack), -1.);
    if (constraint && constraint->a && !constraint->b) {
      AstTerminal * t = ast_left_terminal (call);
      System * s = interpreter_get_data (stack);
      const char * msg = s->warn ? "warning" : "error";
      fprintf (stderr, "%s:%d: %s: the arguments of 'atan2' must have the same dimensions\n",
	       t->file, t->line, msg);
      print_params (call, params, msg);
      dimension_error (stack);
      return value;
    }
    add_constraint (interpreter_get_data (stack), constraint, stack);
    value_dimension (value) = dimension_zero (stack_static_alloc (stack), call);
  }
  else if (!strcmp (name, "pow")) {
    bool error;
    value_dimension (value) = dimension_pow (call, stack, get_dimension (params[0], stack),
					     value_double (params[1], &error, stack));
  }
  else if (!strcmp (name, "fabs"))
    value_dimension (value) = get_dimension (params[0], stack);
  else if (!strcmp (name, "reset_field_value")) {
    char * field = value_data (params[0], char *);
    Dimension * d = *((Dimension **)(field + ast_base_type_size (params[2]->type)));
    if (value_data (params[3], int)) // block != 0
      d->c = NULL; // block scalars are not used, set their dimension to zero
    else { // block == 0
      const char * name = value_data (params[1], char *);
      if (name) {
	d->c[0]->field = allocate (stack_static_alloc (stack), (strlen(name) + 1)*sizeof (char));
	strcpy (d->c[0]->field, name);
      }
      else
	d->c[0]->field = NULL;
    }
  }
  else {
    static char * funcs[] = {
      "exp", "log", "log10",
      "sin", "cos", "tan",
      "asin", "acos", "atan",
      "sinh", "cosh", "tanh",
      "asinh", "acosh", "atanh",
      NULL
    };
    for (char ** i = funcs; *i; i++)
      if (!strcmp (name, *i)) {
	Dimension * constraint = get_dimension (params[0], stack);
	if (constraint && constraint->a && !constraint->b) {
	  AstTerminal * t = ast_left_terminal (call);
	  System * s = interpreter_get_data (stack);
	  const char * msg = s->warn ? "warning" : "error";
	  fprintf (stderr, "%s:%d: %s: the argument of '%s' must be dimensionless\n",
		   t->file, t->line, msg, ast_crop_before (ast_str_append (call, NULL)));
	  print_params (call, params, msg);
	  dimension_error (stack);
	  return value;
	}
	constraint->origin = call;
	add_constraint (interpreter_get_data (stack), constraint, stack);
	value_dimension (value) = dimension_zero (stack_static_alloc (stack), call);
	return value;
      }
  }
  return value;
}

static
Value * dimension_choose (Ast * n, Stack * stack, Value * a, Value * b)
{
  if (!can_have_dimension (a)) return b;
  if (!can_have_dimension (b)) return a;
  if ((value_flags (a) & unset) && !value_dimension (a))
    return b;
  if ((value_flags (b) & unset) && !value_dimension (b))
    return a;
  homogeneous_dimensions (n, a, b, stack);
  return value_dimension (a) == dimension_any ? b :
    unknowns (get_dimension (a, stack)) < unknowns (get_dimension (b, stack)) ? a : b;
}

static
bool is_expression (Ast * n)
{
  AstTerminal * t = ast_right_terminal (n);
  return (t->start && !strcmp (t->start, EXPR_MAGIC));
}

static
int compare_dimensions (const void * pa, const void * pb)
{
  const Dimension * a = *((Dimension * const *) pa);
  const Dimension * b = *((Dimension * const *) pb);

  if (unknowns (a) != 1)
    return 1;
  if (unknowns (b) != 1)
    return -1;
  
  /**
  Expressions last. */

  const Key * ac = a->c[0], * bc = b->c[0];
  if (is_expression (ac->parent)) {
    if (!is_expression (bc->parent))
      return 1;
  }
  else if (is_expression (bc->parent))
    return -1;
    
  /**
  Sort first according to the number of non-zero dimensions. */
  
  int na = 0, nb = 0;
  if (a->a)
    for (double * i = a->a; *i != END; i++)
      if (*i) na++;
  if (b->a)
    for (double * i = b->a; *i != END; i++)
      if (*i) nb++;
  if (na != nb)
    return na > nb ? 1 : -1;

  /**
  Then sort according to the sum of the absolute values of the
  dimensions. */

  double sa = 0., sb = 0.;
  if (a->a)
    for (double * i = a->a; *i != END; i++)
      sa += fabs (*i);
  if (b->a)
    for (double * i = b->a; *i != END; i++)
      sb += fabs (*i);
  sa /= fabs (*a->b), sb /= fabs (*b->b);
  if (fabs (sa - sb) > ROUNDOFF)
    return sa > sb ? 1 : -1;

  /**
  Then sort according to the number of dimensions. */
  
  na = dimensions (a), nb = dimensions (b);
  if (na != nb)
    return na > nb ? 1 : -1;

  /**
  Then sort according to the sum of the dimensions. */

  sa = 0., sb = 0.;
  if (a->a)
    for (double * i = a->a; *i != END; i++)
      sa += *i;
  if (b->a)
    for (double * i = b->a; *i != END; i++)
      sb += *i;
  sa /= *a->b, sb /= *b->b;
  if (fabs (sa - sb) > ROUNDOFF)
    return sa > sb ? 1 : -1;
  
  /**
  "Internal" parser values are displayed first. */
  
  AstTerminal * ta = ast_left_terminal (ac->parent), * tb = ast_left_terminal (bc->parent);
  na = strncmp ("ast/", ta->file, 4), nb = strncmp ("ast/", tb->file, 4);
  if (!na && nb) return -1;
  if (na && !nb) return  1;
  
  /**
  Values which are "deeper" in the file system are displayed first. */

  na = 0; for (const char * s = ta->file; *s != '\0'; s++) if (*s == '/') na++;
  nb = 0; for (const char * s = tb->file; *s != '\0'; s++) if (*s == '/') nb++;
  if (na != nb)
    return na > nb ? -1 : 1;
  
  /**
  .c files are displayed last. */
  
  na = strcmp (".c", ta->file + strlen (ta->file) - 2);
  nb = strcmp (".c", tb->file + strlen (tb->file) - 2);
  if (!na && nb) return  1;
  if (na && !nb) return -1;

  /**
  Alphabetic order of file names. */
  
  na = strcmp (ta->file, tb->file);
  if (na)
    return na;

  /**
  Order by line number. */

  if (ta->line != tb->line)
    return ta->line > tb->line ? 1 : -1;

  /**
  Order by label name. */

  return strcmp (ac->label, bc->label);
}

static
bool is_finite_constant (const Key * c)
{
  double val = ast_evaluate_constant_expression (c->parent);
  return val > 1e-30 && val < 1e30 && val != 1234567890;
}

static
void dimension_after_run (Ast * n, Stack * stack)
{
  if (((StackData *)stack_get_data (stack))->maxcalls < 0)
    return;
  System * s = interpreter_get_data (stack);
  s->alloc = stack_alloc (stack);
  if (!s->output) {
    if (!system_pivot (s, stack))
      ((StackData *)stack_get_data (stack))->maxcalls = -1;
    return;
  }
  if (!system_solve (s, stack)) {
    ((StackData *)stack_get_data (stack))->maxcalls = -1;
    return;
  }
  int u = 0;
  foreach_constraint (s, r) {
    u++;
    if (r->c)
      r->c[0]->label = key_label (r->c[0], ' ', LINENO);
  }
  qsort (s->r, u, sizeof (Dimension *), compare_dimensions);
  FILE * fp = s->output;
  Key ** unconstrained = system_unconstrained (s);
  if (unconstrained) {
    int nfinite = 0;
    for (Key ** u = unconstrained; u[0]; u++)
      if (!s->finite || is_finite_constant (u[0]))
	nfinite++;
    if (nfinite == 1 && s->n - s->m == 1)
      fputs ("The following constant is unconstrained\n", fp);
    else {
      int n = 0;
      for (Key ** u = unconstrained; u[0]; u++)
	if (!s->finite || nfinite < s->n - s->m || is_finite_constant (u[0]))
	  n++;
      fprintf (fp, "There %s %d unconstrained constant%s within the following %d\n",
	       s->n - s->m > 1 ? "are" : "is",
	       s->n - s->m, s->n - s->m > 1 ? "s" : "",
	       n);
    }
    for (Key ** u = unconstrained; u[0]; u++)
      if (!s->finite || nfinite < s->n - s->m || is_finite_constant (u[0])) {
	fputs ("  ", fp); print_key_label (u[0], ' ', fp, LINENO /*| INDEX*/);
	fputc ('\n', fp);
      }
    free (unconstrained);
  }
  double * a = NULL, ca = 0.;
  bool firstconst = true, firstexpr = true;
  foreach_constraint (s, r) {
    if (unknowns (r) > 1) {
      // constraint_print (r, fp, LINENO | INDEX | NORIGIN);
      continue;
    }
    if (r->c && r->c[0]->used &&
	(!s->finite ||
	 (is_finite_constant (r->c[0]) // only lists finite constants
	  && (r->a || in_c_file (r->c[0]->parent)) // list dimensionless constants only in C files
	  ))) {
      assert (r->b[0]);
      double coef = - 1./r->b[0];
      bool identical;
      if (a && r->a) {
	double * i, * j;
	for (i = r->a, j = a; fabs ((*i)*coef - (*j)*ca) <= ROUNDOFF && *i != END; i++, j++);
	identical = (*i == END && *j == END);
      }
      else
	identical = (fabs (ca) > ROUNDOFF && a == NULL && r->a == NULL);
      if (!identical) {
	a = r->a, ca = coef;
	if (s->lineno) {
	  if (is_expression (r->c[0]->parent)) {
	    if (firstexpr) {
	      fputs ("Dimensions of expressions\n", fp);
	      firstexpr = false;
	    }
	  }
	  else if (firstconst) {
	    fprintf (fp, "Dimensions of %sconstants\n", s->finite ? "finite " : "");
	    firstconst = false;
	  }
	  double * b = r->b;
	  r->b = NULL;
	  fputs ("  ", fp), dimension_print (r, fp, coef, true);
	  r->b = b;
	  fputc ('\n', fp);
	}
      }

      if (s->redundant || !identical || _p == s->r ||
	  strcmp (r->c[0]->label, (_p - 1)[0]->c[0]->label)) {
	if (!s->lineno) {
	  double * b = r->b;
	  r->b = NULL;
	  dimension_print (r, fp, coef, true);
	  r->b = b;
	}
	fputs (s->lineno ? "    " : "  ", fp);
	if (s->lineno)
	  fputs (r->c[0]->label, fp);
	else {
	  char * s = r->c[0]->label;
	  while (*s != ':')
	    fputc (*s++, fp);
	  s++;
	  while (*s != ':') s++;
	  while (*s != '\0')
	    fputc (*s++, fp);
	}
	fputc ('\n', fp);
      }
    }
  }
  foreach_constraint (s, r)
    if (r->c) {
      free (r->c[0]->label);
      r->c[0]->label = NULL;
    }
}

/**
# Interface

This function performs dimensional analysis and returns `true` if the
code defined by `root` and `n` is dimensionally correct and `false`
otherwise. The parameters are:

* `verbosity` : verbosity level for the interpreter (1 to 4).
* `maxcalls` : maximum number of calls or instructions interpreted.
* `dimensions` : where to output a summary of the dimensions of constants.
* `finite` : if `true` only "finite" constants (i.e. constants
*            different from "zero" or "infinity") are displayed.
* `redundant` : if `true` all constants are listed, even repeated ones.
* `lineno` : if `true` a different output format is used.
* `warn` : if `true` the code will only warn on dimensional errors.

See also [interpreter.c]().
*/

bool ast_check_dimensions (AstRoot * root, Ast * n, int verbosity, int maxcalls,
			   FILE * dimensions, int finite, int redundant, int lineno,
			   int warn)
{
  ast_type_size = dimension_type_size;
  run = dimension_run;
  after_run = dimension_after_run;
  ast_value_print_hook = value_print_dimension_hook;
  ast_assign_hook = dimension_assign;
  ast_binary_operation_hook = dimension_binary_operation;
  ast_internal_functions_hook = dimension_internal_functions;
  ast_choose_hook = dimension_choose;
  System * system = system_new();
  system->output = dimensions;
  system->finite = finite;
  system->redundant = redundant;
  system->lineno = lineno;
  system->warn = warn;
  int remainingcalls = ast_run (root, n, verbosity, maxcalls, system);
  if (!remainingcalls && dimensions) {
    AstTerminal * t = ast_left_terminal (n);
    fprintf (dimensions, "%s:%d: warning (dimensions): run stopped after %d instructions (see also -maxcalls)\n",
	     t->file, t->line, maxcalls);
  }
  system_destroy (system);
  return remainingcalls >= 0;
}
