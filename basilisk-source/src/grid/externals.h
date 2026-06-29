typedef struct _External External;

struct _External {
  char * name;    // the name of the variable
  void * pointer; // a pointer to the data
  int type;       // the type of the variable
  int nd;         // the number of pointer dereferences or attribute offset or enum constant
  char reduct;    // the reduction operation
  char global;    // is it a global variable?
  char constant;  // is it a constant?
  void * data;    // the dimensions (int *) for arrays or the code (char *) for functions
  scalar s;       // used for reductions on GPUs
  External * externals, * next;
  int used;
};

#if _GPU
static size_t external_size (const External * g)
{
  switch (g->type) {
  case sym_INT: return sizeof (int);
  case sym_LONG: return sizeof (long);
  case sym_FLOAT: return sizeof (float);
  case sym_DOUBLE: return sizeof (double);
  case sym__COORD: return sizeof (_coord);
  case sym_COORD: return sizeof (coord);
  case sym_BOOL: return sizeof (bool);
  case sym_VEC4: return sizeof (vec4);
  default: return 0;
  }
}

static bool is_external_constant (const External * g) {
  return g->constant && !g->data && external_size (g);
}

static bool is_normal_variable (const External * g)
{
  if (g->name[0] == '.') return false;
  if (g->type == sym_function_declaration || g->type == sym_function_definition) return false;
  if (g->type == sym_INT && g->global &&
      (!strcmp (g->name, "N") ||
       !strcmp (g->name, "nl") ||
       !strcmp (g->name, "bc_period_x") ||
       !strcmp (g->name, "bc_period_y")))
    return false;
  return true;
}

static bool is_external_variable (const External * g)
{
  if (is_external_constant (g) || !is_normal_variable (g))
    return false;
  if (g->type == sym_INT ||
      g->type == sym_LONG ||
      g->type == sym_FLOAT ||
      g->type == sym_DOUBLE ||
      g->type == sym__COORD ||
      g->type == sym_COORD ||
      g->type == sym_BOOL ||
      g->type == sym_VEC4) {
    assert (!g->nd);
    assert (!g->data || ((int *)g->data)[1] == 0);
    return true;
  }
  return false;
}
#endif // _GPU
