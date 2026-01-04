/**
# External references 

Returns a list of external references. */

#include <stdlib.h>
#include <string.h>
#include "ast.h"
#include "symbols.h"

typedef struct {
  Ast * scope;
  Stack * nonlocals, * attributes;
  int n, dimension;
} Accelerator;

static
bool is_local_declaration (Ast * n, Stack * stack, Ast * scope)
{
  if (!strcmp (ast_terminal (n)->start, "point"))
    return true;
  Ast ** d;
  for (int i = 0; (d = stack_index (stack, i)); i++)
    if (*d == n)
      return true;
    else if (*d == scope)
      return false;
  return false;
}

static
void external_references (Ast * n, Stack * stack, Accelerator * a);

static
void add_external_reference (const char * name, Stack * stack, Accelerator * a)
{
  Ast * ref = ast_identifier_declaration (stack, name);
  if (!ref) {
#if 0      
    fprintf (stderr, "%s:%d: warning: '%s' undeclared\n", ast_terminal (n)->file, ast_terminal (n)->line,
	     ast_terminal (n)->start);
#endif
    return; // assumes this is OK i.e. this corresponds mostly with macros and undeclared library functions
  }

  if (!strcmp (ast_terminal (ref)->file, "ast/defaults.h")) // ignore "internal" variables and macros
    return;

  if (!is_local_declaration (ref, stack, a->scope) &&
      !fast_stack_find (a->nonlocals, ast_terminal (ref)->start)) {

    /**
    Function call */

    Ast * definition = ast_parent (ref, sym_function_definition), * def;
    if (definition && (def = ast_find (definition, sym_direct_declarator,
				       0, sym_direct_declarator,
				       0, sym_generic_identifier,
				       0, sym_IDENTIFIER)) &&
	!strcmp (ast_terminal (def)->start, name)) {
      Accelerator b = *a;
      b.scope = definition;
      stack_push (stack, &definition);
      external_references (definition, stack, &b);
      ast_pop_scope (stack, definition);
      a->nonlocals = b.nonlocals;
      a->attributes = b.attributes;
    }

    stack_push (a->nonlocals, &ref);
  }
}

static
void external_references (Ast * n, Stack * stack, Accelerator * a)
{
  if (!n || n == ast_placeholder ||
      (n->sym == sym_initializer && n->parent->sym == sym_parameter_declaration))
    return;

  Ast * scope = ast_push_declarations (n, stack);

  if (n->child)
    for (Ast ** c = n->child; *c; c++)
      external_references (*c, stack, a);
  
  if (ast_schema (n->parent, sym_primary_expression,
		  0, sym_IDENTIFIER) &&
      strcmp (ast_terminal (n)->start, "_attribute"))
    add_external_reference (ast_terminal (n)->start, stack, a);
  else if ((n->sym == sym_IDENTIFIER && ast_attribute_access (ast_ancestor (n, 3), stack)) ||
	   (n = ast_attribute_array_access (ast_ancestor (n, 3)))) {
      
    /**
    Scalar attribute */

    Ast * found = fast_stack_find (a->attributes, ast_terminal (n)->start);
    if (!found) {
      Ast * attributes = ast_find (ast_ancestor (ast_identifier_declaration (stack, "_Attributes"), 6),
				   sym_struct_declaration_list);
      assert (attributes);
      found = NULL;
      foreach_item (attributes, 1, decl) {
	Ast * list = ast_schema (decl, sym_struct_declaration,
				 1, sym_struct_declarator_list);
	foreach_item (list, 2, j) {
	  Ast * identifier = ast_find (j, sym_IDENTIFIER);
	  if (identifier && !strcmp (ast_terminal (identifier)->start, ast_terminal (n)->start)) {
	    found = identifier; break;
	  }
	}
	if (found)
	  break;
      }
      if (found)
	stack_push (a->attributes, &found);
    }
  }

  ast_pop_scope (stack, scope);
}

static
char * add_reference (Ast * ref, char * references, Ast * scope, Stack * stack, Stack * functions)
{
  const char * start = ast_terminal (ref)->start;
  if (!strcmp (start, "NULL"))
    return references;

  AstDimensions dim = {0};
  Ast * type = ast_identifier_type (ref, &dim, stack);
  Ast * attributes = ast_parent (ref, sym_struct_or_union_specifier);
  
  if (type == (Ast *) &ast_function) {
    if (!strcmp (start, "qassert"))
      return references;
    
    /**
    Function pointers */
  
    if (ast_schema (ast_ancestor (ref, 4), sym_direct_declarator,
		    1, sym_declarator,
		    0, sym_pointer)) {
      str_append (references, "{.name=\"", attributes ? "." : "", start,
		  "\",.type=sym_function_declaration");
      if (attributes)
	str_append (references, ",.nd=attroffset(", start, ")},");
      else
	str_append (references, ",.pointer=(void *)(long)", start, "},");
    }

    /**
    Function definitions */
    
    else if (ast_ancestor (ref, 6)->sym == sym_function_definition) {
      str_append (references, "{.name=\"", attributes ? "." : "", start,
		  "\",.type=sym_function_definition,.pointer=(void *)(long)", start, "},");
      if (!fast_stack_find (functions, ast_terminal (ref)->start))
	stack_push (functions, &ref);
    }
    
    return references;
  }
  
  str_append (references, "{.name=\"", attributes ? "." : "", !strcmp (start, "val") ? "_val" : start, "\"");
      
  /**
  Type */

  Ast * def;
  if (ast_schema (ast_ancestor (type, 5), sym_declaration,
		  0, sym_declaration_specifiers,
		  0, sym_storage_class_specifier,
		  0, sym_TYPEDEF) &&
      (def = ast_schema (ast_ancestor (type, 5), sym_declaration,
			 1, sym_init_declarator_list,
			 0, sym_init_declarator,
			 0, sym_declarator,
			 0, sym_direct_declarator,
			 0, sym_generic_identifier,
			 0, sym_IDENTIFIER))) {
    // typedef
    if (!strcmp (ast_terminal (def)->start, "scalar"))
      str_append (references, ",.type=sym_SCALAR");
    else if (!strcmp (ast_terminal (def)->start, "vector"))
      str_append (references, ",.type=sym_VECTOR");
    else if (!strcmp (ast_terminal (def)->start, "tensor"))
      str_append (references, ",.type=sym_TENSOR");
    else if (!strcmp (ast_terminal (def)->start, "coord"))
      str_append (references, ",.type=sym_COORD");
    else if (!strcmp (ast_terminal (def)->start, "_coord"))
      str_append (references, ",.type=sym__COORD");
    else if (!strcmp (ast_terminal (def)->start, "vec4"))
      str_append (references, ",.type=sym_VEC4");
    else if (!strcmp (ast_terminal (def)->start, "ivec"))
      str_append (references, ",.type=sym_IVEC");
    else if (!strcmp (ast_terminal (def)->start, "bool"))
      str_append (references, ",.type=sym_BOOL");
    else
      str_append (references, ",.type=sym_TYPEDEF");
  }
  else if (ref->parent->sym == sym_enumeration_constant)
    str_append (references, ",.type=sym_enumeration_constant");
  else {
    char s[20]; snprintf (s, 19, "%d", type->sym);
    str_append (references, ",.type=", s);
  }

  /**
  Is this a global variable? */

  if (!attributes && !ast_parent (ref, sym_compound_statement) && !ast_parent (ref, sym_parameter_declaration))
    str_append (references, ",.global=1");
  
  /**
  Assumes 'double *' are references to arrays with 'nl'
  elements. Fixme: this is very specific and should be made more
  general e.g. systematically using 'fat pointers' to get array
  sizes. */

  Ast * nl;
  if (type->sym == sym_DOUBLE && dim.pointer == 1 && !dim.dimension && strcmp (ast_terminal (ref)->start, "_constant") &&
      (nl = ast_identifier_declaration (stack, "nl"))) {
    dim.pointer = 0;
    dim.dimension = malloc (2*sizeof (Ast *));
    dim.dimension[0] = nl;
    dim.dimension[1] = NULL;
  }
    
  /**
  Pointer */

  if (!(attributes || ref->parent->sym == sym_enumeration_constant))
    str_append (references, ",.pointer=(void *)", dim.pointer || (dim.dimension && (*dim.dimension)->sym != sym_VOID) ?
		"" : "&", start);

  /**
  Attribute offset or enumeration constant or number of pointer dereferences */

  if (attributes)
    str_append (references, ",.nd=attroffset(", start, ")");
  else if (ref->parent->sym == sym_enumeration_constant)
    str_append (references, ",.nd=", start);
  else if (dim.pointer) {
    char s[10];
    snprintf (s, 10, "%d", dim.pointer);
    str_append (references, ",.nd=", s);
  }
    
  /**
  Reduction */
    
  Ast * parameters = ast_child (scope, sym_argument_expression_list);
  if (parameters)
    foreach_item (parameters, 2, item) {
      Ast * reductions = ast_find (item, sym_reduction_list);
      foreach_item (reductions, 1, reduction) {
	Ast * identifier = ast_schema (reduction, sym_reduction,
				       4, sym_reduction_array,
				       0, sym_generic_identifier,
				       0, sym_IDENTIFIER);
	if (!strcmp (ast_terminal (identifier)->start, start)) {
	  char * operator = ast_left_terminal (reduction->child[2])->start;
	  Ast * array = ast_schema (reduction, sym_reduction,
				    4, sym_reduction_array,
				    3, sym_expression);
	  if (array) {
	    // fixme: not implemented yet
	  }
	  else
	    str_append (references, ",.reduct=",
			!strcmp(operator, "min") ? "'m'" :
			!strcmp(operator, "max") ? "'M'" :
			!strcmp(operator, "+")   ? "'+'" :
			"'?'");
	}
      }
    }
  
  /**
  Array dimensions */
  
  if (dim.dimension) {
    if ((*dim.dimension)->sym != sym_VOID) {
      str_append (references, ",.data=(int[]){");    
      for (Ast ** d = dim.dimension; *d; d++)
	references = ast_str_append (*d, references);
      str_append (references, ",0}");
    }
    else if (*(dim.dimension + 1) == NULL) { // a single undefined dimension
      Ast * initializer = ast_schema (ast_parent (ref, sym_init_declarator), sym_init_declarator,
				      2, sym_initializer,
				      1, sym_initializer_list);
      if (initializer) {
	int n = 0;
	foreach_item (initializer, 2, item)
	  n++;
	char s[20];
	snprintf (s, 19, "%d", n);
	str_append (references, ",.data=(int[]){", s, ",0}");
      }
    }
  }
  free (dim.dimension);
    
  str_append (references, "},");
  
  return references;
}

char * ast_external_references (Ast * scope, char * references, Stack * functions)
{
  AstRoot * root = ast_get_root (scope);
  Stack * stack = root->stack;

  stack_push (stack, &scope);
  Accelerator a = { scope };
  a.nonlocals = stack_new (sizeof (Ast *));
  a.attributes = stack_new (sizeof (Ast *));
  external_references (scope, stack, &a);
  ast_pop_scope (stack, scope);

  Ast ** n;
  for (int i = 0; (n = stack_indexi (a.attributes, i)) && (!references || !strstr (references, "@error ")); i++)
    references = add_reference (*n, references, scope, stack, functions);
  for (int i = 0; (n = stack_indexi (a.nonlocals, i)) && (!references || !strstr (references, "@error ")); i++)
    references = add_reference (*n, references, scope, stack, functions);

  stack_destroy (a.nonlocals);
  stack_destroy (a.attributes);
  
  return references;
}
