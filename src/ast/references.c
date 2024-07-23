/**
# Non-local references 

Returns a list of non-local references. */

#include <stdlib.h>
#include <string.h>
#include "ast.h"
#include "symbols.h"

typedef struct {
  Ast * scope;
  Stack * nonlocals;
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
void non_local_references (Ast * n, Stack * stack, void * data);

static
void add_non_local_reference (const char * name, Stack * stack, Accelerator * a)
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
    stack_push (a->nonlocals, &ref);

    /**
    Function call */

    Ast * definition = ast_parent (ref, sym_function_definition);
    if (definition && (ref = ast_find (definition, sym_direct_declarator,
				       0, sym_direct_declarator,
				       0, sym_generic_identifier,
				       0, sym_IDENTIFIER)) &&
	!strcmp (ast_terminal (ref)->start, name)) {
      Accelerator b = *a;
      b.scope = definition;
      stack_push (stack, &definition);
      ast_traverse (definition, stack, non_local_references, &b);
      ast_pop_scope (stack, definition);
    }
  }
}

static
void non_local_references (Ast * n, Stack * stack, void * data)
{
  if (n->sym == sym_IDENTIFIER &&
      n->parent->sym == sym_primary_expression)  // do not consider structure members
    add_non_local_reference (ast_terminal (n)->start, stack, data);
}

void ast_non_local_references (Ast * scope, Ast * argument)
{
  AstRoot * root = ast_get_root (scope);
  Stack * stack = root->stack;
  
  stack_push (stack, &scope);
  Accelerator a = { scope };
  a.nonlocals = stack_new (sizeof (Ast *));
  add_non_local_reference ("X0", stack, &a);
  add_non_local_reference ("Y0", stack, &a);
  add_non_local_reference ("Z0", stack, &a);
  add_non_local_reference ("L0", stack, &a);
  add_non_local_reference ("N", stack, &a);
  ast_traverse (scope, stack, non_local_references, &a);
  ast_pop_scope (stack, scope);

  ast_after (argument, "{");
  Ast ** n;
  while ((n = stack_pop (a.nonlocals))) {
    AstDimensions dim = {0};
    Ast * type = ast_identifier_type (*n, &dim, stack), * def;
    
    if (type == (Ast *) &ast_function)
      continue; // fixme: external functions not handled yet

    const char * start = ast_terminal (*n)->start;
    if (!strcmp (start, "NULL"))
      continue;
    
    ast_after (argument, "{\"", start, "\""); // name
    
    /**
    Type */

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
			   0, sym_IDENTIFIER)))
      // typedef
      ast_after (argument, ",\"", ast_terminal (def)->start, "\"");
    else switch (type->sym) {
      case sym_INT: case sym_DOUBLE: case sym_FLOAT:
	ast_after (argument, ",\"", ast_terminal (type)->start, "\""); break;
	
      default:
	ast_after (argument, ",\"not implemented yet\""); break;
    }
    

    /**
    Pointer */

    ast_after (argument, ",(void *)", dim.pointer || dim.dimension ? "" : "&", start);

    /**
    Array dimensions */

    if (dim.dimension) {
      ast_after (argument, ",(int[]){");
      for (Ast ** d = dim.dimension; *d; d++)
	ast_right_terminal (argument)->after = ast_str_append (*d, ast_right_terminal (argument)->after);
      ast_after (argument, ",0}");
    }
    else
      ast_after (argument, ",NULL");
    free (dim.dimension);

    /**
    Number of pointer dereferences */
    
    char s[10];
    snprintf (s, 10, "%d", dim.pointer);
    ast_after (argument, ",", s);
    
    /**
    Reduction */
    
    Ast * parameters = ast_child (scope, sym_foreach_parameters);
    if (parameters)
      foreach_item (parameters, 2, item) {
	if (item->child[0]->sym == sym_reduction_list) {
	  Ast * reductions = item->child[0];
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
		ast_after (argument, ",",
			   !strcmp(operator, "min") ? "'m'" :
			   !strcmp(operator, "max") ? "'M'" :
			   !strcmp(operator, "+")   ? "'+'" :
			   "'?'");
	    }
	  }
	}
      }
    
    ast_after (argument, "},");
  }
  ast_after (argument, "{0}}");
  stack_destroy (a.nonlocals);
}
