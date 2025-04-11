/**
# Computation kernels */

#include <stdlib.h>
#include <string.h>
#include "ast.h"
#include "symbols.h"

typedef struct {
  char * error;
  bool nolineno;
  int return_macro_index;
  Ast * macroscope;
} KernelData;

/**
## Implicit type casting

GLSL does not support implicit type casting, so we insert the
necessary explicit casts using the function below. */

static
Ast * type_cast (Ast * n, const char * type)
{
  Ast * parent = n->parent;
  int child = ast_child_index (n);
  Ast * call = NN(n, sym_function_call,
		  NN(n, sym_postfix_expression,
		     NN(n, sym_primary_expression,
			NA(n, sym_IDENTIFIER, type))),
		  NCA(n, "("),
		  n,
		  NCA(n, ")"));
  ast_set_child (parent, child, call);
  return call;
}

static
Ast * implicit_type_cast (Ast * n, Stack * stack)
{
  if (!n) return NULL;
  
  Ast * type = NULL;
  
  switch (n->sym) {

  case sym_I_CONSTANT: case sym_ENUMERATION_CONSTANT:
    return (Ast *) &ast_int;

  case sym_F_CONSTANT:
    return (Ast *) &ast_double;

  case sym_types:
    if (ast_terminal (n->child[0]))
      return n->child[0];
    else
      return NULL;

  case sym_IDENTIFIER:
    if (n->parent->sym == sym_primary_expression) {
      Ast * ref = ast_identifier_declaration (stack, ast_terminal (n)->start);
      if (ref) {
	AstDimensions dim = {0};
	type = ast_identifier_type (ref, &dim, stack);
	if (type && type->sym == sym_INT &&
	    (ref = ast_schema (ast_ancestor (type, 5), sym_declaration,
			       1, sym_init_declarator_list,
			       0, sym_init_declarator,
			       0, sym_declarator,
			       0, sym_direct_declarator,
			       0, sym_generic_identifier,
			       0, sym_IDENTIFIER)) && !strcmp (ast_terminal (ref)->start, "bool"))
	  type = (Ast *) &ast_bool;
	if (dim.pointer)
	  type = NULL;
      }
      return type;
    }
    if (ast_schema (ast_ancestor (n, 2), sym_member_identifier,
		    0, sym_generic_identifier)) {
      n = ast_expression_type (n, stack, false);
      if (ast_schema (ast_ancestor (n, 2), sym_direct_declarator))
	return implicit_type_cast (ast_schema (ast_parent (n, sym_struct_declaration), sym_struct_declaration,
					       0, sym_specifier_qualifier_list,
					       0, sym_type_specifier,
					       0, sym_types),
				   stack);
      return NULL;
    }
    if (ast_schema (ast_ancestor (n, 2), sym_direct_declarator))
      return implicit_type_cast (ast_schema (ast_parent (n, sym_declaration), sym_declaration,
					     0, sym_declaration_specifiers,
					     0, sym_type_specifier,
					     0, sym_types),
				 stack);
    break;

  case sym_unary_expression:
    switch (n->child[0]->sym) {
    case sym_SIZEOF: case sym_ALIGNOF:
      return (Ast *) &ast_int;
    case sym_unary_operator:
      if (n->child[0]->child[0]->sym == token_symbol('!')) {
	Ast * a = implicit_type_cast (n->child[1], stack);
	if (!a || a->sym != sym_BOOL)
	  type_cast (n->child[1], "bool");
	return (Ast *) &ast_bool;
      }
    }
    break;
    
  case sym_additive_expression:
  case sym_multiplicative_expression:
    if (n->child[1]) {
      Ast * a = implicit_type_cast (n->child[0], stack);
      if (a && a->sym == sym_BOOL)
	type_cast (n->child[0], "int");
      Ast * b = implicit_type_cast (n->child[2], stack);
      if (b && b->sym == sym_BOOL)
	type_cast (n->child[2], "int");
      if (b && (b->sym == sym_DOUBLE || b->sym == sym_FLOAT))
	return b;
      return a;
    }
    break;

  case sym_assignment_expression:
  case sym_init_declarator:
    if (n->child[1]) {
      Ast * a = implicit_type_cast (n->child[0], stack);
      Ast * b = implicit_type_cast (n->child[2], stack);
      if (a && b && a->sym != b->sym && (a->sym == sym_INT || b->sym == sym_BOOL))
	type_cast (n->child[2], "int");
      return a;
    }
    break;

  case sym_cast_expression: case sym_primary_expression:
    if (n->child[1])
      return implicit_type_cast (n->child[1], stack);
    break;

  case sym_relational_expression:
  case sym_equality_expression:
    if (n->child[1])
      return (Ast *)&ast_bool;
    break;    
    
  case sym_logical_and_expression:
  case sym_logical_or_expression:
    if (n->child[1]) {
      Ast * a = implicit_type_cast (n->child[0], stack);
      if (!a || a->sym != sym_BOOL)
	type_cast (n->child[0], "bool");
      Ast * b = implicit_type_cast (n->child[2], stack);
      if (!b || b->sym != sym_BOOL)
	type_cast (n->child[2], "bool");
      return (Ast *)&ast_bool;
    }
    break;

  case sym_conditional_expression:
    if (n->child[1]) {
      Ast * type = implicit_type_cast (n->child[0], stack);
      if (!type || type->sym != sym_BOOL) {
	str_prepend (ast_left_terminal (n->child[0])->before, "bool(");
	ast_after (n->child[0], ")");
      }
      return implicit_type_cast (n->child[2], stack);
    }
    break;

  case sym_selection_statement: {
    Ast * type = implicit_type_cast (n->child[2], stack);
    if (!type || type->sym != sym_BOOL)
      type_cast (n->child[2], "bool");
    return type;
  }

  case sym_function_call: {
    Ast * identifier = ast_function_call_identifier (n);
    if (identifier) {
      if (!strcmp (ast_terminal (identifier)->start, "val") ||
	  !strcmp (ast_terminal (identifier)->start, "val_out_"))
      return (Ast *) &ast_double;
      if (!strcmp (ast_terminal (identifier)->start, "_attr"))
	return implicit_type_cast (n->child[2], stack);
      identifier = ast_identifier_declaration (stack, ast_terminal (identifier)->start);
      Ast * declaration = ast_parent (identifier, sym_function_declaration);
      if (!declaration)
	declaration = ast_parent (identifier, sym_declaration);
      if (declaration)
	return implicit_type_cast (ast_child (declaration, sym_declaration_specifiers), stack);
      else
	return NULL;
    }
    break;
  }

  }
  if (n->child)
    for (Ast ** c = n->child; *c; c++)
      type = implicit_type_cast (*c, stack);
    
  return type;
}

static
void kernel (Ast * n, Stack * stack, void * data)
{
  KernelData * d = data;

  if (d->error)
    return;
  
  switch (n->sym) {
    
  case sym_IDENTIFIER: {

    /**
    ## Function pointers */

    if (ast_is_function_pointer (n, stack))
      str_prepend (ast_terminal (n)->start, "_p");

    /**
    ## 'val'

    Buggy GLSL preprocessors do not make the difference between 'val'
    as a variable identifier and 'val()' as a macro call. */

    else if (!ast_schema (ast_ancestor (n, 3), sym_function_call,
			  0, sym_postfix_expression,
			  0, sym_primary_expression) &&
	     !strcmp (ast_terminal (n)->start, "val"))
      free (ast_terminal (n)->start), ast_terminal (n)->start = strdup ("_val");
    
    break;
  }
    
  /**
  ## Assumes that pointers to structures are used through "inout" parameter */

  case sym_PTR_OP:
    ast_terminal(n)->start[0] = '.';
    ast_terminal(n)->start[1] = '\0';
    break;

  /**
  ## Remove some reserved GLSL keywords */
    
  case sym_STATIC: case sym_INLINE:
    ast_terminal (n)->start[0] = '\0';
    break;

  /**
  ## 'unsigned' or 'unsigned int' is replaced by 'uint' */

  case sym_UNSIGNED: {
    Ast * identifier;
    if ((identifier = ast_schema (ast_ancestor (n, 3), sym_declaration_specifiers,
				  1, sym_declaration_specifiers,
				  0, sym_type_specifier,
				  0, sym_types,
				  0, sym_INT))) {
      ast_terminal (n)->start[0] = '\0';
      free (ast_terminal (identifier)->start);
      ast_terminal (identifier)->start = strdup ("uint");
    }
    else if (ast_ancestor (n, 4)->sym == sym_declaration)
      strcpy (ast_terminal (n)->start, "uint");
    break;
  }

  /** 
  ## Initialization of `coord` (must include three components) */

  case sym_TYPEDEF_NAME: {
    if (n->parent->sym == sym_types && !strcmp (ast_terminal(n)->start, "coord") &&
	(n = ast_schema (ast_parent (n, sym_declaration), sym_declaration,
			 1, sym_init_declarator_list)))
      foreach_item (n, 2, init_declarator) {
	Ast * list;
	if ((list = ast_schema (init_declarator, sym_init_declarator,
				2, sym_initializer,
				1, sym_initializer_list))) {
	  int i = 0;
	  foreach_item (list, 2, init)
	    i++;
	  Ast * last = ast_schema (list->parent, sym_initializer,
				   2, token_symbol('}'));
	  if (last)
	    while (i++ < 3)
	      str_prepend (ast_terminal (last)->start, ",0");
	}
      }
    break;
  }

  /**
  ## Implicit type casts */

  case sym_jump_statement:
    implicit_type_cast (ast_child (n, sym_expression), stack);
    break;
    
  case sym_assignment_expression:
    if (n->child[1])
      implicit_type_cast (n, stack);
    break;
    
  case sym_init_declarator:
  case sym_selection_statement:
    implicit_type_cast (n, stack);
    break;

  case sym_for_declaration_statement:
  case sym_iteration_statement:
    if (n->child[0]->sym == sym_WHILE || n->child[0]->sym == sym_DO) {
      Ast * expr = ast_child (n, sym_expression);
      Ast * type = implicit_type_cast (expr, stack);
      if (!type || type->sym != sym_BOOL)
	type_cast (expr, "bool");
    }
    else if (n->child[0]->sym == sym_for_scope) {
      Ast * expr = ast_schema (n->child[3], sym_expression_statement,
			       0, sym_expression);
      if (expr) {
	Ast * type = implicit_type_cast (expr, stack);
	if (!type || type->sym != sym_BOOL)
	  type_cast (expr, "bool");
      }
    }
    break;

  /**
  ## Typedef struct */

  case sym_TYPEDEF: {
    Ast * struct1, * identifier;
    if ((struct1 = ast_schema (ast_ancestor (n, 2), sym_declaration_specifiers,
			       1, sym_declaration_specifiers,
			       0, sym_type_specifier,
			       0, sym_types,
			       0, sym_struct_or_union_specifier,
			       0, sym_struct_or_union,
			       0, sym_STRUCT)) &&
	(identifier = ast_schema (ast_ancestor (n, 3), sym_declaration,
				  1, sym_init_declarator_list,
				  0, sym_init_declarator,
				  0, sym_declarator,
				  0, sym_direct_declarator,
				  0, sym_generic_identifier,
				  0, sym_IDENTIFIER))) {
      char * s = ast_terminal (n)->start; s[0] = '\0';
      ast_terminal (n)->start = ast_terminal (struct1)->start;
      ast_terminal (struct1)->start = ast_terminal (identifier)->start;
      ast_terminal (identifier)->start = s;
    }
    break;
  }
    
  case sym_postfix_expression: {
    Ast * list;
    if ((list = ast_schema (n, sym_postfix_expression,
			    3, sym_postfix_initializer,
			    1, sym_initializer_list))) {
      
      /**
      ## Postfix initializers */
      
      Ast * a = n->child[0];
      ast_set_child (n, 0, n->child[1]);
      ast_set_child (n, 1, a);
      a = n->child[3];
      ast_set_child (n, 3, n->child[2]);
      ast_set_child (n, 2, list);
      ast_destroy (a);
    }
    else if (ast_attribute_access (n, stack) || ast_attribute_array_access (n)) {

      /**
      ## Attribute access */

      if (n->parent->sym == sym_function_call) {
	Ast * identifier = ast_schema (n, sym_postfix_expression,
				       2, sym_member_identifier,
				       0, sym_generic_identifier,
				       0, sym_IDENTIFIER);
	ast_before (n->parent, "_attr_", ast_terminal (identifier)->start, "(");
	ast_terminal (ast_schema (n, sym_postfix_expression,
				  1, token_symbol('.')))->start[0] = ',';
	ast_terminal (identifier)->start[0] = '\0';
	ast_after (n->parent, ")");
      }
      else {
	if (ast_schema (n, sym_postfix_expression,
			0, sym_postfix_expression,
			0, sym_array_access)) {
	  Ast * scalar = ast_find (n, sym_unary_expression,
				   0, sym_postfix_expression);
	  ast_destroy (scalar->child[1]);
	  ast_destroy (scalar->child[2]);
	  scalar->child[1] = NULL;
	  Ast * array = ast_find (n, sym_array_access);
	  scalar = ast_find (array, sym_expression);
	  ast_set_child (n, 0, scalar);
	  ast_destroy (array);
	}
	ast_terminal (ast_schema (n, sym_postfix_expression,
				  1, token_symbol('.')))->start[0] = ',';
	type_cast (n, "_attr");
      }
    }
    
    break;
  }

  /**
  ## Arrays as parameters 

  This forces arrays passed as parameters to functions to behave like
  in C99 i.e. passing by reference (inout) rather than by value.  */
    
  case sym_parameter_declaration:
    if (ast_schema (n, sym_parameter_declaration,
		    1, sym_declarator,
		    0, sym_direct_declarator,
		    2, sym_assignment_expression))
      ast_before (n, "inout ");
    break;
    
  /**
  ## Cast expressions */

  case sym_cast_expression:
    if (ast_schema (n, sym_cast_expression,
		    1, sym_type_name)) {
      Ast * a = n->child[0];
      ast_set_child (n, 0, n->child[1]);
      ast_set_child (n, 1, a);
      a = n->child[3];
      ast_set_child (n, 3, n->child[2]);
      ast_set_child (n, 2, a);
    }
    break;

  case sym_pointer: {
    Ast * p, * type, * identifier;
    if ((p = ast_schema (n, sym_pointer,
			 0, token_symbol('*'))) &&
	(type = ast_schema (ast_ancestor (n, 2), sym_parameter_declaration,
			    0, sym_declaration_specifiers,
			    0, sym_type_specifier,
			    0, sym_types,
			    0, sym_TYPEDEF_NAME)) &&
	(!strcmp (ast_terminal (type)->start, "scalar") ||
	 !strcmp (ast_terminal (type)->start, "vector") ||
	 !strcmp (ast_terminal (type)->start, "tensor")) &&
	(identifier = ast_schema (n->parent, sym_declarator,
				  1, sym_direct_declarator,
				  0, sym_generic_identifier,
				  0, sym_IDENTIFIER))) {

      /**
      ## Scalar, vector and tensor lists parameters */

      ast_terminal (p)->start[0] = '\0';
      ast_after (identifier, "[2]"); // fixme: need to set the correct fixed size
      //      ast_print_tree (ast_ancestor (n, 2), stderr, 0, 0, -1);
    }
    else if ((identifier = ast_schema (ast_ancestor (n, 2), sym_parameter_declaration,
				       1, sym_declarator,
				       0, sym_pointer,
				       0, token_symbol ('*'))) &&
	     (type = ast_schema (ast_ancestor (n, 2), sym_parameter_declaration,
				 0, sym_declaration_specifiers,
				 0, sym_type_specifier,
				 0, sym_types)) &&
	     (type->child[0]->sym == sym_DOUBLE ||
	      type->child[0]->sym == sym_FLOAT ||
	      type->child[0]->sym == sym_INT ||
	      type->child[0]->sym == sym_TYPEDEF_NAME)) {

      /**
      ## "inout" function parameter */
      
      ast_before (type->child[0], "inout ");
      ast_terminal (identifier)->start[0] = '\0';
    }
    
    break;
  }

  /**
  ## forin_declaration_statement */

  case sym_forin_declaration_statement: {
    Ast * declarator = n->child[3];
    Ast * identifier = ast_schema (declarator, sym_declarator,
				   0, sym_direct_declarator,
				   0, sym_generic_identifier,
				   0, sym_IDENTIFIER);
    if (!identifier) {
      AstTerminal * t = ast_left_terminal (n);
      fprintf (stderr,
	       "%s:%d: error: incorrect declaration\n",
	       t->file, t->line);
      exit (1);
    }
    ast_before (n, "{");
    ast_after ((Ast *)ast_right_terminal (n->child[0]), "in");
    ast_after ((Ast *)ast_right_terminal (n->child[2]), ",");
    free (ast_terminal (n->child[4])->start);
    ast_terminal (n->child[4])->start = strdup (",");
    ast_after (n, " endforin()}");
    break;
  }

  /**
  ## forin_statement */

  case sym_forin_statement: {
    int narg = 0;
    foreach_item (n->child[2], 2, expr)
      narg++;
    ast_before (n, "{");
    char suffix[20]; snprintf (suffix, 19, "%d", narg);
    ast_after ((Ast *)ast_right_terminal (ast_child (n, sym_for_scope)), "in", suffix);
    free (ast_terminal (ast_child (n, sym_IDENTIFIER))->start);
    ast_terminal (ast_child (n, sym_IDENTIFIER))->start = strdup (",");
    ast_after (n, " endforin", suffix, "()}");
    break;
  }

  case sym_unary_operator: {
    
    /**
    ## Dereference of "inout" parameters */

    Ast * identifier, * ref, * type;
    if (n->child[0]->sym == token_symbol ('*') &&
	(identifier = ast_schema (n->parent, sym_unary_expression,
				  1, sym_cast_expression,
				  0, sym_unary_expression,
				  0, sym_postfix_expression,
				  0, sym_primary_expression,
				  0, sym_IDENTIFIER)) &&
	(ref = ast_identifier_declaration (stack, ast_terminal (identifier)->start)) &&
	(type = ast_schema (ast_ancestor (ref, 4), sym_parameter_declaration,
			    0, sym_declaration_specifiers,
			    0, sym_type_specifier,
			    0, sym_types)) &&
	ast_terminal (type->child[0])->before &&
	!strcmp (ast_terminal (type->child[0])->before + strlen (ast_terminal (type->child[0])->before) - 6, "inout ")) {
      free (ast_terminal (n->child[0])->start);
      ast_terminal (n->child[0])->start = strdup("");
    }

    /**
    ## References 

    We replace the '&' with the 'ast_pointer()' macro. */
    
    Ast * ampersand;
    if ((ampersand = ast_schema (n, sym_unary_operator,
				 0, token_symbol ('&')))) {
      free (ast_terminal (ampersand)->start);
      ast_terminal (ampersand)->start = strdup ("ast_pointer(");
      Ast * cast = ast_schema (n->parent, sym_unary_expression,
			       1, sym_cast_expression);
      ast_after (cast, ")");
    }
    break;
  }
    
  case sym_function_call: {
    Ast * identifier = ast_function_call_identifier (n);
    if (!identifier) break;
    AstTerminal * t = ast_terminal (identifier);
    if (!t) break;
    
    /**
    ## Field assignments 
  
    Kernels often need to know the type of access to fields (i.e. read
    or write). Here we append "_out" to stencil access functions linked
    to assignments (i.e. "write" operations). */
    
    if (!strcmp (t->start, "val")) {
      if (ast_child (ast_parent (n, sym_assignment_expression), sym_assignment_operator))
	str_append (t->start, "_out_");
      break;
    }

    if (!strcmp (t->start, "coarse") || !strcmp (t->start, "fine"))
      break;
    
    /**
    ## Undeclared or unsupported functions */
    
    if (!(identifier = ast_identifier_declaration (stack, t->start))) {
      char s[1000];
      snprintf (s, 999, "\\n@error %s:%d: GLSL: error: unknown function '%s'\\n",
		t->file, d->nolineno ? 0 : t->line, t->start);
      d->error = strdup (s);
      return;
    }
    
    /**
    ## Function pointers */

    if (ast_schema (ast_ancestor (identifier, 3), sym_declarator,
		    0, sym_pointer,
		    0, token_symbol('*'))) {
      char * s = NULL;
      str_append (s, "_f", t->start);
      free (t->start);
      t->start = s;
      AstTerminal * o = ast_terminal (ast_child (n, token_symbol('(')));
      free (o->start); o->start = strdup ("((");
      AstTerminal * c = ast_terminal (ast_child (n, token_symbol(')')));
      free (c->start); c->start = strdup ("))");
      break;
    }
    
    break;
  }
    
  case sym_macro_statement: {
    Ast * identifier = ast_schema (n, sym_macro_statement,
				   0, sym_MACRO);
    
    /**
    ## Diagonalize */

    if (identifier && !strcmp (ast_terminal (identifier)->start, "diagonalize")) {
      Ast * field = ast_schema (n, sym_macro_statement,
				2, sym_argument_expression_list,
				0, sym_argument_expression_list_item,
				0, sym_assignment_expression);
      if (field && (field = ast_is_identifier_expression (field))) {
	stack_push (stack, &n);
	ast_traverse (n, stack, ast_diagonalize, field);
	ast_pop_scope (stack, n);
      }
    }

    /**
    ## Macro call 

    These are the remaining "postmacros" which have not been expanded
    yet. */

    else {
      ast_before (n, "{");
      ast_after (n, "end_", ast_terminal (identifier)->start, "()}");
    }

    break; 
  }
    
  }
}

static
char * stringify (Ast * n, char * output, bool nolineno)
{
  AstTerminal * t = ast_left_terminal (n);
  char * before = t->before;
  if (n->sym == sym_function_definition)
    t->before = NULL;
  char * s = ast_str_print (n, NULL, 0, 1);
  t->before = before;
  for (char * i = s; *i; i++) {
    char a[] = "1";
    switch (*i) {
    case '\n': str_append (output, "\\n"); break;
    case '\\': str_append (output, "\\\\"); break;
    case '"':  str_append (output, "\\\""); break;
    case '#':
      if (i[-1] == '\n') {
	str_append (output, "// #");
	if (nolineno && !strncmp (i, "#line ", 6)) {
	  str_append (output, "line 0");
	  i += 6;
	  while (*i >= '0' && *i <= '9') i++;
	  i--;
	}
	break;
      }
      // fall through
    default:  a[0] = *i; str_append (output, a); break;
    }
  }
  free (s);
  ast_destroy (n);
  return output;
}

static void postmacros (Ast * n, Stack * stack, void * data)
{
  if (n->sym == sym_statement || n->sym == sym_function_call) {
    KernelData * d = data;
    ast_macro_replacement (n, n, stack, d->nolineno, 1, false, &d->return_macro_index, d->macroscope);
  }
}

char * ast_kernel (Ast * n, char * argument, bool nolineno, Ast * macroscope)
{
  AstRoot * root = ast_get_root (n);
  Stack * stack = root->stack;
  stack_push (stack, &n);
  KernelData d = {0, nolineno, 0, macroscope};
  Ast * statement = n->sym == sym_function_definition ?
    ast_copy (n) : ast_copy (ast_child (n, sym_statement));
  ast_traverse (statement, stack, postmacros, &d);
  ast_traverse (statement, stack, kernel, &d);

  if (d.error)
    str_append (argument, "\"", d.error, "\"");
  else {
    str_append (argument, "\"");
    argument = stringify (statement, argument, nolineno);
    str_append (argument, "\"");
  }
  free (d.error);
  ast_pop_scope (stack, n);
  return argument;
}
