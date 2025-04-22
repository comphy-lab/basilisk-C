#include <stdlib.h>
#include <string.h>
#include "ast.h"
#include "symbols.h"

AstTerminal
  ast_int = { {sym_INT}, .start = "ast_int", .file = __FILE__, .line = __LINE__ },
  ast_long = { {sym_LONG}, .start = "ast_long", .file = __FILE__, .line = __LINE__ },
  ast_enum = { {sym_INT}, .start = "ast_int", .file = __FILE__, .line = __LINE__ },
  ast_char = { {sym_CHAR}, .start = "ast_char", .file = __FILE__, .line = __LINE__ },
  ast_void = { {sym_VOID}, .start = "ast_void", .file = __FILE__, .line = __LINE__ },
  ast_float = { {sym_FLOAT}, .start = "ast_float", .file = __FILE__, .line = __LINE__ },
  ast_double = { {sym_DOUBLE}, .start = "ast_double", .file = __FILE__, .line = __LINE__ },
  ast_function = { {sym_function_definition}, .start = "ast_function", .file = __FILE__, .line = __LINE__ },
  ast_bool = { {sym_BOOL}, .start = "ast_bool", .file = __FILE__, .line = __LINE__ };

Ast * ast_base_type (Ast * type, AstDimensions * d, Stack * stack)
{
  if (type && (type->sym == sym_TYPEDEF_NAME || type->sym == sym_IDENTIFIER)) {
    Ast * type_identifier = ast_identifier_declaration (stack, ast_terminal (type)->start);
    if (!type_identifier)
      return NULL; // fixme: message (NULL, type, "unknown typedef '%s'\n", warning_verbosity, stack);
    return ast_identifier_type (type_identifier, d, stack);
  }
  else
    return type;
}

Ast * ast_get_array_dimensions (Ast * direct_declarator, int symbol, AstDimensions * d, int nd, Stack * stack)
{
  assert (d->dimension == NULL);
  d->dimension = calloc (1, (nd + 1)*sizeof (Ast *));
  d->dimension[nd] = NULL;
  nd = 0;
  while (ast_schema (direct_declarator->parent, symbol,
		     0, token_symbol ('[')) ||
	 ast_schema (direct_declarator->parent, symbol,
		     1, token_symbol ('['))) {
    Ast * nelem = ast_schema (direct_declarator->parent, symbol,
			      0, token_symbol ('[')) ?
      ast_schema (direct_declarator->parent, symbol,
		  1, sym_assignment_expression) :
      ast_schema (direct_declarator->parent, symbol,
		  2, sym_assignment_expression);
    if (nelem)
      d->dimension[nd] = nelem;
    else if (ast_schema (direct_declarator->parent, symbol,
			 1, token_symbol (']')) ||
	     ast_schema (direct_declarator->parent, symbol,
			 2, token_symbol (']'))) {
      if (nd != 0) {
	AstTerminal * t = ast_left_terminal (direct_declarator->parent);
	fprintf (stderr, "%s:%d: error: only the first dimension of a multidimensional array can be undefined\n",
		 t->file, t->line);
	exit (1);
      }
      d->dimension[0] = (Ast *)&ast_void; // empty dimension
    }
    nd++;
    direct_declarator = direct_declarator->parent;
  }
  return direct_declarator;
}

static
Ast * direct_declarator_type (Ast * direct_declarator, AstDimensions * d, Stack * stack)
{
  if (direct_declarator->sym == sym_enumerator)
    return (Ast *) &ast_enum;
  
  if (direct_declarator->sym == sym_struct_or_union_specifier)
    return direct_declarator;
  
  /**
  Array sizes. */
    
  if (ast_schema (direct_declarator, sym_direct_declarator,
		  0, sym_generic_identifier)) {
    int nd = 0;
    Ast * array = direct_declarator;
    while (ast_schema (array->parent, sym_direct_declarator,
		       1, token_symbol ('['))) {
      nd++;
      array = array->parent;
    }
    if (nd && !(direct_declarator = ast_get_array_dimensions (direct_declarator, sym_direct_declarator, d, nd, stack)))
      return NULL;
  }

  Ast * pointer = ast_schema (direct_declarator->parent, sym_declarator);
  while ((pointer = ast_child (pointer, sym_pointer)))
    d->pointer++;

  if (ast_schema (ast_ancestor (direct_declarator, 2), sym_direct_declarator,
		  1, sym_declarator))
    direct_declarator = ast_ancestor (direct_declarator, 2);

  if (direct_declarator->sym == sym_struct_or_union_specifier)
    return direct_declarator;
  else if (ast_schema (direct_declarator->parent, sym_direct_declarator,
		       1, token_symbol ('('))) {    
    // Function declaration
    if (!d->pointer)
      d->pointer = 1;
    return (Ast *) &ast_function;
  }
  else {
    // Standard declaration
    Ast * specifiers = NULL;
    if (ast_schema (ast_ancestor (direct_declarator, 2), sym_init_declarator,
		    0, sym_declarator))
      specifiers = ast_schema (ast_parent (direct_declarator, sym_declaration), sym_declaration,
			       0, sym_declaration_specifiers);
    else if (ast_schema (ast_ancestor (direct_declarator, 2), sym_struct_declarator,
			 0, sym_declarator))
      specifiers = ast_schema (ast_parent (direct_declarator, sym_struct_declaration),
			       sym_struct_declaration,
			       0, sym_specifier_qualifier_list);
    else if ((specifiers = ast_schema (ast_ancestor (direct_declarator, 2),
				       sym_forin_declaration_statement,
				       2, sym_declaration_specifiers)))
      ;
    else if ((specifiers = ast_schema (ast_ancestor (direct_declarator, 2),
				       sym_parameter_declaration,
				       0, sym_declaration_specifiers)))
      ;
    assert (specifiers);
    Ast * type = ast_find (specifiers, sym_types);
    assert (type);
    return type->child[0];
  }
  ast_print_tree (direct_declarator, stderr, 0, 0, -1);
  assert (false);
  return NULL;
}

Ast * ast_identifier_type (Ast * n, AstDimensions * d, Stack * stack)
{
  Ast * type = direct_declarator_type (ast_ancestor (n, 2), d, stack);
  if (!type)
    return NULL;

  if (type->sym == sym_struct_or_union_specifier &&
      !ast_child (type, sym_struct_declaration_list)) {
    if (ast_child (type, sym_YYerror))
      return NULL;
    type = ast_child (type, sym_generic_identifier)->child[0];
  }
  
  if (n->sym == sym_IDENTIFIER && type->sym == sym_IDENTIFIER &&
      !strcmp (ast_terminal (n)->start, ast_terminal (type)->start))
    // type and identifier are identical e.g. 'typedef struct Foo Foo;'
    return type;

  return ast_base_type (type, d, stack);
}
