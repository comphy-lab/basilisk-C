/**
# The Basilisk C to C99 translator

Uses the [AST](README) library to transform the AST obtained when
parsing code with the [Basilisk C grammar](basilisk.yacc) into an AST
respecting the C99 grammar (with added macros).

## Utility functions */

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "ast.h"
#include "symbols.h"
#include "einstein_sum.h"

/**
By default grammar checks are turned off. */

#if 0
# define CHECK(x, recursive) ast_check_grammar(x, recursive)
#else
# define CHECK(x, recursive) ((void) x)
#endif

Ast * ast_is_typedef (const Ast * identifier)
{
  const Ast * declaration = identifier;
  while (declaration && declaration->sym != sym_declaration)
    declaration = declaration->parent;
  if (declaration)
    return ast_schema (declaration, sym_declaration,
		       0, sym_declaration_specifiers,
		       0, sym_storage_class_specifier,
		       0, sym_TYPEDEF);
  return NULL;
}

Ast * ast_find_function (Ast * n, const char * name)
{
  Ast * found = NULL;
  if (n->sym == sym_function_definition) {
    Ast * identifier = ast_find (n, sym_direct_declarator,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER);
    if (!strcmp (ast_terminal(identifier)->start, name))
      found = n;
  }
  if (n->child)
    for (Ast ** c = n->child; *c && !found; c++)
      found = ast_find_function (*c, name);
  return found;
}

Ast * ast_function_identifier (const Ast * function_definition)
{
  return ast_schema (function_definition, sym_function_definition,
		     0, sym_function_declaration,
		     1, sym_declarator,
		     0, sym_direct_declarator,
		     0, sym_direct_declarator,
		     0, sym_generic_identifier,
		     0, sym_IDENTIFIER);
}

Ast * ast_function_call_identifier (const Ast * n)
{
  return ast_schema (n, sym_function_call,
		     0, sym_postfix_expression,
		     0, sym_primary_expression,
		     0, sym_IDENTIFIER);
}

/**
Appends (block) list `list1` to (block) list `list`. */

Ast * ast_list_append_list (Ast * list, Ast * list1)
{
  assert (list->sym == list1->sym);
  Ast * oldparent = list->parent;
  int index = ast_child_index (list);
  Ast * parent = list1;
  while (parent->child[1])
    parent = parent->child[0];
  Ast * item = parent->child[0];
  ast_new_children (parent, list, item);
  ast_set_child (oldparent, index, list1);
  return list1;
}

/**
Appends `item` to (block) `list`. The list item symbol is `item_sym`. */

Ast * ast_block_list_append (Ast * list, int item_sym, Ast * item)
{
  ast_set_line (item, ast_right_terminal (list));
  Ast * parent = list->parent;
  int index = ast_child_index (list);
  Ast * l = ast_new_children (ast_new (parent, list->sym),
			      list, 
			      ast_attach (ast_new (list, item_sym), item));
  ast_set_child (parent, index, l);
  return l;
}

/**
Appends `item` to (comma-separated) `list`. The list item symbol is
`item_sym`. */

Ast * ast_list_append (Ast * list, int item_sym, Ast * item)
{
  ast_set_line (item, ast_right_terminal (list));
  Ast * parent = list->parent;
  int index = ast_child_index (list);
  Ast * l;
  if (item->sym == item_sym)
    l = ast_new_children (ast_new (parent, list->sym),
			  list, 
			  ast_terminal_new_char (item, ","),
			  item);
  else {
    l =  ast_new_children (ast_new (parent, list->sym),
			   list, 
			   ast_terminal_new_char (item, ","),
			   ast_new (item, item_sym));
    ast_attach (l->child[2], item);
  }
  ast_set_child (parent, index, l);
  return l;
}

/**
Prepends `item` to `list`. The list item symbol is `item_sym`. */

Ast * ast_block_list_prepend (Ast * list, int item_sym, Ast * item)
{
  Ast * r = list;
  while (r->child[0]->sym != item_sym)
    r = r->child[0];
  Ast * l = ast_block_list_append (r, item_sym, item), * tmp = r->child[0];
  ast_set_child (r, 0, l->child[1]);
  ast_set_child (l, 1, tmp);
  return r != list ? list : l;
}

/**
Prepends `item` to (comma-separated) `list`. The list item symbol is
`item_sym`. */

Ast * ast_list_prepend (Ast * list, int item_sym, Ast * item)
{
  Ast * r = list;
  while (r->child[0]->sym != item_sym)
    r = r->child[0];
  Ast * l = ast_list_append (r, item_sym, item), * tmp = r->child[0];
  ast_set_child (r, 0, l->child[2]);
  ast_set_child (l, 2, tmp);
  return r != list ? list : l;
}

/**
Removes `item` from the (comma-separated) `list` and returns the new
list or NULL if the list contains only *item*. */

Ast * ast_list_remove (Ast * list, Ast * item)
{
  Ast * grand_parent = item->parent->parent;
  if (ast_child_index (item) == 0) {
    if (grand_parent->sym == list->sym) {
      ast_replace_child (grand_parent, 0, grand_parent->child[2]);
      ast_destroy (grand_parent->child[1]);
      grand_parent->child[1] = NULL;
    }
    else
      return NULL;
  }
  else {
    Ast * parent = item->parent;
    list = parent->child[0];
    ast_replace_child (grand_parent, ast_child_index (parent), list);
  }
  return list;
}

/**
Removes `item` from the (block) `list` and returns the new
list or NULL if the list contains only *item*. */

Ast * ast_block_list_remove (Ast * list, Ast * item)
{
  Ast * grand_parent = item->parent->parent;
  if (ast_child_index (item) == 0) {
    if (grand_parent->sym == list->sym) {
      ast_replace_child (grand_parent, 0, grand_parent->child[1]);
      grand_parent->child[1] = NULL;
    }
    else
      return NULL;
  }
  else {
    Ast * parent = item->parent;
    list = parent->child[0];
    ast_replace_child (grand_parent, ast_child_index (parent), list);
  }
  return list;
}

/**
Transforms a list of expressions into a list of arguments. */

void ast_argument_list (Ast * expression)
{
  while (expression->sym == sym_expression) {
    int child = expression->child[1] ? 2 : 0;
    expression->sym = sym_argument_expression_list;
    Ast * item = ast_new (expression, sym_argument_expression_list_item);
    ast_new_children (item, expression->child[child]);
    ast_set_child (expression, child, item);
    expression = expression->child[0];
  }
}

/**
Transforms a list of arguments ('argument_expression_list') into a
list of initializers ('initializer_list'). */

Ast * ast_initializer_list (Ast * list)
{
  Ast * start = list;
  while (list->sym == sym_argument_expression_list) {
    list->sym = sym_initializer_list;
    Ast * initializer = list->child[1] ? list->child[2] : list->child[0];
    if (initializer) {
      initializer->sym = sym_initializer;
      Ast * equals = ast_schema (initializer, sym_initializer,
				 0, sym_assignment_expression,
				 1, sym_assignment_operator,
				 0, token_symbol('='));
      if (equals) {
	Ast * name = ast_schema (initializer, sym_initializer,
				 0, sym_assignment_expression,
				 0, sym_unary_expression,
				 0, sym_postfix_expression,
				 0, sym_primary_expression,
				 0, sym_IDENTIFIER);
	if (!name)
	  name = ast_schema (initializer, sym_initializer,
			     0, sym_assignment_expression,
			     0, sym_TYPEDEF_NAME);
	if (name) {
	  Ast * designator = ast_new (initializer, sym_designator);
	  Ast * identifier = ast_new (initializer, sym_generic_identifier);
	  Ast * dot = ast_terminal_new_char (initializer, ".");
	  ast_new_children (designator, dot, identifier);
	  ast_new_children (identifier, name);
	  AstTerminal * left = ast_left_terminal (identifier);
	  ast_terminal (dot)->line = left->line;
	  ast_terminal (dot)->before = left->before; left->before = NULL;      
	  Ast * designator_list =
	    ast_new_children (ast_new (initializer, sym_designator_list),
			      designator);
	  Ast * designation =
	    ast_new_children (ast_new (initializer, sym_designation),
			      designator_list, equals);
	  if (initializer->child[0]->child[2]->sym == sym_assignment_expression)
	    ast_set_child (initializer, 0, initializer->child[0]->child[2]);
	  else {
	    assert (initializer->child[0]->child[2]->sym ==
		    sym_postfix_initializer);
	    initializer = initializer->child[0]->child[2];
	    initializer->sym = sym_initializer;
	  }
	  if (list->child[1])
	    ast_new_children (list,
			      list->child[0], list->child[1], designation,
			      initializer);
	  else
	    ast_new_children (list, designation, initializer);
	}
      }
      else if (ast_schema (initializer, sym_initializer,
			   0, sym_postfix_initializer)) {	
	initializer = initializer->child[0];
	initializer->sym = sym_initializer;
	ast_set_child (list, list->child[1] ? 2 : 0, initializer);
      }
    }
    list = list->child[0];    
  }
  return start;
}

Ast * ast_new_cast_expression (Ast * parent)
{
  return ast_new (parent,
		  sym_assignment_expression,
		  sym_conditional_expression,
		  sym_logical_or_expression,
		  sym_logical_and_expression,
		  sym_inclusive_or_expression,
		  sym_exclusive_or_expression,
		  sym_and_expression,
		  sym_equality_expression,
		  sym_relational_expression,
		  sym_shift_expression,
		  sym_additive_expression,
		  sym_multiplicative_expression,
		  sym_cast_expression);
}

Ast * ast_new_unary_expression (Ast * parent)
{
  return ast_attach (ast_new_cast_expression (parent),
		     ast_new (parent, sym_unary_expression));
}

Ast * ast_is_unary_expression (const Ast * n)
{
  if (!n)
    return NULL;
  int sym[] = {
    sym_assignment_expression,
    sym_conditional_expression,
    sym_logical_or_expression,
    sym_logical_and_expression,
    sym_inclusive_or_expression,
    sym_exclusive_or_expression,
    sym_and_expression,
    sym_equality_expression,
    sym_relational_expression,
    sym_shift_expression,
    sym_additive_expression,
    sym_multiplicative_expression,
    sym_cast_expression,
    sym_unary_expression,
    -1
  }, * i;
  for (i = sym; *i >= 0 && *i != n->sym; i++);
  for (; n != ast_placeholder && *i == n->sym && n->child; i++, n = n->child[0])
    if (n->sym == sym_unary_expression)
      return (Ast *) n;
  return NULL;
}

Ast * ast_new_assignment_function_call (Ast * parent, const char * func)
{
  return ast_attach (ast_new_unary_expression (parent),
		     NN(parent, sym_postfix_expression,
			NN(parent, sym_function_call,
			   NN(parent, sym_postfix_expression,
			      NN(parent, sym_primary_expression,
				 NA(parent, sym_IDENTIFIER, func))),
			   NCA(parent, "("),
			   NCA(parent, ")"))));
}

Ast * ast_new_function_call (Ast * parent, const char * func)
{
  return NN(parent, sym_statement,
	    NN(parent, sym_expression_statement,
	       NN(parent, sym_expression,
		  ast_new_assignment_function_call (parent, func)),
	       NCA(parent, ";")));
}

Ast * ast_is_identifier_expression (const Ast * n)
{
  n = ast_is_unary_expression (n);
  if (n)
    n = ast_schema (n, sym_unary_expression,
		    0, sym_postfix_expression,
		    0, sym_primary_expression,
		    0, sym_IDENTIFIER);
  return (Ast *) n;
}

Ast * ast_is_simple_expression (const Ast * n)
{
  n = ast_schema (ast_is_unary_expression (n), sym_unary_expression,
		  0, sym_postfix_expression,
		  0, sym_primary_expression);
  if (n) {
    n = n->child[0];
    if (n->sym == sym_IDENTIFIER ||
	n->sym == sym_constant ||
	n->sym == sym_string)
      return (Ast *) n;
  }
  return NULL;
}

Ast * ast_is_iteration_statement (const Ast * n)
{
  if (n && (n->sym == sym_iteration_statement ||
	    n->sym == sym_foreach_statement ||
	    n->sym == sym_foreach_inner_statement ||
	    n->sym == sym_forin_declaration_statement ||
	    n->sym == sym_forin_statement))
    return (Ast *) n;
  return NULL;
}

Ast * ast_new_constant (Ast * parent, int symbol, const char * value)
{
  return ast_attach (ast_new_unary_expression (parent),
		     ast_new (parent,
			      sym_postfix_expression,
			      sym_primary_expression,
			      sym_constant),
		     ast_terminal_new (parent, symbol, value));
}

Ast * ast_new_empty_scalar (Ast * parent)
{
  return NN(parent, sym_initializer,
	    NCA(parent, "{"),
	    NN(parent, sym_initializer_list,
	       NN(parent, sym_initializer,
		  ast_attach(ast_new_cast_expression (parent),
			     NN(parent, sym_unary_expression,
				NN(parent, sym_unary_operator,
				   NCA(parent, "-")),
				NN(parent, sym_cast_expression,
				   NN(parent, sym_unary_expression,
				      NN(parent, sym_postfix_expression,
					 NN(parent, sym_primary_expression,
					    NN(parent, sym_constant,
					       NA(parent, sym_I_CONSTANT, "1")))))))))),
	    NCA(parent, "}"));
}

Ast * ast_new_empty_vector (Ast * parent, int dimension)
{
  Ast * list = NN(parent, sym_initializer_list,
		  ast_new_empty_scalar (parent));
  Ast * ret = NN(parent, sym_initializer,
		 NCA(parent, "{"), list, NCA(parent, "}"));
  for (int i = 1; i < dimension; i++)
    ast_list_append (list, sym_initializer, ast_new_empty_scalar (parent));
  return ret;
}

Ast * ast_new_empty_tensor (Ast * parent, int dimension)
{
  Ast * list = NN(parent, sym_initializer_list,
		  ast_new_empty_vector (parent, dimension));
  Ast * ret = NN(parent, sym_initializer,
		 NCA(parent, "{"), list, NCA(parent, "}"));
  for (int i = 1; i < dimension; i++)
    ast_list_append (list, sym_initializer, ast_new_empty_vector (parent, dimension));
  return ret;
}

Ast * ast_new_identifier (Ast * parent, const char * name)
{
  return ast_attach (ast_new (parent,
			      sym_postfix_expression,
			      sym_primary_expression),
		     ast_terminal_new (parent, sym_IDENTIFIER, name)); 
}

Ast * ast_new_member_identifier (Ast * parent, const char * name)
{
  return ast_attach (ast_new (parent,
			      sym_member_identifier,
			      sym_generic_identifier),
		     ast_terminal_new (parent, sym_IDENTIFIER, name));
}

Ast * ast_get_struct_name (Ast * declaration_specifiers)
{
  return ast_schema (declaration_specifiers, sym_declaration_specifiers,
		     0, sym_type_specifier,
		     0, sym_types,
		     0, sym_struct_or_union_specifier,
		     1, sym_generic_identifier,
		     0, sym_IDENTIFIER);
}

static Ast * find_struct_member (Ast * n, const char * member)
{
  if (!n)
    return NULL;
  Ast * identifier = ast_schema (n, sym_struct_declarator,
				 0, sym_declarator,
				 0, sym_direct_declarator,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER);
  if (identifier && !strcmp (ast_terminal (identifier)->start, member))
    return identifier;
  if (n->child)
    for (Ast ** c = n->child; *c; c++) {
      Ast * found = find_struct_member (*c, member);
      if (found)
	return found;
    }
  return NULL;
}

Ast * ast_declaration_from_type (const Ast * type)
{
  if (!type)
    return NULL;
  while (type->sym != sym_declaration &&
	 type->sym != sym_function_declaration &&
	 type->sym != sym_parameter_declaration &&
	 type->sym != sym_struct_declaration &&
	 type->sym != sym_forin_declaration_statement)
    type = type->parent;
  assert (type);
  return (Ast *) type;
}

Ast * ast_expression_type (Ast * expr, Stack * stack, bool higher_dimension)
{
  if (!expr || expr == ast_placeholder)
    return NULL;
  switch (expr->sym) {

  case sym_IDENTIFIER:
    if (ast_ancestor (expr, 2)->sym == sym_member_identifier)
      return ast_expression_type (ast_ancestor (expr, 3), stack,
				  higher_dimension);
    else
      return ast_identifier_declaration (stack, ast_terminal (expr)->start);

  case sym_unary_expression:
  case sym_primary_expression:
  case sym_argument_expression_list_item:
    return ast_expression_type (expr->child[0], stack, higher_dimension);
    
  case sym_initializer:
  case sym_assignment_expression:
    while (expr->child && expr->sym != sym_postfix_expression)
      expr = expr->child[0];
    return expr->sym == sym_postfix_expression ?
      ast_expression_type (expr, stack, higher_dimension) : NULL;
    
  case sym_postfix_expression:
    assert (expr->child && expr->child[0]);
    if (expr->child[1] == NULL || expr->child[2] == NULL)
      return ast_expression_type (expr->child[0], stack, higher_dimension);
    if (expr->child[1]->sym == token_symbol('.') || expr->child[1]->sym == sym_PTR_OP) {
      // struct member access
      Ast * str = ast_expression_type (expr->child[0], stack, higher_dimension);
      if (str) {
	Ast * member = ast_find (expr->child[2], sym_member_identifier,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER);
	Ast * declaration = ast_find (ast_declaration_from_type (str), sym_types);
	assert (declaration);
	AstTerminal * typename = (AstTerminal *)
	  ast_schema (declaration, sym_types,
		      0, sym_TYPEDEF_NAME);
	if (!typename)
	  typename = (AstTerminal *)
	    ast_schema (declaration, sym_types,
			0, sym_struct_or_union_specifier,
			1, sym_generic_identifier,
			0, sym_IDENTIFIER);
	if (typename) {
	  Ast * type = ast_identifier_declaration (stack, typename->start);
	  if (!type) {
#if 0	    
	    fprintf (stderr, "%s:%d: warning: unknown type name '%s'\n",
		     typename->file, typename->line, typename->start);
#endif
	    return NULL;
	  }

	  if (!member)
	    return NULL;
	  
	  /**
	  Special treatment of vector and tensor fields, to deal with
	  possibly undefined components in lower dimensions. */

	  const char * mname = 
	    (higher_dimension &&
	     ast_terminal (member)->start[1] == '\0' &&
	     strchr ("xyz", ast_terminal (member)->start[0]) &&
	     (!strcmp (ast_terminal (type)->start, "vector") ||
	      !strcmp (ast_terminal (type)->start, "tensor"))) ? "x" :
	    ast_terminal (member)->start;

	  if (!strcmp (ast_terminal (type)->start, "scalar"))
	    type = ast_identifier_declaration (stack, "_Attributes");

	  while (type->sym != sym_declaration)
	    type = type->parent;
	  return
	    find_struct_member (ast_find (type, sym_struct_declaration_list),
				mname);
	}
	else if ((str = ast_schema (declaration, sym_types,
				    0, sym_struct_or_union_specifier,
				    2, sym_struct_declaration_list)))
	  return member ? find_struct_member (str, ast_terminal (member)->start) : NULL;
      }
    }
    break;
    
  }  
  return NULL;
}

static char * typedef_name_from_declaration (Ast * declaration)
{
  Ast * types = ast_find (declaration, sym_types), * n;
  if ((n = ast_schema (types, sym_types, 0, sym_TYPEDEF_NAME)))
    return ast_terminal(n)->start;
  return NULL;
}

AstTerminal * ast_type (const Ast * identifier)
{
  if (!identifier)
    return NULL;
  const Ast * declarator = ast_parent (identifier, sym_declarator);
  if (!ast_schema (declarator, sym_declarator,
		   0, sym_direct_declarator,
		   0, sym_generic_identifier,
		   0, sym_IDENTIFIER) &&
      !ast_schema (declarator, sym_declarator,
		   0, sym_direct_declarator,
		   0, sym_direct_declarator,
		   0, sym_generic_identifier,
		   0, sym_IDENTIFIER))
    return NULL; // this is a pointer
  return ast_terminal (ast_find (ast_declaration_from_type (identifier),
				 sym_types)->child[0]);
}

char * ast_typedef_name (const Ast * identifier)
{
  AstTerminal * type = ast_type (identifier);
  if (!type || ((Ast *)type)->sym != sym_TYPEDEF_NAME)
    return NULL;
  return type->start;
}

static Ast * inforeach (Ast * n)
{
  Ast * parent = n->parent;
  while (parent) {
    if (parent->sym == sym_foreach_statement)
      return parent;
    parent = parent->parent;
  }
  return NULL;
}

static bool point_declaration (Stack * stack)
{
  const char * typename =
    ast_typedef_name (ast_identifier_declaration (stack, "point"));
  return typename && !strcmp (typename, "Point");
}

/**
Add arguments ('0') to `function_call` so that the call has exactly
`n` arguments. */

static void complete_arguments (Ast * function_call, int n)
{
  Ast * args = ast_child (function_call, sym_argument_expression_list);
  if (!args) { // function call without arguments
    ast_new_children (function_call,
		      function_call->child[0],
		      function_call->child[1],
		      ast_attach (ast_new (function_call,
					   sym_argument_expression_list,
					   sym_argument_expression_list_item),
				  ast_new_constant (function_call->child[1],
						    sym_I_CONSTANT, "0")),
		      function_call->child[2]);
    args = ast_child (function_call, sym_argument_expression_list);
  }
  
  int i = 0;
  foreach_item (args, 2, item)
    i++;
  for (; i < n; i++) {
    args = ast_list_append (args,
			    sym_argument_expression_list_item,
			    ast_new_constant (function_call->child[3],
					      sym_I_CONSTANT, "0"));
    ast_set_child (function_call, 2, args);
  }
}

static Ast * rotate_arguments (Ast * list, int dimension)
{
  for (int i = 0; i < 3 - dimension; i++) {
    assert (list->child[1]);
    list = list->child[0];
  }
  if (!list->child[1])
    ast_print (list, stderr, 0);
  assert (list->child[1]);
  Ast * next = list->child[0], * item = list->child[2];
  for (int i = 1; i < dimension && next; i++) {
    if (next->child[1]) {
      ast_set_child (list, 2, next->child[2]);
      list = next;
      next = list->child[0];
    }
    else {
      ast_set_child (list, 2, next->child[0]);
      list = next;
      next = NULL;
    }	    
  }
  if (list->child[1])
    ast_set_child (list, 2, item);
  else
    ast_set_child (list, 0, item);
  return list;
}

typedef struct {
  Ast * identifier;
  int type, index, dimension, symmetric;
} Field;

static char * field_value (Field * c, const char * prefix, int type)
{
  bool constant = false;
  int cindex = c->index;
  if (cindex >= 65535)
    cindex -= 65535, constant = true;
  char * src = NULL;
  if (c->type == 3) { // tensor
    int index = cindex, m[c->dimension][c->dimension];    
    for (int j = 0; j < c->dimension; j++) {
      if (type > 1)
	str_append (src, "{");
      for (int i = 0; i < c->dimension; i++) {
	char s[20];
	if (c->symmetric) {
	  m[i][j] = i >= j ? index++ : m[j][i];
	  snprintf (s, 19, "%s%d", constant ? "_NVARMAX+" : "", m[i][j]);
	}
	else
	  snprintf (s, 19, "%s%d", constant ? "_NVARMAX+" : "", index++);
	str_append (src, "{", prefix, s, "}",
		    i < c->dimension - 1 ? "," : "");
      }
      str_append (src, type > 1 ? "}" : "", j < c->dimension - 1 ? "," : "");
    }
  }
  else if (c->type == 2) // vector
    for (int i = 0; i < c->dimension; i++) {
      char s[20];
      snprintf (s, 19, "%s%d", constant ? "_NVARMAX+" : "", cindex + i);
      str_append (src, "{", prefix, s, "}",
		  i < c->dimension - 1 ? "," : "");
    }
  else if (c->type == 1) { // scalar
    char s[30];
    snprintf (s, 29, "%s%d", constant ? "_NVARMAX+" : "", cindex);
    str_append (src, prefix, s);
  }
  if (type >= c->type) {
    str_prepend (src, "{");
    str_append (src, "}");
  }
  return src;
}

static void field_init (Field * c, const char * typename,
			int dimension, int * index)
{
  c->index = *index;
  if (!strcmp (typename, "scalar") ||
      !strcmp (typename, "vertex scalar"))
    c->type = 1, c->dimension = 1, *index += 1;
  else if (!strcmp (typename, "vector") ||
	   !strcmp (typename, "face vector"))
    c->type = 2, c->dimension = dimension, *index += dimension;
  else if (!strcmp (typename, "tensor")) {
    c->type = 3, c->dimension = dimension;
    if (c->symmetric)
      *index += dimension*(dimension + 1)/2;
    else
      *index += dimension*dimension;
  }
  else if (!strcmp (typename, "symmetric tensor"))
    c->type = 3, c->dimension = dimension, c->symmetric = 1,
      *index += dimension*(dimension + 1)/2;
}

static Field * field_append (Field ** fields, Ast * identifier,
			     const char * typename, int dimension, int * index)
{
  int len = 0;
  for (Field * c = *fields; c->identifier; c++, len++);
  *fields = realloc (*fields, (len + 2)*sizeof (Field));
  (*fields)[len + 1] = (Field){0};
  Field * c = &(*fields)[len];
  c->identifier = identifier;
  c->symmetric = 0;
  field_init (c, typename, dimension, index);
  return c;
}

typedef struct {
  int dimension;
  bool nolineno, parallel, cpu;
  Field * constants;
  int constants_index, fields_index, nboundary;
  Ast * init_solver, * init_events, * init_fields, * last_events;
  Ast * boundary;
  char * swigname, * swigdecl, * swiginit;
} TranslateData;

static Ast * in_stencil_point_function (Ast * n)
{
  n = ast_parent (n, sym_function_definition);
  if (!n)
    return NULL;
  if (ast_is_stencil_function (n))
    return n;
  return NULL;
}

static int stencil_access_function (const char * name)
{
  if (!strcmp (name, "val") ||
      !strcmp (name, "val_a") || !strcmp (name, "val_r") || !strcmp (name, "val_o") ||
      !strcmp (name, "fine") ||
      !strcmp (name, "coarse"))
    return 4;
  else if (!strcmp (name, "allocated") ||
	   !strcmp (name, "allocated_child") ||
	   !strcmp (name, "neighbor") ||
	   !strcmp (name, "neighborp") ||
	   !strcmp (name, "aparent") ||
	   !strcmp (name, "aparent_a") || !strcmp (name, "aparent_r") || !strcmp (name, "aparent_o") ||
	   !strcmp (name, "child"))
    return 3;
  return 0;
}

static void rotate (Ast * n, Stack * stack, void * data)
{
  TranslateData * d = data;
  switch (n->sym) {
    
  case sym_IDENTIFIER: case sym_FOREACH: {
    AstTerminal * t = ast_terminal (n);
    int len = strlen (t->start);
    if (len >= 2 && t->start[len - 2] == '_' &&
	strchr ("xyz", t->start[len - 1]))
      t->start[len - 1] = 'x' + (t->start[len - 1] + 1 - 'x') % d->dimension;
    else if (d->dimension > 1) {
      if (!strcmp (t->start, "right"))
	free (t->start), t->start = strdup ("top");
      else if (!strcmp (t->start, "left"))
	free (t->start), t->start = strdup ("bottom");
      else if (!strcmp (t->start, "top"))
	free (t->start), t->start = strdup ("front");
      else if (!strcmp (t->start, "bottom"))
	free (t->start), t->start = strdup ("back");
      else if (!strcmp (t->start, "front"))
	free (t->start), t->start = strdup ("right");
      else if (!strcmp (t->start, "back"))
	free (t->start), t->start = strdup ("left");
    }
    break;
  }

  case sym_member_identifier: {
    AstTerminal * t = ast_terminal (ast_schema (n, sym_member_identifier,
						0, sym_generic_identifier,
						0, sym_IDENTIFIER));
    if (t->start[1] == '\0' && strchr ("xyz", *t->start))
      *t->start = 'x' + (*t->start + 1 - 'x') % d->dimension;
    break;
  }

  case sym_function_call: {
    if (d->dimension > 1) {
      Ast * identifier = ast_function_call_identifier (n);
      if (identifier) {
	const char * name = ast_terminal (identifier)->start;
	if (strcmp (name, "child") && stencil_access_function (name) &&
	    (inforeach (n) || point_declaration (stack) ||
	     in_stencil_point_function (n)))
	  rotate_arguments (n->child[2], d->dimension);
      }
    }
    break;
  }
    
  }
}

static void rotate_list_item (Ast * item, Ast * n,
			      Stack * stack, TranslateData * d)
{
  int dimension = d->dimension;
  if (n->child[4]) {
    d->dimension = atoi (ast_terminal (n->child[2])->start);
    if (d->dimension > dimension)
      d->dimension = dimension;
  }
  
  Ast * list = item->parent;
  Ast * body = ast_last_child (n), * copy = body;
  if (d->dimension == 1) {
    stack_push (stack, &copy);
    ast_traverse (copy, stack, rotate, d);
    ast_pop_scope (stack, copy);    
  }
  else
    for (int i = 1; i < d->dimension; i++) {
      copy = ast_copy (copy);
      stack_push (stack, &copy);
      ast_traverse (copy, stack, rotate, d);
      ast_pop_scope (stack, copy);
      list = ast_block_list_append (list, item->sym, copy);
    }
  ast_set_child (item, 0, body);
  ast_remove (n, ast_left_terminal (body));

  d->dimension = dimension;
}

/**
This function returns a block_item containing *statement*. */

Ast * ast_block_list_get_item (Ast * statement)
{
  assert (statement->sym == sym_statement ||
	  statement->sym == sym_declaration);
  Ast * item = statement->parent;

  /**
  if *item* is not already a block item we need to replace it with a
  compound statement containing a new block_item_list. */
  
  if (item->sym != sym_block_item) {
    AstTerminal * l = ast_left_terminal (statement);
    Ast * left = ast_terminal_new_char ((Ast *) l, "{"),
      * right =
      ast_terminal_new_char ((Ast *) ast_right_terminal (statement), "}");
    ast_terminal (left)->before = l->before, l->before = NULL;
    Ast * parent = item;
    int index = ast_child_index (statement);
    item = ast_new_children (ast_new (parent, sym_block_item), statement);
    Ast * list = ast_new_children (ast_new (parent, sym_block_item_list),
				   item);
    Ast * compound =
      ast_new_children (ast_new (parent, sym_statement),
			ast_new_children (ast_new (parent,
						   sym_compound_statement),
					  left, list, right));
    ast_replace_child (parent, index, compound);
  }
  
  return item;
}

static
void maybeconstfield (Ast * n, Stack * stack,
		      void func (Ast * n, Ast * type, void * data),
		      void * data)
{
  Ast * identifier = ast_schema (n, sym_primary_expression,
				 0, sym_IDENTIFIER);
  if (identifier) {
    Ast * type = ast_identifier_declaration (stack,
					     ast_terminal (identifier)->start);
    if (type) {
      Ast * declaration = type;
      while (declaration &&
	     declaration->sym != sym_declaration &&
	     declaration->sym != sym_parameter_declaration &&
	     declaration->sym != sym_forin_declaration_statement)
	declaration = declaration->parent;
      if (ast_schema (ast_child (declaration, sym_declaration_specifiers),
		      sym_declaration_specifiers,
		      0, sym_type_qualifier,
		      0, sym_MAYBECONST)) {
	if (!strcmp (ast_terminal (identifier)->start, "fs") &&
	    ast_terminal (identifier)->line == 233) {
	  ast_stack_print (stack, stderr);
	  ast_print (identifier, stderr, 0);
	  abort();
	}
	func (n, type, data);
      }
    }
  }
  if (n->child)
    for (Ast ** c = n->child; *c; c++)
      maybeconstfield (*c, stack, func, data);  
}

Ast * ast_is_point_point (const Ast * n)
{
  Ast * identifier = ast_schema (n, sym_IDENTIFIER);
  if (!identifier)
    identifier = ast_schema (n, sym_generic_identifier,
			     0, sym_IDENTIFIER);
  if (identifier && identifier->parent->parent->sym == sym_direct_declarator &&
      !strcmp (ast_terminal (identifier)->start, "point")) {    
    const Ast * decl = n;
    while (decl->sym != sym_declaration &&
	   decl->sym != sym_parameter_declaration)
      decl = decl->parent;
    Ast * type = ast_schema (decl->child[0],
			     sym_declaration_specifiers,
			     0, sym_type_specifier,
			     0, sym_types,
			     0, sym_TYPEDEF_NAME);
    if (type && !strcmp (ast_terminal (type)->start, "Point")) {
      if (decl->sym == sym_declaration)
	return (Ast *) decl;
      else if (decl->sym == sym_parameter_declaration) {
	while (decl->sym != sym_parameter_type_list)
	  decl = decl->parent;
	if ((decl = decl->parent)->sym != sym_direct_declarator ||
	    (decl = decl->parent)->sym != sym_declarator ||
	    (decl = decl->parent)->sym != sym_function_declaration ||
	    (decl = decl->parent)->sym != sym_function_definition)
	  return NULL;
	return ast_last_child (decl)->child[0];
      }
    }
  }
  return NULL;
}

Ast * ast_is_point_function (const Ast * declarator)
{
  Ast * parameters = ast_find (declarator, sym_parameter_type_list,
			       0, sym_parameter_list);
  if (parameters)
    foreach_item (parameters, 2, param) {
      Ast * identifier = ast_find (param, sym_IDENTIFIER);
      if (identifier &&
	  identifier->parent->parent->sym == sym_direct_declarator &&
	  !strcmp (ast_terminal (identifier)->start, "point")) {
	const Ast * decl = identifier;
	while (decl->sym != sym_declaration &&
	       decl->sym != sym_parameter_declaration)
	  decl = decl->parent;
	Ast * type = ast_schema (decl->child[0],
				 sym_declaration_specifiers,
				 0, sym_type_specifier,
				 0, sym_types,
				 0, sym_TYPEDEF_NAME);
	if (type && !strcmp (ast_terminal (type)->start, "Point"))
	  return identifier;
      }
    }
  return NULL;
}

typedef struct {
  void (* func) (Ast * n, Ast * type, void * data);
  void * data;
} ConstData;

static
void maybeconst_traverse (Ast * n, Stack * stack, void * vd)
{
  Ast * identifier = ast_function_call_identifier (n);
  if (identifier) {
    ConstData * d = vd;
    const char * name = ast_terminal (identifier)->start;
    if (!strcmp (name, "val") || !strcmp (name, "fine") || !strcmp (name, "coarse"))
      maybeconstfield (ast_find (n->child[2], sym_argument_expression_list_item),
		       stack, d->func, d->data);
  }
}

static
void maybeconst (Ast * n, Stack * stack,
		 void func (Ast * n, Ast * type, void * data),
		 void * data)
{
  stack_push (stack, &n);
  ConstData d = { func, data };
  ast_traverse (n, stack, maybeconst_traverse, &d);
  ast_pop_scope (stack, n);
}

static
void append_const (Ast * n, Ast * type, void * data)
{
  Ast *** m = data;
  if (!*m) {
    *m = malloc (2*sizeof (Ast *));
    (*m)[0] = type;
    (*m)[1] = NULL;
  }
  else {
    int size = 0;
    Ast ** c;
    for (c = *m; *c && *c != type; c++, size++);
    if (*c != type) {
      *m = realloc (*m, (size + 2)*sizeof (Ast *));
      (*m)[size] = type;
      (*m)[size + 1] = NULL;
    }
  }
}

/**
Replaces child at `index` of `parent` with `replacement` or with a
parent of `replacement` of the same symbol as the child. */

void ast_replace_child_same_symbol (Ast * parent, int index, Ast * replacement)
{
  while (replacement && replacement->sym != parent->child[index]->sym)
    replacement = replacement->parent;
  assert (replacement);
  ast_replace_child (parent, index, replacement);
}

static double sq (double x) { return x*x; }
static double cube (double x) { return x*x*x; }

/**
Evaluates a constant (numerical) expression. Return DBL_MAX if the
expression is not a constant. */

double ast_evaluate_constant_expression (const Ast * n)
{
  if (!n)
    return DBL_MAX;
  
  switch (n->sym) {

  case sym_I_CONSTANT: case sym_F_CONSTANT: case sym_ENUMERATION_CONSTANT:
    return ast_terminal (n)->start ? atof (ast_terminal (n)->start) : 0.;

  case sym_constant:
    return ast_evaluate_constant_expression (n->child[0]);

  case sym_expression:
    return ast_evaluate_constant_expression (ast_child (n, sym_assignment_expression));
    
  case sym_expression_error:
    return ast_evaluate_constant_expression (n->child[0]);
    
  case sym_primary_expression:
    if (n->child[0]->sym == sym_constant)
      return ast_evaluate_constant_expression (n->child[0]);
    else if (n->child[1])
      return ast_evaluate_constant_expression (n->child[1]);
    break;
    
  case sym_assignment_expression:
    if (!n->child[1])
      return ast_evaluate_constant_expression (n->child[0]);
    break;

  case sym_postfix_expression:
    if (n->child[0]->sym == sym_primary_expression ||
	n->child[0]->sym == sym_function_call ||
	n->child[0]->sym == sym_array_access)
      return ast_evaluate_constant_expression (n->child[0]);
    break;
    
  case sym_cast_expression:
    if (n->child[0]->sym == sym_unary_expression)
      return ast_evaluate_constant_expression (n->child[0]);
    break;
    
  case sym_unary_expression:
    if (n->child[0]->sym == sym_postfix_expression)
      return ast_evaluate_constant_expression (n->child[0]);
    if (n->child[0]->sym == sym_unary_operator &&
	strchr ("+-", ast_terminal (n->child[0]->child[0])->start[0])) {
      double v = ast_evaluate_constant_expression (n->child[1]);
      return v < DBL_MAX ? (ast_terminal (n->child[0]->child[0])->start[0] == '+' ? 1. : - 1.)*v : DBL_MAX;
    }
    break;

  case sym_conditional_expression: { // fixme not sure that this is correct
    double cond = ast_evaluate_constant_expression (n->child[0]);
    if (!n->child[1])
      return cond;
    if (cond == DBL_MAX)
      return DBL_MAX;
    if (cond)
      return ast_evaluate_constant_expression (n->child[2]);
    else
      return ast_evaluate_constant_expression (n->child[4]);
  }

  case sym_logical_or_expression: {
    double v = ast_evaluate_constant_expression (n->child[0]);
    if (!n->child[1])
      return v;
    else {
      double v1 = ast_evaluate_constant_expression (n->child[2]);
      if (v < DBL_MAX && v1 < DBL_MAX)
	return v || v1;
    }
    break;
  }
    
  case sym_logical_and_expression: {
    double v = ast_evaluate_constant_expression (n->child[0]);
    if (!n->child[1])
      return v;
    else {
      double v1 = ast_evaluate_constant_expression (n->child[2]);
      if (v < DBL_MAX && v1 < DBL_MAX)
	return v && v1;
    }
    break;
  }
    
  case sym_inclusive_or_expression: {
    double v = ast_evaluate_constant_expression (n->child[0]);
    if (!n->child[1])
      return v;
    else {
      double v1 = ast_evaluate_constant_expression (n->child[2]);
      if (v < DBL_MAX && v1 < DBL_MAX)
	return ((int) v) | ((int) v1);
    }
    break;
  }
    
  case sym_exclusive_or_expression: {
    double v = ast_evaluate_constant_expression (n->child[0]);
    if (!n->child[1])
      return v;
    else {
      double v1 = ast_evaluate_constant_expression (n->child[2]);
      if (v < DBL_MAX && v1 < DBL_MAX)
	return ((int) v) ^ ((int) v1);
    }
    break;
  }
    
  case sym_and_expression: {
    double v = ast_evaluate_constant_expression (n->child[0]);
    if (!n->child[1])
      return v;
    else {
      double v1 = ast_evaluate_constant_expression (n->child[2]);
      if (v < DBL_MAX && v1 < DBL_MAX)
	return ((int) v) & ((int) v1);
    }
    break;
  }

  case sym_equality_expression: {
    double v = ast_evaluate_constant_expression (n->child[0]);
    if (!n->child[1])
      return v;
    else {
      double v1 = ast_evaluate_constant_expression (n->child[2]);
      if (v < DBL_MAX && v1 < DBL_MAX) {
	if (n->child[1]->sym == sym_EQ_OP)
	  return v == v1;
	else
	  return v != v1;
      }
    }
    break;
  }
    
  case sym_relational_expression: {
    double v = ast_evaluate_constant_expression (n->child[0]);
    if (!n->child[1])
      return v;
    else if (v < DBL_MAX) {
      double v1 = ast_evaluate_constant_expression (n->child[2]);
      if (v1 < DBL_MAX) {
	if (n->child[1]->sym == sym_LE_OP)
	  return v <= v1;
	else if (n->child[1]->sym == sym_GE_OP)
	  return v >= v1;
	else if (n->child[1]->sym == token_symbol('>'))
	  return v > v1;
	else if (n->child[1]->sym == token_symbol('<'))
	  return v < v1;
      }
    }
    break;
  }

  case sym_shift_expression: {
    double v = ast_evaluate_constant_expression (n->child[0]);
    if (!n->child[1])
      return v;
    else {
      double v1 = ast_evaluate_constant_expression (n->child[2]);
      if (v < DBL_MAX && v1 < DBL_MAX) {
	if (n->child[1]->sym == sym_LEFT_OP)
	  return ((int) v) << ((int) v1);
	else
	  return ((int) v) >> ((int) v1);
      }
    }
    break;
  }
      
  case sym_additive_expression: {
    double v = ast_evaluate_constant_expression (n->child[0]);
    if (!n->child[1])
      return v;
    else {
      double v1 = ast_evaluate_constant_expression (n->child[2]);
      if (v < DBL_MAX && v1 < DBL_MAX) {
	if (n->child[1]->sym == token_symbol('+'))
	  return v + v1;
	else
	  return v - v1;
      }
    }
    break;
  }
    
  case sym_multiplicative_expression: {
    double v = ast_evaluate_constant_expression (n->child[0]);
    if (!n->child[1])
      return v;
    else {
      double v1 = ast_evaluate_constant_expression (n->child[2]);
      if (v < DBL_MAX && v1 < DBL_MAX) {
	if (n->child[1]->sym == token_symbol('*'))
	  return v*v1;
	if (n->child[1]->sym == token_symbol('/'))
	  return v/v1;
	else
	  return ((int) v) % ((int) v1);
      }
    }
    break;
  }

  case sym_array_access:
    return ast_evaluate_constant_expression (n->child[0]);


  case sym_function_call: {
    Ast * name = ast_schema (n, sym_function_call,
			     0, sym_postfix_expression,
			     0, sym_primary_expression,
			     0, sym_IDENTIFIER);
    if (name && n->child[3]) {
      struct {
	const char * name;
	double (* func) (double);
      } funcs[] = {
	{"fabs", fabs}, {"sqrt", sqrt}, {"exp", exp}, {"log", log}, {"log10", log10},
	{"sin", sin}, {"cos", cos}, {"tan", tan},
	{"asin", asin}, {"acos", acos}, {"atan", atan},
	{"sinh", sinh}, {"cosh", cosh}, {"tanh", tanh},
	{"asinh", sinh}, {"acosh", cosh}, {"atanh", tanh},
	{"sq", sq}, {"cube", cube},
	{NULL}
      }, * i = funcs;
      for (; i->name; i++)
	if (!strcmp (ast_terminal (name)->start, i->name)) {
	  double arg = ast_evaluate_constant_expression (n->child[2]);
	  if (arg == DBL_MAX)
	    return arg;
	  return i->func (arg);
	}
    }    
    break;
  }

  case sym_argument_expression_list:
  case sym_argument_expression_list_item:
    return ast_evaluate_constant_expression (n->child[0]);    
    
  }

  return DBL_MAX;
}

/**
### (const) fields combinations */

typedef struct {
  Ast ** consts;
  int bits;
} ReplaceConst;

static
void replace_const (Ast * n, Ast * type, void * data)
{
  ReplaceConst * r = data;
  int index = 0;
  Ast ** c;
  for (c = r->consts; *c && *c != type; c++, index++);
  assert (*c == type);
  if (r->bits & (1 << index)) {
    Ast * unary = n;
    while (unary->sym != sym_unary_expression)
      unary = unary->parent;
    unary = unary->parent;
    while (unary->sym != sym_unary_expression)
      unary = unary->parent;
    unary = unary->parent;

    Ast * p = unary;
    while (p->sym != sym_expression && !p->child[1])
      p = p->parent;
    if (!p->child[1] && p->parent->sym == sym_expression_statement) {
      /**
      Remove statement with no effect (to avoid compiler warnings) */
      p = p->parent;
      ast_destroy (p->child[0]);
      p->child[0] = p->child[1]; p->child[1] = NULL;
    }
    else {
      str_prepend (ast_terminal (n->child[0])->start, "_const_");
      ast_replace_child_same_symbol (unary, 0, n);
    }
  }
}

static
char * combination_constants (TranslateData * d, Ast ** consts, int bits,
			      char * constants)
{
  int nmaybeconst = 0;
  for (Ast ** c = consts; *c; c++, nmaybeconst++);
  for (int i = 0; i < nmaybeconst; i++)
    if (bits & (1 << i)) {
      const char * name = ast_terminal (consts[i])->start;
      const char * typename = ast_typedef_name (consts[i]);
      if (!strcmp (typename, "vector") ||
	  !strcmp (typename, "face vector")) {
	str_append (constants, "_coord _const_", name, "={_constant[",
		    name,  ".x.i-_NVARMAX]");
	for (int j = 1; j < d->dimension; j++) {
	  char s[] = ".x.i-_NVARMAX]"; s[1] = 'x' + j;
	  str_append (constants, ",_constant[", name, s);
	}
	str_append (constants, "};");
      }
      else
	str_append (constants,
		    "double _const_", name, "=_constant[",
		    name, ".i-_NVARMAX];");
      str_append (constants, "NOT_UNUSED(_const_", name, ");");
    }
  return constants;
}

static void combinations (Ast * n, Stack * stack, TranslateData * d,
			  Ast ** consts,
			  Ast * list, Ast * item, const char * key)
{
  int nmaybeconst = 0;
  for (Ast ** c = consts; *c; c++, nmaybeconst++);
  int n2 = 1 << nmaybeconst;
  char * condition = NULL;
  for (int bits = 0; bits < n2; bits++) {
    if (bits > 0)
      str_append (condition, "else ");
    if (bits == n2 - 1)
      str_append (condition, "{");
    else {
      str_append (condition, "if(");
      for (int i = 0; i < nmaybeconst; i++) {
	const char * name = ast_terminal (consts[i])->start;
	const char * typename = ast_typedef_name (consts[i]);
	str_append (condition,
		    (bits & (1 << i)) ? "" : "!",
		    "is_constant(", name,
		    !strcmp (typename, "vector") ||
		    !strcmp (typename, "face vector") ? ".x" : "",
		    ")");
	if (i < nmaybeconst - 1)
	  str_append (condition, " && ");
      }
      str_append (condition, "){");
    }
    condition = combination_constants (d, consts, bits, condition);
    char index[20];
    snprintf (index, 19, "%d", bits);
    str_append (condition, key, "{_statement", index, "_;}}");
  }
  Ast * conditional = ast_parse_expression (condition, ast_get_root (n));
  free (condition);
  for (int bits = 1; bits < n2; bits++) {
    Ast * copy = ast_copy (n);
    maybeconst (copy, stack, replace_const, &(ReplaceConst){consts, bits});
    char statement[100];
    snprintf (statement, 99, "_statement%d_", bits);
    assert (ast_replace (conditional, statement, copy));
  }
  assert (ast_replace (conditional, "_statement0_", n));
  ast_replace_child (item, 0, ast_new_children (ast_new (list, sym_statement),
						conditional));
}

static int field_list_type (Ast * list, Stack * stack, bool mustbe)
{
  int type = 4; // tensor
  foreach_item (list, 2, expr) {
    const char * typename =
      ast_typedef_name (ast_expression_type (expr, stack, false));
    if (!typename ||
	(strcmp (typename, "scalar") &&
	 strcmp (typename, "vector") &&
	 strcmp (typename, "tensor"))) {
      if (mustbe) {
	AstTerminal * t = ast_left_terminal (expr);
	fprintf (stderr,
		 "%s:%d: error: '%s' is not a scalar, vector or tensor\n",
		 t->file, t->line, ast_str_append (expr, NULL));
	exit (1);
      }
      return -1;
    }
    if (type > 1 && !strcmp (typename, "scalar")) type = 1;
    else if (type > 2 && !strcmp (typename, "vector")) type = 2;
    else if (type > 3 && !strcmp (typename, "tensor")) type = 3;
  }
  return type > 3 ? -1 : type;
}

bool ast_is_field (const char * typename)
{
  return typename && (!strcmp (typename, "scalar") ||
		      !strcmp (typename, "vertex scalar") ||
		      !strcmp (typename, "vector") ||
		      !strcmp (typename, "face vector") ||
		      !strcmp (typename, "tensor") ||
		      !strcmp (typename, "symmetric tensor"));
}

static Ast * declarator_is_allocator (Ast * declarator)
{
  Ast * allocator;
  if ((allocator = ast_schema (declarator, sym_declarator,
			       0, sym_direct_declarator)) &&
      allocator->child[0]->sym == sym_direct_declarator &&
      allocator->child[1]->sym == token_symbol('[') &&
      !allocator->child[3] &&
      (allocator = allocator->child[0]->child[0])->sym
      == sym_generic_identifier)
    return allocator;
  return NULL;
}

static Ast * automatic_argument (const Ast * init_declarator)
{
  Ast
    * initializer = ast_child (init_declarator, sym_initializer),
    * unary = ast_is_unary_expression (ast_child (initializer,
						  sym_assignment_expression)),
    * function_call = ast_schema (unary, sym_unary_expression,
				  0, sym_postfix_expression,
				  0, sym_function_call),
    * function_name = ast_function_call_identifier (function_call),
    * argument = ast_schema (function_call, sym_function_call,
			     2, sym_argument_expression_list,
			     0, sym_argument_expression_list_item,
			     0, sym_assignment_expression);
  if (function_name && argument &&
      !strcmp (ast_terminal (function_name)->start, "automatic"))
    return argument;
  return NULL;
}

static Ast * declarator_is_automatic (const Ast * declarator)
{
  Ast * allocator = ast_schema (declarator, sym_declarator,
				0, sym_direct_declarator,
				0, sym_generic_identifier);
  if (!allocator)
    return NULL;
  if (automatic_argument (declarator->parent))
    return allocator;
  return NULL;
}

static
void foreach_field_allocator (Stack * stack, TranslateData * t, Ast * scope,
			      void func (Stack *, TranslateData *,
					 const char *,
					 Ast *, Ast *, Ast *, void *),
			      void * data)
{
  Ast ** d;
  for (int i = 0; (d = stack_index (stack, i)) && *d != scope; i++)
    if (*d && (*d)->sym == sym_IDENTIFIER) {
      Ast * declarator, * init_declarator, * allocator;
      if (((declarator = ast_ancestor (*d, 4)) &&
	   (init_declarator = declarator->parent)->sym == sym_init_declarator &&
	   (allocator = declarator_is_allocator (declarator))) ||
	  ((declarator = ast_ancestor (*d, 3)) &&
	   (init_declarator = declarator->parent)->sym == sym_init_declarator &&
	   (allocator = declarator_is_automatic (declarator)))) {
	Ast * declaration = ast_declaration_from_type (allocator);
	const char * typename = typedef_name_from_declaration (declaration);
	if (ast_is_field (typename))
	  func (stack, t, typename, init_declarator, declarator, allocator,
		data);
      }
    }
}

static void
field_deallocation (Stack * stack, TranslateData * d,
		    const char * typename,
		    Ast * init_declarator, Ast * declarator, Ast * allocator,
		    void * data)
{
  char ** delete = data;
  Ast * argument = automatic_argument (init_declarator);
  if (argument) {
    char * arg = ast_str_append (argument, NULL);
    if (strchr (typename, ' '))
      typename = strchr (typename, ' ') + 1;
    str_append (delete[1], "if((", arg, ")",
		!strcmp (typename, "scalar") ? ".i" :
		!strcmp (typename, "vector") ? ".x.i" : ".x.x.i",
		"<=0)delete((scalar*){",
		ast_terminal (allocator->child[0])->start,
		"});");
    free (arg);
  }
  else
    str_append (delete[0], ast_terminal (allocator->child[0])->start, ",");
}

static char * delete_fields (char ** delete)
{
  char * fields = delete[0], * automatics = delete[1];
  if (fields || automatics) {
    if (fields) {
      str_prepend (fields, "delete((scalar*){");
      fields[strlen (fields) - 1] = '\0';
      str_append (fields, "});");
    }
    if (automatics)
      str_append (fields, automatics);
    delete[0] = fields;
    return fields;
  }
  return NULL;
}

static void
field_allocation (Stack * stack, TranslateData * d,
		  const char * typename,
		  Ast * init_declarator, Ast * declarator, Ast * allocator,
		  void * data)
{
  field_deallocation (stack, d, typename,
		      init_declarator, declarator, allocator, data);
  
  char * src = strdup (typename);
  for (char * s = src; *s != '\0'; s++)
    if (*s == ' ')
      *s = '_';
  const char * name = ast_terminal (allocator->child[0])->start;
  if (strchr (typename, ' '))
    typename = strchr (typename, ' ') + 1;

  Ast * argument = automatic_argument (init_declarator);
  if (argument) {
    char * arg = ast_str_append (argument, NULL);
    str_prepend (src, typename, " _field_=(", arg, ")",
		 !strcmp (typename, "scalar") ? ".i" :
		 !strcmp (typename, "vector") ? ".x.i" : ".x.x.i",
		 ">0?(", arg, "):new_"); // fixme: should be >= 0
    free (arg);
  }
  else
    str_prepend (src, typename, " _field_=new_");
  str_append (src, "(\"", name, "\");");
	    
  Ast * expr = ast_parse_expression (src, ast_get_root (init_declarator));
  free (src);  
  ast_set_line (expr, ast_right_terminal (declarator));
  declarator = ast_find (expr, sym_init_declarator);
  ast_replace_child (declarator, 0, init_declarator->child[0]);
  ast_replace_child (init_declarator->parent,
		     ast_child_index (init_declarator), declarator);
  ast_destroy (expr);

  /**
  Remove '[]' from declarator if necessary. */

  declarator = declarator->child[0];
  Ast * direct = ast_schema (declarator, sym_declarator,
			     0, sym_direct_declarator,
			     0, sym_direct_declarator);
  if (direct)
    ast_replace_child (declarator, 0, direct);
}

static Ast * compound_jump (Ast * return_statement, Ast * function_definition,
			    const char * expression)
{
  assert (return_statement->sym == sym_jump_statement);
  Ast * ret = ast_child (return_statement, sym_RETURN);  
  if (ret && return_statement->child[2] &&
      !ast_is_simple_expression (return_statement->child[1]->child[0])) {
    // return sthg (complicated);
    char * src = NULL;
    str_append (src, "{int ");
    Ast * pointer = ast_schema (function_definition, sym_function_definition,
				0, sym_function_declaration,
				1, sym_declarator,
				0, sym_pointer);
    if (pointer)
      src = ast_str_append (pointer, src);
    str_append (src, "_ret=val;", expression, "return _ret;}");
    Ast * compound =
      ast_parse_expression (src, ast_get_root (function_definition));
    free (src);
    ast_replace (compound, "val", ast_find (return_statement,
					    sym_assignment_expression));
    if (function_definition->sym == sym_function_definition) {
      Ast * func = ast_find (function_definition, sym_direct_declarator);
      while (func->child[0]->sym == sym_direct_declarator)
	func = func->child[0];
      Ast * declarator = ast_flatten (ast_copy (func, sym_IDENTIFIER),
				      ast_left_terminal (return_statement));
      AstTerminal * t = ast_terminal (ast_find (declarator, sym_IDENTIFIER));
      free (t->start); t->start = strdup ("_ret");
      ast_replace (compound, "_ret", declarator);
      
      Ast * type_specifier =
	ast_flatten (ast_copy (ast_find (function_definition,
					 sym_declaration_specifiers,
					 0, sym_type_specifier)),
		     ast_left_terminal (return_statement));
      ast_replace (compound, "int", type_specifier);
    }
    else
      assert (function_definition->sym == sym_event_definition);
    
    ast_replace_child (return_statement->parent, 0, compound);
    return compound;
  }
  else {
    // return;
    char * src = NULL;
    str_append (src, "{", expression, "return _ret;}");
    Ast * compound =
      ast_parse_expression (src, ast_get_root (function_definition));
    free (src);
    Ast * parent = return_statement->parent;
    ast_replace (compound, "_ret", return_statement);
    ast_replace_child (parent, 0, compound);
    return compound;
  }
  return NULL;
}

/**
### Boundary conditions 

This function replaces neumann/dirichlet(...) with
neumann/dirichlet(0) and returns the number of replacements. */

static int homogeneize (Ast * n)
{
  int nh = 0;
  if (n->sym == sym_function_call) {
    Ast * identifier = ast_function_call_identifier (n);
    if (identifier &&
	(!strcmp (ast_terminal (identifier)->start, "_neumann") ||
	 !strcmp (ast_terminal (identifier)->start, "_dirichlet") ||
	 !strcmp (ast_terminal (identifier)->start, "_dirichlet_face"))) {
      str_append (ast_terminal (identifier)->start, "_homogeneous");
      nh = 1;
    }
  }
  if (n->child)
    for (Ast ** c = n->child; *c; c++)
      nh += homogeneize (*c);
  return nh;
}

static void boundary_expr (Ast * n, Stack * stack, void * data)
{
  switch (n->sym) {

  case sym_postfix_expression: {
    
    /**
    Replaces `.n`, `.t` and `.r` relative vector components with the
    corresponding absolute `.x`, `.y` or `.z` absolute vector
    components. */
    
    if (n->child[1] && n->child[1]->sym == token_symbol('.')) {
      const char * typename =
	ast_typedef_name (ast_expression_type (n->child[0], stack, false));
      if (typename && (!strcmp (typename, "vector") ||
		       !strcmp (typename, "face vector"))) {
	Ast * member = ast_find (n->child[2], sym_member_identifier,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER);
	TranslateData * d = data;
	char * name = ast_terminal(member)->start,
	  * dir = ast_left_terminal (d->boundary->child[2])->start;
	if (!strcmp (dir, "left") || !strcmp (dir, "right")) {
	  if (!strcmp (name, "n"))
	    name[0] = 'x';
	  else if (!strcmp (name, "t"))
	    name[0] = 'y';
	  else if (!strcmp (name, "r"))
	    name[0] = 'z';
	}
	else if (!strcmp (dir, "top") || !strcmp (dir, "bottom")) {
	  if (!strcmp (name, "n"))
	    name[0] = 'y';
	  else if (!strcmp (name, "t"))
	    name[0] = d->dimension > 2 ? 'z' : 'x';
	  else if (!strcmp (name, "r"))
	    name[0] = 'x';
	}
	else if (!strcmp (dir, "front") || !strcmp (dir, "back")) {
	  if (!strcmp (name, "n"))
	    name[0] = 'z';
	  else if (!strcmp (name, "t"))
	    name[0] = 'x';
	  else if (!strcmp (name, "r"))
	    name[0] = 'y';
	}
      }
    }

    /**
    Replaces a boundary field with its local value `_s`. */
    
    TranslateData * d = data;
    if (ast_are_identical (n, d->boundary->child[0]))
      ast_replace_child (n->parent, ast_child_index (n),
			 ast_new_identifier (d->boundary, "_s"));
    
    break;
  }
    
  /**
  Replaces `ghost` with the corresponding indices. */

  case sym_array_access: {
    Ast * identifier;
    if (n->child[3] &&
	(identifier = ast_is_identifier_expression (n->child[2]->child[0])) &&
	!strcmp (ast_terminal(identifier)->start, "ghost")) {
      TranslateData * d = data;
      char * dir = ast_left_terminal (d->boundary->child[2])->start,
	* index = (!strcmp (dir, "left") ? "a[-1,0,0];" :
		   !strcmp (dir, "right") ? "a[1,0,0];" :
		   !strcmp (dir, "bottom") ? "a[0,-1,0];" :
		   !strcmp (dir, "top") ? "a[0,1,0];" :
		   !strcmp (dir, "back") ? "a[0,0,-1];" :
		   !strcmp (dir, "front") ? "a[0,0,1];" : NULL);
      assert (index);
      Ast * expr = ast_parse_expression (index, ast_get_root (d->boundary));
      ast_replace_child (n, 2, ast_find (expr, sym_array_access,
					 2, sym_expression));
      ast_destroy (expr);
    }
    break;
  }

  /**
  Dirichlet boundary conditions for normal components of face fields. */

  case sym_function_call: {
    TranslateData * d = data;
    Ast * member = ast_schema (d->boundary->child[0], sym_postfix_expression,
			       2, sym_member_identifier,
			       0, sym_generic_identifier,
			       0, sym_IDENTIFIER);
    if (member && !strcmp (ast_terminal(member)->start, "x")) {
      Ast * identifier = ast_function_call_identifier (n);
      if (identifier &&
	  !strcmp (ast_terminal (identifier)->start, "_dirichlet")) {
	const char * typename =
	  ast_typedef_name (ast_expression_type
			    (d->boundary->child[0]->child[0],
			     stack, false));
	if (!strcmp (typename, "face vector"))
	  str_append (ast_terminal (identifier)->start, "_face");
      }
    }
    break;
  }
    
  }
}

#define foreach_map(map)					\
  Ast ** _m, * map;						\
  for (int _i = 0; (_m = stack_index (stack, _i)); _i++)	\
    if ((map = ast_schema (*_m, sym_macro_statement,		\
			   1, sym_compound_statement,		\
			   1, sym_block_item_list)))

static Ast * boundary_function (Ast * expr, Stack * stack, TranslateData * d,
				char * before, char * ind)
{
  char * src = NULL;
  snprintf (ind, 19, "%d", d->nboundary++);
  str_append (src,
	      "static double _boundary", ind,
	      "(Point point,Point neighbor,scalar _s,void *data){{"); // The double brackets are important
      
  char * index[] = {"i","j","k"}, * dir[] = {"x","y","z"};
  for (int i = 0; i < d->dimension; i++)
    str_append (src, "int ",
		index[i], "g=neighbor.", index[i], "-point.", index[i], ";"
		"if(", index[i], "g==0)", index[i], "g=_attribute[_s.i].d.",
		dir[i], ";",
		"NOT_UNUSED(", index[i], "g);");
  assert (before);
  str_append (src, "POINT_VARIABLES;");
  str_append (src, "{return(", before, "_expr_);}}}");
  free (before);
  Ast * boundary =
    ast_child (ast_parse_external_declaration (src, ast_get_root (expr)),
	       sym_function_definition);
  ast_get_root (boundary)->alloc = ast_get_root (expr)->alloc;
  free (src);
  assert (expr->sym == sym_assignment_expression);
  ast_replace (boundary, "_expr_", expr);
  ast_set_line (boundary, ast_left_terminal (expr));
  Ast * compound = ast_find (ast_find (boundary, sym_compound_statement)->child[1], sym_compound_statement);
  Ast * list = ast_schema (compound, sym_compound_statement,
			   1, sym_block_item_list);
  foreach_map (m)
    foreach_item (m, 1, item)
      ast_block_list_prepend (list, sym_block_item, ast_copy (item->child[0]));
  stack_push (stack, &expr);
  ast_traverse (expr, stack, boundary_expr, d);
  ast_pop_scope (stack, expr);
  return boundary;
}

static void set_boundary_component (Ast * member_identifier)
{
  Ast * member = ast_schema (member_identifier, sym_member_identifier,
			     0, sym_generic_identifier,
			     0, sym_IDENTIFIER);
  if (member) {
    if (!strcmp (ast_terminal(member)->start, "n"))
      ast_terminal(member)->start[0] = 'x';
    else if (!strcmp (ast_terminal(member)->start, "t"))
      ast_terminal(member)->start[0] = 'y';
    else if (!strcmp (ast_terminal(member)->start, "r"))
      ast_terminal(member)->start[0] = 'z';
  }  
}

static char * set_boundary (Ast * array, char * ind)
{
  assert (array->sym == sym_array_access);
  char * bc = ast_str_append (array->child[2], NULL);
  char * scalar = ast_str_append (array->child[0], NULL);
  char * set = NULL;
  str_append (set,
	      "_attribute[", scalar, ".i].dirty=1,",
	      "_attribute[", scalar, ".i].boundary[", bc,
	      "]=_boundary", ind, ",",
	      "_attribute[", scalar, ".i].boundary_homogeneous[", bc,
	      "]=_boundary", ind);
  free (scalar);
  free (bc);
  return set;
}

static Ast * function_scope (Ast * n, Stack * stack)
{
  if (point_declaration (stack))
    return NULL;
  while (n) {
    if (n->sym == sym_foreach_statement)
      return NULL;
    if (n->sym == sym_function_definition ||
	n->sym == sym_event_definition)
      return n;
    n = n->parent;
  }
  return NULL;
}

/**
Inserts `item` after `insert` in the list containing `insert`. */

Ast * ast_list_insert_after (Ast * insert, Ast * item)
{
  Ast * list_item = insert->parent, * list = list_item->parent,
    * parent = list->parent;
  int item_sym = list_item->sym;
  int index = ast_child_index (list);
  Ast * nlist = NN(list, list->sym,
		   list,
		   NCB (list, ","),
		   NN (list, item_sym,
		       item));
  if (parent->sym != list->sym)
    ast_set_child (parent, index, nlist);
  else 
    ast_set_child (parent, 0, nlist);
  return list;
}

/**
Inserts `item` after `insert` in the (block) list containing `insert`. */

Ast * ast_block_list_insert_after (Ast * insert, Ast * item)
{
  Ast * list_item = insert->parent, * list = list_item->parent,
    * parent = list->parent;
  int item_sym = list_item->sym;
  assert (parent->sym == list->sym);	
  ast_set_child (parent, 0,
		 ast_new_children (ast_new (list, list->sym),
				   list,
				   ast_new_children (ast_new (list, item_sym),
						     item)));
  return list;
}

/**
Inserts `item` before `insert` in the (block) list containing `insert`. */

Ast * ast_block_list_insert_before (Ast * insert, Ast * item)
{
  return ast_block_list_insert_after
    (insert->parent->parent->child[0]->child[1]->child[0], item);
}

Ast * ast_block_list_insert_before2 (Ast * insert, Ast * item)
{
  // fixme: merge with above
  Ast * parent = insert->child[0];
  Ast * list = ast_block_list_append (insert->parent, insert->sym, item);
  ast_set_child (insert, 0, list->child[1]->child[0]);
  ast_set_child (list->child[1], 0, parent);
  return list;
}

const char * ast_crop_before (const char * s)
{
  while (strchr (" \t\n\r", *s)) s++;
  while (*s == '#' || *s == '@') {
    s++;
    while (*s != '\0' && *s != '\n') s++;
    while (strchr (" \t\n\r", *s)) s++;
  }
  return s;
}

static
void compound_init (Ast * compound, Ast * n)
{
  Ast * list = NN(compound, sym_block_item_list,
		  NN(compound, sym_block_item,
		     n));
  ast_new_children (compound, compound->child[0], list, compound->child[1]);
}

static
void compound_append (Ast * compound, Ast * n)
{
  Ast * list = ast_schema (compound, sym_compound_statement,
			   1, sym_block_item_list);
  if (!list)
    compound_init (compound, n);
  else
    ast_block_list_append (list, sym_block_item, n);
}
  
static
void compound_prepend (Ast * compound, Ast * n)
{
  Ast * list = ast_schema (compound, sym_compound_statement,
			   1, sym_block_item_list);
  if (!list)
    compound_init (compound, n);
  else
    ast_block_list_prepend (list, sym_block_item, n);
}
  
/**
# First pass: Global boundaries and stencils */

static void global_boundaries_and_stencils (Ast * n, Stack * stack, void * data)
{
  switch (n->sym) {

  /**
  ## Warnings for Basilisk C parse errors */
    
  case sym_YYerror: {
    AstTerminal * t = ast_left_terminal (n);
    char * s = NULL;
    s = ast_str_append (n, s);
    fprintf (stderr, "%s:%d: warning: Basilisk C parse error near `%s'\n",
	     t->file, t->line, ast_crop_before (s));
    free (s);
    break;
  }

  /**
  ## Local boundary conditions */
    
  case sym_array_access: {

    Ast * assign = ast_ancestor (n, 3), * scope;
    if (assign->sym == sym_assignment_expression &&
	(scope = function_scope (n, stack))) {
      const char * typename =
	ast_typedef_name (ast_expression_type (n->child[0], stack, false));
      Ast * member = NULL;
      if ((typename &&
	   (!strcmp (typename, "scalar") ||
	    !strcmp (typename, "vertex scalar"))) ||
	  ((member = ast_schema (n->child[0], sym_postfix_expression,
				 2, sym_member_identifier,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER)) &&
	   (!strcmp (ast_terminal (member)->start, "n") ||
	    !strcmp (ast_terminal (member)->start, "t") ||
	    !strcmp (ast_terminal (member)->start, "r")) &&
	   (typename =
	    ast_typedef_name (ast_expression_type (n->child[0]->child[0],
						   stack, false))) &&
	   (!strcmp (typename, "vector") ||
	    !strcmp (typename, "face vector")))) {
	AstTerminal * t = ast_left_terminal (assign);
	char * before = t->before;
	t->before = NULL;
	char ind[20];
	TranslateData * d = data;
	d->boundary = n;
	Ast * boundary =
	  boundary_function (ast_child (assign, sym_assignment_expression),
			     stack, data, before, ind);
	char * set = set_boundary (n, ind);
	ast_block_list_insert_before (scope, boundary);
	
	Ast * homogeneous = ast_copy (boundary);
	if (!homogeneize (homogeneous)) {
	  ast_destroy (homogeneous);
	  str_append (set, ";\n");
	}
	else {
	  Ast * func = ast_find (homogeneous, sym_IDENTIFIER);
	  str_append (ast_terminal (func)->start, "_homogeneous");
	  str_append (set, "_homogeneous;\n");
	  ast_block_list_insert_before (scope, homogeneous);
	}	
	Ast * expr = ast_parse_expression (set, ast_get_root (n));
	free (set);
	Ast * parent = ast_ancestor (assign, 2);
	assert (parent->sym == sym_expression_statement);
	ast_replace_child (parent, 0, ast_child (expr, sym_expression));
	ast_destroy (expr);
      }
    }
    break;
  }

  /**
  ## Global boundary conditions */
    
  case sym_boundary_definition: {
    Ast * expr = ast_schema (n, sym_boundary_definition,
			     0, sym_assignment_expression,
			     2, sym_assignment_expression);
    Ast * array = ast_find (n, sym_array_access);
    if (expr && array) {
      AstTerminal * t = ast_left_terminal (n);
      char * before = t->before;
      t->before = NULL;
      set_boundary_component (ast_schema (array->child[0],
					  sym_postfix_expression,
					  2, sym_member_identifier));
      char ind[20];
      TranslateData * d = data;
      d->boundary = array;      
      Ast * boundary = boundary_function (expr, stack, data, before, ind);
      char * set = set_boundary (array, ind);
      ast_replace_child (n->parent, 0, boundary);
      
      Ast * homogeneous = ast_copy (boundary);
      if (!homogeneize (homogeneous)) {
	ast_destroy (homogeneous);
	str_append (set, ";\n");
      }
      else {
	Ast * func = ast_find (homogeneous, sym_IDENTIFIER);
	str_append (ast_terminal (func)->start, "_homogeneous");
	str_append (set, "_homogeneous;\n");
	ast_block_list_insert_after (n, homogeneous);
      }
      Ast * expr = ast_parse_expression (set, ast_get_root (n));
      compound_append (d->last_events, NN(n, sym_statement, expr));
      free (set);
    }
    break;
  }

  /**
  ## Stencils */

  case sym_foreach_statement: {
    if (!strcmp (ast_terminal (n->child[0])->start, "foreach") ||
	!strcmp (ast_terminal (n->child[0])->start, "foreach_vertex") ||
	!strcmp (ast_terminal (n->child[0])->start, "foreach_face") ||
	!strcmp (ast_terminal (n->child[0])->start, "foreach_visible") ||
	!strcmp (ast_terminal (n->child[0])->start, "foreach_point") ||
	!strcmp (ast_terminal (n->child[0])->start, "foreach_region")) {
      bool overflow = false, nowarning = false, gpu = false;
      Ast * parameters = ast_child (n, sym_foreach_parameters);
      foreach_item (parameters, 2, item) {
	Ast * identifier = ast_is_identifier_expression (item->child[0]);
	bool noauto;
	if (identifier) {
	  const char * start = ast_terminal (identifier)->start;
	  if (!strcmp (start, "gpu"))
	    gpu = true;
	  else if ((noauto = !strcmp (start, "noauto")) ||
		   !strcmp (start, "overflow") ||
		   !strcmp (start, "nowarning")) {
	    if (!strcmp (start, "overflow"))
	      overflow = true;
	    else if (!strcmp (start, "nowarning"))
	      nowarning = true;
	    parameters = ast_list_remove (parameters, item);
	    if (parameters == NULL) {
	      ast_destroy (n->child[2]);
	      for (Ast ** c = n->child + 2; *c; c++)
		*c = *(c + 1);
	    }
	    if (noauto)
	      return;
	  }
	}
      }
      TranslateData * d = data;
      bool parallel = d->parallel &&
	strcmp (ast_terminal (n->child[0])->start, "foreach_visible");
      Ast * stencil = ast_copy (n);
      if (!ast_stencil (stencil, parallel, overflow, nowarning)) {
	ast_destroy (stencil);
	if (!gpu)
	  break;
	else
	  stencil = NN(n, sym_foreach_statement,
		       NB(n, sym_FOREACH, ast_terminal (n->child[0])->start),
		       NCB(n, "("),
		       NCB(n, ")"),
		       NN(n, sym_statement,
			  NN(n, sym_expression_statement,
			     NCB(n, ";"))));
      }
      str_append (ast_terminal (ast_child (stencil, sym_FOREACH))->start, "_stencil");
      if (n->child[4])
	ast_new_children (n, n->child[0], n->child[1], n->child[2],
			  n->child[3], n->child[4], stencil);
      else
	ast_new_children (n, n->child[0], n->child[1], n->child[2], n->child[3], stencil);
    }
    break;
  }

  /**
  ## Einstein summation */

  case sym_macro_statement: {
    Ast * identifier = ast_schema (n, sym_macro_statement,
				   0, sym_function_call,
				   0, sym_postfix_expression,
				   0, sym_primary_expression,
				   0, sym_IDENTIFIER);
    if (!strcmp (ast_terminal (identifier)->start, "einstein_sum"))
      einstein_sum_global (n, stack, ((TranslateData *)data)->dimension);
  }

  }
}

/**
# Second pass: Most transformations */

static void diagonalize (Ast * n, Stack * stack, void * field)
{
  if (n->sym == sym_function_call) {
    Ast * identifier = ast_function_call_identifier (n);
    if (identifier) {
      Ast * arg;
      if (!strcmp (ast_terminal (identifier)->start, "val") &&
	  (inforeach (n) || point_declaration (stack)) &&
	  (arg = ast_is_identifier_expression
	   (ast_find (n, sym_assignment_expression))) &&
	  !strcmp (ast_terminal (arg)->start,
		   ast_terminal ((Ast *)field)->start))
	str_append (ast_terminal (identifier)->start, "_diagonal");
    }
  }
}

bool ast_is_foreach_stencil (const Ast * n)
{
  Ast * foreach = ast_schema (n, sym_foreach_statement,
			      0, sym_FOREACH);
  if (!foreach)
    return false;
  int len = strlen (ast_terminal (foreach)->start) - 8;
  return len > 0 && !strcmp (ast_terminal (foreach)->start + len, "_stencil");
}

static Ast * higher_dimension (Ast * n)
{
  char * s = in_stencil_point_function (n) ? "_stencil_val_higher_dimension" : "_val_higher_dimension";
  return ast_attach (ast_new (n, sym_primary_expression),
		     ast_terminal_new (n, sym_IDENTIFIER, s));
}

static char * append_initializer (char * init, Ast * initializer, const char * typename)
{
  if (!strcmp (typename, "vector")) {
    str_append (init, "(double[])");
    Ast * list = ast_schema (initializer, sym_initializer,
			     1, sym_initializer_list);
    if (list) {
      char * initialize = ast_str_append (list, NULL);
      str_append (init, "{", initialize);
      free (initialize);
      int nr = 3;
      foreach_item (list, 2, item) nr--;
      while (nr--)
	str_append (init, ",0.");
      str_append (init, "}");
      return init;
    }
  }
  char * initialize = ast_str_append (initializer, NULL);
  str_append (init, initialize);
  free (initialize);
  return init;
}

/**
## Attribute declaration */

static void attribute (Ast * n, Stack * stack, void * data)
{
  if (n->sym != sym_attribute)
    return;
  Ast * identifier = ast_schema (n, sym_attribute,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER);
  if (identifier &&
      !strcmp (ast_terminal (identifier)->start, "attribute")) {

    /**
    Remove 'attribute' from external declarations. */

    Ast * translation = n->parent->parent;
    assert (translation->child[1]);
    Ast * next = translation->child[0];
    assert (translation->parent->child[1]);
    ast_set_child (translation->parent, 0, next);
    assert (next->child[1]);
    str_prepend (ast_left_terminal (translation->parent->child[1])->before,
		 ast_left_terminal (n)->before);

    /**
    Add attributes to typedef '_Attributes'. */
      
    Ast * attr = ast_identifier_declaration (stack, "_Attributes");
    while (attr->sym != sym_declaration)
      attr = attr->parent;
    ast_list_append_list (ast_find (attr, sym_struct_declaration_list),
			  n->child[2]);

    /**
    Cleanup. */
      
    ast_destroy (translation);
  }
}

static Ast * abstract_declarator_from_declarator (Ast * n)
{
  switch (n->sym) {
  case sym_generic_identifier: {
    ast_destroy (n);
    return NULL;
  }
  case sym_declarator: n->sym = sym_abstract_declarator; break;
  case sym_direct_declarator: n->sym = sym_direct_abstract_declarator; break;
  }
  if (!n->child)
    return n;
  for (Ast ** c = n->child; *c; c++)
    if (!abstract_declarator_from_declarator (*c)) {
      for (Ast ** d = c; *d; d++)
	*d = *(d + 1);
      c--;
    }
  if (!*n->child) {
    ast_destroy (n);
    return NULL;
  }
  return n;
}

static Ast * obsolete_function_declaration (const Ast * type)
{
  Ast * parameters = ast_find (type, sym_parameter_list);
	
  /**
  Obsolete optional arguments syntax using 'struct ...' parameters. */

  Ast * struct_name;
  if (parameters && !parameters->child[1] &&
      (struct_name = ast_get_struct_name (ast_schema (parameters, sym_parameter_list,
						      0, sym_parameter_declaration,
						      0, sym_declaration_specifiers))) &&
      ast_schema (parameters, sym_parameter_list,
		  0, sym_parameter_declaration,
		  1, sym_declarator,
		  0, sym_direct_declarator,
		  0, sym_generic_identifier,
		  0, sym_IDENTIFIER))
    return struct_name;
  else
    return NULL;
}

Ast * ast_constant_postfix_expression (const Ast * n, Stack * stack)
{
  Ast * identifier = ast_schema (n, sym_postfix_expression,
				 0, sym_primary_expression,
				 0, sym_IDENTIFIER);
  if (!identifier)
    identifier = ast_schema (n, sym_postfix_expression,
			     0, sym_postfix_expression,
			     0, sym_primary_expression,
			     0, sym_IDENTIFIER);
  if (!identifier)
    identifier = ast_schema (n, sym_postfix_expression,
			     0, sym_postfix_expression,
			     0, sym_postfix_expression,
			     0, sym_primary_expression,
			     0, sym_IDENTIFIER);
  if (identifier) {
    Ast * type = ast_identifier_declaration (stack, ast_terminal (identifier)->start);
    while (type && type->sym != sym_declaration)
      type = type->parent;
    if (type)
      return ast_schema (type->child[0], sym_declaration_specifiers,
			 0, sym_type_qualifier,
			 0, sym_CONST);
  }
  return NULL;
}

/**
## Stencil access 

This function transforms stencil accesses of the form `s[i,j]` into the
function call `val(s,i,j,0)`. */

void ast_stencil_access (Ast * n, Stack * stack, int dimension)
{
  const char * typename =
    ast_typedef_name (ast_expression_type (n->child[0], stack, false));
  Ast * member, * foreach = NULL;
  if (typename &&
      (!strcmp (typename, "scalar") ||
       !strcmp (typename, "vertex scalar")) &&
      ((foreach = inforeach (n)) || point_declaration (stack))) {
    n->sym = sym_function_call;
    ast_set_char (ast_child (n, token_symbol('[')), '(');

    Ast * list = ast_child (n, sym_expression);
    if (list)
      ast_argument_list (list);
    complete_arguments (n, 3);
    list = ast_child (n, sym_argument_expression_list);
    char * before = ast_left_terminal (n)->before;
    ast_left_terminal (n)->before = NULL;
    Ast * func;
    if (ast_constant_postfix_expression (n->child[0], stack))
      func = ast_new_identifier (n, "_val_constant");
    else
      func = ast_new_identifier (n, "val");
    ast_set_child (n, 2,
		   ast_list_prepend (list,
				     sym_argument_expression_list_item,
				     ast_attach (ast_new_unary_expression (n),
						 n->child[0])));
    ast_set_char (n->child[3], ')');
    ast_set_child (n, 0, func);
    ast_left_terminal (n)->before = before;
  }

  /**
  Check whether we are trying to access (undeclared) 'y' or 'z'
  members of a vector or tensor field (i.e. higher dimension members). */
    
  else if ((member = ast_schema (n->child[0], sym_postfix_expression,
				 2, sym_member_identifier,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER)) &&
	   ((dimension < 2 &&
	     (!strcmp (ast_terminal (member)->start, "y") ||
	      !strcmp (ast_terminal (member)->start, "t"))) ||
	    (dimension < 3 &&
	     (!strcmp (ast_terminal (member)->start, "z") ||
	      !strcmp (ast_terminal (member)->start, "r")))) &&
	   (typename =
	    ast_typedef_name (ast_expression_type (n->child[0]->child[0],
						   stack, false))) &&
	   (!strcmp (typename, "vector") ||
	    !strcmp (typename, "face vector")) &&
	   ((foreach = inforeach (n)) || point_declaration (stack)))
    ast_replace_child (n->parent, 0, higher_dimension (n));
  else if ((member = ast_schema (n->child[0], sym_postfix_expression,
				 0, sym_postfix_expression,
				 2, sym_member_identifier,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER)) &&
	   ((dimension < 2 &&
	     !strcmp (ast_terminal (member)->start, "y")) ||
	    (dimension < 3 &&
	     !strcmp (ast_terminal (member)->start, "z"))) &&
	   (typename =
	    ast_typedef_name (ast_expression_type
			      (n->child[0]->child[0]->child[0],
			       stack, false))) &&
	   !strcmp (typename, "tensor") &&
	   (inforeach (n) || point_declaration (stack)))
    ast_replace_child (n->parent, 0, higher_dimension (n));
}

static void translate (Ast * n, Stack * stack, void * data)
{
  typedef struct {
    char * target, * replacement;
  } Replacement;
  
  switch (n->sym) {

  /**
  ## foreach_dimension() */

  case sym_foreach_dimension_statement: {
    Ast * item = ast_block_list_get_item (n->parent->parent);
    rotate_list_item (item, n, stack, data);
    break;
  }

  /**
  ## External foreach_dimension() */

  case sym_external_foreach_dimension: {
    rotate_list_item (n->parent, n, stack, data);
    break;
  }
    
  /**
  ## Diagonalize */

  case sym_macro_statement: {
    Ast * identifier = ast_schema (n, sym_macro_statement,
				   0, sym_function_call,
				   0, sym_postfix_expression,
				   0, sym_primary_expression,
				   0, sym_IDENTIFIER);
    if (!strcmp (ast_terminal (identifier)->start, "diagonalize")) {
      Ast * field = ast_schema (n, sym_macro_statement,
				0, sym_function_call,
				2, sym_argument_expression_list,
				0, sym_argument_expression_list_item,
				0, sym_assignment_expression);
      if (field && (field = ast_is_identifier_expression (field))) {
	stack_push (stack, &n);
	ast_traverse (n, stack, diagonalize, field);
	ast_pop_scope (stack, n);
      }
    }
    break; 
  }
    
  /**
  ## Foreach statements */

  case sym_foreach_statement: {

    /**
    ### foreach_face() statements */

    bool is_face_stencil = !strcmp (ast_terminal (n->child[0])->start,
				    "foreach_face_stencil");
    if (is_face_stencil ||
	!strcmp (ast_terminal (n->child[0])->start, "foreach_face")) {
      char order[] = "xyz";

      /**
      The complicated stuff below is just to read each (optional) x, y
      and z arguments, in the correct order, and update the *order*
      string. */

      Ast * parameters = ast_schema (n, sym_foreach_statement, 2, sym_foreach_parameters);
      if (parameters) {
	char * s = order + 2;
	foreach_item (parameters, 2, param)
	  if (param->child[0]->sym == sym_assignment_expression) {
	    Ast * identifier = ast_find (param, sym_postfix_expression,
					 0, sym_primary_expression,
					 0, sym_IDENTIFIER);
	    if (identifier && ast_terminal (identifier)->start[1] == '\0' &&
		strchr ("xyz", ast_terminal (identifier)->start[0])) {
	      *s-- = ast_terminal (identifier)->start[0];
	      parameters = ast_list_remove (parameters, param);
	    }
	  }
	if (s != order + 2 && s >= order)
	  memmove (order, s + 1, strlen(s));
	if (parameters == NULL) {
	  ast_destroy (n->child[2]);
	  for (Ast ** c = n->child + 2; *c; c++)
	    *c = *(c + 1);
	}
      }

      /**
      Here we add the `is_face_x()` condition to the loop statement. */

      Ast * expr = ast_parse_expression
	(is_face_stencil ? "_stencil_is_face_x(){;}" : "is_face_x(){;}",
	 ast_get_root (n));
      Ast * cond = ast_find (expr, sym_IDENTIFIER);
      ast_terminal (cond)->start[strlen(ast_terminal (cond)->start) - 1] = order[0];
      Ast * statement = ast_child (n, sym_statement);
      int index = ast_child_index (statement);
      ast_replace (expr, ";", statement);
      ast_set_line (expr, ast_left_terminal (n));
      ast_replace_child (n, index, ast_new_children (ast_new (n, sym_statement), expr));
	
      /**
      Finally, we "dimension-rotate" the statement. */

      if (strlen (order) > 1) {
	Ast * statement = ast_child (n, sym_statement);
	Ast * item = ast_block_list_get_item (statement);
	TranslateData * d = data;
	int dimension = d->dimension;
	d->dimension = strlen (order);
	if (d->dimension > dimension) d->dimension = dimension;
	
	Ast * list = item->parent, * copy = statement;
	for (int i = 1; i < d->dimension; i++) {
	  copy = ast_copy (copy);
	  stack_push (stack, &copy);
	  ast_traverse (copy, stack, rotate, d);
	  ast_pop_scope (stack, copy);
	  Ast * cond = ast_find (copy, sym_IDENTIFIER);
	  ast_terminal (cond)->start[strlen(ast_terminal (cond)->start) - 1] =
	    order[i];
	  list = ast_block_list_append (list, item->sym, copy);
	}
	if (statement->sym != sym_statement)
	  statement = ast_new_children (ast_new (n, sym_statement), statement);
	ast_set_child (item, 0, statement);

	d->dimension = dimension;
      }
    }

    /**
    ### (const) fields combinations (except for stencils) */

    if (!ast_is_foreach_stencil (n)) {
      Ast ** consts = NULL;
      maybeconst (n, stack, append_const, &consts);
      if (consts) {
	Ast * item = ast_block_list_get_item (n->parent->parent);
	Ast * list = item->parent;
	combinations (n, stack, data, consts, list, item, "foreach()");
	free (consts);
      }
    }
    
    break;
  }

  /**
  ## (const) fields combinations for Point functions */

  case sym_function_definition: {
    if (obsolete_function_declaration (n)) {
      AstTerminal * t = ast_left_terminal (n);
      fprintf (stderr, "%s:%d: warning: obsolete optional/named arguments syntax\n", t->file, t->line);      
    }
    if (ast_is_point_function (ast_schema (n, sym_function_definition,
					   0, sym_function_declaration,
					   1, sym_declarator)) &&
	!ast_is_stencil_function (n)) {
      Ast ** consts = NULL;
      maybeconst (n, stack, append_const, &consts);
      if (consts) {
	Ast * compoundi = ast_schema (n, sym_function_definition,
				      1, sym_compound_statement);
	Ast * compound = ast_copy (compoundi);
	Ast * list = ast_child (compoundi, sym_block_item_list);
	Ast * item = list->child[0];
	if (list->child[1]) {
	  ast_destroy (list->child[1]);
	  list->child[1] = NULL;
	}
	item->sym = sym_block_item;
	ast_destroy (item->child[0]);
	if (item->child[1]) {
	  ast_destroy (item->child[1]);
	  item->child[1] = NULL;
	}
	combinations (compound, stack, data, consts, list, item, "");
	free (consts);
      }
    }
    break;
  }

  case sym_array_access:
    ast_stencil_access (n, stack, ((TranslateData *)data)->dimension);
    break;

  case sym_attribute:
    attribute (n, stack, data);
    break;
    
  /**
  ## Replacement of some identifiers */
  
  case sym_IDENTIFIER: {
    if (n->parent->sym == sym_primary_expression) {
      static Replacement replacements[] = {
	{ "stderr", "ferr" },
	{ "stdout", "fout" },
	{ "qerr", "qstderr()" },
	{ "qout", "qstdout()" },
	{ NULL, NULL }
      };
      Replacement * i = replacements;
      AstTerminal * identifier = ast_terminal (n);
      while (i->target) {
	if (identifier->start && !strcmp (identifier->start, i->target)) {
	  free (identifier->start);
	  identifier->start = strdup (i->replacement);
	}
	i++;
      }
    }
    break;
  }

  /**
  ## Breaks within foreach_inner loops */

  case sym_BREAK: {
    Ast * loop = n->parent;
    while (loop &&
	   loop->sym != sym_foreach_inner_statement &&
	   loop->sym != sym_foreach_statement &&
	   loop->sym != sym_forin_declaration_statement &&
	   loop->sym != sym_forin_statement &&
	   loop->sym != sym_iteration_statement &&
	   (loop->sym != sym_selection_statement ||
	    loop->child[0]->sym != sym_SWITCH))
      loop = loop->parent;
    if (loop && loop->sym == sym_foreach_inner_statement) {
      ast_before (n, ast_terminal (loop->child[0])->start, "_");
      ast_after (n, "()");
    }
    break;
  }

  /**
  ## Constant field and global field allocations */

  case sym_init_declarator: {
    Ast * declarator = declarator_is_allocator (n->child[0]);
    if (declarator) {
      Ast * declaration = ast_declaration_from_type (declarator);
      const char * typename = typedef_name_from_declaration (declaration);
      if (ast_is_field (typename)) {
	char * func = strdup (typename);
	for (char * s = func; *s != '\0'; s++)
	  if (*s == ' ')
	    *s = '_';
	AstTerminal * field = ast_terminal (declarator->child[0]);
	const char * name = field->start;
	if (strchr (typename, ' '))
	  typename = strchr (typename, ' ') + 1;
	
	/**
	### Constant fields initialization */

	if (ast_schema (declaration, sym_declaration,
			0, sym_declaration_specifiers,
			0, sym_type_qualifier,
			0, sym_CONST)) {
	  const char * const_func = strchr (func, '_');
	  const_func = const_func ? const_func + 1 : func;
	  TranslateData * d = data;
	  
	  if (!n->child[1]) {
	    AstTerminal * t = ast_left_terminal (n);
	    fprintf (stderr,
		     "%s:%d: error: constant fields must be initialized\n",
		     t->file, t->line);
	    exit (1);
	  }

	  if (declaration->parent->sym == sym_external_declaration) {

	    /**
	    #### Global constant field declaration */

	    Field * c = field_append (&d->constants, declarator->child[0],
				      typename, d->dimension,
				      &d->constants_index);
	    field->value = (void *)(long) c->index + 65536;
	    char * src = field_value (c, "_NVARMAX+", c->type);
	    char * init = NULL;
	    str_append (init,
			"init_const_", const_func, "((", typename, ")",
			src, ",\"", name, "\",");
	    init = append_initializer (init, n->child[2], typename);
	    str_append (init, ");");
	    Ast * finit = ast_parse_expression (init, ast_get_root (n));
	    free (init);
	    ast_set_line (finit, ast_left_terminal (n->child[2]));
	    compound_append (d->init_fields, NN(n, sym_statement, finit));
	    
	    str_prepend (src, typename, " s=");
	    str_append (src, ";");
	    Ast * expr = ast_parse_expression (src, ast_get_root (n));
	    free (src);
	    ast_replace_child (n, 2, ast_find (expr, sym_initializer));
	    ast_destroy (expr);
	  }
	  else {
	    
	    /**
	    #### Local constant field declaration */

	    char * src = NULL, ind[10];
	    snprintf (ind, 9, "%d", d->constants_index);
	    d->constants_index +=
	      !strcmp (typename, "scalar") ? 1 : d->dimension;
	    str_append (src, "double a = new_const_",
			const_func, "(\"", name, "\",",
			ind, ",");
	    src = append_initializer (src, n->child[2], typename);
	    str_append (src, ");");
	    Ast * expr = ast_parse_expression (src, ast_get_root (n));
	    free (src);
	    ast_replace_child (n, 2, ast_find (expr, sym_initializer));	    
	    ast_destroy (expr);
	  }
	}

	/**
	### Global field allocation */
	
	else if (declaration->parent->sym == sym_external_declaration) {
	  TranslateData * d = data;
	  Field c;
	  c.symmetric = !strcmp (func, "symmetric_tensor");
	  field_init (&c, typename, d->dimension, &d->fields_index);
	  field->value = (void *)(long) c.index + 1;
	  char * src = field_value (&c, "", c.type);
	  char * init = NULL;
	  str_append (init,
		      "  init_", func, "((", typename, ")", src, ",\"",
		      name, "\");");
	  Ast * finit = ast_parse_expression (init, ast_get_root (n));
	  free (init);
	  compound_append (d->init_fields, NN(n, sym_statement, finit));
	  str_prepend (src, typename, " _field_=");
	  str_append (src, ";");

	  Ast * expr = ast_parse_expression (src, ast_get_root (n));
	  free (src);
	  ast_set_line (expr, ast_right_terminal (n->child[0]));
	  declarator = ast_find (expr, sym_init_declarator);
	  ast_replace_child (declarator, 0, n->child[0]);
	  ast_replace_child (n->parent, ast_child_index (n), declarator);
	  n = declarator;
	  ast_destroy (expr);

	  /**
	  #### SWIG interface */

	  if (d->swigname) {
	    str_append (d->swigdecl, "extern ", typename, " ", name, ";\n");
	    str_append (d->swiginit, name, "=", typename,
			 "(_", d->swigname, ".cvar.", name, ")\n");
	  }
	}

	/**
	This is a an automatic (local) field allocations, which is
	treated [at the end of the
	scope](#automatic-field-allocation-and-deallocation) (together
	with deallocation). */
	
	else {
	  free (func);
	  break;
	}
	
	/**
	Remove '[]' from declarator. */
	
	declarator = n->child[0];
	Ast * direct = declarator->child[0];
	ast_replace_child (declarator, 0, direct->child[0]);
	free (func);
      }
    }
    else if (n->child[1] && ast_declaration_from_type (n)->parent->sym
	     == sym_external_declaration) {
	    
      /**
      ### Global constant field initialization */

      Ast * identifier = ast_is_identifier_expression (n->child[2]->child[0]);
      if (identifier) {
	TranslateData * d = data;
	for (Field * c = d->constants; c->identifier; c++)
	  if (!strcmp (ast_terminal (c->identifier)->start,
		       ast_terminal (identifier)->start)) {
	    char * src = field_value (c, "_NVARMAX+", c->type);
	    str_prepend (src, "double s=");
	    str_append (src, ";");
	    Ast * expr = ast_parse_expression (src, ast_get_root (n));
	    free (src);
	    ast_replace_child (n, 2, ast_find (expr, sym_initializer));
	    ast_destroy (expr);
	    break;
	  }
      }
    }
        
    break;
  }

  /**
  ## Function calls */

  case sym_function_call: {
    Ast * identifier = ast_function_call_identifier (n);
    if (identifier) {
      AstTerminal * t = ast_terminal (identifier);
      TranslateData * d = data;

      /**
      ### Memory allocation tracing */

      static Replacement replacements[] = {
	{ "malloc",  "pmalloc" },
	{ "calloc",  "pcalloc" },
	{ "realloc", "prealloc" },
	{ "free",    "pfree" },
	{ "strdup",  "pstrdup" },
	{ NULL, NULL }
      };
      Replacement * i = replacements;
      while (i->target) {
	if (!strcmp (t->start, i->target)) {
	  free (t->start);
	  t->start = strdup (i->replacement);
	  assert (n->child[3]);
	  ast_before (n->child[3], ",__func__,__FILE__,",
		      d->nolineno ? "0" : "__LINE__");
	  return;
	}
	i++;
      }

      /**
      ### Stencil functions */

      int args = stencil_access_function (t->start);
      if (args && (inforeach (n) || point_declaration (stack)))
	complete_arguments (n, args);
      if (!strcmp (ast_terminal (identifier)->start, "_overflow") ||
	  !strcmp (ast_terminal (identifier)->start, "_assign") ||
	  !strcmp (ast_terminal (identifier)->start, "r_assign")) {
	Ast * val = ast_find (n->child[2], sym_function_call);
	if (val) {
	  Ast * name = ast_function_call_identifier (val);
	  str_append (ast_terminal (name)->start,
		      !strcmp (ast_terminal (identifier)->start, "_overflow") ?
		      "_o" :
		      ast_terminal (identifier)->start[0] == '_' ? "_a" : "_r");
	  ast_replace_child (n->parent, ast_child_index (n), val);
	  return;
	}
      }
      
      /**
      ### Macro statement */

      if (n->parent->sym == sym_macro_statement) {
	char * name = NULL;
	str_append (name, "begin_", t->start);
	Ast * type = ast_identifier_declaration (stack, name);
	if (type &&
	    ast_declaration_from_type (type)->sym == sym_function_declaration) {
	  ast_before (identifier, "{");
	  ast_after (n, ";");
	  ast_after (n->parent, "end_", t->start, "();}");
	  free (t->start);
	  t->start = name;
	}
	else // fixme: should this be an error?
	  free (name);
      }
      
      /**
      ### Functions with optional arguments */

      Ast * type = ast_identifier_declaration (stack, t->start);
      if (type) {
	while (type->sym != sym_declarator)
	  type = type->parent;
	if (!ast_schema (type, sym_declarator,
			 0, sym_pointer)) { // exclude function pointers
	  while (type->sym != sym_declaration &&
		 type->sym != sym_function_definition)
	    type = type->parent;
	  Ast * parameters = ast_find (type, sym_parameter_list);
	
	  /**
	  Obsolete optional arguments syntax using 'struct ...' parameters. */

	  Ast * struct_name = obsolete_function_declaration (type);
	  if (struct_name) {
	    Ast * arguments = ast_find (n, sym_argument_expression_list);
	    if (!arguments) {
	      Ast * expr = ast_parse_expression ("func((struct Name){0});",
						 ast_get_root (n));
	      Ast * list = ast_find (expr, sym_argument_expression_list);
	      AstTerminal * t = ast_terminal (ast_find (list, sym_IDENTIFIER));
	      free (t->start);
	      t->start = strdup (ast_terminal (struct_name)->start);
	      ast_set_line (list, ast_terminal (n->child[1]));
	      ast_new_children (n, n->child[0], n->child[1],
				ast_placeholder,
				n->child[2]);
	      ast_replace_child (n, 2, list);
	      ast_destroy (expr);
	    }
	    else {
	      Ast * struct_arg = arguments->child[1] ? NULL :
		ast_is_identifier_expression (arguments->child[0]->child[0]);
	      if (struct_arg) {
		Ast * type =
		  ast_identifier_declaration (stack,
					      ast_terminal (struct_arg)->start);
		while (type &&
		       type->sym != sym_declaration &&
		       type->sym != sym_parameter_declaration)
		  type = type->parent;
		Ast * struct_namep = 
		  ast_get_struct_name (ast_child (type,
						  sym_declaration_specifiers));
		if (!struct_namep ||
		    strcmp (ast_terminal (struct_namep)->start,
			    ast_terminal (struct_name)->start))
		  struct_arg = NULL;
	      }
	      if (!struct_arg) {
		Ast * expr = ast_parse_expression ("func((struct Name){a});",
						   ast_get_root (n));
		Ast * list = ast_find (expr, sym_argument_expression_list);
		AstTerminal * t = ast_terminal (ast_find (list, sym_IDENTIFIER));
		free (t->start);
		t->start = strdup (ast_terminal (struct_name)->start);
		Ast * initializer_list = ast_initializer_list (arguments);
		ast_replace (list, "a", initializer_list);
		ast_replace_child (n, 2, list);
		if (initializer_list->child[1] &&
		    initializer_list->child[1]->sym == token_symbol (',') &&
		    !initializer_list->child[2]) {
		  Ast * postfix = initializer_list->parent;
		  assert (postfix->sym == sym_postfix_initializer &&
			  postfix->child[2]->sym == token_symbol ('}'));
		  ast_new_children (postfix,
				    postfix->child[0],
				    initializer_list->child[0],
				    initializer_list->child[1],
				    postfix->child[2]);
		}
		ast_destroy (expr);
	      }
	    }
	  }
	
	  /**
	  Check for optional or named function call arguments. */
	  
	  else if (ast_find (parameters, sym_parameter_declaration,
			     3, sym_initializer) ||
		   ast_find (n, sym_argument_expression_list_item,
			     0, sym_assignment_expression,
			     1, sym_assignment_operator)) {
	    parameters = ast_copy (parameters); // fixme: memory is leaking from here
	    Ast * parameters1 = parameters;
	    while (parameters && parameters->child[0]->sym == parameters->sym)
	      parameters = parameters->child[0];
	    Ast * arguments = ast_schema (n, sym_function_call,
					  2, sym_argument_expression_list);
	    if (arguments) {
	      foreach_item_r (arguments, sym_argument_expression_list_item, argument) {
		Ast * identifier = ast_schema (argument, sym_argument_expression_list_item,
					       0, sym_assignment_expression,
					       1, sym_assignment_operator) ?
		  ast_schema (argument, sym_argument_expression_list_item,
			      0, sym_assignment_expression,
			      0, sym_unary_expression,
			      0, sym_postfix_expression,
			      0, sym_primary_expression,
			      0, sym_IDENTIFIER) : NULL;
		Ast * parameter;
		if (identifier) {
		  parameter = NULL;
		  foreach_item (parameters1, 2, i) {
		    Ast * id = ast_find (i, sym_direct_declarator,
					 0, sym_generic_identifier,
					 0, sym_IDENTIFIER);
		    if (!strcmp (ast_terminal (identifier)->start, ast_terminal (id)->start)) {
		      parameter = i;
		      break;
		    }
		  }
		  if (!parameter) {
		    AstTerminal * t = ast_terminal (identifier);
		    fprintf (stderr, "%s:%d: error: unknown function parameter '%s'\n",
			     t->file, t->line, t->start);
		    exit (1);
		  }
		  argument = ast_schema (argument, sym_argument_expression_list_item,
					 0, sym_assignment_expression)->child[2];
		}
		else {
		  parameter = ast_child (parameters, sym_parameter_declaration);
		  parameters = parameters->parent;
		  argument = argument->child[0];
		}
		assert (parameter);
		if (ast_schema (parameter, sym_parameter_declaration,
				3, sym_initializer))
		  ast_set_child (parameter->child[3], 0, argument);
		else if (ast_schema (parameter, sym_parameter_declaration,
				     1, sym_declarator))
		  ast_new_children (parameter,
				    parameter->child[0],
				    parameter->child[1],
				    NCA(n, "="),
				    NN(n, sym_initializer,
				       argument));
		else
		  assert (false); // not implemented
		Ast * comma = ast_schema (parameter->parent, sym_parameter_list,
					  1, token_symbol (','));
		if (comma) {
		  AstTerminal * t = ast_terminal (comma), * ta = ast_left_terminal (argument);
		  t->file = ta->file, t->line = ta->line;
		}
	      }
	    }
	    foreach_item (parameters1, 2, parameter) {
	      Ast * initializer = ast_schema (parameter, sym_parameter_declaration,
					      3, sym_initializer);
	      if (!initializer) {
		Ast * id = ast_find (parameter, sym_direct_declarator,
				     0, sym_generic_identifier,
				     0, sym_IDENTIFIER);
		AstTerminal * t = ast_left_terminal (n);
		fprintf (stderr, "%s:%d: error: missing compulsory parameter '%s' in function call\n",
			 t->file, t->line, ast_terminal (id)->start);
		exit (1);
	      }
	      Ast * assign = ast_schema (initializer, sym_initializer,
					 0, sym_assignment_expression);
	      if (assign)
		ast_new_children (parameter, assign);
	      else {
		if (ast_schema (initializer, sym_initializer,
				0, sym_postfix_initializer))
		  initializer = initializer->child[0];
		else
		  assert (ast_schema (initializer, sym_initializer,
				      1, sym_initializer_list));
		initializer->sym = sym_postfix_initializer;
		Ast * type_specifier = ast_find (parameter, sym_declaration_specifiers,
						 0, sym_type_specifier);
		Ast * declarator = ast_schema (parameter, sym_parameter_declaration,
					       1, sym_declarator);
		Ast * abstract = abstract_declarator_from_declarator (declarator);
		assert (type_specifier);
		AstTerminal * ob = NCA(parameter, "("), * cb = NCA(parameter, ")");
		Ast * type_name = abstract ?
		  NN(n, sym_type_name,
		     NN(n, sym_specifier_qualifier_list,
			type_specifier),
		     abstract) :
		  NN(n, sym_type_name,
		     NN(n, sym_specifier_qualifier_list,
			type_specifier));
		ast_new_children (parameter, ast_attach
				  (ast_new_unary_expression (parameter),				 
				   NN(n, sym_postfix_expression,
				      ob, type_name, cb,
				      initializer)));
	      }
	      parameter->sym = sym_argument_expression_list_item;
	      parameter->parent->sym = sym_argument_expression_list;
	    }
	    if (n->child[3])
	      ast_set_child (n, 2, parameters1);
	    else
	      ast_new_children (n, n->child[0], n->child[1], parameters1, n->child[2]);
	  }	  
	}
      }
    }
    break;
    
    if (!identifier || strcmp (ast_terminal (identifier)->start, "automatic"))
      break;
    else {

      /**
      This is a call to automatic() which will be treated with
      sym_NEW_FIELD below. */

      n = identifier;
    }
  }

  /**
  ## `New' and `automatic' fields */

  case sym_NEW_FIELD: {
    Ast * parent = n;
    while (parent &&
	   parent->sym != sym_init_declarator &&
	   (parent->sym != sym_assignment_expression || !parent->child[1]))
      parent = parent->parent;
    if (!parent) {
      AstTerminal * t = ast_terminal (n);
      fprintf (stderr,
	       "%s:%d: error: '%s' must be used within a declarator "
	       "or an assignment expression\n", t->file, t->line, t->start);
      exit (1);
    }
    Ast * identifier = NULL, * declaration = NULL;
    if ((identifier = ast_schema (parent, sym_init_declarator,
				  0, sym_declarator,
				  0, sym_direct_declarator,
				  0, sym_generic_identifier,
				  0, sym_IDENTIFIER)))
      declaration = ast_declaration_from_type (identifier);
    else if ((identifier = ast_schema (parent, sym_assignment_expression,
				       0, sym_unary_expression,
				       0, sym_postfix_expression,
				       0, sym_primary_expression,
				       0, sym_IDENTIFIER))) {
      AstTerminal * t = ast_terminal (identifier);
      declaration = ast_identifier_declaration (stack, t->start);
      if (!declaration) {
	fprintf (stderr,
		 "%s:%d: error: undeclared variable '%s'\n",
		 t->file, t->line, t->start);
	exit (1);
      }
      declaration = ast_declaration_from_type (declaration);
    }
    else {
      AstTerminal * t = ast_terminal (n);
      fprintf (stderr,
	       "%s:%d: error: '%s' must be used to initialize a named field\n",
	       t->file, t->line, t->start);
      exit (1);
    }
    const char * typename = typedef_name_from_declaration (declaration);
    if (ast_is_field (typename)) {
      if (!strstr (ast_terminal (n)->start, typename)) {
	AstTerminal * t = ast_terminal (n);
	fprintf (stderr,
		 "%s:%d: error: type mismatch for `new', "
		 "expected '%s' got '%s'\n",
		 t->file, t->line, typename, t->start);
	exit (1);	
      }

      char * src = strdup (typename);
      for (char * s = src; *s; s++)
	if (*s == ' ')
	  *s = '_';
      Ast * layers = ast_schema (n->parent, sym_new_field,
				 2, sym_postfix_expression);
      if (layers) {
	str_append (src, "(\"", ast_terminal (identifier)->start,
		    !strcmp (src, "scalar") ? "\",\"\"," : "\",");
	src = ast_str_append (layers, src);
	str_append (src, ");");
	str_prepend (src, "new_block_");
      }
      else {
	str_prepend (src, "new_");
	str_append (src, "(\"", ast_terminal (identifier)->start, "\");");
      }
      Ast * expr = ast_parse_expression (src, ast_get_root (n));
      free (src);
      ast_set_line (expr, ast_terminal (n));

      Ast * r = ast_find (expr, sym_assignment_expression);
      ast_remove (n, ast_left_terminal (r));
      if (parent->sym == sym_init_declarator) {
	parent = ast_schema (parent, sym_init_declarator,
			     2, sym_initializer);
	ast_replace_child (parent, 0, r);
      }
      else
	ast_replace_child (parent, 2, r);
      ast_destroy (expr);
    }
    else {
      AstTerminal * t = ast_terminal (n);
      fprintf (stderr,
	       "%s:%d: error: '%s' must be used to initialize a "
	       "scalar, vector or tensor field\n",
	       t->file, t->line, t->start);
      exit (1);      
    }
    break;
  }

  /**
  ## Static FILE * */

  case sym_declaration: {
    Ast * type, * pointer, * identifier, * equal;
    if (ast_schema (n, sym_declaration,
		    0, sym_declaration_specifiers,
		    0, sym_storage_class_specifier,
		    0, sym_STATIC) &&
	(type = ast_schema (n, sym_declaration,
			    0, sym_declaration_specifiers,
			    1, sym_declaration_specifiers,
			    0, sym_type_specifier,
			    0, sym_types,
			    0, sym_TYPEDEF_NAME)) &&
	!strcmp (ast_terminal(type)->start, "FILE") &&
	(pointer = ast_schema (n, sym_declaration,
			       1, sym_init_declarator_list,
			       0, sym_init_declarator,
			       0, sym_declarator,
			       0, sym_pointer)) &&
	!pointer->child[1] &&
	ast_parent (n, sym_event_definition) &&
	(identifier = ast_schema (n, sym_declaration,
				  1, sym_init_declarator_list,
				  0, sym_init_declarator,
				  0, sym_declarator,
				  1, sym_direct_declarator,
				  0, sym_generic_identifier,
				  0, sym_IDENTIFIER)) &&
	(equal = ast_schema (n, sym_declaration,
			     1, sym_init_declarator_list,
			     0, sym_init_declarator,
			     1, token_symbol ('='))))
      ast_after (equal, "NULL;if(!",
		 ast_terminal (identifier)->start,
		 "||i==0)",
		 ast_terminal (identifier)->start,
		 "=pid()>0?fopen(\"/dev/null\",\"w\"):");
    break;
  }

  /**
  ## Events */

  case sym_event_definition: {
    if (!strcmp (ast_left_terminal (n)->start, "event")) {

      /**
      Make the name unique. */
      
      AstTerminal * t = ast_left_terminal (n->child[1]);
      char * name = malloc (strlen (t->start) + 20),
	* suffix = name + strlen(t->start);
      strcpy (name, t->start);
      long last = 0;
      Ast * parent = ast_identifier_declaration (stack, name);
      if (parent)
	last = (long) ast_terminal (parent)->value;
      int i = 0;
      while (ast_identifier_declaration (stack, name))
	snprintf (suffix, 19, "_%d", i++);

      /**
      Define the event expressions. */

      char * iarray = NULL, * tarray = NULL, anexpr[20];
      AstRoot * root = ast_get_root (n);
      int nexpr = 0;
      foreach_item (n->child[3], 2, event_parameter) {
	Ast * initializer = ast_child (event_parameter, sym_postfix_initializer);
	if (initializer) {
	  Ast * identifier = ast_is_identifier_expression
	    (ast_child (event_parameter, sym_unary_expression));
	  if (!identifier || (strcmp (ast_terminal (identifier)->start, "t") &&
			      strcmp (ast_terminal (identifier)->start, "i"))) {
	    AstTerminal * t = ast_left_terminal (event_parameter);
	    fprintf (stderr,
		     "%s:%d: error: an event list can only be used "
		     "to set 't' or 'i'\n", t->file, t->line);
	    exit (1);
	  }
	  snprintf (anexpr, 19, "%d", nexpr++);
	  {
	    char * expr = NULL;
	    str_append (expr, "static int ", name, "_expr", anexpr,
			"(int *ip,double *tp,Event *_ev)"
			"{int i=*ip;double t=*tp;"
			"int ret=(1);*ip=i;*tp=t;return ret;}");
	    Ast * expr0 = ast_parse_external_declaration (expr, root);
	    ast_set_line (expr0, ast_left_terminal (n));
	    free (expr);
	    ast_block_list_insert_before (n->parent->child[0], expr0->child[0]);
	  }
	  {
	    char * expr = NULL;
	    if (!strcmp (ast_terminal (identifier)->start, "t")) {
	      str_append (tarray, name, "_array");
	      str_append (expr, "static double ", tarray, "[]=");
	    }
	    else {
	      str_append (iarray, name, "_array");
	      str_append (expr, "static int ", iarray, "[]=");
	    }
	    ast_before (ast_last_child(initializer), ",-1");
	    expr = ast_str_append (initializer, expr);
	    str_append (expr, ";");
	    Ast * expr0 = ast_parse_external_declaration (expr, root);
	    ast_set_line (expr0, ast_left_terminal (n));
	    free (expr);
	    ast_block_list_insert_before (n->parent->child[0], expr0->child[0]);
	  }
	  break;
	}
	else {
	  Ast * identifier;
	  if (!ast_child (event_parameter, sym_assignment_operator) &&
	      (identifier = ast_is_identifier_expression
	       (ast_child (event_parameter, sym_conditional_expression))) &&
	      (!strcmp (ast_terminal (identifier)->start, "last") ||
	       !strcmp (ast_terminal (identifier)->start, "first"))) {
	    if (!strcmp (ast_terminal (identifier)->start, "last"))
	      last = 1;
	    else
	      last = 0;
	  }
	  else {
	    snprintf (anexpr, 19, "%d", nexpr++);
	    char * expr = NULL;
	    str_append (expr, "static int ", name, "_expr", anexpr,
			"(int *ip,double *tp,Event *_ev)"
			"{int i=*ip;double t=*tp;"
			"int ret=(");
	    Ast * rhs = ast_child (event_parameter,
				   sym_conditional_expression), * identifier;
	    if (rhs && (identifier = ast_is_identifier_expression (rhs)) &&
		!strcmp (ast_terminal (identifier)->start, "end")) {
	      free (ast_terminal (identifier)->start);
	      ast_terminal (identifier)->start = strdup ("TEND_EVENT");
	    }
	    expr = ast_str_append (event_parameter, expr);
	    str_append (expr, ")!=0;*ip=i;*tp=t;return ret;}");
	    Ast * expr0 = ast_parse_external_declaration (expr, root);
	    ast_set_line (expr0, ast_left_terminal (n));
	    free (expr);
	    ast_block_list_insert_before (n->parent->child[0], expr0->child[0]);
	  }
	}
      }
      
      /**
      Register the event. */

      char * reg = NULL;
      snprintf (anexpr, 19, "%d", nexpr);
      str_append (reg, "  event_register((Event){0,", anexpr, ",", name, ",{");
      for (int i = 0; i < nexpr; i++) {
	snprintf (anexpr, 19, "%d", i);
	str_append (reg, name, "_expr", anexpr, i < nexpr - 1 ? "," : "");
      }
      TranslateData * d = data;
      str_append (reg, "},",
		  iarray ? iarray : "((int *)0)",
		  ",",
		  tarray ? tarray : "((double *)0)",
		  ",",
		  ast_file_line (t, d->nolineno), ",\"", t->start, "\"});\n");
      Ast * registration = NN(n, sym_statement,
			      ast_parse_expression (reg, root));
      ast_set_line (registration, t);
      if (last)
	compound_append (d->last_events, registration);
      else
	compound_append (d->init_events, registration);
      free (reg);
      free (iarray);
      free (tarray);
      
      /**
      Define the action fonction. */
      
      char * src = NULL;
      Ast * statement = ast_child (n, sym_statement);
      str_append (src,
		  ast_schema (statement, sym_statement,
			      0, sym_compound_statement,
			      1, token_symbol ('}')) ||
		  ast_schema (statement, sym_statement,
			      0, sym_expression_statement,
			      0, token_symbol (';'))
		  ? "" : "trace ",
		  "static int ", name,
		  "(const int i,const double t,Event *_ev)"
		  "{_statement_;return 0;}");
      Ast * def = ast_parse_external_declaration (src, root);
      Ast * identifier = ast_schema (def, sym_external_declaration,
				     0, sym_function_definition,
				     0, sym_function_declaration,
				     1, sym_declarator,
				     0, sym_direct_declarator,
				     0, sym_direct_declarator,
				     0, sym_generic_identifier,
				     0, sym_IDENTIFIER);
      ast_terminal (identifier)->value = (void *) last;
      free (src);
      ast_replace (def, "_statement_", statement);
      ast_replace_child (n->parent, 0, def->child[0]);
      ast_destroy (def);

      free (name);      
    }
    break;
  }
    
  /**
  ## Automatic field deallocation before jump statements */

  case sym_jump_statement: {
    if (n->child[0]->sym == sym_GOTO) {
      AstTerminal * t = ast_terminal (n->child[0]);
      fprintf (stderr, "%s:%d: warning: goto statements are unsafe in Basilisk "
	       "(and are bad programming style)\n",
	       t->file, t->line);
      break;
    }

    int jump_sym = n->child[0]->sym;
    Ast * parent = n;
    while (parent &&
	   ((jump_sym == sym_BREAK &&
	     parent->child[0]->sym != sym_SWITCH &&
	     !ast_is_iteration_statement (parent)) ||
	    (jump_sym == sym_CONTINUE &&
	     !ast_is_iteration_statement (parent)) ||
	    (jump_sym == sym_RETURN &&
	     parent->sym != sym_function_definition &&
	     parent->sym != sym_event_definition)))
      parent = parent->parent;
    Ast * scope = ast_find (parent, sym_compound_statement);
    if (scope) {
      char * delete[2] = {NULL};
      foreach_field_allocator (stack, data, scope, field_deallocation, delete);
      char * fields = delete_fields (delete);
      if (fields)
	compound_jump (n, parent, fields);
      free (delete[0]);
      free (delete[1]);
    }
    break;
  }

  /**
  Boundary ids */

  case sym_external_declaration: {
    Ast * identifier = ast_schema (n, sym_external_declaration,
				   0, sym_declaration,
				   0, sym_declaration_specifiers,
				   0, sym_type_specifier,
				   0, sym_types,
				   0, sym_TYPEDEF_NAME);
    if (identifier && !strcmp (ast_terminal (identifier)->start, "bid")) {
      Ast * list = ast_schema (n, sym_external_declaration,
			       0, sym_declaration,
			       1, sym_init_declarator_list);
      if (list)
	foreach_item (list, 2, item)
	  if ((identifier = ast_schema (item, sym_init_declarator,
					0, sym_declarator,
					0, sym_direct_declarator,
					0, sym_generic_identifier,
					0, sym_IDENTIFIER))) {
	    Ast * init =
	      NN(n, sym_statement,
		 NN(n, sym_expression_statement,
		    NN(n, sym_expression,
		       NN(n, sym_assignment_expression,
			  NN(n, sym_unary_expression,
			     NN(n, sym_postfix_expression,
				NN(n, sym_primary_expression,
				   ast_copy (identifier)))),
			  NN(n, sym_assignment_operator,
			     NCA(n, "=")),
			  ast_new_assignment_function_call (n, "new_bid"))),
		    NCA(n, ";")));
	    TranslateData * d = data;
	    compound_append (d->init_fields, init);
	  }
    }
    break;
  }

  }

  /**
  ## Automatic field allocation and deallocation */

  if (n->sym == token_symbol('}') && n->parent->sym == sym_compound_statement) {
    char * delete[2] = {NULL};
    foreach_field_allocator (stack, data, n->parent, field_allocation, delete);
    
    /**
    ### Field deallocation */

    char * fields = delete_fields (delete);
    if (fields) {
      Ast * expr = ast_parse_expression (fields, ast_get_root (n));
      ast_block_list_append (ast_child (n->parent, sym_block_item_list),
			     sym_block_item,
			     ast_new_children (ast_new (n, sym_statement),
					       expr));
    }
    free (delete[0]);
    free (delete[1]);
  }
}

static void trace_return (Ast * n, Stack * stack, void * data)
{
  Ast * function_definition = ((void **)data)[0];
  AstTerminal * function_identifier = ((void **)data)[1];
  if (ast_schema (n, sym_jump_statement, 0, sym_RETURN)) {
    char * end_tracing = NULL;
    TranslateData * d = ((void **)data)[2];
    str_append (end_tracing,
		"end_tracing(\"", function_identifier->start, "\",",
		ast_file_line (n->child[0], d->nolineno), ");");
    compound_jump (n, function_definition, end_tracing);
    free (end_tracing);
  }
}

static const char * get_field_type (Ast * declaration, AstTerminal * t)
{
  if (declaration)
    declaration = ast_declaration_from_type (declaration);
  const char * typename = NULL;
  if (!declaration ||
      !(typename = typedef_name_from_declaration (declaration)) ||
      (strcmp (typename, "scalar") &&
       strcmp (typename, "vector") &&
       strcmp (typename, "tensor"))) {
    fprintf (stderr,
	     "%s:%d: error: '%s' is not a scalar, vector or tensor\n",
	     t->file, t->line, t->start);
    exit (1);
  }
  return typename;
}

static void mpi_operator (Ast * n, Ast * op)
{
  char * operator = ast_left_terminal (op)->start;
  ast_after (n,
	     !strcmp(operator, "min") ? "MPI_MIN" :
	     !strcmp(operator, "max") ? "MPI_MAX" :
	     !strcmp(operator, "+")   ? "MPI_SUM" :
	     !strcmp(operator, "||")  ? "MPI_LOR" :
	     "Unknown", ",");
}

static void maps (Ast * n, Stack * stack, void * data)
{
  if (n->sym == sym_foreach_statement && !ast_is_foreach_stencil (n)) {
    if (!strcmp (ast_terminal (n->child[0])->start, "foreach_face")) {
      free (ast_terminal (n->child[0])->start);
      ast_terminal (n->child[0])->start = strdup ("foreach_face_generic");
    }
    else { // maps for !foreach_face() loops
      foreach_map (m) {
	Ast * list = ast_block_list_get_item (ast_child (n, sym_statement))->parent;
	foreach_item (m, 1, item)
	  ast_block_list_prepend (list, sym_block_item, ast_copy (item->child[0]));
      }
    }    
  }
}

/**
# Third pass: foreach stencils */

static void stencils (Ast * n, Stack * stack, void * data)
{
  switch (n->sym) {
    
  case sym_foreach_statement:
    assert (ast_last_child(n)->sym == sym_statement); // make sure all stencils have been dealt with
    if (ast_is_foreach_stencil (n)) {
      Ast * foreach = n->parent;
      if (foreach->sym == sym_foreach_statement && ast_last_child (foreach) == n) {
	Ast ** c;
	for (c = foreach->child; *c; c++);
	*(--c) = NULL;
	Ast * statement = foreach->parent->parent;
	Ast * item = ast_block_list_get_item (statement), * list = item->parent;
	list = ast_block_list_append
	  (list, item->sym,
	   ast_new_children (ast_new (foreach, sym_statement),
			     ast_new_children (ast_new (foreach,
							sym_basilisk_statements),
					       n)));
	ast_set_child (item, 0, list->child[1]->child[0]);
	ast_set_child (list->child[1], 0, statement);
      }

      Ast * parameters = ast_child (n, sym_foreach_parameters);
      if (parameters) {
	ast_destroy (parameters);
	for (Ast ** c = n->child + 2; *c; c++)
	  *c = *(c + 1);
      }
      
      /**
      ## Kernel for GPUs */

      Ast * params = ast_schema (n, sym_foreach_statement,
				 1, token_symbol('('));
      assert (params);
      parameters = ast_child (foreach, sym_foreach_parameters);
      bool gpu = !((TranslateData *)data)->cpu;
      int parallel = (gpu ?
		      1 :  // parallel on CPU || GPU
		      2 ); // parallel on CPU
      if (parameters)
	foreach_item (parameters, 2, item) {
	  Ast * identifier = ast_is_identifier_expression (item->child[0]);
	  if (identifier) {
	    AstTerminal * t = ast_terminal (identifier);
	    if (!strcmp (t->start, "serial"))
	      parallel = 0, gpu = false;
	    else if (!strcmp (t->start, "cpu"))
	      parallel = 2, gpu = false;
	    else if (!strcmp (t->start, "gpu"))
	      parallel = 3, gpu = true;
	  }
	}

      ast_after (params,
		 parallel == 0 ? "0" : parallel == 1 ? "1" : parallel == 2 ? "2" : "3",
		 ",{(NonLocal[])");
      ast_non_local_references (foreach, params);
      if (gpu) {
	ast_after (params, ",");
	ast_kernel (foreach, params);
	ast_after (params, "}");
      }
      else
	ast_after (params, ",NULL}");
    }

    break;
    
  /**
  ## Foreach inner statements */

  case sym_foreach_inner_statement: {
    AstTerminal * t = ast_left_terminal(n);
    if (!strcmp (t->start, "foreach_block") &&
	(inforeach (n) || point_declaration (stack)))
      str_append (t->start, "_inner");
    ast_before (n, "{");
    ast_after (n, "end_", t->start, "()}");
    break;
  }
    
  case sym_macro_statement: {
    Ast * identifier = ast_schema (n, sym_macro_statement,
				   0, sym_function_call,
				   0, sym_postfix_expression,
				   0, sym_primary_expression,
				   0, sym_IDENTIFIER);
    if (identifier) {
      AstTerminal * t = ast_terminal (identifier);

      /**
      ## is_face_... statements */

      if (!strcmp (t->start, "is_face_x") ||
	  !strcmp (t->start, "is_face_y") ||
	  !strcmp (t->start, "is_face_z")) {
	foreach_map (m) {
	  Ast * list = ast_schema (n, sym_macro_statement,
				   1, sym_compound_statement,
				   1, sym_block_item_list);
	  foreach_item (m, 1, item)
	    ast_block_list_prepend (list, sym_block_item, ast_copy (item->child[0]));
	}
	ast_after (n, "end_", t->start, "()");
      }

      /**
      ## _stencil_is_face_... statements */

      else if (!strcmp (t->start, "_stencil_is_face_x") ||
	       !strcmp (t->start, "_stencil_is_face_y") ||
	       !strcmp (t->start, "_stencil_is_face_z"))
	ast_after (n, "end_", t->start, "()");

      /**
      ## Map */

      else if (!strcmp (ast_terminal (identifier)->start, "map")) {
	Ast * item = n->parent;
	if (item->sym == sym_external_declaration) {
	  assert (ast_child_index (item) == 1);
	  Ast * parent = item->parent, * grand_parent = parent->parent;
	  ast_set_child (grand_parent, ast_child_index (parent),
			 parent->child[0]);
	}
      }
    }
    break;
  }
    
  /**
  ## Point point */
  
  case sym_IDENTIFIER: {
    Ast * decl = ast_is_point_point (n);
    if (decl) {
      AstTerminal * t = ast_terminal (n);
      if (ast_parent (n, sym_declaration) &&
	  strncmp (t->file, BASILISK "/grid/", strlen (BASILISK "/grid/")) &&
	  strncmp (t->file, "./grid/", strlen ("./grid/")))
	fprintf (stderr,
		 "%s:%d: warning: 'Point point' is obsolete, use 'foreach_point/region' instead\n",
		 t->file, t->line);
      TranslateData * d = data;
      static const char * name[3] = {"ig", "jg", "kg"};
      for (int i = 0; i < d->dimension; i++)
	ast_after (decl, "int ", name[i], "=0;"
		   "NOT_UNUSED(", name[i], ");");
      str_append (ast_right_terminal (decl)->after, "POINT_VARIABLES;");
      if (decl->sym == sym_declaration) {
	foreach_map (m) {
	  Ast * list = ast_block_list_get_item (decl)->parent;
	  foreach_item (m, 1, item)
	    ast_block_list_append (list, sym_block_item, ast_copy (item->child[0]));
	}
      }
      else {
	Ast * list = ast_schema (decl->parent, sym_compound_statement,
				 1, sym_block_item_list);
	foreach_map (m)
	  foreach_item (m, 1, item)
	    ast_block_list_prepend (list, sym_block_item, ast_copy (item->child[0]));
      }
    }
    break;
  }
    
  /**
  ## Hide Basilisk C keywords */

  case sym_MAYBECONST: ast_hide (ast_terminal (n)); break;
  case sym_TYPEDEF_NAME: {
    AstTerminal * t = ast_terminal (n);
    if (!strcmp (t->start, "face vector") ||
	!strcmp (t->start, "vertex scalar") ||
	!strcmp (t->start, "symmetric tensor")) {
      char * s = strchr (t->start, ' ') + 1;
      memmove (t->start, s, strlen (s) + 1);
    }
    break;
  }

  /**
  ## Remove '_val_higher_dimension' statements with no effect (to avoid compiler warnings) */

  case sym_expression_statement: {
    Ast * id;
    if ((id = ast_is_identifier_expression (ast_schema (n, sym_expression_statement,
							0, sym_expression,
							0, sym_assignment_expression))) &&      
	!strcmp (ast_terminal (id)->start, "_val_higher_dimension")) {
      ast_destroy (n->child[0]);
      n->child[0] = n->child[1]; n->child[1] = NULL;
    }
    break;
  }
    
  }
}

static char * get_type (const char * name, Stack * stack)
{
  Ast * decl = ast_find (ast_declaration_from_type (ast_identifier_declaration (stack, name)),
			 sym_declaration_specifiers);
  if (!decl) return NULL;
  AstTerminal * t = ast_left_terminal (decl);
  char * before = t->before;
  t->before = NULL;
  char * type = ast_str_append (decl, NULL);
  t->before = before;
  return type;
}

/**
# Fourth pass: "macro" expressions 

This pass should regroup all transformations which require the use of
macros which are not included in the Basilisk C grammar. */

static void macros (Ast * n, Stack * stack, void * data)
{
  switch (n->sym) {

  /**
  ## Stencil access function calls */

  case sym_function_call: {
    Ast * identifier = ast_function_call_identifier (n);
    if (identifier) {
      AstTerminal * t = ast_terminal (identifier);
      Ast * foreach = NULL;
      if (stencil_access_function (t->start) &&
	  (((foreach = inforeach (n)) && ast_is_foreach_stencil (foreach)) ||
	   in_stencil_point_function (n)))
	str_prepend (t->start, "_stencil_");
    }
    break;
  }

  case sym_foreach_statement: {

    assert (ast_last_child(n)->sym == sym_statement); // make sure all stencils have been dealt with

    /**
    ## Foreach stencils */
    
    if (ast_is_foreach_stencil (n)) {
      ast_after (n, "end_", ast_left_terminal(n)->start, "();");
      break;
    }

    /**
    ## Foreach statements */

    Ast * foreach = inforeach (n);
    if (foreach) {
      AstTerminal * t = ast_terminal (n->child[0]);
      AstTerminal * p = ast_terminal (foreach->child[0]);
      fprintf (stderr,
	       "%s:%d: error: this %s cannot include\n", p->file, p->line,
	       foreach->sym == sym_foreach_statement ?
	       ast_is_foreach_stencil (foreach) ?
	       "'Point point' scope" :
	       "foreach*() iterator" :
	       "point function");
      fprintf (stderr,
	       "%s:%d: error: this %s\n", t->file, t->line,
	       n->sym == sym_foreach_statement ?
	       ast_is_foreach_stencil (n) ?
	       "'Point point' scope" :
	       "foreach*() iterator" :
	       "point function");
      exit (1);
    }
        
    ast_after (n, "end_", ast_left_terminal(n)->start, "();");

    /**
    ### Reductions */

    Ast * parameters = ast_child (n, sym_foreach_parameters);
    bool serial = false;
    char * sreductions = NULL;
    if (parameters) {
      foreach_item (parameters, 2, item) {
	Ast * identifier = ast_is_identifier_expression (item->child[0]);
	if (identifier && !strcmp (ast_terminal (identifier)->start, "serial")) {
	  serial = true;
	  parameters = ast_list_remove (parameters, item);
	}
	else if (identifier && (!strcmp (ast_terminal (identifier)->start, "cpu") ||
				!strcmp (ast_terminal (identifier)->start, "gpu")))
	  parameters = ast_list_remove (parameters, item);
	else if (item->child[0]->sym == sym_reduction_list) {
	  Ast * reductions = item->child[0];
	  foreach_item (reductions, 1, reduction) {
	    Ast * identifier = ast_schema (reduction, sym_reduction,
					   4, sym_reduction_array,
					   0, sym_generic_identifier,
					   0, sym_IDENTIFIER);
	    AstTerminal * t = ast_terminal (identifier);
	    Ast * array = ast_schema (reduction, sym_reduction,
				      4, sym_reduction_array,
				      3, sym_expression);
	    char * type = get_type (t->start, stack);
	    if (!type) {
	      fprintf (stderr,
		       "%s:%d: error: cannot determine type of '%s'\n",
		       t->file, t->line, t->start);
	      exit (1);
	    }
	    if (strcmp (type, "coord") &&
		strcmp (type, "mat3") &&
		strcmp (type, "double") &&
		strcmp (type, "int") &&
		strcmp (type, "long") &&
		strcmp (type, "bool") &&
		strcmp (type, "unsigned char")) {
	      fprintf (stderr,
		       "%s:%d: error: does not know how to reduce "
		       "type '%s' of '%s'\n",
		       t->file, t->line, type, t->start);
	      exit (1);
	    }
	    if (array) {
	      if (strcmp (type, "coord") && strcmp (type, "mat3")) {
		ast_after (n, "mpi_all_reduce_array(", t->start, ",", type, ",");
		mpi_operator (n, reduction->child[2]);
		ast_right_terminal (n)->after = ast_str_append (array, ast_right_terminal (n)->after);
		ast_after (n, ");");
	      } else {
		ast_after (n, "mpi_all_reduce_array((double *)", t->start, ",double,");
		mpi_operator (n, reduction->child[2]);
		char s[100];
		snprintf (s, 99, "sizeof(%s)/sizeof(double)", t->start);
		ast_after (n, s, ");");
	      }
	    }
	    else {
	      char s[100] = "1";
	      if (strcmp (type, "coord") && strcmp (type, "mat3"))
		ast_after (n, "mpi_all_reduce_array(&", t->start,",", type);
	      else {
		// cast the adress of the first member into a double for coord and mat3
		ast_after (n, "mpi_all_reduce_array((double *)&", t->start,",double");
		snprintf (s, 99, "sizeof(%s)/sizeof(double)", t->start);
	      }
	      ast_after (n, ",");
	      mpi_operator (n, reduction->child[2]);
	      ast_after (n, s, ");");
	    }
	    sreductions = ast_str_append (reduction, sreductions);
	    free (type);
	  }
	  parameters = ast_list_remove (parameters, item);
	}
      }
      if (parameters == NULL) {
	ast_destroy (n->child[2]);
	for (Ast ** c = n->child + 2; *c; c++)
	  *c = *(c + 1);
      }
    }

    if (serial)
      ast_before (n, "\n"
		  "#if _OPENMP\n"
		  "  #undef OMP\n"
		  "  #define OMP(x)\n"
		  "#endif\n");
    if (sreductions) {
      ast_before (n, "\n"
		  "#undef OMP_PARALLEL\n"
		  "#define OMP_PARALLEL()\n"
		  "OMP(omp parallel ", sreductions, "){");
      ast_after (n,
		 "\n"
		 "#undef OMP_PARALLEL\n"
		 "#define OMP_PARALLEL() OMP(omp parallel)\n}");
      free (sreductions);
    }
    else {
      ast_before (n, "{");
      ast_after (n, "}");
    }
    if (serial)
      ast_after (n, "\n"
		 "#if _OPENMP\n"
		 "  #undef OMP\n"
		 "  #define OMP(x) _Pragma(#x)\n"
		 "#endif\n");
    
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
    const char * typename =
      get_field_type (n->child[2], ast_terminal(identifier));
    char * src = NULL, * name = ast_terminal(identifier)->start;
    str_append (src, "{", typename, "*_i=(", typename, "*)(list);if(_i)"
		"for(", typename, " ", name, "=*_i;(&",
		name,
		!strcmp (typename, "scalar") ? ")->i" :
		!strcmp (typename, "vector") ? ")->x.i" :
		")->x.x.i",
		">=0;", name, "=*++_i){_statement_;}}");
    Ast * expr = ast_parse_expression (src, ast_get_root (n));
    free (src);
    Ast * parent = n->parent;
    Ast * arg = ast_child (n, sym_forin_arguments)->child[0];
    if (arg->sym == sym_expression) {
      Ast * initializer = ast_find (expr, sym_expression_error);
      ast_replace_child (initializer, 0, arg);
    }
    else {
      arg = ast_find (arg, sym_postfix_initializer);
      Ast * initializer = ast_find (expr, sym_cast_expression);
      ast_replace_child (initializer, 3, arg);
      initializer->sym = sym_postfix_expression;
      Ast * parent = initializer->parent;
      int index = ast_child_index (initializer);
      Ast * unary = ast_new_children (ast_new (n, sym_unary_expression),
				      initializer);
      Ast * cast = ast_new_children (ast_new (n, sym_cast_expression),
				     unary);
      char * before = ast_left_terminal (arg)->before;
      ast_replace_child (parent, index, cast);
      ast_left_terminal (arg)->before = before;
      ast_left_terminal (parent->child[index])->before = NULL;
    }
    assert (ast_replace (expr, "_statement_", ast_child (n, sym_statement)));
    ast_replace_child (parent->parent, ast_child_index (parent), expr);
    break;
  }

  /**
  ## forin_statement */

  case sym_forin_statement: {
    Ast * arg = ast_child (n, sym_forin_arguments)->child[0];
    char * decl = strdup ("{"), * fors = strdup ("if(_i0)for("), * fore = NULL;
    int index = 0;
    foreach_item (n->child[2], 2, expr) {
      Ast * identifier = ast_is_identifier_expression (expr);
      if (!identifier) {
	AstTerminal * t = ast_left_terminal (expr);
	fprintf (stderr,
		 "%s:%d: error: not a scalar, vector or tensor\n",
		 t->file, t->line);
	exit (1);
      }
      AstTerminal * t = ast_terminal (identifier);
      const char * typename =
	get_field_type (ast_identifier_declaration (stack, t->start), t);
      if (!arg) {
	fprintf (stderr,
		 "%s:%d: error: lists must have the same size\n",
		 t->file, t->line);
	exit (1);
      }
      Ast * l;
      if (arg->sym == sym_postfix_initializer || !arg->child[1]) {
	l = arg;
	arg = NULL;
      }
      else {
	l = arg->child[2];
	arg = arg->child[0];
      }
      char ind[20];
      snprintf (ind, 19, "%d", index);
      str_append (decl, typename, "*_i", ind, "=");
      decl = ast_str_append (l, decl);
      str_append (decl, ";");
      str_append (fors, index > 0 ? "," : "", t->start, "=*_i", ind);
      if (!fore)
	str_append (fore, "_i", ind,
		    !strcmp (typename, "scalar") ? "->i" :
		    !strcmp (typename, "vector") ? "->x.i" :
		    "->x.x.i",
		    ">= 0;");
      str_append (fore, index > 0 ? "," : "", t->start, "=*++_i", ind);
      index++;
    }
    str_append (decl, fors, ";", fore, "){_statement_;}}");
    Ast * expr = ast_parse_expression (decl, ast_get_root (n));
    free (decl); free (fors); free (fore);
    assert (ast_replace (expr, "_statement_", ast_child (n, sym_statement)));
    ast_replace_child (n->parent->parent, ast_child_index (n->parent), expr);
    break;
  }

  /**
  ## Attribute access */

  case sym_postfix_expression: {
    if (n->child[1] && n->child[1]->sym == token_symbol('.')) {
      const char * typename =
	ast_typedef_name (ast_expression_type (n->child[0], stack, false));
      if (typename && (!strcmp (typename, "scalar") ||
		       !strcmp (typename, "vertex scalar"))) {
	Ast * member = ast_find (n->child[2], sym_member_identifier,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER);
	Ast * type = ast_identifier_declaration (stack, "scalar");
	assert (type);
	while (type->sym != sym_declaration)
	  type = type->parent;
	if (!find_struct_member (ast_find (type, sym_struct_declaration_list),
				 ast_terminal (member)->start)) {
	  Ast * expr = ast_parse_expression ("_attribute[_field_.i];",
					     ast_get_root (n));
	  ast_replace (expr, "_field_", n->child[0]);
	  ast_replace_child (n, 0, ast_find (expr, sym_postfix_expression));
	  ast_destroy (expr);
	}
      }

      /**
      ## Boundary vector component access */
      
      else if (typename && (!strcmp (typename, "vector") ||
			    !strcmp (typename, "face vector")))
	set_boundary_component (ast_find (n->child[2], sym_member_identifier));
    }
    break;
  }

  /**
  ## Field lists */

  case sym_postfix_initializer: {
      
    /**
    Do not consider lists explicitly cast as structures. */

    if (n->parent->sym == sym_postfix_expression &&
	ast_child_index (n) == 3 &&
	ast_schema (n->parent->child[1], sym_type_name,
		    0, sym_specifier_qualifier_list,
		    0, sym_type_specifier,
		    0, sym_types,
		    0, sym_struct_or_union_specifier,
		    0, sym_struct_or_union,
		    0, sym_STRUCT))
      break;
    // else fall through
  }

  case sym_initializer: {
    if (n->child[1] && n->child[2]) {
      Ast * list = n->child[1];
      int type = field_list_type (list, stack, false);
      if (type > 0) {

	/**
	### External/global lists */
	
	bool external = true;
	Ast * scope = n->parent;
	while (external && scope) {
	  if (scope->sym == sym_compound_statement)
	    external = false;
	  scope = scope->parent;
	}
	if (external) {
	  char * src = NULL;
	  foreach_item (list, 2, expr) {
	    const char * typename =
	      ast_typedef_name (ast_expression_type (expr, stack, false));
	    assert (typename);
	    Ast * unary = ast_is_unary_expression (expr->child[0]);
	    if (!unary) {
	      AstTerminal * t = ast_terminal (expr);
	      fprintf (stderr,
		       "%s:%d: error: global lists can only be initialized "
		       "with simple expressions\n", t->file, t->line);
	      exit (1);
	    }
	    int stype = 1; // identifier
	    Ast * identifier = ast_schema (unary, sym_unary_expression,
					   0, sym_postfix_expression,
					   0, sym_primary_expression,
					   0, sym_IDENTIFIER);
	    if (!identifier) {
	      stype = 2;  // identifier.x
	      identifier = ast_schema (unary, sym_unary_expression,
				       0, sym_postfix_expression,
				       0, sym_postfix_expression,
				       0, sym_primary_expression,
				       0, sym_IDENTIFIER);
	    }
	    if (!identifier) {
	      stype = 3;  // identifier.x.x
	      identifier = ast_schema (unary, sym_unary_expression,
				       0, sym_postfix_expression,
				       0, sym_postfix_expression,
				       0, sym_postfix_expression,
				       0, sym_primary_expression,
				       0, sym_IDENTIFIER);
	    }
	    if (identifier) {
	      AstTerminal * t = ast_terminal (identifier);
	      Ast * declaration = ast_identifier_declaration (stack, t->start);
	      const char * typename1 = get_field_type (declaration, t);
	      if (!ast_terminal (declaration)->value) {
		fprintf (stderr,
			 "%s:%d: error: variable `%s' is not initialized\n",
			 t->file, t->line, t->start);
		exit (1);			
	      }
	      Field c = {
		.identifier = NULL,
		.type = (!strcmp (typename1, "scalar") ? 1 :
			 !strcmp (typename1, "vector") ? 2 : 3),
		.index = ((long) ast_terminal (declaration)->value) - 1,
		.dimension = ((TranslateData *)data)->dimension
	      };
	      if (stype > 1) { // .x or .x.x
		Ast * member = ast_schema (unary, sym_unary_expression,
					   0, sym_postfix_expression,
					   2, sym_member_identifier,
					   0, sym_generic_identifier,
					   0, sym_IDENTIFIER);
		if (member) {
		  if (stype == 2) { // .x
		    c.index += (ast_terminal(member)->start[0] - 'x')*
		      (c.type == 3 ? c.dimension : 1);
		    if (type == 1 && c.type == 3)
		      c.type = 2;
		    else
		      c.type = type;
		  }
		  else if (stype == 3) { // .x.x
		    Ast * member1 = ast_schema (unary, sym_unary_expression,
						0, sym_postfix_expression,
						0, sym_postfix_expression,
						2, sym_member_identifier,
						0, sym_generic_identifier,
						0, sym_IDENTIFIER);
		    if (member1) {
		      c.index += (ast_terminal(member)->start[0] - 'x') +
			(ast_terminal(member1)->start[0] - 'x')*c.dimension;
		      c.type = 1;
		      ast_terminal(member1)->start[0] = '\0';
		      Ast * dot = ast_schema (unary, sym_unary_expression,
					      0, sym_postfix_expression,
					      0, sym_postfix_expression,
					      1, token_symbol ('.'));
		      ast_terminal(dot)->start[0] = '\0';
		    }
		  }
		  ast_terminal(member)->start[0] = '\0';
		  Ast * dot = ast_schema (unary, sym_unary_expression,
					  0, sym_postfix_expression,
					  1, token_symbol ('.'));
		  ast_terminal(dot)->start[0] = '\0';
		}
	      }
	      str_prepend (src, field_value (&c, "", type), src ? "," : "");
	    }
	  }
	  if (src) {
	    str_prepend (src, "int a = {");
	    str_append (src, "};");
	    Ast * expr = ast_parse_expression (src, ast_get_root (n));
	    free (src);
	    Ast * initializer = ast_find (expr, sym_initializer);
	    ast_replace_child (n->parent, ast_child_index (n), initializer);
	    ast_destroy (expr);
	    n = initializer;
	  }
	}

	/**
	### Local lists */

	else
	  foreach_item (list, 2, expr) {
	    const char * typename =
	      ast_typedef_name (ast_expression_type (expr, stack, false));
	    if ((type == 1 && !strcmp (typename, "vector")) ||
		(type == 2 && !strcmp (typename, "tensor"))) {
	      Ast * unary = ast_find (expr, sym_unary_expression);
	      ast_set_child (unary, 0,
			     NN(expr, sym_postfix_expression,
				unary->child[0],
				NCA(expr, "."),
				NN(expr, sym_member_identifier,
				   NN(expr, sym_generic_identifier,
				      NA(expr, sym_IDENTIFIER, "x")))));
	      Ast * a = expr->child[0];
	      TranslateData * d = data;
	      for (int i = 1; i < d->dimension; i++) {
		Ast * b = ast_copy (a);
		Ast * id = ast_find (b, sym_member_identifier,
				     0, sym_generic_identifier,
				     0, sym_IDENTIFIER);
		ast_terminal (id)->start[0] = 'x' + i;
		ast_list_insert_after (a, b);
		a = b;
	      }
	    }
	    else if (type == 1 && !strcmp (typename, "tensor")) {
	      Ast * unary = ast_find (expr, sym_unary_expression);
	      ast_set_child (unary, 0,
			     NN(expr, sym_postfix_expression,
				NN(expr, sym_postfix_expression,
				   unary->child[0],
				   NCA(expr, "."),
				   NN(expr, sym_member_identifier,
				      NN(expr, sym_generic_identifier,
					 NA(expr, sym_IDENTIFIER, "x")))),
				NCA(expr, "."),
				NN(expr, sym_member_identifier,
				   NN(expr, sym_generic_identifier,
				      NA(expr, sym_IDENTIFIER, "x")))));
	      Ast * a = expr->child[0];
	      TranslateData * d = data;
	      for (int i = 0; i < d->dimension; i++)
		for (int j = 0; j < d->dimension; j++)
		  if (i || j) {
		    Ast * b = ast_copy (a);
		    Ast * i1 = ast_find (b, sym_postfix_expression,
					 0, sym_postfix_expression,
					 2, sym_member_identifier,
					 0, sym_generic_identifier,
					 0, sym_IDENTIFIER);
		    ast_terminal (i1)->start[0] = 'x' + i;
		    Ast * i2 = ast_find (b, sym_postfix_expression,
					 2, sym_member_identifier,
					 0, sym_generic_identifier,
					 0, sym_IDENTIFIER);
		    ast_terminal (i2)->start[0] = 'x' + j;
		    ast_list_insert_after (a, b);
		    a = b;
		  }
	    }
	  }

	/**
	Finalize both global and local lists */
	
	if (type == 1) // scalar
	  ast_list_append (n->child[1], sym_initializer, ast_new_empty_scalar (n));
	else if (type == 2) { // vector
	  TranslateData * d = data;
	  ast_list_append (n->child[1], sym_initializer, ast_new_empty_vector (n, d->dimension));
	}
	else { // tensor
	  TranslateData * d = data;
	  ast_list_append (n->child[1], sym_initializer, ast_new_empty_tensor (n, d->dimension));
	}

	Ast * type_name = NN(n, sym_type_name,
			     NN(n, sym_specifier_qualifier_list,
				NN(n, sym_type_specifier,
				   NN(n, sym_types,
				      NA(n, sym_TYPEDEF_NAME,
					 (type == 1 ? "scalar" : type == 2 ? "vector" : "tensor"))))),
			     NN(n, sym_abstract_declarator,
				NN(n, sym_direct_abstract_declarator,
				   NCA(n, "["), NCA(n, "]"))));
	int index = ast_child_index (n);
	Ast * parent = n->parent;
	Ast * unary = NN(parent, sym_unary_expression,
			 NN(parent, sym_postfix_expression,
			    NN(parent, sym_primary_expression,
			       NCA(n, "("),
			       NN(parent, sym_expression_error,
				  NN(parent, sym_expression,
				     ast_attach (ast_new_unary_expression (parent),
						 NN(parent, sym_postfix_expression,
						    NCA(n, "("), type_name, NCA(n, ")"),
						    n)))),
			       NCA(n, ")"))));
	n->sym = sym_postfix_initializer;
	if (parent->sym == sym_forin_arguments)
	  ast_set_child (parent, index,
			 NN(parent, sym_expression,
			    ast_attach (ast_new_cast_expression (parent), unary)));
	else if (parent->sym == sym_initializer_list ||
		 parent->sym == sym_init_declarator)
	  ast_set_child (parent, index,
			 NN(parent, sym_initializer,
			    ast_attach (ast_new_cast_expression (parent), unary)));
	else if (parent->sym == sym_postfix_expression) {
	  assert (index == 3);
	  ast_set_child (parent, 3,
			 NN(parent, sym_cast_expression, unary));
	  Ast * ancestor = ast_ancestor (parent, 2);
	  if (ancestor->sym == sym_cast_expression)
	    ast_set_child (ancestor->parent, ast_child_index (ancestor), parent);
	  else
	    ast_set_child (ast_ancestor (parent, 3), 0, parent);
	  parent->sym = sym_cast_expression;
	}
	else
	  ast_set_child (parent, index,
			 ast_attach (ast_new_cast_expression (parent), unary));
      }
    }
    break;
  }

  case sym_function_definition: {
    
    /**
    ## Function profiling with `trace` */
    
    Ast * trace = ast_schema (n, sym_function_definition,
			      0, sym_function_declaration,
			      0, sym_declaration_specifiers,
			      0, sym_storage_class_specifier,
			      0, sym_TRACE);
    if (trace) {
      TranslateData * d = data;
      ast_hide (ast_terminal (trace));      
      Ast * identifier = ast_find (n, sym_direct_declarator,
				   0, sym_generic_identifier,
				   0, sym_IDENTIFIER);
      Ast * compound_statement = ast_child (n, sym_compound_statement);
      ast_after (compound_statement->child[0],
		 "tracing(\"", ast_terminal (identifier)->start, "\",",
		 ast_file_line (identifier, d->nolineno), ");");
      Ast * end = ast_child (compound_statement, token_symbol ('}'));
      ast_before (end,
		  "end_tracing(\"", ast_terminal (identifier)->start, "\",",
		  ast_file_line (end, d->nolineno), ");");
      if (compound_statement->child[1]->sym == sym_block_item_list) {
	void * adata[] = { n, identifier, data };
	ast_traverse (compound_statement, stack, trace_return, adata);
      }
    }

    /**
    ## Solver initialization and termination. */

    Ast * identifier = ast_function_identifier (n);
    char * init;
    if (identifier && !strcmp (ast_terminal (identifier)->start, "main")) {
      Ast * compound_statement = ast_child (n, sym_compound_statement);
      compound_prepend (compound_statement, ast_new_function_call (n, "_init_solver"));
      compound_append (compound_statement, ast_new_function_call (n, "free_solver"));
    }
    else if (identifier && ast_left_terminal (n)->before &&
	     (init = strstr (ast_left_terminal (n)->before, "@init_solver"))) {
      for (int i = 0; i < 12; i++)
	init[i] = ' ';
      TranslateData * d = data;
      compound_prepend (d->last_events, ast_new_function_call (n, ast_terminal (identifier)->start));
    }
    break;
  }

  }
}

/**
# Traversal functions 

These functions traverse the tree while maintaining a stack of
declared symbols. */

void ast_push_declaration (Stack * stack, Ast * n)
{
  if (n == ast_placeholder)
    return;
  if (n->sym == sym_parameter_type_list ||
      n->sym == sym_struct_declaration_list)
    return; // skip function arguments and struct members
  Ast * identifier = ast_schema (n, sym_direct_declarator,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER);
  if (!identifier)
    identifier = ast_schema (n, sym_enumeration_constant,
			     0, sym_IDENTIFIER);
  if (!identifier && n->sym == sym_struct_or_union_specifier &&
      n->child[2])
    identifier = ast_schema (n, sym_struct_or_union_specifier,
			     1, sym_generic_identifier,
			     0, sym_IDENTIFIER);
  if (identifier)
    stack_push (stack, &identifier);
  if (n->child)
    for (Ast ** c = n->child; *c; c++)
      ast_push_declaration (stack, *c);
}

void ast_pop_scope (Stack * stack, Ast * scope)
{
  if (!scope)
    return;
  while (*((Ast **)stack_pop (stack)) != scope);
}

Ast * ast_push_function_definition (Stack * stack, Ast * declarator)
{
  Ast * identifier = ast_find (declarator, sym_direct_declarator,
			       0, sym_generic_identifier,
			       0, sym_IDENTIFIER);
  stack_push (stack, &identifier);
  stack_push (stack, &declarator);
  Ast * parameters = ast_find (declarator, sym_parameter_list);
  if (parameters)
    ast_push_declaration (stack, parameters);
  return identifier;
}

static void declare_point_variables (Stack * stack)
{
  Ast * list = ast_find (ast_parent (ast_identifier_declaration (stack, "_Variables"),
				     sym_function_definition),
			 sym_block_item_list);
  assert (list);
  foreach_item (list, 1, item) {
    Ast * declaration = ast_schema (item, sym_block_item,
				    0, sym_declaration);
    ast_push_declaration (stack, declaration);
  }
}

Ast * ast_push_declarations (Ast * n, Stack * stack)
{
  switch (n->sym) {

  /**
  These should match the corresponding action/mid-action rules in
  [basilisk.yacc](). */
    
  case sym_function_definition: {
    Ast * declarator = ast_find (n, sym_direct_declarator);
    ast_push_function_definition (stack, declarator);
    if (ast_is_point_function (ast_schema (n, sym_function_definition,
					   0, sym_function_declaration,
					   1, sym_declarator)) &&
	!ast_is_stencil_function (n))
      declare_point_variables (stack);
    return declarator;
  }
    
  case sym_compound_statement:
  case sym_for_declaration_statement:
    stack_push (stack, &n);
    return n;
    
  case sym_forin_declaration_statement:
    stack_push (stack, &n);
    ast_push_declaration (stack, n->child[3]);
    return n;

  case sym_foreach_statement:
    stack_push (stack, &n);
    declare_point_variables (stack);
    return n;

  case sym_macro_statement: {
    Ast * identifier = ast_schema (n, sym_macro_statement,
				   0, sym_function_call,
				   0, sym_postfix_expression,
				   0, sym_primary_expression,
				   0, sym_IDENTIFIER);
    if (!strcmp (ast_terminal (identifier)->start, "map"))
      stack_push (stack, &n);
    return NULL;
  }

  case sym_declaration:
    ast_push_declaration (stack, n);
    return NULL;

  /**
  ## Local boundary conditions */

  case sym_assignment_expression: {
    Ast * array;
    if (n->child[1] && function_scope (n, stack) &&
	(array = ast_schema (n, sym_assignment_expression,
			     0, sym_unary_expression,
			     0, sym_postfix_expression,
			     0, sym_array_access))) {
      const char * typename =
	ast_typedef_name (ast_expression_type (array->child[0], stack, false));
      Ast * member = NULL;
      if ((typename &&
	   (!strcmp (typename, "scalar") ||
	    !strcmp (typename, "vertex scalar"))) ||
	  ((member = ast_schema (array->child[0], sym_postfix_expression,
				 2, sym_member_identifier,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER)) &&
	   (!strcmp (ast_terminal (member)->start, "n") ||
	    !strcmp (ast_terminal (member)->start, "t") ||
	    !strcmp (ast_terminal (member)->start, "r")) &&
	   (typename =
	    ast_typedef_name (ast_expression_type (array->child[0]->child[0],
						   stack, false))) &&
	   (!strcmp (typename, "vector") ||
	    !strcmp (typename, "face vector")))) {
	stack_push (stack, &n);
	declare_point_variables (stack);
	return n;
      }
    }
    return NULL;
  }

  /**
  ## Global boundary conditions */
    
  case sym_boundary_definition: {
    Ast * expr = ast_schema (n, sym_boundary_definition,
			     0, sym_assignment_expression,
			     2, sym_assignment_expression);
    Ast * array = ast_find (n, sym_array_access);
    if (expr && array) {
      stack_push (stack, &n);
      declare_point_variables (stack);
      return n;
    }
    return NULL;
  }

  /**
  ## Point point */

  case sym_direct_declarator:
    if (ast_schema (ast_is_point_point (ast_schema (n, sym_direct_declarator,
						    0, sym_generic_identifier,
						    0, sym_IDENTIFIER)),
		    sym_declaration))
      declare_point_variables (stack);
    return NULL;
    
  }

  return NULL;
}

void ast_traverse (Ast * n, Stack * stack,
		   void func (Ast *, Stack *, void *),
		   void * data)
{
  if (!n || n == ast_placeholder)
    return;

  Ast * scope = ast_push_declarations (n, stack);

  if (n->child)
    for (Ast ** c = n->child; *c; c++)
      ast_traverse (*c, stack, func, data);
  func (n, stack, data);

  ast_pop_scope (stack, scope);
}

/**
# The entry function

Called by [qcc](/src/qcc.c) to trigger the translation. */

static void checks (AstRoot * root, AstRoot * d, TranslateData * data)
{
  CHECK ((Ast *) root, true);
  CHECK ((Ast *) d, true);
  CHECK (data->init_solver, true);
}

void * endfor (FILE * fin, FILE * fout,
	       const char * grid, int dimension,
	       bool nolineno, bool progress, bool catch, bool parallel, bool cpu,
	       FILE * swigfp, char * swigname)
{
  char * buffer = NULL;
  size_t len = 0, maxlen = 0;
  int c;
  while ((c = fgetc (fin)) != EOF) {
    if (len >= maxlen) {
      maxlen += 4096;
      buffer = realloc (buffer, maxlen);      
    }
    buffer[len++] = c;
  }
  if (len >= maxlen) {
    maxlen++;
    buffer = realloc (buffer, maxlen);      
  }
  buffer[len++] = '\0';

  FILE * fp = fopen (BASILISK "/ast/defaults.h", "r");
  assert (fp);
  AstRoot * d = ast_parse_file (fp, NULL);
  fclose (fp);
  
  AstRoot * root = ast_parse (buffer, d);
  free (buffer);
  if (!root) {
    fprintf (stderr, "qcc: error: cannot parse input (missing closing braces?)\n");
    exit (1);
  }
  root->stack = d->stack; d->stack = NULL;
  root->alloc = d->alloc; d->alloc = NULL;

  TranslateData data = {
    .dimension = dimension, .nolineno = nolineno, .parallel = parallel, .cpu = cpu,
    .constants_index = 0, .fields_index = 0, .nboundary = 0,
    // fixme: splitting of events and fields is not used yet
    .init_solver = NULL, .init_events = NULL, .init_fields = NULL,
    .swigname = NULL, .swigdecl = NULL, .swiginit = NULL
  };
  data.constants = calloc (1, sizeof (Field));
  data.swigname = swigfp ? swigname : NULL;

  fp = fopen (BASILISK "/ast/init_solver.h", "r");
  AstRoot * init = ast_parse_file (fp, root);
  fclose (fp);
  data.init_solver = ast_find ((Ast *) init, sym_function_definition);
  assert (data.init_solver);
  str_prepend (ast_left_terminal (data.init_solver)->before, "\n");  
  ast_block_list_append (ast_find ((Ast *)root, sym_translation_unit),
			 sym_external_declaration, data.init_solver);
  data.last_events = 
    ast_find (ast_find (data.init_solver, sym_compound_statement)->child[1],
	      sym_compound_statement);
  assert (data.last_events);
  data.init_fields =
    ast_find (ast_find (data.last_events, sym_compound_statement)->child[1],
	      sym_compound_statement);
  assert (data.init_fields);
  data.init_events =
    ast_find (ast_find (data.init_fields, sym_compound_statement)->child[1],
	      sym_compound_statement);
  assert (data.init_events);
  ast_destroy ((Ast *) init);

  typedef void (* TraverseFunc) (Ast *, Stack *, void *);
  for (TraverseFunc * pass = (TraverseFunc[]){ global_boundaries_and_stencils, translate, maps, stencils, macros, NULL }; *pass; pass++) {
    stack_push (root->stack, &root);
    ast_traverse ((Ast *) root, root->stack, *pass, &data);
    ast_pop_scope (root->stack, (Ast *) root);
    checks (root, d, &data);    
  }
  
  if (data.fields_index) {
    Ast * call_init_solver = ast_find (data.init_solver, sym_function_call);
    char n[10];
    snprintf (n, 9, "%d", data.fields_index);
    char * src = NULL;
    str_append (src, "datasize=", n, "*sizeof(double);");
    Ast * expr = ast_parse_expression (src, root);
    free (src);
    ast_block_list_insert_before2 (ast_parent (call_init_solver, sym_block_item),
				   NN(call_init_solver, sym_statement, expr));
  }

  char methods[strlen(grid) + strlen("_methods") + 1];
  strcpy (methods, grid);
  strcat (methods, "_methods");
  compound_prepend (data.last_events, ast_new_function_call (data.last_events, methods));
  if (catch)
    compound_append (data.last_events, 
		     ast_new_function_call (data.last_events, "catch_fpe"));
  if (progress)
    compound_append (data.last_events, 
		     ast_new_function_call (data.last_events, "last_events"));


  /* SWIG interface */
  if (data.swigname) {
    if (data.swigdecl) {
      fprintf (swigfp,
	       "\n%%{\n"
	       "%s"
	       "%%}\n"
	       "\n"
	       "%s",
	       data.swigdecl,
	       data.swigdecl);
      free (data.swigdecl);
    }
    if (data.swiginit) {
      fprintf (swigfp,
	       "\n"
	       "%%pythoncode %%{\n"
	       "%s"
	       "%%}\n",
	       data.swiginit);
      free (data.swiginit);
    }
    fclose (swigfp);
  }

  checks (root, d, &data);
  
  free (data.constants);
  
  ast_print ((Ast *) root, fout, 0);

  ((Ast *)root)->parent = (Ast *) d;
  return root;
}

bool check_dimensions (AstRoot * root,
		       bool nolineno,
		       int run, FILE * dimensions,
		       int finite, int redundant, int warn, int maxcalls)
{
  Ast * d = ((Ast *)root)->parent;
  ((Ast *)root)->parent = NULL;
  Ast * main = ast_parent (ast_identifier_declaration (root->stack, "main"),
			   sym_function_definition);
  bool ret = true;
  if (main) {
    if (dimensions != stdout)
      ret = ast_check_dimensions (root, main, run >= 0 ? run : 0,
				  maxcalls, dimensions, finite, redundant, !nolineno, warn);  
    else if (run >= 0)
      ast_run (root, main, run, maxcalls, NULL);
  }

  ast_destroy (d);
  ast_destroy ((Ast *) root);
  return ret;
}

/**
# Kernels */

static
void kernel (Ast * n, Stack * stack, void * data)
{
  switch (n->sym) {
    
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

  /**
  ## Dereference of "inout" parameters */

  case sym_unary_operator: {
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
	!strcmp (ast_terminal (type->child[0])->before + strlen (ast_terminal (type->child[0])->before) - 6, "inout ")) {
      free (ast_terminal (n->child[0])->start);
      ast_terminal (n->child[0])->start = strdup("");
    }
    break;
  }
    
  /**
  ## Field assignments 
  
  Kernels often need to know the type of access to fields (i.e. read
  or write). Here we append "_out" to stencil access functions linked
  to assignments (i.e. "write" operations ). */

  case sym_function_call: {
    Ast * identifier = ast_function_call_identifier (n);
    if (!identifier) break;
    AstTerminal * t = ast_terminal (identifier);
    if (!t) break;
    if (!strcmp (t->start, "val")) {
      if (ast_child (ast_parent (n, sym_assignment_expression), sym_assignment_operator))
	str_append (t->start, "_out_");
      break;
    }

    /**
    ## Function calls */
    
    if (!(identifier = ast_identifier_declaration (stack, t->start))) break;
    Ast * definition = ast_parent (identifier, sym_function_definition);
    if (!(identifier = ast_find (definition, sym_direct_declarator,
				 0, sym_direct_declarator,
				 0, sym_generic_identifier,
				 0, sym_IDENTIFIER)) ||
	strcmp (ast_terminal (identifier)->start, t->start)) break;

    Ast * previous = fast_stack_find (data, t->start);
    if (!previous) {
      definition = ast_copy (definition);
      identifier = ast_find (definition, sym_direct_declarator,
			     0, sym_direct_declarator,
			     0, sym_generic_identifier,
			     0, sym_IDENTIFIER);
      stack_push (data, &identifier);
    }

    /**
    Here we check for "inout" arguments (i.e. pointers taken using '&'). */
    
    Ast * parameters = ast_find (ast_schema (definition, sym_function_definition,
					     0, sym_function_declaration,
					     1, sym_declarator),
				 sym_direct_declarator,
				 2, sym_parameter_type_list,
				 0, sym_parameter_list);
    Ast * arguments = ast_schema (n, sym_function_call,
				  2, sym_argument_expression_list);
    if (arguments) {
      while (parameters->child[0]->sym == parameters->sym)
	parameters = parameters->child[0];
      Ast * parameter = ast_child (parameters, sym_parameter_declaration);
      foreach_item_r (arguments, sym_argument_expression_list_item, argument) {
	Ast * ampersand, * type, * pointer;
	if ((ampersand = ast_schema (ast_is_unary_expression (argument->child[0]), sym_unary_expression,
				     0, sym_unary_operator,
				     0, token_symbol ('&'))) &&
	    (pointer = ast_schema (parameter, sym_parameter_declaration,
				   1, sym_declarator,
				   0, sym_pointer,
				   0, token_symbol ('*'))) &&
	    (type = ast_schema (parameter, sym_parameter_declaration,
				0, sym_declaration_specifiers,
				0, sym_type_specifier,
				0, sym_types)) &&
	    (type->child[0]->sym == sym_DOUBLE ||
	     type->child[0]->sym == sym_FLOAT ||
	     type->child[0]->sym == sym_INT)) {
	  if (!previous) {
	    ast_before (type->child[0], "inout ");
	    free (ast_terminal (pointer)->start);
	    ast_terminal (pointer)->start = strdup ("");
	  }
	  free (ast_terminal (ampersand)->start);
	  ast_terminal (ampersand)->start = strdup ("");
	}
	parameters = parameters->parent, parameter = ast_child (parameters, sym_parameter_declaration);
	if (!parameter) break;
      }
    }

    if (!previous)
      ast_traverse (definition, stack, kernel, data);
    break;
  }
    
  }
}

static
void stringify (Ast * n, Ast * argument)
{
  char * s = ast_str_print (n, NULL, 0, 1);
  for (char * i = s; *i; i++) {
    char a[] = "1";
    switch (*i) {
    case '\n': ast_after (argument, "\\n"); break;
    case '\\': ast_after (argument, "\\\\"); break;
    case '"': ast_after (argument, "\\\""); break;
    case '#':
      if (i[-1] == '\n') {
	ast_after (argument, "// #"); break;
      }
      // fall through
    default:  a[0] = *i; ast_after (argument, a); break;
    }
  }
  free (s);
  ast_destroy (n);
}

void ast_kernel (Ast * n, Ast * argument)
{
  AstRoot * root = ast_get_root (n);
  Stack * stack = root->stack;
  stack_push (stack, &n);
  Stack * funcs = stack_new (sizeof (Ast *));
  Ast * statement = ast_copy (ast_child (n, sym_statement));
  ast_traverse (statement, stack, kernel, funcs);
  ast_after (argument, "_(\"");
  Ast ** func;
  while ((func = stack_pop (funcs)))
    stringify (ast_parent (*func, sym_function_definition), argument);
  stack_destroy (funcs);
  ast_after (argument, "\"),_(\"");
  stringify (statement, argument);
  ast_after (argument, "\")");
  ast_pop_scope (stack, n);
}
