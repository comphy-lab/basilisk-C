/**
# Macros

The [C preprocessor](https://gcc.gnu.org/onlinedocs/cpp/Macros.html)
(also known as CPP) macros are familiar to C programmers (and are
briefly introduced in the [Tutorial](/Tutorial#using-macros)). While
they are very useful, their shortcomings are also well-known, for
example:

* multiline macros are hard to write, read and debug,
* macros only deal with text i.e. they do not enforce the C syntax/grammar 
  (which also makes them hard to write and debug),
* macro statements cannot 'return' a value,
* while macros can be used to implement simple iterators, the syntax is cumbersome.

As an illustration, consider the following simple iteration macro

~~~c
#define iterator(start, end, index, expr) do {    \
  for (int index = start; index <= end; index++)  \
    expr					  \
  printf ("do something after the loop\n");	  \
} while(0)

int main() {
  iterator (0, 10, i, { printf ("%d\n", i); });
}
~~~

Using Basilisk macros, this would be written

~~~literatec
macro iterator (int start, int end, int index) {
  for (int index = start; index <= end; index++)
    {...}
  printf ("do something after the loop\n");
}

int main() {
  iterator (0, 10, i)
    printf ("%d\n", i);
}
~~~

The syntax/grammar of the macro definition is almost identical to the
definition of a standard C function and is checked as such by the
compiler. Using the macro also respects the C grammar (i.e. that of
standard C iterators: `for`, `while` etc.).

Two new reserved keywords have been introduced: 

* `macro` which indicates that what follows is a macro definition, and 
* `{...}` which is expanded to the statement used when the macro is
  called (i.e. to `printf ("%d\n", i);` in this example).

## Macro expansion

The expansion of the arguments passed to Basilisk macros is what
makes them different from functions. To illustrate this, consider the
following example

~~~literatec
macro initialize (scalar s, double expr)
{
  foreach()
    s[] = expr;
}

int main() {
  init_grid (16);
  scalar a[];
  initialize (a, sqrt (x*x + y*y));
}
~~~

It is clear that if `initialize` was a function (i.e. `macro
initialize...` would just be replaced by `void initialize...`), this
would not compile, because `x` and `y` are not defined when
`initialize` is called in `main`. This works because Basilisk macros
are expanded before compilation (like their CPP counterpart). After
macro expansion the code would read

~~~literatec
...
int main() {
  init_grid (16);
  scalar a[];
  foreach()
    a[] = sqrt (x*x + y*y);
}
~~~

Note also that in this example the macro does not use the "ellipsis
macro" operator `{...}`.

## Macros returning a value

As mentioned above, standard C macros cannot 'return' a value like
functions do. Basilisk macros can. There are two types of "return
macros" in Basilisk.

### Inlined return macros

They are simple macros which behave very much like their standard CPP
counterpart. The body of the macro must consist of a single return
statement e.g.

~~~literatec
macro number sq (number x) {
  return x*x;
}
~~~

The returned expression (`x*x` here) is expanded directly where the
macro is called (i.e. it is "inlined"). Note that, unlike what the
equivalent C function (not macro) would do, the arguments and the
returned value are not implicitly cast to the specified type (`number`
here). So the sq() macro will preserve the type of its arguments
i.e. it will work with any numeric type (`int`, `double`, `float`
etc.).

### Complex return macros

Complex return macros return a value which is more complex than a
simple expression. They cannot be inlined, because they are composed
of several statements. Instead they are expanded as additional
statements, which are inserted before the calling statement. In
contrast with inlined macros, the type of the return value is
automatically cast to the type given in the macro definition.

## Optional arguments

Since macros share their syntax with that of functions, they can also
define [optional arguments](/Basilisk
C#namedoptional-arguments-in-function-calls). The default value `None`
is specific to macros and will be expanded to nothing.

## Breaking out of macro iterators

By default trying to `break` out of a macro iterator (like the
`iterator` example above) will cause an error. Breaking must
explicitly be set when defining the macro, for example

~~~literatec
macro iterator (int start, int end, int index, break = break) {
  ...
}
~~~

where breaking uses the standard `break` statement, or

~~~literatec
macro iterator (int start, int end, int index, break = (index = end + 1)) {
  ...
}
~~~

where the given expression will be expanded in replacement of the
`break` statement.

## Overloading

Unlike for C functions, one is allowed to define the same macro
several times. The latest definition (in order of appearance in the
source file) "overloads" the previous ones.

It is sometimes useful to break this rule to specify a "default"
definition which will only apply if no other definition is
present. This can be done by adding the `auto` "storage class
specifier" to the macro definition.

## Inheritance

Macros implement a simple mechanism for "inheritance". A "recursive"
definition of macros is possible, which allows to include the
expansion of a previous definition of a macro into the new
(overloaded) definition.

## Non-casting of arguments

The comment made above regarding the non-type-casting of the return
values of inlined macros, also applies to the arguments of any
macro. This is different from what happens for standard C functions.

## Postmacros

Macro definitions can use an "order" parameter which define their
expansion priority when preprocessing by [qcc](/src/qcc.c). This is
used to control which macros are expanded and thus seen or not by the
[interpreter](interpreter/interpreter.c) (only `macro` are expanded)
or by [computation kernels](kernels.c) (`macro` and `macro1` are expanded).

*This is a low-level functionality which must not be used for "user"
macros.*

## See also

* [Macro tests](/src/test/macro.c)

# Implementation */

typedef struct {
  Ast * statement, * arguments, * parameters, * call, * label, * initial, * scope;
  Stack * sparameters;
  int * returnindex;  
  bool nolineno, complex_call;
  int postmacros;
} MacroReplacement;

static Ast * argument_value (Ast * identifier, Stack * stack, const MacroReplacement * r)
{
  Ast * decl = ast_identifier_declaration (r->sparameters, ast_terminal (identifier)->start);
  if (!decl)
    return NULL;
  Ast * value = NULL;
  if (r->arguments) {
    Ast * parent = ast_parent (decl, sym_parameter_declaration);
    Ast * parameters = r->parameters->child[0];
    while (parameters && parameters->child[0]->sym == parameters->sym)
      parameters = parameters->child[0];
    Ast * arguments = r->arguments;
    foreach_item_r (arguments, sym_argument_expression_list_item, argument) {
      Ast * parameter = ast_child (parameters, sym_parameter_declaration);
      parameters = parameters->parent;
      assert (parameter);
      if (parameter == parent) {
	value = argument;
	break;
      }
    }
  }
  if (!value) {
    AstTerminal * t = ast_left_terminal (r->call);
    fprintf (stderr, "%s:%d: error: missing '%s' macro parameter\n",
	     t->file, t->line, ast_terminal (identifier)->start);
    exit (1);
  }
  return value;
}

static Ast * remove_leading_spaces (Ast * n)
{
  AstTerminal * t = ast_left_terminal (n);
  for (char * s = t->before; s && *s != '\0'; s++)
    if (!strchr (" \t\n", *s))
      return n;
  free (t->before); t->before = NULL;
  return n;
}

static void replace_arguments (Ast * n, Stack * stack, void * data)
{
  const MacroReplacement * r = data;
  Ast * identifier;
  if ((identifier = ast_schema (n, sym_postfix_expression,
				0, sym_primary_expression,
				0, sym_IDENTIFIER))) {
    Ast * value = argument_value (identifier, stack, r);
    if (value) {
      if (value->child[0]->sym == sym_assignment_expression) {
	Ast * primary = ast_is_simple_expression (value->child[0]);
	if (primary) {
	  primary = remove_leading_spaces (ast_copy (primary->parent));
	  ast_set_line (primary, ast_terminal (identifier), true);
	  Ast * identifier = ast_schema (primary, sym_primary_expression,
					 0, sym_IDENTIFIER);
	  if (identifier && !strcmp (ast_terminal (identifier)->start, "None"))
	    ast_terminal (identifier)->start[0] = '\0';
	}
	else
	  primary = NN(n, sym_primary_expression,
		       NCA(n, "("),
		       NN(n, sym_expression_error,
			  NN(n, sym_expression,
			     remove_leading_spaces (ast_copy (value->child[0])))),
		       NCA(n, ")"));
	ast_replace_child (n, 0, primary);
      }
      else { // postfix_initializer
	assert (value->child[0]->sym == sym_postfix_initializer);
	Ast * parent = ast_ancestor (n, 3);
	assert (ast_schema (parent, sym_cast_expression,
			    1, sym_type_name));
	parent->sym = sym_postfix_expression;
	ast_replace_child (parent, 3, value->child[0]);
	int index = ast_child_index (parent);
	Ast * grandparent = parent->parent;
	ast_set_child (grandparent, index,
		       NN(grandparent, sym_cast_expression,
			  NN(grandparent, sym_unary_expression,
			     parent)));
      }
    }
    else if (r->initial) {
      if (!strcmp (ast_terminal (identifier)->start, "S__FILE__")) {
	char * filename = NULL;
	str_append (filename, "\"", ast_left_terminal (r->initial)->file, "\"");
	ast_replace_child (identifier->parent, 0, NN(n, sym_string,
						     NB(n, sym_STRING_LITERAL, filename)));
	free (filename);
      }
      else if (!strcmp (ast_terminal (identifier)->start, "S_LINENO")) {
	int line = r->nolineno ? 0 : ast_left_terminal (r->initial)->line;
	char sline[20];
	snprintf (sline, 19, "%d", line);
	ast_replace_child (identifier->parent, 0, NN(n, sym_constant,
						     NB(n, sym_I_CONSTANT, sline)));
      }
    }
  }
  else if ((identifier = ast_schema (n, sym_direct_declarator,
				     0, sym_generic_identifier,
				     0, sym_IDENTIFIER))) {
    Ast * value = argument_value (identifier, stack, r);
    if (value) {
      assert (value->child[0]->sym == sym_assignment_expression); // does not know yet how to deal with postfix_initializer
      Ast * replacement = ast_is_identifier_expression (ast_schema (value, sym_argument_expression_list_item,
								    0, sym_assignment_expression));
      if (!replacement) {
	AstTerminal * t = ast_left_terminal (value);
	fprintf (stderr, "%s:%d: error: macro argument '%s' must be a simple identifier\n",
		 t->file, t->line, ast_terminal (identifier)->start);
	t = ast_terminal (identifier);
	fprintf (stderr, "%s:%d: error: because it is used here as a declarator\n",
		 t->file, t->line);
	exit (1);
      }
      free (ast_terminal (identifier)->start);
      ast_terminal (identifier)->start = strdup (ast_terminal (replacement)->start);
    }
  }
}

static void replace_ellipsis (Ast * n, Stack * stack, void * data)
{
  if (ast_schema (ast_child (n, sym_statement), sym_statement,
		  0, sym_basilisk_statements,
		  0, sym_ELLIPSIS_MACRO)) {
    Ast * child = ast_child (n, sym_statement);
    MacroReplacement * r = data;
    ast_set_child (n, ast_child_index (child), ast_copy (r->statement));
    ast_destroy (child);
  }
}

Ast * ast_is_macro_declaration (const Ast * function_declaration)
{
  return ast_find (ast_schema (function_declaration, sym_function_declaration,
			       0, sym_declaration_specifiers),
		   sym_declaration_specifiers,
		   0, sym_storage_class_specifier,
		   0, sym_MACRODEF) ?
    ast_schema (function_declaration, sym_function_declaration,
		1, sym_declarator,
		0, sym_direct_declarator,
		0, sym_direct_declarator,
		0, sym_generic_identifier,
		0, sym_IDENTIFIER) : NULL;
}

static void replace_break (Ast * n, Ast * breaking, Ast * parent)
{
  if (ast_schema (n, sym_statement,
		  0, sym_jump_statement,
		  0, sym_BREAK)) {
    Ast * loop = n->parent;
    while (loop != parent &&
	   loop->sym != sym_forin_declaration_statement &&
	   loop->sym != sym_forin_statement &&
	   loop->sym != sym_iteration_statement &&
	   (loop->sym != sym_selection_statement ||
	    loop->child[0]->sym != sym_SWITCH) &&
	   !ast_is_foreach_statement (loop))
      loop = loop->parent;
    if (loop == parent) {
      if (breaking) {
	if (breaking->sym == sym_BREAK)
	  ; // do nothing
	else {
	  assert (breaking->sym == sym_assignment_expression);
	  n->child[0]->sym = sym_expression_statement;
	  ast_replace_child (n->child[0], 0,
			     NN(n, sym_expression,
				ast_copy (breaking)));
	}
      }
      else {
	Ast * breaking = ast_schema (n, sym_statement,
				     0, sym_jump_statement,
				     0, sym_BREAK);
	Ast * identifier = ast_schema (ast_parent (breaking, sym_macro_statement), sym_macro_statement,
				       0, sym_MACRO);
	fprintf (stderr, "%s:%d: error: cannot break out of macro '%s'\n",
		 ast_terminal (breaking)->file, ast_terminal (breaking)->line,
		 ast_terminal (identifier)->start);
	exit (1);
      }
    }
  }
  else if (n->child)
    for (Ast ** c = n->child; *c; c++)
      replace_break (*c, breaking, parent);
}

static void replace_return (Ast * n, Stack * stack, void * data)
{
  if (n->sym == sym_RETURN) {
    Ast * expr = ast_schema (n->parent, sym_jump_statement,
			     1, sym_expression);

    n->sym = sym_GOTO;
    free (ast_terminal (n)->start);
    ast_terminal (n)->start = strdup ("goto");
    MacroReplacement * r = data;
    if (!r->label) {
      char rindex[25]; snprintf (rindex, 24, "_return_%d", (*r->returnindex)++);
      r->label = NN(n, sym_generic_identifier,
		    NA(n, sym_IDENTIFIER, rindex));
    }
    Ast * label = ast_copy (r->label);
    ast_before (label, " ");
    ast_set_line (label, ast_terminal (n), true);
    if (expr) // 'return expr;' replaced with 'goto label;'
      ast_set_child (n->parent, 1, label);
    else // 'return;' replaced with 'goto label;'
      ast_new_children (n->parent, n->parent->child[0], label, n->parent->child[1]);
      
    Ast * statement = ast_parent (n, sym_statement);
    Ast * parent = statement->parent;
    int index = ast_child_index (statement);

    Ast * compound;
    if (r->call->sym == sym_function_call) {
      if (!expr) {
	fprintf (stderr, "%s:%d: error: 'void' return value in a macro returning non-void\n",
		 ast_terminal (n)->file, ast_terminal (n)->line);
	exit (1);
      }
      char * val = strdup (ast_terminal (ast_child(label, sym_IDENTIFIER))->start);
      str_append (val, "val");
      AstTerminal * open = NCB(parent, "{"), * close = NCB(parent, "}");
      Ast * assign = NN(n, sym_statement,
			NN(n, sym_expression_statement,
			   NN(n, sym_expression,
			      NN(n, sym_assignment_expression,
				 NN(n, sym_unary_expression,
				    NN(n, sym_postfix_expression,
				       NN(n, sym_primary_expression,
					  NB(n, sym_IDENTIFIER, val)))),
				 NN(n, sym_assignment_operator,
				    NCB(n, "=")),
				 ast_child (expr, sym_assignment_expression))),	  
			   NCB(n, ";")));
      free (val);

      compound = NN(parent, sym_statement,
		    NN(parent, sym_compound_statement,
		       open,
		       NN(parent, sym_block_item_list,
			  NN(parent, sym_block_item_list,
			     NN(parent, sym_block_item,
				assign)),
			  NN(parent, sym_block_item,
			     statement)),
		       close));
    }
    else
      compound = statement;

    ast_set_child (parent, index, compound);
  }
}

static
Ast * get_macro_definition (Stack * stack, const Ast * identifier, Ast * scope)
{
  Ast * decl = scope ?
    ast_identifier_declaration_from_to (stack, ast_terminal (identifier)->start, scope, NULL) :
    ast_identifier_declaration (stack, ast_terminal (identifier)->start);
  Ast * macro_definition = ast_schema (ast_ancestor (decl, 6), sym_function_definition);

  if (!macro_definition) {
    fprintf (stderr, "%s:%d: error: undefined macro '%s'\n",
	     ast_terminal (identifier)->file, ast_terminal (identifier)->line, ast_terminal (identifier)->start);
    exit (1);
  }

  Ast * non_auto_macro = macro_definition;
  while (ast_find (ast_schema (non_auto_macro, sym_function_definition,
			       0, sym_function_declaration),
		   sym_declaration_specifiers,
		   0, sym_storage_class_specifier,
		   0, sym_AUTO)) {
    Ast * other = ast_identifier_declaration_from_to (stack, ast_terminal (identifier)->start, decl, NULL);
    if (other == decl)
      break;
    non_auto_macro = ast_schema (ast_ancestor (other, 6), sym_function_definition);
    decl = other;
  }
  if (non_auto_macro)
    macro_definition = non_auto_macro;
  
  if (macro_definition == ast_parent (identifier, sym_function_definition)) { // This is a "recursive" macro
    Ast * other = ast_identifier_declaration_from_to (stack, ast_terminal (identifier)->start, decl, NULL);
    if (other == decl) {
      fprintf (stderr, "%s:%d: error: cannot find ancestor of recursive macro '%s'\n",
	       ast_terminal (decl)->file, ast_terminal (decl)->line, ast_terminal (decl)->start);
      exit (1);
    }
    macro_definition = ast_schema (ast_ancestor (other, 6), sym_function_definition);
  }  

  return macro_definition;
}

static void replace_macros (Ast * n, Stack * stack, void * data)
{
  if (n->sym == sym_statement || n->sym == sym_function_call) {
    MacroReplacement * r = data;
    ast_macro_replacement (n, r->initial, stack, r->nolineno, r->postmacros, true, r->returnindex, r->scope);
  }
}

void ast_macro_replacement (Ast * statement, Ast * initial, Stack * stack,
			    bool nolineno, int postmacros, bool expand_definitions,
			    int * return_macro_index, Ast * scope)
{
  Ast * identifier, * macro_statement = ast_schema (statement, sym_statement,
						    0, sym_basilisk_statements,
						    0, sym_macro_statement);
  if (macro_statement)
    identifier = ast_child (macro_statement, sym_MACRO);
  else if ((identifier = ast_schema (statement, sym_function_call,
				     0, sym_postfix_expression,
				     0, sym_primary_expression,
				     0, sym_MACRO)))
    macro_statement = statement;
  if (!macro_statement)
    return;

  if (!strcmp (ast_terminal (identifier)->start, "OMP_PARALLEL"))
    return;

  if (!expand_definitions) {
    Ast * definition = ast_parent (statement, sym_function_definition);
    if (definition && ast_is_macro_declaration (definition->child[0]))
      return;
  }
    
  if (!strcmp (ast_terminal (identifier)->start, "foreach_block") &&
      (inforeach (statement) || point_declaration (stack)))
    str_append (ast_terminal (identifier)->start, "_inner");

  Ast * macro_definition = get_macro_definition (stack, identifier, scope);
  if (postmacros >= 0) {
    Ast * macrodef = ast_find (ast_schema (macro_definition, sym_function_definition,
					   0, sym_function_declaration),
			       sym_declaration_specifiers,
			       0, sym_storage_class_specifier,
			       0, sym_MACRODEF);
    assert (macrodef);
    const char * suffix = ast_terminal (macrodef)->start + 5;
    if (atoi (suffix) > postmacros)
      return;
  }
  
  optional_arguments (macro_statement, stack);

  MacroReplacement r = {
    .statement = ast_child (macro_statement, sym_statement),
    .arguments = ast_child (macro_statement, sym_argument_expression_list),
    .parameters = ast_schema (macro_definition, sym_function_definition,
			      0, sym_function_declaration,
			      1, sym_declarator,
			      0, sym_direct_declarator,
			      2, sym_parameter_type_list),
    .call = macro_statement,
    .nolineno = nolineno,
    .complex_call = !ast_schema (ast_find (macro_definition, sym_compound_statement),
				 sym_compound_statement,
				 1, sym_block_item_list,
				 0, sym_block_item,
				 0, sym_statement,
				 0, sym_jump_statement,
				 0, sym_RETURN),
    .returnindex = return_macro_index,
    .postmacros = postmacros,
    .scope = scope
  };

  if (r.parameters) {
    r.sparameters = stack_new (sizeof (Ast *));
    ast_push_declaration (r.sparameters, ast_child (r.parameters, sym_parameter_list));
  }
  
  if (initial) {
    Ast * definition = ast_parent (initial, sym_function_definition);
    if (definition && ast_is_macro_declaration (definition->child[0]))
      r.initial = NULL;
    else
      r.initial = initial;
  }

  int na = 0;
  if (r.arguments)
    foreach_item (r.arguments, 2, argument) na++;
  int np = 0;
  if (r.parameters) {
    Ast * list = r.parameters->child[0];
    foreach_item (list, 2, parameter)
      if (!ast_schema (parameter, sym_parameter_declaration,
		       0, sym_BREAK))
	np++;
  }

  if (na != np) {
    AstTerminal * t = ast_terminal (identifier);
    fprintf (stderr, "%s:%d: error: too %s arguments for macro '%s'\n",
	     t->file, t->line, na > np ? "many" : "few", ast_terminal (identifier)->start);
    t = ast_left_terminal (macro_definition);
    fprintf (stderr, "%s:%d: error: macro '%s' is defined here\n", t->file, t->line, ast_terminal (identifier)->start);
#if 0
    ast_print_tree (macro_definition->child[0], stderr, 0, 0, -1);
    ast_print (r.call, stderr, 0);
#endif
    exit (1);
  }
  
#if 0
  fprintf (stderr, "***** replacing macro: \n");
  ast_print (macro_statement, stderr, 0);
  fprintf (stderr, "\n===== with: \n");
  ast_print (ast_child (macro_definition, sym_function_declaration), stderr, 0);
  ast_print (ast_child (macro_definition, sym_compound_statement), stderr, 0);
  fputc ('\n', stderr);
#endif
  
  /**
  Replace 'break' with its macro definition (if it exists). */

  if (r.statement) {
    Ast * breaking = NULL;
    if (r.parameters && !(breaking = ast_find (r.parameters, sym_parameter_declaration,
					       2, sym_initializer,
					       0, sym_assignment_expression)))
      breaking = ast_find (r.parameters, sym_parameter_declaration,
			   2, sym_BREAK);
    Ast * copy = breaking;
    if (breaking) {
      if (r.parameters) {
	copy = ast_copy (breaking);
	stack_push (stack, &copy);
	ast_traverse (copy, stack, replace_arguments, &r);
	ast_pop_scope (stack, copy);
      }
    }
    replace_break (r.statement, copy, r.statement);
    if (copy != breaking)
      ast_destroy (copy);
  }

  Ast * copy = ast_copy (ast_find (macro_definition, sym_compound_statement));
  stack_push (stack, &copy);
  if (r.parameters) {
    ast_traverse (copy, stack, replace_arguments, &r);
    stack_destroy (r.sparameters);
  }
  if (r.complex_call)
    ast_traverse (copy, stack, replace_return, &r);
  if (r.statement)
    ast_traverse (copy, stack, replace_ellipsis, &r);
  ast_pop_scope (stack, copy);
  str_prepend (ast_left_terminal (copy)->before, ast_left_terminal (macro_statement)->before);

  stack_push (stack, &copy);
  ast_traverse (copy, stack, replace_macros, &r);
  ast_pop_scope (stack, copy);

  /**
  Return label */
  
  if (r.label) {
    ast_set_line (r.label, ast_right_terminal (copy), true);
    ast_block_list_append (ast_schema (copy, sym_compound_statement,
				       1, sym_block_item_list), sym_block_item,
			   NN(copy, sym_statement,
			      NN(copy, sym_labeled_statement,
				 r.label,
				 NCA(copy, ":"),
				 NN(copy, sym_statement,
				    NN(copy, sym_expression_statement,
				       NCA(copy, ";"))))));
  }
  
  /**
  Statement */
  
  if (statement->sym == sym_statement) {
    if (statement->parent->sym == sym_block_item) {
      Ast * list = ast_schema (copy, sym_compound_statement,
			       1, sym_block_item_list);
      Ast * parent = ast_ancestor (statement, 2);
      assert (parent->sym == sym_block_item_list);
      foreach_item (list, 1, item) {
	if (ast_child_index (item))
	  ast_block_list_append (parent, sym_block_item, item->child[0]);
	else {
	  ast_before (item, ast_left_terminal (statement)->before);
	  ast_set_child (parent, parent->child[1] ? 1 : 0, item);
	}
      }
      ast_destroy (copy);
    }
    else
      ast_set_child (statement, 0, copy);
    ast_destroy (macro_statement);
  }
  
  /**
  Function call */
  
  else {
    assert (statement->sym == sym_function_call);
    Ast * jump = ast_schema (copy, sym_compound_statement,
			     1, sym_block_item_list,
			     0, sym_block_item,
			     0, sym_statement,
			     0, sym_jump_statement);
    if (ast_schema (jump, sym_jump_statement,
		    0, sym_RETURN)) {
      assert (!r.complex_call);

      /**
      This is a simple return macro i.e. 'macro int func(...){ return ...; }'. */
      
      Ast * expr = ast_schema (jump, sym_jump_statement,
			       1, sym_expression);
      if (!expr) {
	AstTerminal * t = ast_terminal (identifier);
	fprintf (stderr, "%s:%d: error: using value of macro '%s' returning void\n",
		 t->file, t->line, ast_terminal (identifier)->start);
	t = ast_terminal (ast_schema (jump, sym_jump_statement,
				      0, sym_RETURN));
	fprintf (stderr, "%s:%d: error: return value is defined here\n", t->file, t->line);
	exit (1);
      }
      AstTerminal * o = NCB(statement, "("), * c = NCB(statement, ")");
      statement->sym = sym_primary_expression;
      for (Ast ** c = statement->child; *c; c++)
	ast_destroy (*c);
      ast_set_line (expr, o, true);
      ast_new_children (statement, o,
			NN(expr, sym_expression_error,
			   expr), c);
      ast_destroy (copy);
    }
    else {

      /**
      This is a complex return macro. */
      
      assert (r.complex_call);
      assert (r.label);
      
      Ast * returntype = ast_schema (macro_definition, sym_function_definition,
				     0, sym_function_declaration,
				     0, sym_declaration_specifiers,
				     1, sym_declaration_specifiers);
      if (!returntype) {
	fprintf (stderr, "%s:%d: error: the type of a macro returning a value must be defined\n",
		 ast_left_terminal (macro_definition)->file,
		 ast_left_terminal (macro_definition)->line);
	exit (1);
      }

      char * label = strdup (ast_terminal (ast_schema (r.label, sym_generic_identifier,
						       0, sym_IDENTIFIER))->start);
      str_append (label, "val");
      Ast * declaration = NN(copy, sym_declaration,
			     ast_copy (returntype),
			     NN(copy, sym_init_declarator_list,
				NN(copy, sym_init_declarator,
				   NN(copy, sym_declarator,
				      NN(copy, sym_direct_declarator,
					 NN(copy, sym_generic_identifier,
					    NB(copy, sym_IDENTIFIER, label)))))),
			     NCB(copy, ";"));

      AstTerminal * val =  NB(statement, sym_IDENTIFIER, label);
      for (Ast ** c = statement->child; *c; c++)
	ast_destroy (*c), *c = NULL;
      statement->sym = sym_primary_expression;
      ast_set_child (statement, 0, (Ast *) val);
      
      Ast * parent = statement->parent;
      while (parent->sym != sym_statement && parent->sym != sym_declaration)
	parent = parent->parent;
      
      Ast * gp = parent->parent;
      if (gp->sym == sym_block_item) {
	Ast * s = NN(gp, sym_statement, copy);
	ast_block_list_insert_before2 (gp, s);
	assert (s->parent->sym == sym_block_item);
	ast_block_list_insert_before2 (s->parent, declaration);
      }
      else {
	int index = ast_child_index (parent);
	AstTerminal * open = NCB(gp, "{"), * close = NCA(parent, "}");
	Ast * compound = NN(gp, sym_statement,
			    NN(gp, sym_compound_statement,
			       open,
			       NN(gp, sym_block_item_list,
				  NN(gp, sym_block_item_list,
				     NN(gp, sym_block_item_list,
					NN(gp, sym_block_item,
					   declaration)),
				     NN(gp, sym_block_item,
					NN(gp, sym_statement,
					   copy))),
				  NN(gp, sym_block_item,
				     parent)),
			       close));
	ast_set_child (gp, index, compound);
      }

      free (label);
    }
  }
}
