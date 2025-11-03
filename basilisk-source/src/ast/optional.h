/**
## Optional arguments for function calls

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

static void optional_arguments (Ast * call, Stack * stack)
{
  AstTerminal * t = ast_terminal (call->sym == sym_function_call ? ast_function_call_identifier (call) :
				  ast_schema (call, sym_macro_statement, 0, sym_MACRO));
  Ast * type = ast_identifier_declaration (stack, t->start);
  if (type && !ast_schema (ast_ancestor (type, 3), sym_declarator,
			   0, sym_pointer)) {
    AstTerminal * tname = t;
    type = ast_parent (type, sym_declarator);
    while (type->sym != sym_declaration &&
	   type->sym != sym_function_definition)
      type = type->parent;
    Ast * parameters = ast_find (type, sym_parameter_list);
	
    /**
    Obsolete optional arguments syntax using 'struct ...' parameters. */

    Ast * struct_name = obsolete_function_declaration (type);
    if (struct_name) {
      Ast * arguments = ast_find (call, sym_argument_expression_list);
      if (!arguments) {
	Ast * expr = ast_parse_expression ("func((struct Name){0});",
					   ast_get_root (call));
	Ast * list = ast_find (expr, sym_argument_expression_list);
	AstTerminal * t = ast_terminal (ast_find (list, sym_IDENTIFIER));
	free (t->start);
	t->start = strdup (ast_terminal (struct_name)->start);
	ast_set_line (list, ast_terminal (call->child[1]), false);
	ast_new_children (call, call->child[0], call->child[1],
			  ast_placeholder,
			  call->child[2]);
	ast_replace_child (call, 2, list);
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
					     ast_get_root (call));
	  Ast * list = ast_find (expr, sym_argument_expression_list);
	  AstTerminal * t = ast_terminal (ast_find (list, sym_IDENTIFIER));
	  free (t->start);
	  t->start = strdup (ast_terminal (struct_name)->start);
	  Ast * initializer_list = ast_initializer_list (arguments);
	  ast_replace (list, "a", initializer_list);
	  ast_replace_child (call, 2, list);
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

    else if (parameters &&
	     (ast_find (parameters, sym_parameter_declaration,
			3, sym_initializer) ||
	      ast_find (call, sym_argument_expression_list_item,
			0, sym_assignment_expression,
			1, sym_assignment_operator))) {
      parameters = ast_copy (parameters); // fixme: memory is leaking from here
      Ast * breaking = ast_find (parameters, sym_parameter_declaration,
				 0, sym_BREAK);
      if (breaking) {
	if (ast_ancestor (breaking, 2) != parameters) {
	  fprintf (stderr, "%s:%d: error: 'break' must be the last parameter\n",
		   ast_terminal (breaking)->file, ast_terminal (breaking)->line);
	  exit (1);
	}
	if (ast_child_index (breaking->parent) == 0) {
	  ast_destroy (parameters);
	  parameters = NULL;
	}
	else {
	  ast_destroy (parameters->child[1]);
	  ast_destroy (parameters->child[2]);
	  parameters->child[1] = parameters->child[2] = NULL;
	  parameters = parameters->child[0];
	}
      }
      Ast * parameters1 = parameters;
      while (parameters && parameters->child[0]->sym == parameters->sym)
	parameters = parameters->child[0];
      Ast * arguments = ast_child (call, sym_argument_expression_list);
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
	      fprintf (stderr, "%s:%d: error: unknown %s parameter '%s'\n",
		       t->file, t->line,
		       call->sym == sym_macro_statement ? "macro" : "function",
		       ast_terminal (identifier)->start);
	      t = ast_left_terminal (type);
	      fprintf (stderr, "%s:%d: error: %s '%s' is defined here\n",
		       t->file, t->line,
		       call->sym == sym_macro_statement ? "macro" : "function",
		       tname->start);
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
	  if (!parameter) {
	    fprintf (stderr, "%s:%d: error: too many arguments when calling '%s'\n",
		     t->file, t->line,
		     t->start);
	    AstTerminal * tt = ast_left_terminal (type);
	    fprintf (stderr, "%s:%d: error: %s '%s' is defined here\n",
		     tt->file, tt->line,
		     call->sym == sym_macro_statement ? "macro" : "function",
		     tname->start);
	    exit (1);
	  }
	  if (ast_schema (parameter, sym_parameter_declaration,
			  3, sym_initializer))
	    ast_set_child (parameter->child[3], 0, argument);
	  else if (ast_schema (parameter, sym_parameter_declaration,
			       1, sym_declarator))
	    ast_new_children (parameter,
			      parameter->child[0],
			      parameter->child[1],
			      NCA(call, "="),
			      NN(call, sym_initializer,
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
	  AstTerminal * t = ast_left_terminal (call);
	  fprintf (stderr, "%s:%d: error: missing compulsory parameter '%s' in %s call\n",
		   t->file, t->line, ast_terminal (id)->start,
		   call->sym == sym_macro_statement ? "macro" : "function");
	  exit (1);
	}
	Ast * assign = ast_schema (initializer, sym_initializer,
				   0, sym_assignment_expression);
	if (assign)
	  ast_new_children (parameter, assign); // fixme: possible memory leak here
	else if (ast_schema (initializer, sym_initializer,
			     0, sym_reduction_list))
	  ast_new_children (parameter, ast_schema (initializer, sym_initializer,
						   0, sym_reduction_list)); // fixme: possible memory leak here
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
	    NN(call, sym_type_name,
	       NN(call, sym_specifier_qualifier_list,
		  type_specifier),
	       abstract) :
	    NN(call, sym_type_name,
	       NN(call, sym_specifier_qualifier_list,
		  type_specifier));
	  ast_new_children (parameter, ast_attach
			    (ast_new_unary_expression (parameter),
			     NN(call, sym_postfix_expression,
				ob, type_name, cb,
				initializer)));
	}
	parameter->sym = sym_argument_expression_list_item;
	parameter->parent->sym = sym_argument_expression_list;
      }
	
      /**
      Here we replace the arguments of the called function or macro with the completed arguments. */

      if (call->child[2]->sym == sym_argument_expression_list)
	ast_set_child (call, 2, parameters1);
      else if (call->sym == sym_macro_statement)
	ast_new_children (call, call->child[0], call->child[1], parameters1, call->child[2], call->child[3]);
      else // function call
	ast_new_children (call, call->child[0], call->child[1], parameters1, call->child[2]);
    }	  
  }
}
