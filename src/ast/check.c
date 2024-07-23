#include <stdio.h>
#include <stdlib.h>
#include "ast.h"
#include "symbols.h"

static bool check (Ast * n)
{
  int sym_error = sym_YYerror;
  #include "grammar.h"
  return false;
}

/**
If `stencils` is true the "hack" grammatical expression:

~~~
foreach_statement
├─FOREACH
├─'('
.
├─')'
├─statement
└─foreach_statement
~~~

used for constructing stencils in [/src/ast/translate.c](), is authorized. */

Ast * ast_check_grammar (Ast * n, bool recursive, bool stencils)
{
  if (!n)
    return NULL;
  if (n->child) {
    if (!check (n) &&
	!(stencils && n->sym == sym_foreach_statement &&
	  ((n->child[0] && n->child[0]->sym == sym_FOREACH &&
	    n->child[1] && n->child[1]->sym == token_symbol('(') &&
	    n->child[2] && n->child[2]->sym == token_symbol(')') &&
	    n->child[3] && n->child[3]->sym == sym_statement &&
	    n->child[4] && n->child[4]->sym == sym_foreach_statement && !n->child[5]) ||
	   (n->child[0] && n->child[0]->sym == sym_FOREACH &&
	    n->child[1] && n->child[1]->sym == token_symbol('(') &&
	    n->child[2] && n->child[2]->sym == sym_foreach_parameters &&
	    n->child[3] && n->child[3]->sym == token_symbol(')') &&
	    n->child[4] && n->child[4]->sym == sym_statement &&
	    n->child[5] && n->child[5]->sym == sym_foreach_statement && !n->child[6])))) {	
      fprintf (stderr, "\ngrammatical error:\n");
      ast_print_tree (n, stderr, 0, 0, 10);
      abort ();
    }
    if (recursive)
      for (Ast ** c = n->child; *c; c++) {
	assert ((*c)->parent == n);
	ast_check_grammar (*c, true, stencils);
      }
  }
  else {
    AstTerminal * t = ast_terminal (n);
    assert (t->file && t->line);
  }
  return n;
}
