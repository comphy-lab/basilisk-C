#include "ast.h"
#include "symbols.h"

int main (int argc, char * argv[])
{
  if (argc < 2) {
    fprintf (stderr, "usage: %s 'code' [OPTIONS]\n", argv[0]);
    return 1;
  }
  Ast * n = (Ast *) ast_parse (argv[1], NULL);
  if (!n)
    n = ast_parse_expression (argv[1], NULL);
  if (!n) {
    fprintf (stderr, "%s: error: could not parse code\n", argv[0]);
    return 1;
  }
  if (argc > 2) {
    ast_print_constructor (n, stderr, 0);
    fputc ('\n', stderr);
  }
  else
    ast_print_tree (n, stderr, 0, false, -1);
  ast_destroy (n);
}
