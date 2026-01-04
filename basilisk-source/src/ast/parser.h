#include <string.h>
#include <assert.h>
#include <stdlib.h>
  
#include "ast.h"

int  yylex (Ast ** lvalp, AstRoot * parse);
void yyerror (AstRoot * parse, Ast * root, char const *);
int  yylex_destroy();
int  token_symbol (int token);
int  sym_type (const char * name);
