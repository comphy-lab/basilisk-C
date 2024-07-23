# 2 "ast/interpreter/declarations.h"

/**
# Declarations for the interpreter
*/

enum AstBoolean { false = 0, true = 1 };
static const void * NULL = 0;
static const int _NVARMAX = 65536, INT_MAX = 2147483647;

FILE * stderr, * stdout, * systderr, * systdout;
