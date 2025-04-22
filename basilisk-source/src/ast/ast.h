/**
# Interface for the Abstract Syntax Tree (AST) library

It is divided in two parts: a first part for functions which are meant
to be independent from the grammar, and a second part for (Basilisk C)
grammar-dependent functions. 

## Grammar-independent functions */

#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include "allocator.h"
#include "stack.h"

typedef struct _Ast Ast;

struct _Ast {
  int sym;
  Ast ** child, * parent;
};

typedef struct {
  Ast ast;
  char * before, * start, * after;
  const char * file;
  int line;
  void * value;
} AstTerminal;

typedef struct {
  Ast ast;
  char * before, * after;
  const char * file;
  Allocator * alloc;
  Stack * stack;
  bool type_already_specified;
} AstRoot;

AstRoot * ast_parse            (const char * code, AstRoot * parent);
Ast *     ast_parse_expression (const char * expr, AstRoot * parent);
Ast *     ast_parse_external_declaration (const char * decl, AstRoot * parent);
void      ast_destroy          (Ast * n);
AstRoot * ast_parse_file       (FILE * fp, AstRoot * parent);
char *    ast_str_print        (const Ast * n, char * s, int kind, int real);
void      ast_print            (const Ast * n, FILE * fp, int kind);
void      ast_print_tree       (Ast * n, FILE * fp, const char * indent,
				bool compress, int maxdepth);
void      ast_print_constructor (Ast * n, FILE * fp, const char * indent);
void      ast_print_file_line  (Ast * n, FILE * fp);
AstRoot * ast_get_root         (const Ast * n);
void      ast_identifier_print (Ast * identifier, FILE * fp);
void      ast_stack_print      (Stack * stack, FILE * fp);
const
char *    ast_crop_before (const char * s);
double    ast_evaluate_constant_expression (const Ast * n);

static inline Ast * ast_last_child (const Ast * n)
{
  if (!n)
    return NULL;
  Ast ** c = n->child;
  if (!c)
    return NULL;
  while (*(c + 1)) c++;
  return *c;
}

extern Ast * const ast_placeholder;

static inline Ast * ast_child (const Ast * n, int sym)
{
  if (!n)
    return NULL;
  Ast ** c = n->child;
  if (c == NULL)
    return NULL;
  while (*c && (*c == ast_placeholder || (*c)->sym != sym)) c++;
  return *c;
}

static inline Ast * ast_parent (const Ast * n, int sym)
{
  if (!n || n == ast_placeholder)
    return NULL;
  Ast * parent = n->parent;
  while (parent && parent->sym != sym)
    parent = parent->parent;
  return parent;
}

static inline AstTerminal * ast_left_terminal (const Ast * n)
{
  while (n && n != ast_placeholder && n->child) {
    Ast ** c = n->child;
    while (*c && *c == ast_placeholder) c++;
    n = *c;
  }
  return (AstTerminal *) (n != ast_placeholder ? n : NULL);
}

static inline AstTerminal * ast_right_terminal (const Ast * n)
{
  while (n && n != ast_placeholder && n->child)
    n = ast_last_child (n);
  return (AstTerminal *) (n != ast_placeholder ? n : NULL);
}
AstTerminal * ast_next_terminal (const Ast * n);

AstTerminal * ast_terminal_new (Ast * parent, int symbol, const char * start);
#define ast_terminal_new_char(p,s)			\
  (Ast *) ast_terminal_new (p, token_symbol((s)[0]), s)

Ast * ast_new_internal (Ast * parent, ...);
#define ast_new(parent,...) ast_new_internal (parent, __VA_ARGS__, -1)
Ast * ast_schema_internal (const Ast * n, ...);
#define ast_schema(n,...) ast_schema_internal (n, __VA_ARGS__, -1)
Ast * ast_find_internal (const Ast * n, const char * identifier, ...);
#define ast_find(n,...) ast_find_internal (n, NULL, __VA_ARGS__, -1)
#define ast_find_identifier(id,n,...) ast_find_internal (n, id, __VA_ARGS__, -1)
Ast * ast_copy_single (const Ast * n,
		       AstRoot ** dst_root, AstRoot ** src_root);
Ast * ast_copy_internal (const Ast * n, ...);
#define ast_copy(n,...) ast_copy_internal (n, __VA_ARGS__ + 0, -1)
Ast * ast_attach_internal (Ast * parent, ...);
#define ast_attach(p,...) ast_attach_internal (p, __VA_ARGS__, NULL)
Ast * ast_new_children_internal (Ast * parent, ...);
#define ast_new_children(p,...) \
  ast_new_children_internal (p, __VA_ARGS__, NULL)

void ast_set_child (Ast * parent, int index, Ast * child);
void ast_replace_child (Ast * parent, int index, Ast * child);

static inline int ast_child_index (const Ast * n)
{
  assert (n->parent);
  int index = 0;
  Ast ** c;
  for (c = n->parent->child; *c && *c != n; c++, index++);
  return *c == n ? index : - 1;
}

static inline Ast * ast_ancestor (const Ast * n, int i)
{
  while (n && i)
    n = n->parent, i--;
  return (Ast *) n;
}

char * str_append_realloc (char * dst, ...);
#define str_append(dst, ...)						\
  do { dst = str_append_realloc (dst, __VA_ARGS__, NULL); } while(0)
char * str_prepend_realloc (char * dst, ...);
#define str_prepend(dst, ...)						\
  do { dst = str_prepend_realloc (dst, __VA_ARGS__, NULL); } while(0)

#define ast_before(n,...) str_append(ast_left_terminal (n)->before, __VA_ARGS__)
#define ast_after(n,...)  str_append(ast_right_terminal (n)->after, __VA_ARGS__)
#define ast_terminal(n) ((n)->child ? NULL : (AstTerminal *)(n))
#define ast_root(n) ((n)->parent ? NULL : (AstRoot *)(n))
char *  ast_str_append (const Ast * n, char * s);

static inline void ast_hide (AstTerminal * n)
{
  for (char * s = n->start; *s != '\0'; s++)
    *s = ' ';
}

char * ast_line (AstTerminal * t);
#define ast_file_line(n, nolineno)					\
  "\"", ast_terminal((Ast *)n)->file, "\",",				\
    nolineno ? "0" : ast_line(ast_terminal((Ast *)n))
void ast_set_line (Ast * n, AstTerminal * l, bool overwrite);
Ast * ast_flatten (Ast * n, AstTerminal * t);
AstTerminal * ast_replace (Ast * n, const char * terminal, Ast * with);

#define NN(parent,sym,...) ast_new_children (ast_new (parent, sym), __VA_ARGS__)

static inline AstTerminal * NB (Ast * parent, int sym, const char * name)
{
  AstTerminal * t = ast_terminal_new (parent, sym, name);
  AstTerminal * r = ast_left_terminal (parent);
  t->before = r->before, r->before = NULL;
  return t;
}

static inline AstTerminal * NA (Ast * parent, int sym, const char * name)
{
  AstTerminal * t = ast_terminal_new (parent, sym, name);
  AstTerminal * r = ast_right_terminal (parent);
  t->line = r->line;
  return t;
}

#define NCB(parent,sym) NB(parent, token_symbol((sym)[0]), sym)
#define NCA(parent,sym) NA(parent, token_symbol((sym)[0]), sym)

/**
## Grammar-specific functions */

void  ast_push_declaration         (Stack * stack, Ast * declaration);
Ast * ast_push_function_definition (Stack * stack, Ast * declarator);
void  ast_pop_scope                (Stack * stack, Ast * scope);
Ast * ast_push_declarations        (Ast * n, Stack * stack);
void  ast_traverse                 (Ast * n, Stack * stack,
				    void func (Ast *, Stack *, void *),
				    void * data);

#define foreach_item(list, index, item)					\
  for (Ast * _list = list, * item = _list && _list != ast_placeholder ?	\
	 (_list->child[1] ? _list->child[index] : _list->child[0]) : NULL; \
       (_list = _list && _list != ast_placeholder &&			\
	_list->child[1] ? _list->child[0] : NULL), item;		\
       item = _list && _list != ast_placeholder ?			\
	 (_list->child[1] ? _list->child[index] :			\
	  _list->child[0]) : NULL					\
       )

#define foreach_item_r(list, symbol, arg)				\
  while (list->child[0]->sym == list->sym)				\
    list = list->child[0];						\
  for (Ast * arg = ast_child (list, symbol); arg;			\
       list = list->parent, arg = ast_child (list, symbol))

Ast * ast_identifier_declaration (Stack * stack, const char * identifier);
int   ast_identifier_parse_type (Stack * stack, const char * identifier, bool call,
				 const char * file, int line);
Ast * ast_identifier_declaration_from_to (Stack * stack, const char * identifier,
					  const Ast * start, const Ast * end);
Ast * ast_function_identifier (const Ast * function_definition);
Ast * ast_function_call_identifier (const Ast * n);

void ast_set_char (Ast * n, int c);
void ast_remove_internal (Ast * n, AstTerminal * before);
void ast_remove (Ast * n, AstTerminal * before);
void ast_erase (Ast * n);
void ast_check (Ast * n);
Ast * ast_is_typedef (const Ast * identifier);
Ast * ast_find_function (Ast * n, const char * name);
Ast * ast_list_append_list (Ast * list, Ast * list1);
Ast * ast_block_list_append (Ast * list, int item_sym, Ast * item);
Ast * ast_block_list_prepend (Ast * list, int item_sym, Ast * item);
Ast * ast_block_list_insert_after (Ast * insert, Ast * item);
Ast * ast_block_list_insert_before (Ast * insert, Ast * item);
Ast * ast_block_list_insert_before2 (Ast * insert, Ast * item);
Ast * ast_block_list_get_item (Ast * statement);
Ast * ast_list_append (Ast * list, int item_sym, Ast * item, const char * separator);
Ast * ast_list_prepend (Ast * list, int item_sym, Ast * item);
Ast * ast_list_remove (Ast * list, Ast * item);
Ast * ast_list_insert_after (Ast * insert, Ast * item);
void  ast_argument_list (Ast * expression);
Ast * ast_new_cast_expression (Ast * parent);
Ast * ast_new_unary_expression (Ast * parent);
Ast * ast_new_constant (Ast * parent, int symbol, const char * value);
Ast * ast_new_identifier (Ast * parent, const char * name);
Ast * ast_new_member_identifier (Ast * parent, const char * name);
Ast * ast_new_function_call (Ast * parent, const char * func);
Ast * ast_new_empty_scalar (Ast * parent);
Ast * ast_new_empty_vector (Ast * parent, int dimension);
Ast * ast_new_empty_tensor (Ast * parent, int dimension);
Ast * ast_is_unary_expression (const Ast * n);
Ast * ast_is_identifier_expression (const Ast * n);
Ast * ast_is_simple_expression (const Ast * n);
Ast * ast_get_struct_name (Ast * declaration_specifiers);
bool  ast_are_identical (const Ast * a, const Ast * b);
Ast * ast_declaration_from_type (const Ast * type);
Ast * ast_expression_type (Ast * expr, Stack * stack, bool higher_dimension);
AstTerminal * ast_type (const Ast * identifier);
char * ast_typedef_name (const Ast * identifier);
bool ast_is_field (const char * typename);

Ast * ast_check_grammar (Ast * n, bool recursive, bool stencils);

Ast * ast_get_function_definition (Stack * stack, Ast * identifier, Ast * declaration);
bool  ast_is_foreach_stencil (const Ast * n);
bool  ast_is_stencil_function (Ast * n);
Ast * ast_is_point_function (const Ast * declarator);
Ast * ast_stencil (Ast * n, bool parallel, bool overflow, bool nowarning);
void  ast_stencil_access (Ast * n, Stack * stack, int dimension);
const Ast * ast_attribute_access (const Ast * n, Stack * stack);
Ast * ast_attribute_array_access (Ast * n);
Ast * ast_constant_postfix_expression (const Ast * n, Stack * stack);
Ast * ast_is_function_pointer (const Ast * n, Stack * stack);
bool  ast_is_foreach_statement (const Ast * n);

/**
## Types */

extern AstTerminal ast_int, ast_long, ast_enum, ast_char, ast_void, ast_float, ast_double, ast_function, ast_bool;

typedef struct {
  int pointer;
  Ast ** dimension;
} AstDimensions;

Ast * ast_identifier_type      (Ast * n, AstDimensions * d, Stack * stack);
Ast * ast_base_type            (Ast * type, AstDimensions * d, Stack * stack);
Ast * ast_get_array_dimensions (Ast * direct_declarator, int symbol, AstDimensions * d, int nd, Stack * stack);

/**
## Kernels */

void ast_diagonalize (Ast * n, Stack * stack, void * field);
char * ast_external_references (Ast * n, char * references, Stack * functions);
char * ast_kernel              (Ast * n, char * argument, bool nolineno, Ast * macroscope);

/**
## Macros */

Ast * ast_is_macro_declaration (const Ast * function_declaration);
void ast_macro_replacement (Ast * statement, Ast * initial, Stack * stack,
			    bool nolineno, int postmacros, bool expand_definitions,
			    int * return_macro_index, Ast * scope);

/**
## Interface for the generic C interpreter */

int ast_run (AstRoot * root, Ast * n, int verbosity, int maxcalls, void * data);

/**
## Interface for dimensional analysis */

bool ast_check_dimensions (AstRoot * root, Ast * n, int verbosity, int maxcalls,
			   FILE * dimensions, int finite, int redundant, int lineno, int warn);

/**
## The entry function

Called by [qcc](/src/qcc.c) to trigger the translation. */

AstRoot * endfor (FILE * fin, FILE * fout,
		  const char * grid, int dimension,
		  bool nolineno, bool progress, bool catch,
		  bool parallel, bool cpu, bool gpu,
		  bool prepost,
		  FILE * swigfp, char * swigname);
