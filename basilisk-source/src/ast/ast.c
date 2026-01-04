#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <stdarg.h>

#include "parser.h"
#include "basilisk.h"
#include "symbols.h"

Ast * const ast_placeholder = (Ast *) 128;

Ast * ast_new_children_internal (Ast * parent, ...)
{
  va_list ap;
  va_start (ap, parent);
  int n = 1;
  Ast * c = va_arg (ap, Ast *);
  while (c) {
    n++;
    c = va_arg (ap, Ast *);
  }
  va_end (ap);

  if (n > 1) {
    parent->child = allocate (ast_get_root (parent)->alloc, n*sizeof (Ast *));
    parent->child[n - 1] = NULL;
    va_start (ap, parent);
    n = 0;
    c = va_arg (ap, Ast *);
    while (c) {
      ast_set_child (parent, n++, c);
      c = va_arg (ap, Ast *);
    }
    va_end (ap);
  }
  
  return parent;
}

static void cancel_file_line (Ast * n)
{
  AstTerminal * t = ast_terminal(n);
  if (t)
    t->file = NULL, t->line = 0;
  else
    for (Ast ** c = n->child; *c; c++)
      cancel_file_line (*c);
}

static Ast * ast_attach_single (Ast * parent, Ast * n)
{
  AstRoot * root = ast_get_root (parent);
  while (parent->child)
    parent = parent->child[0];
  parent->child = allocate (root->alloc, 2*sizeof (Ast *));
  parent->child[0] = NULL;
  parent->child[1] = NULL;
  ast_set_child (parent, 0, n);
  return parent;
}

Ast * ast_attach_internal (Ast * parent, ...)
{
  va_list ap;
  va_start (ap, parent);
  Ast * p = parent, * n = va_arg (ap, Ast *);
  while (n) {
    ast_attach_single (p, n);
    p = n;
    n = va_arg (ap, Ast *);
  }
  va_end (ap);
  return parent;
}

static void ast_destroy_internal (Ast * n)
{
  if (n->child) {
    for (Ast ** c = n->child; *c; c++)
      if (*c != ast_placeholder)
	ast_destroy_internal (*c);
  }
  else {
    AstTerminal * t = ast_terminal (n);
    free (t->before);
    free (t->start);
    free (t->after);
    t->before = t->start = t->after = NULL;
  }
  AstRoot * r = ast_root (n);
  if (r) {
    free (r->before);
    free (r->after);
    if (r->stack)
      stack_destroy (r->stack);
    if (r->alloc)
      free_allocator (r->alloc);
  }
}

void ast_destroy (Ast * n)
{
  if (!n || n == ast_placeholder)
    return;
  if (n->parent && n->parent->child) {
    Ast ** c;
    for (c = n->parent->child; *c && *c != n; c++);
    if (*c == n)
      *c = ast_placeholder;
  }
  ast_destroy_internal (n);
}

char * ast_line (AstTerminal * t)
{
  static char s[20];
  snprintf (s, 19, "%d", t->line);
  return s;
}

AstRoot * ast_get_root (const Ast * n)
{
  const Ast * root = n;
  while (root->parent) {
    assert (root != root->parent);
    root = root->parent;
  }
  return ast_root (root);
}

static int count_lines (const char * s)
{
  int line = 0;
  while (*s != '\0') {
    if (*s == '\n')
      line++;
    s++;
  }
  return line;
}

void ast_set_child (Ast * parent, int index, Ast * child)
{
  if (!child)
    return;
  if (child == ast_placeholder) {
    parent->child[index] = child;
    return;
  }
  if (child->parent) {
    Ast * oldparent = child->parent;
    if (oldparent->child) {
      Ast ** c;
      for (c = oldparent->child; *c && *c != child; c++);
      if (*c == child)
	*c = ast_placeholder;
    }
  }
  parent->child[index] = child;
  child->parent = parent;
}

static AstTerminal * ast_left_line (Ast * n)
{
  AstTerminal * t = ast_terminal (n);
  if (t)
    return t->file ? t : NULL;
  else
    for (Ast ** c = n->child; *c; c++) {
      t = ast_left_line (*c);
      if (t)
	return t;
    }
  return NULL;
}

void ast_replace_child (Ast * parent, int index, Ast * child)
{
  Ast * oldchild = parent->child[index];
  ast_set_child (parent, index, child);
  AstTerminal * left = ast_left_terminal (child);
  if (oldchild != ast_placeholder) {
    AstTerminal * oldleft = ast_left_terminal (oldchild);
    if (oldleft && left) {
      if (!left->file) {
	if (left->before)
	  free (left->before), left->before = NULL;
	assert (!left->before);
	left->file = oldleft->file;
	left->line = oldleft->line;
	ast_set_line (child, left, false);
      }
      if (left->before)
	str_append (left->before, oldleft->before);
      else
	left->before = oldleft->before, oldleft->before = NULL;
    }
    ast_destroy (oldchild);
  }
  if (left && !left->file) {
    AstTerminal * line = ast_left_line (child);
    if (line) {
      if (left->before)
	str_append (left->before, line->before);
      else
	left->before = line->before, line->before = NULL;      
      left->file = line->file;
      left->line = line->line;
      ast_set_line (child, left, false);
    }      
  }
}

AstTerminal * ast_replace (Ast * n, const char * terminal, Ast * with)
{
  AstTerminal * t = ast_terminal (n);
  if (t && !strcmp (t->start, terminal)) {
    while (n && n->sym != with->sym)
      n = n->parent;
    if (n) {
      ast_set_child (n->parent, ast_child_index (n), with);
      ast_destroy (n);
      return ast_right_terminal (with);
    }
    return NULL;
  }
  if (n->child)
    for (Ast ** c = n->child; *c; c++) {
      AstTerminal * right = ast_replace (*c, terminal, with);
      if (right) {
	int after = count_lines (right->start);
	right->line += after;
	for (c++; *c; c++)
	  ast_set_line (*c, right, false);
	right->line -= after;
	return right;
      }
    }
  return NULL;
}

typedef struct {
  const char * name;
  int line, macro;
} File;

static void update_file_line (const char * preproc, File * file)
{
  const char * s = preproc;
  while (*s != '\0') {
    if (*s == '#') {
      const char * u = s > preproc ? s - 1 : preproc;
      while (u != preproc && !strchr ("\n\r", *u) && strchr (" \t", *u))
	u--;
      if (u == preproc || strchr ("\n\r", *u)) {
	s++;
	while (*s && strchr (" \t", *s)) s++;
	if (strchr ("0123456789", *s)) {
	  int line = atoi (s) - 1;
	  while (*s && !strchr (" \t", *s)) s++;
	  while (*s && strchr (" \t", *s)) s++;
	  if (*s == '"') {
	    file->line = line;
	    file->name = s + 1;
	  }
	}
      }
    }
    else if (*s == '\n')
      file->line++;
    s++;
  }  
}

typedef struct {
  char * s;
  int len, size;
} String;

void string_append (void * a, ...)
{
  va_list ap;
  va_start (ap, a);
  char * s = va_arg(ap, char *);
  while (s) {
    String * b = a;
    b->len += strlen (s);
    if (b->len >= b->size) {
      b->size = b->len + 1024;
      char * n = realloc (b->s, b->size);
      if (!b->s) n[0] = '\0';
      b->s = n;
    }
    strcat (b->s, s);
    s = va_arg(ap, char *);
  }
  va_end (ap);
}

static const char * only_spaces (const char * s, const File * file, const AstTerminal * t)
{
  const char * spaces = NULL;
  int line = file->line;
  while (*s != '\0') {
    if (!strchr (" \t\n\r", *s))
      return NULL;
    if (strchr ("\n\r", *s)) {
      spaces = s + 1;
      line++;
    }
    s++;
  }
  return line > t->line ? spaces : NULL;
}

static void str_print_internal (const Ast * n, int sym, int real, File * file,
				void (* output) (void *, ...), void * data)
{
  if (n == ast_placeholder) {
    output (data, "$", NULL);
    return;
  }

  //  ast_print_file_line (n, stderr);
  AstTerminal * t = ast_terminal (n);
  if (t) {
    if (t->before) {
      const char * spaces = only_spaces (t->before, file, t);
      if (spaces)
	output (data, spaces, NULL);
      else {
	output (data, t->before, NULL);
	update_file_line (t->before, file);
      }
    }
    if (t->file) {
      int len = strlen (t->file);
      if (!file->name || strncmp (t->file, file->name, len) ||
	  (file->name[len] != '"' && file->name[len] != '\0')) {
	char line[20]; snprintf (line, 19, "%d", t->line);
	if (sym == 3)
	  output (data, "\n", t->file, ":", line, " ", NULL);
	else
	  output (data, "\n#line ", line, " \"", t->file, "\"\n", NULL);
	file->line = t->line;
	file->name = t->file;      
      }
      else if (t->line != file->line) {
	if (file->line > t->line || t->line - file->line > 10) {
	  char line[20]; snprintf (line, 19, "%d", t->line);
	  if (sym == 3)
	    output (data, "\n", t->file, ":", line, " ", NULL);
	  else
	    output (data, "\n#line ", line, "\n", NULL);
	  file->line = t->line;
	}
	else
	  while (file->line < t->line) {
	    output (data, "\n", NULL);
	    file->line++;
	  }
      }
    }
    if (sym == 1)
      output (data, "|", symbol_name (n->sym), "|", NULL);
    else if (sym == 2) {
      if (t->value)
	output (data, "^", NULL);
    }
    if (real && (n->sym == sym_DOUBLE || n->sym == sym_FLOAT))
      output (data, "real", NULL);
    else if (n->sym == sym_MACRODEF) {
      if (ast_schema (ast_ancestor (n, 2), sym_declaration_specifiers,
		      1, sym_declaration_specifiers))
	output (data, "", NULL);
      else
	output (data, "void", NULL);
    }
    else if (n->sym == sym_AUTO) {
      if (ast_schema (ast_ancestor (n, 2), sym_declaration_specifiers,
		      1, sym_declaration_specifiers,
		      0, sym_storage_class_specifier,
		      0, sym_MACRODEF) ||
	  ast_schema (ast_ancestor (n, 3), sym_parameter_declaration))
	output (data, "    ", NULL);
      else
	output (data, t->start, NULL);
    }
    else if (n->sym == sym_ELLIPSIS_MACRO)
      output (data, "{}", NULL);
    else if (n->sym == sym_IDENTIFIER) {
      if (ast_find (ast_schema (ast_ancestor (n, 5), sym_function_declaration,
				0, sym_declaration_specifiers), sym_MACRODEF)) {
	char s[20]; snprintf (s, 19, "%d", file->macro++);
	output (data, "macro", s, "_", NULL);
      }
      output (data, t->start, NULL);
    }
    else
      output (data, t->start, NULL);
    file->line += count_lines (t->start);
    if (sym == 1)
      output (data, "/", NULL);
    if (t->after) {
      output (data, t->after, NULL);
      file->line += count_lines (t->after);
    }
  }
  else { // !terminal

    /**
    Do not output macro definitions. */

    if (ast_is_macro_declaration (ast_schema (n, sym_function_definition,
					      0, sym_function_declaration))) {
      AstTerminal * t = ast_left_terminal (n);
      if (t->before) {
	const char * spaces = only_spaces (t->before, file, t);
	if (spaces)
	  output (data, spaces, NULL);
	else {
	  output (data, t->before, NULL);
	  update_file_line (t->before, file);
	}
      }
      return;
    }
    
    /**
    Ignore 'break =' macro parameters. */
    
    if (ast_schema (ast_child (n, sym_parameter_declaration), sym_parameter_declaration,
		    0, sym_BREAK)) {
      Ast * list = ast_child (n, sym_parameter_list);
      if (list)
	str_print_internal (list, sym, real, file, output, data);
      return;
    }

    AstRoot * r = ast_root (n);
    if (r && r->before) {
      output (data, r->before, NULL);
      update_file_line (r->before, file);
    }

    /** 
    Ignore dimensioned constant expressions e.g. `1.23
    [1,-1]`, `(2.*3 + 4.)[0,1]` etc. */
    
    if (n->sym == sym_array_access && n->child[2] &&
	ast_evaluate_constant_expression (n->child[0]) < DBL_MAX)
      str_print_internal (n->child[0], sym, real, file, output, data);

    /**
    Ignore the values of optional parameters. */

    else if (ast_schema (n, sym_parameter_declaration,
			 3, sym_initializer))
      for (int i = 0; i < 2; i++)
	str_print_internal (n->child[i], sym, real, file, output, data);
    
    else
      for (Ast ** c = n->child; *c; c++)
	str_print_internal (*c, sym, real, file, output, data);
    
    if (r && r->after) {
      output (data, r->after, NULL);
      file->line += count_lines (r->after);
    }
  }
}

char * ast_str_print (const Ast * n, char * s, int sym, int real)
{
  String a = {0};
  str_print_internal (n, sym, real, &(File){0}, string_append, &a);
  if (s) {
    str_append (s, a.s);
    free (a.s);
    return s;
  }
  return a.s;
}

static void file_append (void * a, ...)
{
  va_list ap;
  va_start (ap, a);
  char * s = va_arg(ap, char *);
  while (s) {
    fputs (s, a);
    s = va_arg(ap, char *);
  }
  va_end (ap);
}

void ast_print (const Ast * n, FILE * fp, int sym)
{
  str_print_internal (n, sym, 0, &(File){0}, file_append, fp);
}

char * ast_str_append (const Ast * n, char * s)
{
  AstTerminal * t = ast_terminal (n);
  if (t) {
    if (t->before)
      str_append (s, t->before);
    str_append (s, t->start);
    if (t->after)
      str_append (s, t->after);
  }
  else {
    AstRoot * r = ast_root (n);
    if (r && r->before)
      str_append (s, r->before);
    for (Ast ** c = n->child; *c; c++)
      s = ast_str_append (*c, s);
    if (r && r->after)
      str_append (s, r->after);
  }
  return s;
}

void ast_print_file_line (Ast * n, FILE * fp)
{
  assert (n);
  AstTerminal * t = ast_left_terminal (n);
  assert (t);
  fprintf (fp, "%s:%d: %s\n", t->file, t->line,
	   ast_terminal(n) ? t->start : symbol_name (n->sym));
}

static void print_child_tree (Ast * n, FILE * fp,
			      const char * indent, bool isLast,
			      bool compress, int maxdepth)
{
  char * ind;
  if (indent) {
    fputs (indent, fp);
    ind = strdup (indent);
  }
  else
    ind = strdup ("");
    
  if (isLast) {
    fputs ("└─", fp);
    str_append (ind, "  ");
  }
  else {
    fputs ("├─", fp);
    str_append (ind, "│ ");
  }
  ast_print_tree (n, fp, ind, compress, maxdepth);
  free (ind);
}

void ast_print_tree (Ast * n, FILE * fp, const char * indent,
		     bool compress, int maxdepth)
{
  if (n == ast_placeholder) {
    fputs ("_placeholder_\n", fp);
    return;
  }
  if (!maxdepth) {
    fputs ("...\n", fp);
    return;
  }
  if (compress)
    while (n->child && !n->child[1] && n->child[0] != ast_placeholder)
      n = n->child[0];
  fprintf (fp, "%s", symbol_name (n->sym));
  AstTerminal * t = ast_terminal (n);
  if (t)
    fprintf (fp, " %s %s:%d\n", t->start, t->file, t->line);
  else {
    fputc ('\n', fp);
    for (Ast **c = n->child; *c; c++)
      print_child_tree (*c, fp, indent, *(c + 1) == NULL,
			compress, maxdepth - 1);
  }
}

void ast_print_constructor (Ast * n, FILE * fp, const char * indent)
{
  if (indent)
    fputs (indent, fp);
  AstTerminal * t = ast_terminal (n);
  if (t) {
    if (t->start[1] == '\0')
      fprintf (fp, "NCA(n, \"%s\")", t->start);
    else
      fprintf (fp, "NA(n, sym_%s, \"%s\")", symbol_name (n->sym), t->start);
  }
  else {
    fprintf (fp, "NN(n, sym_%s,\n", symbol_name (n->sym));
    char * ind = indent ? strdup (indent) : strdup ("");
    str_append (ind, "   ");
    for (Ast **c = n->child; *c; c++) {
      ast_print_constructor (*c, fp, ind);
      if (*(c + 1))
	fputs (",\n", fp);
    }
    fputs (")", fp);
    free (ind);
  }
}

AstTerminal * ast_terminal_new (Ast * parent, int symbol, const char * start)
{
  AstTerminal * t = allocate (ast_get_root (parent)->alloc,
			      sizeof (AstTerminal));
  memset (t, 0, sizeof (AstTerminal));
  ((Ast *)t)->sym = symbol;
  ((Ast *)t)->parent = parent;
  t->start = strdup (start);
  AstTerminal * r = ast_left_terminal (parent);
  t->file = r->file;
  t->line = r->line;
  return t;
}

static Ast * vast_new_internal (Ast * parent, va_list ap)
{
  int sym = va_arg (ap, int);
  if (sym < 0)
    return NULL;
  AstRoot * root = ast_get_root (parent);
  assert (root->alloc);
  Ast * n = allocate (root->alloc, sizeof (Ast));
  n->sym = sym;
  n->parent = parent;
  n->child = NULL;
  Ast * m = n;
  sym = va_arg (ap, int);
  while (sym >= 0) {
    Ast * c = allocate (root->alloc, sizeof (Ast));
    c->sym = sym;
    c->parent = m;
    m->child = allocate (root->alloc, 2*sizeof (Ast *));
    m->child[0] = c;
    m->child[1] = NULL;
    m = c;
    sym = va_arg (ap, int);
  }
  return n;
}

Ast * ast_new_internal (Ast * parent, ...)
{
  va_list ap;
  va_start (ap, parent);
  Ast * n = vast_new_internal (parent, ap);
  va_end (ap);
  return n;
}

static Ast * vast_schema_internal (const Ast * n, va_list ap)
{
  int sym = va_arg(ap, int);
  if (n == ast_placeholder || n->sym != sym)
    return NULL;
  int c = va_arg(ap, int);
  while (c >= 0) {
    if (!n->child)
      return NULL;
    int i;
    for (i = 0; i <= c && n->child[i]; i++);
    if (i <= c)
      return NULL;
    n = n->child[c];
    int sym = va_arg(ap, int);
    if (n == ast_placeholder || n->sym != sym)
      return NULL;
    c = va_arg(ap, int);
  }
  return (Ast *) n;
}

Ast * ast_schema_internal (const Ast * n, ...)
{
  if (!n)
    return NULL;
  va_list ap;
  va_start (ap, n);
  n = vast_schema_internal (n, ap);
  va_end (ap);
  return (Ast *) n;
}

static Ast * vast_find_internal (const Ast * n, const char * identifier, va_list ap)
{
  if (n == ast_placeholder)
    return NULL;
  va_list bp;
  va_copy (bp, ap);
  Ast * found = vast_schema_internal (n, bp);
  va_end (bp);
  if (found && (!identifier ||
		(ast_terminal (found) && !strcmp (ast_terminal (found)->start, identifier))))
    return found;
  if (!ast_terminal(n))
    for (Ast ** c = n->child; *c; c++)
      if ((found = vast_find_internal (*c, identifier, ap)))
	return found;
  return NULL;
}

Ast * ast_find_internal (const Ast * n, const char * identifier, ...)
{
  if (!n)
    return NULL;
  va_list ap;
  va_start (ap, identifier);
  n = vast_find_internal (n, identifier, ap);
  va_end (ap);
  return (Ast *) n;
}

static Ast * copy_ast (Ast * dst, const Ast * src)
{
  dst->sym = src->sym;
  dst->child = NULL;
  dst->parent = dst;
  return dst;
}

static Ast * terminal_copy (const AstTerminal * src,
			    const AstRoot * dst_root, const AstRoot * src_root)
{
  AstTerminal * dst =
    (AstTerminal *) copy_ast (allocate (dst_root->alloc, sizeof (AstTerminal)),
			      (Ast *)src);
  dst->before = src->before ? strdup (src->before) : NULL;
  dst->start = strdup (src->start);
  dst->after = src->after ? strdup (src->after) : NULL;
  dst->line = src->line;
  dst->file = src->file;
  return (Ast *)dst;
}

static Ast * root_copy (const AstRoot * src)
{
  Allocator * alloc = new_allocator();
  AstRoot * dst = (AstRoot *) copy_ast (allocate (alloc, sizeof (AstRoot)),
					(Ast *)src);
  dst->alloc = alloc;
  dst->before = src->before ? strdup (src->before) : NULL;
  dst->after = src->after ? strdup (src->after) : NULL;
  dst->stack = NULL;
  ((Ast *)dst)->parent = NULL;
  return (Ast *)dst;
}

Ast * ast_copy_single (const Ast * n,
		       AstRoot ** dst_root, AstRoot ** src_root)
{
  Ast * c = NULL;
  AstRoot * r = ast_root (n);
  if (r) {
    c = root_copy (r);
    *src_root = r;
    *dst_root = ast_root (c);
  }
  AstTerminal * t = ast_terminal (n);
  if (t)
    c = terminal_copy (t, *dst_root, *src_root);
  else if (!c)
    c = copy_ast (allocate ((*dst_root)->alloc, sizeof (Ast)), n);
  if (!t) {
    int len = 0;
    for (Ast ** i = n->child; *i; i++, len++);
    c->child = allocate ((*dst_root)->alloc, (len + 1)*sizeof (Ast *));
    c->child[len] = NULL;
  }
  return c;
}

static Ast * vast_copy_internal (const Ast * n, va_list ap, bool * found,
				 AstRoot * dst_root, AstRoot * src_root)
{
  Ast * c = ast_copy_single (n, &dst_root, &src_root);
  va_list cp;
  va_copy (cp, ap);
  if (vast_schema_internal (n, cp))
    *found = true;
  va_end (cp);

  if (!ast_terminal (n)) {
    int len = 0;
    for (Ast ** i = n->child, ** j = c->child;
	 *i && !(*found); i++, j++, len++) {
      if (*i == ast_placeholder)
	*j = ast_placeholder;
      else {
	*j = vast_copy_internal (*i, ap, found, dst_root, src_root);
	(*j)->parent = c;
      }
    }
    c->child[len] = NULL;
  }
  return c;
}

Ast * ast_copy_internal (const Ast * n, ...)
{
  if (!n) return NULL;
  va_list ap;
  va_start (ap, n);
  bool found = false;
  AstRoot * src_root = ast_get_root (n);
  Ast * c = vast_copy_internal (n, ap, &found, src_root, src_root);
  c->parent = (Ast *) n;
  va_end (ap);
  return c;
}

Ast * ast_parse_expression (const char * expr, AstRoot * parent)
{
  char * s = NULL;
  str_append (s, "void main() {", expr, "}");
  Ast * n = (Ast *) ast_parse (s, parent);
  free (s);
  if (n) {
    ast_pop_scope (parent->stack, n);
    Ast * c = ast_find (n, sym_statement);
    if (c)
      c = c->child[0];
    else
      c = ast_find (n, sym_declaration); 
    cancel_file_line (c);
    Ast ** i;
    for (i = c->parent->child; *i && *i != c; i++);
    *i = ast_placeholder;
    ast_destroy (n);
    n = c;
  }
  return n;
}

Ast * ast_parse_external_declaration (const char * decl, AstRoot * parent)
{
  Ast * def = (Ast *) ast_parse (decl, parent);
  Ast * n = ast_find (def, sym_external_declaration);
  if (!n) {
    ast_destroy (def);
    return NULL;
  }
  cancel_file_line (n);
  Ast ** i;
  for (i = n->parent->child; *i && *i != n; i++);
  assert (*i == n);
  *i = ast_placeholder;
  ast_destroy (def);
  return n;
}

AstRoot * ast_parse_file (FILE * fp, AstRoot * parent)
{
  char * buffer = NULL;
  size_t len = 0, maxlen = 0;
  int c;
  while ((c = fgetc (fp)) != EOF) {
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
  AstRoot * root = ast_parse (buffer, parent);
  free (buffer);
  return root;
}

/**
This is to ignore "face ", "vertex " and "symmetric " typedef prefixes. */

static const char * ignore_prefixes (const char * identifier)
{
  if (!strncmp (identifier, "face ", 5))
    return identifier + 5;
  else if (!strncmp (identifier, "vertex ", 7))
    return identifier + 7;
  else if (!strncmp (identifier, "symmetric ", 10))
    return identifier + 10;
  return identifier;
}

Ast * ast_identifier_declaration_from_to (Stack * stack, const char * identifier,
					  const Ast * start, const Ast * end)
{
  if (!identifier) return NULL;
  identifier = ignore_prefixes (identifier);
  
  Ast ** d;
  int i = 0;
  if (start) {
    for (; (d = stack_index (stack, i)) && *d != start; i++);
    if (d) {
      assert (*d == start);
      i++;
    }
    else {
      assert (false); // start not found!!
      return NULL;
    }
  }
  for (; (d = stack_index (stack, i)); i++)
    if (end && *d == end)
      break;
    else if (*d && (*d)->sym == sym_IDENTIFIER) {

      /**
      WARNING: this assumes that the "after" string is never modified
      for the declaration identifiers stored in the stack. */

      if (!ast_terminal(*d)->after) {
	if (ast_terminal(*d)->start &&
	    !strcmp (ast_terminal(*d)->start, identifier))
	  return *d;
      }

      /**
      If 'after' is defined, we assume that the function is called
      from the [lexer](tokens.lex) and that 'after' is the last
      character of the declaration identifier (as set by the
      lexer). */
      
      else {
	char * s = ast_terminal(*d)->start, * end = ast_terminal(*d)->after;
	const char * i = identifier;
	for (; *i != '\0' && s <= end && *s == *i; s++, i++);
	if (*i == '\0' && s == end + 1)
	  return *d;
      }
    }
  return NULL;
}

Ast * ast_identifier_declaration (Stack * stack, const char * identifier)
{
  Ast ** n = fast_stack_find (stack, identifier);
  if (!n)
    n = fast_stack_find (stack, ignore_prefixes (identifier));
  return n ? *n : NULL;
}

char * str_append_realloc (char * src, ...)
{
  va_list ap;
  va_start (ap, src);
  char * i = va_arg(ap, char *);
  int len = src ? strlen(src) : 0;
  while (i) {
    len += strlen (i);
    i = va_arg(ap, char *);
  }
  va_end (ap);

  if (len < 1)
    return NULL;
  
  char * dst = realloc (src, len + 1);
  if (!src)
    dst[0] = '\0';
  
  va_start (ap, src);
  i = va_arg(ap, char *);
  while (i) {
    strcat (dst, i);
    i = va_arg(ap, char *);
  }
  va_end (ap);

  return dst;
}

char * str_prepend_realloc (char * src, ...)
{
  va_list ap;
  va_start (ap, src);
  char * i = va_arg(ap, char *);
  int len = src ? strlen(src) : 0;
  while (i) {
    len += strlen (i);
    i = va_arg(ap, char *);
  }
  va_end (ap);

  if (len < 1)
    return NULL;
  
  char * dst = malloc (len + 1);
  dst[0] = '\0';

  va_start (ap, src);
  i = va_arg(ap, char *);
  while (i) {
    strcat (dst, i);
    i = va_arg(ap, char *);
  }
  va_end (ap);

  if (src)
    strcat (dst, src);

  free (src);
  return dst;
}

void ast_identifier_print (Ast * identifier, FILE * fp)
{
  AstTerminal * t = ast_terminal (identifier);
  fprintf (fp, "%s:%d: %p ", t->file, t->line, t);
  if (!t->after)
    fprintf (fp, "* %s", t->start);
  else {
    char * s = t->start, * end = t->after;
    for (; s <= end; s++)
      fputc (*s, fp);
  }
  fputc ('\n', fp);
}

void ast_stack_print (Stack * stack, FILE * fp)
{
  Ast ** n;
  for (int i = 0; (n = stack_index (stack, i)); i++)
    if (*n) {
      if (ast_terminal (*n))
	ast_identifier_print (*n, fp);
      else if (*n == ast_placeholder)
	fprintf (fp, "_placeholder_\n");
      else
	fprintf (fp, "%s %p\n", symbol_name ((*n)->sym), *n);
    }
}

void ast_set_char (Ast * n, int c)
{
  n->sym = token_symbol(c), ast_terminal (n)->start[0] = c;
}

void ast_remove_internal (Ast * n, AstTerminal * before)
{
  if (n->child) {
    for (Ast ** c = n->child; *c; c++)
      if (*c != ast_placeholder)
	ast_remove_internal (*c, before);
  }
  else {
    AstTerminal * t = ast_terminal (n);
    if (t->before) {
      str_prepend (before->before, t->before);
      free (t->before), t->before = NULL;
    }
    free (t->start), t->start = NULL;
    free (t->after), t->after = NULL;
  }
}

void ast_remove (Ast * n, AstTerminal * before)
{
  if (n->parent) {
    Ast ** c;
    for (c = n->parent->child; *c && *c != n; c++);
    if (*c == n)
      for (; *c; c++)
	*c = *(c + 1);
  }
  ast_remove_internal (n, before);
}

AstTerminal * ast_next_terminal (const Ast * n)
{
  if (!n || n == ast_placeholder)
    return NULL;
  Ast * parent = n->parent;
  while (parent) {
    int index = ast_child_index (n) + 1;
    if (index <= 0)
      return NULL;
    while ((n = parent->child[index])) {
      if (n != ast_placeholder)
	return ast_left_terminal (n);
      index++;
    }
    n = parent;
    parent = n->parent;
  }
  return NULL;
}

void ast_erase (Ast * n)
{
  AstTerminal * t = ast_next_terminal (n);
  if (t)
    ast_remove_internal (n, t);
  ast_destroy (n);
}

static void ast_check_children (Ast * n)
{
  if (n->child)
    for (Ast ** c = n->child; *c; c++) {
      assert ((*c)->parent == n);
      ast_check_children (*c);
    }
}

void ast_check (Ast * n)
{
  ast_check_children (n);
  while (n->parent) {
    assert (n->parent->child);
    Ast ** c;
    for (c = n->parent->child; *c && *c != n; c++);
    assert (*c == n);
    assert (n->parent != n);
    n = n->parent;
  }
}

void ast_set_line (Ast * n, AstTerminal * l, bool overwrite)
{
  if (n == ast_placeholder)
    return;
  AstTerminal * t = ast_terminal (n);
  if (t && (overwrite || !t->file)) {
    t->file = l->file;
    t->line = l->line;
  }
  if (n->child)
    for (Ast ** c = n->child; *c; c++)
      ast_set_line (*c, l, overwrite);
}

Ast * ast_flatten (Ast * n, AstTerminal * t)
{
  AstTerminal * r = ast_terminal (n);
  if (r) {
    if (r->before) {
      free (r->before);
      r->before = strdup (" ");
    }
    if (r->after) {
      free (r->after);
      r->after = NULL;
    }
    r->file = t->file;
    r->line = t->line;
  }
  if (n->child)
    for (Ast ** c = n->child; *c; c++)
      ast_flatten (*c, t);
  return n;
}

bool ast_are_identical (const Ast * a, const Ast * b)
{
  if (a == NULL && b == NULL)
    return true;
  if (a == NULL || b == NULL)
    return false;
  if (a->sym != b->sym)
    return false;
  AstTerminal * ta = ast_terminal (a), * tb = ast_terminal (b);
  if (ta) {
    if (!tb || strcmp (ta->start, tb->start))
      return false;
  }
  else if (tb)
    return false;
  else {
    Ast ** c, ** d;
    for (c = a->child, d = b->child; *c && *d; c++, d++)
      if (!ast_are_identical (*c, *d))
	return false;
    if (*c || *d)
      return false;
  }
  return true;
}
