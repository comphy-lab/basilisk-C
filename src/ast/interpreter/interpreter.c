/**
# A generic C Interpreter

This file defines the function

~~~literatec
void ast_run (AstRoot * root, Ast * n, int verbosity, int maxcalls, 
              void * user_data));
~~~

which
[interprets](https://en.wikipedia.org/wiki/Interpreter_(computing))
the code defined by `n`, an Abstract Syntax Tree (AST) obtained by
parsing a Basilisk C code using the [AST Library](/src/ast/README).

The `root` parameter is the root of the AST associated with `n` (they
can be identical or different).

The `verbosity` parameter defines the level of the messages displayed
by the interpreter ("errors", "warnings", etc. see below).

The `maxcalls` parameter defines the maximum number of (recursive)
calls made by the interpreter. This is a rough measure of the maximum
number of instructions which will be interpreted and can be used to
limit the time taken by the interpreter. A negative value means no
limit.

The `user_data` parameter is a pointer passed to the 'hook functions'.

Two special functions can be used within the interpreted code to
display local information during interpretation.

* `interpreter_verbosity (int verbosity)`: sets the level of verbosity
  within the (local) scope of this call.
* `display_value (expression)`: displays the internal representation
  of `expression`. 
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <limits.h>
#include "ast.h"
#include "symbols.h"

#include "khash.h"
KHASH_MAP_INIT_INT64(INT64, void *)

#if 0 // for debugging memory accesses with valgrind
# define allocate(alloc, size) calloc (1, size)
#endif

static unsigned
  maximum_iterations = 32; // INT_MAX;
  
static const unsigned
  maximum_recursion = 2,
  display_pointers = false,
  
  error_verbosity = 1,
  warning_verbosity = 2,
  function_verbosity = 3,
  expression_verbosity = 4;

enum {
  unset               = 1 << 0,
  constant_expression = 1 << 1
};

typedef struct {
  void * p;
  char * start;
  int size, pscope;
} Pointer;

typedef struct _Value Value;

struct _Value {
  AstTerminal n;
  Ast * type;
  Pointer data; // fixme: should just be a void *
  int size, pointer, vscope;
  int * dimension;
 
  /**
  This could be used in principle to "unset" the values of array
  members accessed with an undefined index, but this breaks some test
  cases in a non-trivial way. The relevant test case is
  [test26.c](). */
  
#if UNSET_ARRAY
  Value * unset_array;
  Ast * unset_member;
#endif
};

static
void display_value (const Value * v);

Value * ast_run_node (Ast * n, Stack * stack);
extern Value * (* run) (Ast *, Stack *);
extern void    (* after_run) (Ast *, Stack *);
extern Value * (* ast_choose_hook) (Ast *, Stack *, Value *, Value *);
extern void    (* ast_value_print_hook) (const Value * v, FILE * fp);
extern Value * (* ast_binary_operation_hook) (Ast *, Stack *, Value *, Value *, Value *);
extern Value * (* ast_internal_functions_hook) (Ast * call, Ast * identifier, Value ** parameters,
						Stack * stack, Value * value);
Value * (* ast_assign_hook) (Ast * n, Value * dst, Value * src, Stack * stack) = NULL;

typedef struct {
  Allocator * alloc, * static_alloc;  
  int verbosity, maxcalls, warnings, scope, call;
  khash_t(INT64) * messages, * nonlocals;
  int conditional;
  Ast * root, * stencil, * end_stencil;
  void * data;
} StackData;

static
Allocator * stack_alloc (const Stack * s)
{
  return ((StackData *)stack_get_data (s))->alloc;
}

static
Allocator * stack_static_alloc (const Stack * s)
{
  return ((StackData *)stack_get_data (s))->static_alloc;
}

static
int stack_verbosity (const Stack * s)
{
  return ((StackData *)stack_get_data (s))->verbosity;
}

static
void * interpreter_get_data (const Stack * s)
{
  return ((StackData *)stack_get_data (s))->data;
}

typedef struct {
  int pointer, size;
  int * dimension;
} Dimensions;

static
Pointer stack_allocate (Stack * stack, int size)
{
  Pointer p;
  p.p = allocate (stack_alloc (stack), size);
  p.start = p.p;
  p.size = size;
  p.pscope = ((StackData *)stack_get_data (stack))->scope;
  return p;
}

static
bool is_new_message (const Ast * n, const Stack * stack)
{
  khash_t(INT64) * messages = ((StackData *)stack_get_data (stack))->messages;
  khiter_t k = kh_get (INT64, messages, (long) n);
  return (k == kh_end (messages));    
}

static
void * message (const Ast * scope, const Ast * n, const char * msg, int level, Stack * stack)
{
  if (stack_verbosity (stack) >= level && is_new_message (n, stack)) {
    AstTerminal * t = ast_left_terminal (n);
    char * s = ast_str_append (n, NULL);
    fprintf (stderr, "%s:%d: %s (interpreter): ", t->file, t->line,
	     level == error_verbosity ? "error" :
	     level == warning_verbosity ? "warning" :
	     "info");
    fprintf (stderr, msg, ast_crop_before (s));
    free (s);
    int ret;
    kh_put (INT64, ((StackData *)stack_get_data (stack))->messages, (long) n, &ret);
    assert (ret > 0);
  }
  assert (!scope);
  if (level <= warning_verbosity)
    ((StackData *)stack_get_data (stack))->warnings++;
  return NULL;  
}

static inline
void * not_implemented_internal (Ast * scope, Ast * n, const char * file, int line, Stack * stack)
{
  if (stack_verbosity (stack) > error_verbosity) {
    char * s = ast_str_append (n, NULL);
    fprintf (stderr, "%s:%d: warning (interpreter): operation '%s' not implemented\n",
	     file, line, ast_crop_before (s));
    free (s);
  }
  assert (!scope);
  return NULL;
}

#define not_implemented(scope, n, stack) not_implemented_internal (scope, n, __FILE__, __LINE__, stack)

#define value_data(v, type) (*((type *)(v)->data.p))

typedef char Flags;

static inline
bool has_flags (const Value * v)
{
  return v && (v->pointer || (v->type->sym != sym_struct_or_union_specifier &&
			      v->type->sym != sym_CHAR));
}

static inline
Flags value_flags (const Value * v)
{
  if (!has_flags (v))
    return 0;
  Flags f = *((Flags *)(((char *)v->data.p) + v->size - sizeof (Flags)));
#if UNSET_ARRAY
  if (v->unset_array)
    f |= unset;
#endif
  return f;
}

static inline
void value_unset_flags (Value * v, Flags flags)
{
  if (has_flags (v))
    *((Flags *)(((char *)v->data.p) + v->size - sizeof (Flags))) &= ~flags;
}

static inline
void value_set_flags (Value * v, Flags flags)
{
  if (has_flags (v)) {
    *((Flags *)(((char *)v->data.p) + v->size - sizeof (Flags))) |= flags;
    if (flags & unset)
      value_unset_flags (v, constant_expression);
  }
}

static inline
bool is_constant_expression (const Value * v)
{
  //  return v->flags & constant_expression;
  return value_flags (v) & constant_expression;
}

static inline
void set_constant_expression (Value * v)
{
  //  v->flags |= constant_expression;
  value_set_flags (v, constant_expression);
}

static
Value * struct_member_value (Ast * n, Value * value, Ast * member, int index, Stack * stack);

static
void unset_value (Value * value, Stack * stack)
{
  if (value->pointer == 0 && value->type->sym == sym_struct_or_union_specifier) {
    int index = 0;
    Value * v;
    while ((v = struct_member_value ((Ast *)value, value, NULL, index++, stack)))
      unset_value (v, stack);
  }
  else
    value_set_flags (value, unset);
}

static inline
void value_set_write (Value * v, Stack * stack)
{
  StackData * d = stack_get_data (stack);
#if 0  
  if (d->nonlocals && stack_verbosity (stack) > 1) {
    fprintf (stderr, "conditional: %d stack->scope: %d considering:\n", d->conditional, d->scope);
    display_value (v);
  }
#endif
  if (d->nonlocals && v->vscope <= d->conditional  // value is not local to the conditional
#if 1 /** 
      Ignore pointers except function pointers : this is necessary
      for test16.c (i.e. "undefined" allocations of scalars). */
      
      && (!v->pointer || (v->pointer == 1 && v->type->sym == sym_function_definition))
#endif
      ) {
    int ret;
    khiter_t k = kh_put (INT64, d->nonlocals, (long) v->data.p, &ret);
    if (ret) { // new key
      char * p = malloc (sizeof (Value) + v->size);
      memcpy (p, v, sizeof (Value));
      memcpy (p + sizeof (Value), v->data.p, v->size);
      kh_value (d->nonlocals, k) = p;
#if 0
      if (stack_verbosity (stack) > 1) {
	// ast_print_tree ((Ast *)v, stderr, 0, 0, -1);
	display_value (v);
	fprintf (stderr, "conditional: %d stack->scope: %d\n", d->conditional, d->scope);
      }
#endif
    }
  }
}

static
void value_set_parent (Value * v, const Ast * parent)
{
  if (!v)
    return;
  memcpy (v, parent, ast_terminal (parent) ? sizeof (AstTerminal) : sizeof (Ast));
  ((AstTerminal *) v)->value = v;
}

void (* ast_value_print_hook) (const Value * v, FILE * fp) = NULL;

static
void value_print (const Value * v, FILE * fp)
{
  if (!v) {
    fputc ('!', fp);
    return;
  }
  if (v->pointer) {
    if (v->type->sym == sym_function_definition && v->pointer == 1) {
      Ast * n = value_data (v, Ast *);
      if (n) {
	AstTerminal * t = ast_left_terminal (n);
	char * s = ast_str_append (n, NULL);
	fprintf (fp, "%s (%s:%d)", ast_crop_before (s), t->file, t->line);
	free (s);
      }
      else
	fputs ("(fnull)", fp);
    }
    else if (display_pointers)
      fprintf (fp, "%p", value_data (v, void *));
    else
      fprintf (fp, "0xaddress");
  }
  else
    switch (v->type->sym) {
    case sym_CHAR :   fprintf (fp, "%d '%c'", value_data (v, char), value_data (v, char)); break;
    case sym_INT :    fprintf (fp, "%d", value_data (v, int));    break;
    case sym_LONG :   fprintf (fp, "%ld", value_data (v, long));    break;
    case sym_FLOAT :  fprintf (fp, "%g", value_data (v, float));  break;
    case sym_DOUBLE : fprintf (fp, "%g", value_data (v, double)); break;
    case sym_SIGNED : fprintf (fp, "%d", value_data (v, signed)); break;
    case sym_UNSIGNED : fprintf (fp, "%u", value_data (v, unsigned)); break;
    case sym_BOOL :   fprintf (fp, "%d", value_data (v, bool));    break;
    case sym_struct_or_union_specifier : fputs ("{...}", fp); break;
    default:
      fputc ('e', fp);
    }
  if (value_flags (v) & unset)
    fputs (" (unset)", fp);
  if (ast_value_print_hook)
    ast_value_print_hook (v, fp);
}

static
void display_value (const Value * v)
{
  if (!v)
    return;
  fputs ("-------------------", stderr);
  ast_print ((Ast *) v, stderr, 0);
  fputc ('\n', stderr);
  ast_print (v->type, stderr, 0);
  fputc ('\n', stderr);
  fprintf (stderr, "size: %d, flags: %d, pointer: %d, scope: %d, value: ",
	   v->size, value_flags (v), v->pointer, v->vscope);
  value_print (v, stderr);
  fputc ('\n', stderr);
  if (v->pointer && (v->pointer > 1 || v->type->sym != sym_function_definition)) {
    Pointer p = value_data (v, Pointer);
    fprintf (stderr, "p: %p, start: %p, size: %d, scope: %d\nchars: '", p.p, p.start, p.size, p.pscope);
    int maxchars = 80;
    for (char * s = p.p; s < p.start + p.size && maxchars; s++, maxchars--)
      fputc (*s, stderr);
    if (!maxchars)
      fputs ("...", stderr);
    fputs ("'\n", stderr);
  }
  if (v->dimension) {
    fputs ("dimensions:", stderr);
    for (int *i = v->dimension; *i >= 0; i++)
      fprintf (stderr, " %d", *i);
    fputc ('\n', stderr);
  }
  fputs ("===================\n", stderr);
}

static void * undefined (Ast * scope, Ast * n, Stack * stack)
{
  message (scope, n, "undefined argument(s) in '%s'\n", warning_verbosity, stack);
  return NULL;
}

int ast_base_type_size (const Ast * type)
{
  switch (type->sym) {
  case sym_VOID: return 0;
  case sym_SHORT: return sizeof (short);
  case sym_INT: return sizeof (int);
  case sym_LONG: return sizeof (long);
  case sym_FLOAT: return sizeof (float);
  case sym_DOUBLE: return sizeof (double);
  case sym_SIGNED: return sizeof (signed);
  case sym_UNSIGNED: return sizeof (unsigned);
  case sym_BOOL: return sizeof (bool);
  case sym_enum_specifier: return sizeof (int);    
  case sym_COMPLEX: return sizeof (complex double);
  case sym_IMAGINARY:
  default:
    assert (false);
  }
}

int (* ast_type_size) (const Ast *) = ast_base_type_size;

static
Ast * identifier_type (Ast * n, Dimensions * d, Stack * stack);

Ast * base_type (Ast * type, Dimensions * d, Stack * stack)
{
  if (type && (type->sym == sym_TYPEDEF_NAME || type->sym == sym_IDENTIFIER)) {
    Ast * type_identifier = ast_identifier_declaration (stack, ast_terminal (type)->start);
    if (!type_identifier)
      return message (NULL, type, "unknown typedef '%s'\n", warning_verbosity, stack);
    return identifier_type (type_identifier, d, stack);
  }
  else
    return type;
}

static
Ast * struct_size (Ast * type,
		   const char * member, int index,
		   Dimensions * d,
		   int * size, Stack * stack);

// fixme: this is a real mess...

static
int type_size (int pointer, Ast * type, Stack * stack)
{
  if (type && pointer <= 1 && type->sym == sym_function_definition)
    return sizeof (void *) + sizeof (Flags);
  if (pointer)
    return sizeof (Pointer) + sizeof (Flags);
  if (!type)
    return 0;
  if (type->sym == sym_CHAR)
    return sizeof (char);
  Dimensions d = { .size = 1 };
  type = base_type (type, &d, stack);
  // fixme: may not always work
  assert ((!d.pointer || type->sym == sym_function_definition) && d.size == 1);
  if (!type) {
    message (NULL, type, "could not find type of '%s'\n", warning_verbosity, stack);
    return 0;
  }
  if (type->sym == sym_struct_or_union_specifier) {
    if (!ast_child (type, sym_struct_declaration_list))
      return type_size (0, ast_child (type, sym_generic_identifier)->child[0], stack);
    int size;
    struct_size (type, NULL, -1, NULL, &size, stack);
    return size;
  }
  else if (type->sym == sym_function_definition)
    return sizeof (void *) + sizeof (Flags);
  else if (type->sym == sym_CHAR)
    return sizeof (char);
  else
    return ast_type_size (type) + sizeof (Flags);
}

static
Value * new_value (Stack * stack, Ast * n, Ast * type, int pointer)
{
  int size = type_size (pointer, type, stack);
  Value * value = allocate (stack_alloc (stack), sizeof (Value) + size);
  value_set_parent (value, n);
  value->type = type, value->size = size, value->pointer = pointer;
  value->data.p = ((char *)value) + sizeof (Value);
  value->data.start = value->data.p;
  value->data.size = size;
  value->vscope = ((StackData *)stack_get_data (stack))->scope;
#if UNSET_ARRAY
  value->unset_array = NULL;
  value->unset_member = NULL;
#endif
  return value;
}

static
Value * value_copy (Value * v, Allocator * alloc, bool data)
{
  if (!v)
    return NULL;
  Value * copy = allocate (alloc, sizeof (Value) + v->size);
  memcpy (copy, v, sizeof (Value));
  copy->data.p = ((char *)copy) + sizeof (Value);
  copy->data.start = copy->data.p;
  copy->data.size = v->size;
  if (data && v->size)
    memcpy (copy->data.p, v->data.p, v->size);
  if (!data || (value_flags (v) & unset))
    value_set_flags (copy, unset);
  else
    value_unset_flags (copy, unset);
  return copy;
}

static
bool value_bool (const Value * v, bool * error, Stack * stack)
{
  if (v->pointer)
    return value_data (v, void *) != NULL;
  switch (v->type->sym) {
  case sym_CHAR: return value_data (v, char) != 0;
  case sym_INT: return value_data (v, int) != 0;
  case sym_UNSIGNED: return value_data (v, unsigned) != 0;
  case sym_LONG: return value_data (v, long) != 0;
  case sym_FLOAT: return value_data (v, float) != 0.;
  case sym_DOUBLE: return value_data (v, double) != 0.;
  default:
    message (NULL, (Ast *)v, "cannot cast '%s' to bool\n", warning_verbosity, stack);
  }
  *error = true;
  return 0;
}

static
int value_int (Value * v, bool * error, Stack * stack)
{
  if (v->pointer)
    return (long) value_data (v, void *);
  switch (v->type->sym) {
  case sym_CHAR: return value_data (v, char);
  case sym_INT: return value_data (v, int);
  case sym_UNSIGNED: return value_data (v, unsigned);
  case sym_LONG: return value_data (v, long);
  case sym_FLOAT: return value_data (v, float);
  case sym_DOUBLE: return value_data (v, double);
  default:
    message (NULL, (Ast *)v, "cannot cast '%s' to int\n", warning_verbosity, stack);
  }
  *error = true;
  return 0;
}

static
long value_long (Value * v, bool * error, Stack * stack)
{
  if (v->pointer)
    return (long) value_data (v, void *);
  switch (v->type->sym) {
  case sym_CHAR: return value_data (v, char);
  case sym_INT: return value_data (v, int);
  case sym_UNSIGNED: return value_data (v, unsigned);
  case sym_LONG: return value_data (v, long);
  case sym_FLOAT: return value_data (v, float);
  case sym_DOUBLE: return value_data (v, double);
  default:
    message (NULL, (Ast *)v, "cannot cast '%s' to long\n", warning_verbosity, stack);
  }
  *error = true;
  return 0;
}

static
double value_double (Value * v, bool * error, Stack * stack)
{
  switch (v->type->sym) {
  case sym_CHAR: return value_data (v, char);
  case sym_INT: return value_data (v, int);
  case sym_UNSIGNED: return value_data (v, unsigned);
  case sym_LONG: return value_data (v, long);
  case sym_FLOAT: return value_data (v, float);
  case sym_DOUBLE: return value_data (v, double);
  default:
    message (NULL, (Ast *)v, "cannot cast '%s' to double\n", warning_verbosity, stack);
  }
  *error = true;
  return 0;
}

static inline
int array_dimension (const Value * array, int * index)
{
  int nd = 0;
  if (array->dimension)
    for (int * i = array->dimension + 1; *i >= 0; i++)
      *index *= *i, nd++;
  return nd;
}

static
Value * array_member_value (Ast * n, Value * array, int index, int unset, Stack * stack)
{
  if (!array) return NULL;
  if (!array->pointer)
    return message (NULL, (Ast *) array, "'%s' is not a pointer\n", error_verbosity, stack);
  Value * value;
  char * data = value_data (array, void *);
  if (!data) {
#if 0    
    if (stack_verbosity (stack) >= warning_verbosity)
      fprintf (stderr, "data: %p \n", data);
#endif
    return message (NULL, n, "unallocated array access in '%s'\n", warning_verbosity, stack);
  }
  int nd = array_dimension (array, &index);
  value = new_value (stack, n, array->type, array->pointer - 1);
#if UNSET_ARRAY
  if (unset == true)
    value->unset_array = array;
#endif
  Pointer p = value_data (array, Pointer);
  value->vscope = p.pscope;
#if 0
  if (stack_verbosity (stack) > 1) {
    display_value (array);
  }
#endif
  if (nd) {
    assert (value->pointer);
    value_data (value, Pointer) = p;
    data += index*type_size (array->pointer - nd - 1, array->type, stack);
    value_data (value, void *) = data;
  }
  else {
    data += index*value->size;
    value->data.p = data;
    value->data.start = p.start;
    value->data.size = p.size;
    if (p.start && // fixme: p.start should always be defined
	(data < p.start || data + value->size > p.start + p.size)) {
#if 0      
      fprintf (stderr, "data: %p, p.start: %p, p.size: %d\n", data, p.start, p.size);
      display_value (array);
      display_value (value);
#endif
      if (unset == 2)
	return NULL;
      else if (!unset)
	return message (NULL, n, data < p.start ?
			"array index underflow in '%s'\n" :
			"array index overflow in '%s'\n", error_verbosity, stack);	
      else {

	/**
	If the index is unset and there is an underflow or overflow, we
	return the first element of the array. */

	data -= index*value->size;
	if (data < p.start || data + value->size > p.start + p.size)
	  return message (NULL, n, data < p.start ?
			  "array index underflow in '%s'\n" :
			  "array index overflow in '%s'\n", error_verbosity, stack);	
	else
	  value->data.p = data;
      }
    }
  }
#if 0
  if (stack_verbosity (stack) >= warning_verbosity) {
    fprintf (stderr, "==== array member value ====\n");
    display_value (array);
    fprintf (stderr, "%p %d %d\n", data, index, value->size);
    display_value (value);
    fprintf (stderr, "==== end array member value ====\n");
  }
#endif
  if (array->dimension)
    value->dimension = *(array->dimension + 1) >= 0 ? array->dimension + 1 : NULL;
  if (value_flags (array) & unset)
    unset_value (value, stack);
  return value;
}

static
Value * assign (Ast * n, Value * dst, Value * src, Stack * stack)
{
  if (!dst)
    return undefined (NULL, n, stack);
  value_set_flags (dst, unset);  
  if (!src)
    return undefined (NULL, n, stack);
  value_set_write (dst, stack);
  bool error = false;
  if (dst->pointer && src->pointer) {
    if (dst->pointer == 1 && dst->type->sym == sym_function_definition) {
      if (src->pointer == 1)
	memcpy (dst->data.p, src->data.p, dst->size);
      else {
	if (src->type->sym != sym_function_definition || src->pointer != 2)
	  return message (NULL, n, "assignment between incompatible types\n", error_verbosity, stack);
	value_data (dst, Ast *) = *value_data (src, Ast **);
      }
    }
    else {
      assert (src->size == dst->size);
#if 0
      if (stack_verbosity (stack) > 1) {
	ast_print (n, stderr, 0);
	display_value (dst);
      }
#endif
      memcpy (dst->data.p, src->data.p, dst->size);
#if 1
      value_data (dst, Pointer).pscope = fmin (value_data (src, Pointer).pscope,
					       ((StackData *)stack_get_data (stack))->scope);
#endif
    }
  }
  else if (dst->type->sym == src->type->sym && dst->pointer == src->pointer) {
    if (dst->size != src->size)
      return message (NULL, n, "assignment between incompatible types\n", error_verbosity, stack);
    memcpy (dst->data.p, src->data.p, src->size);
  }
  else if (dst->type->sym == sym_DOUBLE)
    value_data (dst, double) = value_double (src, &error, stack);
  else if (dst->type->sym == sym_FLOAT)
    value_data (dst, float) = value_double (src, &error, stack);
  else if (dst->type->sym == sym_CHAR)
    value_data (dst, char) = value_int (src, &error, stack);
  else if (dst->type->sym == sym_INT)
    value_data (dst, int) = value_int (src, &error, stack);
  else if (dst->type->sym == sym_UNSIGNED)
    value_data (dst, unsigned) = value_int (src, &error, stack);
  else if (dst->type->sym == sym_LONG)
    value_data (dst, long) = value_long (src, &error, stack);

  /* fixme: what about value_data (dst, Pointer) ? */
  
  else if (dst->pointer && src->type->sym == sym_INT)
    value_data (dst, void *) = (void *)(long)value_data (src, int);
  else if (dst->pointer && src->type->sym == sym_LONG)
    value_data (dst, void *) = (void *)value_data (src, long);
  else if (dst->type->sym == sym_struct_or_union_specifier &&
	   ((src->type == (Ast *)&ast_int && value_data (src, int) == 0) ||
	    (src->type == (Ast *)&ast_long && value_data (src, long) == 0)))
    memset (dst->data.p, 0, dst->size);
  else
    return message (NULL, n, "assignment between incompatible types\n", error_verbosity, stack);
  if (!(value_flags (src) & unset))
    value_unset_flags (dst, unset);
#if UNSET_ARRAY
  int a = 0;
  if (dst->unset_array && !array_dimension (dst->unset_array, &a)) {
    int index = 0;
    Value * v;
    while ((v = array_member_value ((Ast *)dst->unset_array, dst->unset_array, index++, 2, stack))) {
      if (dst->unset_member)
	v = struct_member_value ((Ast *)v, v, dst->unset_member, -1, stack);
      unset_value (v, stack);
    }
  }
#endif
  value_unset_flags (dst, constant_expression);
  return ast_assign_hook ? ast_assign_hook (n, dst, src, stack) : src;
}

static
Value * pointer_unary_operation (Ast * n, Ast * op, Value * a, Stack * stack)
{
  if (!a)
    return undefined (NULL, n, stack);

  if (value_flags (a) & unset)
    return message (NULL, n, "unary operation '%s' on unset pointer\n", warning_verbosity, stack);
    
  char * data = value_data (a, void *);

  switch (op->sym) {
    
  case sym_INC_OP:
    value_set_write (a, stack);
    value_data (a, void *) = data + type_size (a->pointer - 1, a->type, stack);
    return a;
    
  case sym_DEC_OP:
    value_set_write (a, stack);
    value_data (a, void *) = data - type_size (a->pointer - 1, a->type, stack);
    return a;

  case sym_unary_operator:
    if (op->child[0]->sym == token_symbol ('*')) {
      if (a->type->sym == sym_function_definition && a->pointer == 1)
	return a;
      return array_member_value (n, a, 0, false, stack);
    }
    else if (op->child[0]->sym == token_symbol ('!')) {
      Value * v = new_value (stack, n, (Ast *)&ast_int, 0);
      value_data (v, int) = !data;
      return v;
    }
    else if (op->child[0]->sym == token_symbol ('&')) {
      Value * v = new_value (stack, n, a->type, a->pointer + 1);
      v->vscope = a->vscope;
      value_data (v, void *) = a->data.p;
      value_data (v, Pointer).start = a->data.start;
      value_data (v, Pointer).size = a->data.size;
      value_data (v, Pointer).pscope = a->vscope;
      return v;
    }

  default:
    return not_implemented (NULL, n, stack);
    
  }

  return NULL;
}

static
Value * unary_operation (Ast * n, Ast * op, Value * a, Stack * stack)
{  
  if (!a)
    return undefined (NULL, n, stack);

  if (a->pointer)
    return pointer_unary_operation (n, op, a, stack);

  Value * v = NULL;
  switch (op->sym) {
    
  case sym_INC_OP:
    value_set_write (a, stack);
    if (a->type->sym == sym_INT)
      ++value_data (a, int);
    else if (a->type->sym == sym_LONG)
      ++value_data (a, long);
    else if (a->type->sym == sym_FLOAT)
      ++value_data (a, float);
    else if (a->type->sym == sym_DOUBLE)
      ++value_data (a, double);
    else
      return not_implemented (NULL, n, stack);
    return a;
    
  case sym_DEC_OP:
    value_set_write (a, stack);
    if (a->type->sym == sym_INT)
      --value_data (a, int);
    else if (a->type->sym == sym_LONG)
      --value_data (a, long);
    else if (a->type->sym == sym_FLOAT)
      --value_data (a, float);
    else if (a->type->sym == sym_DOUBLE)
      --value_data (a, double);
    else
      return not_implemented (NULL, n, stack);
    return a;
    
  case sym_unary_operator:
    switch (ast_terminal (op->child[0])->start[0]) {
    
    case '&':
      v = new_value (stack, n, a->type, a->pointer + 1);
      v->vscope = a->vscope;
      value_data (v, void *) = a->data.p;
      value_data (v, Pointer).start = a->data.start;
      value_data (v, Pointer).size = a->data.size;
      value_data (v, Pointer).pscope = a->vscope;
      break;
            
    case '+':
      return a;

    case '-': {
      v = value_copy (a, stack_alloc (stack), true);
      if (value_flags (a) & unset)
	value_set_flags (v, unset);
      else
	switch (a->type->sym) {
	case sym_INT:
	  value_data (v, int) = - value_data (a, int);
	  break;
	case sym_LONG:
	  value_data (v, long) = - value_data (a, long);
	  break;
	case sym_FLOAT:
	  value_data (v, float) = - value_data (a, float);
	  break;
	case sym_DOUBLE:
	  value_data (v, double) = - value_data (a, double);
	  break;
	default:
	  return not_implemented (NULL, n, stack);
	}
      break;
    }

    case '!': {
      v = new_value (stack, n, (Ast *)&ast_int, 0);      
      if (value_flags (a) & unset)
	value_set_flags (v, unset);
      else {
	if (is_constant_expression (a))
	  set_constant_expression (v);
	switch (a->type->sym) {
	case sym_INT:
	  value_data (v, int) = !value_data (a, int);
	  break;
	case sym_LONG:
	  value_data (v, int) = !value_data (a, long);
	  break;
	case sym_FLOAT:
	  value_data (v, int) = !value_data (a, float);
	  break;
	case sym_DOUBLE:
	  value_data (v, int) = !value_data (a, double);
	  break;
	default:
	  return not_implemented (NULL, n, stack);
	}
      }
      break;
    }
      
    case '~': {
      v = new_value (stack, n, a->type, 0);
      if (value_flags (a) & unset)
	value_set_flags (v, unset);
      else {     
	if (is_constant_expression (a))
	  set_constant_expression (v);
	if (a->type->sym == sym_INT)
	  value_data (v, int) = ~value_data (a, int);
	else if (a->type->sym == sym_LONG)
	  value_data (v, long) = ~value_data (a, long);
	else
	  return not_implemented (NULL, n, stack);
      }
      break;
    }
      
    default:
      return not_implemented (NULL, n, stack);
      
    }
    break;

  default:
    return not_implemented (NULL, n, stack);
    
  }

  return v;
}

#define GET_VALUE(a, va)			\
  if (a->type->sym == sym_INT)			\
    va = value_data (a, int);			\
  else if (a->type->sym == sym_UNSIGNED)	\
    va = value_data (a, unsigned);		\
  else if (a->type->sym == sym_LONG)		\
    va = value_data (a, long);			\
  else if (a->type->sym == sym_FLOAT)		\
    va = value_data (a, float);			\
  else if (a->type->sym == sym_DOUBLE)		\
    va = value_data (a, double);		\
  else if (a->type->sym == sym_CHAR)		\
    va = value_data (a, char);			\
  else						\
    return not_implemented (NULL, n, stack)

static
Value * no_hook (Ast * n, Stack * stack, Value * a, Value * b, Value * value) { return value; }

/* fixme: change hooking to global */
Value * (* ast_binary_operation_hook) (Ast *, Stack *, Value *, Value *, Value *) = no_hook;

static
Value * binary_operation (Ast * n, Stack * stack)
{
  int op = n->child[1]->sym;
  
  /**
  And/Or operations require special treatment due to partial evaluation. */

  if (op == sym_AND_OP || op == sym_OR_OP) {
    Value * a = run (n->child[0], stack);
    if (!a)
      return undefined (NULL, n->child[0], stack);
    Value * value = new_value (stack, n, (Ast *) &ast_int, 0);
    bool error;
    bool va = value_bool (a, &error, stack);
    if (va && op == sym_OR_OP && !(value_flags (a) & unset)) {
      value_data (value, int) = true;
      if (is_constant_expression (a))
	set_constant_expression (value);
      return value;
    }
    if (!va && op == sym_AND_OP && !(value_flags (a) & unset)) {
      if (is_constant_expression (a))
	set_constant_expression (value);
      return value;
    }
    Value * b = run (n->child[2], stack);
    if (!b)
      return undefined (NULL, n->child[2], stack);
    if (is_constant_expression (a) && is_constant_expression (b))
      set_constant_expression (value);
    bool vb = value_bool (b, &error, stack);
    if (vb && op == sym_OR_OP && !(value_flags (b) & unset)) {
      value_data (value, int) = true;
      return value;
    }
    if (!vb && op == sym_AND_OP && !(value_flags (b) & unset))
      return value;
    if ((value_flags (a) & unset) ||
	(value_flags (b) & unset))
      value_set_flags (value, unset);
    else if (op == sym_AND_OP && va && vb)
      value_data (value, int) = true;
    return value;
  }

  /**
  In other cases we must evaluate both operands. */
  
  Value * a = run (n->child[0], stack), * b = run (n->child[2], stack);
  
  if (!a)
    return undefined (NULL, n->child[0], stack);
  if (!b)
    return undefined (NULL, n->child[2], stack);

  /**
  Pointer arithmetic. */

  if (a->pointer || b->pointer) {

    if (op == sym_EQ_OP || op == sym_NE_OP) {
      Value * value = ast_binary_operation_hook (n, stack, a, b, new_value (stack, n, (Ast *) &ast_int, 0));
      if ((value_flags (a) & unset) || (value_flags (b) & unset)) {
	value_set_flags (value, unset);
	return value;
      }
      bool erra = false, errb = false;
      switch (op) {
      case sym_EQ_OP:
	value_data (value, int) =
	  value_long (a, &erra, stack) == value_long (b, &errb, stack);
	break;
      case sym_NE_OP:
	value_data (value, int) =
	  value_long (a, &erra, stack) != value_long (b, &errb, stack);
	break;
      }
      if (erra || errb)
	return not_implemented (NULL, n, stack);      
      return value;
    }

    if (a->pointer && b->pointer) {
      if (op == token_symbol ('-')) {
	
	Value * value = ast_binary_operation_hook (n, stack, a, b, new_value (stack, n, (Ast *) &ast_long, 0));
	if ((value_flags (a) & unset) || (value_flags (b) & unset)) {
	  value_set_flags (value, unset);
	  return value;
	}
	// fixme: will work only for char * pointers!!!
	char * va = value_data (a, char *);
	char * vb = value_data (b, char *);
	value_data (value, long) = va - vb;
	return value;
      }
    }
    else {
      bool error = false;
      int vb = value_int (b, &error, stack);
      if (a->pointer && !error) {
	Value * value;
	if (n->sym != sym_assignment_expression)
	  value = ast_binary_operation_hook (n, stack, a, b, value_copy (a, stack_alloc (stack), true));
	else {
	  value = ast_binary_operation_hook (n, stack, a, b, a);
	  value_set_write (value, stack);
	}
	if ((value_flags (a) & unset) || (value_flags (b) & unset)) {
	  value_set_flags (value, unset);
	  return value;
	}
	char * data = value_data (a, void *);
	vb *= type_size (a->pointer - 1, a->type, stack);
	if (op == token_symbol('+') ||
	    ast_schema (n->child[1], sym_assignment_operator, 0, sym_ADD_ASSIGN))
	  data += vb;
	else if (op == token_symbol('-') ||
		 ast_schema (n->child[1], sym_assignment_operator, 0, sym_SUB_ASSIGN))
	  data -= vb;
	else
	  return message (NULL, n, "invalid pointer arithmetic operation '%s'\n", error_verbosity, stack);
	value_data (value, void *) = data;
	return value;
      }
    }
    
    return not_implemented (NULL, n, stack);  
  }

  /**
  Float/integer arithmetic. */

#define OPERATION(res, va, vb)						\
  set = true;								\
  if (op == token_symbol('*') ||					\
      ast_schema (n->child[1], sym_assignment_operator, 0, sym_MUL_ASSIGN)) \
    res = va*vb;							\
  else if (op == token_symbol('/') ||					\
	   ast_schema (n->child[1], sym_assignment_operator, 0, sym_DIV_ASSIGN)) \
    res = vb ? va/vb : 0.;						\
  else if (op == token_symbol('+') ||					\
	   ast_schema (n->child[1], sym_assignment_operator, 0, sym_ADD_ASSIGN)) \
    res = va + vb;							\
  else if (op == token_symbol('-') ||					\
	   ast_schema (n->child[1], sym_assignment_operator, 0, sym_SUB_ASSIGN)) \
    res = va - vb;							\
  else									\
    set = false
  
#define LOGICAL_OPERATION(res, va, vb)		\
  set = true;					\
  if (op == token_symbol('<'))			\
    res = va < vb;				\
  else if (op == token_symbol('>'))		\
    res = va > vb;				\
  else if (op == sym_LE_OP)			\
    res = va <= vb;				\
  else if (op == sym_GE_OP)			\
    res = va >= vb;				\
  else if (op == sym_EQ_OP)			\
    res = va == vb;				\
  else if (op == sym_NE_OP)			\
    res = va != vb;				\
  else						\
    set = false

#define BITWISE_OPERATION(res, va, vb)		\
  set = true;					\
  if (op == sym_LEFT_OP)			\
    res = va << vb;				\
  else if (op == sym_RIGHT_OP)			\
    res = va >> vb;				\
  else if (op == token_symbol('&'))		\
    res = va & vb;				\
  else if (op == token_symbol('^'))		\
    res = va ^ vb;				\
  else if (op == token_symbol('|'))		\
    res = va | vb;				\
  else if (op == token_symbol('%'))		\
    res = vb ? va % vb : 0;			\
  else						\
    set = false

#define SET_VALUE(a, va)			\
  if (a->type->sym == sym_INT)			\
    value_data (a, int) = va;			\
  else if (a->type->sym == sym_UNSIGNED)	\
    value_data (a, unsigned) = va;		\
  else if (a->type->sym == sym_LONG)		\
    value_data (a, long) = va;			\
  else if (a->type->sym == sym_FLOAT)		\
    value_data (a, float) = va;			\
  else if (a->type->sym == sym_DOUBLE)		\
    value_data (a, double) = va;		\
  else						\
    return not_implemented (NULL, n, stack)

#define GENERIC_OPERATION(type) do {					\
    type res;								\
    bool set;								\
    GET_VALUE (a, va);							\
    GET_VALUE (b, vb);							\
    OPERATION (res, va, vb);						\
    if (set) {								\
      Value * value;							\
      if (n->sym == sym_assignment_expression)				\
	value_set_write (a, stack), value = ast_binary_operation_hook (n, stack, a, b, a); \
      else {								\
	value = new_value (stack, n, (Ast *) &ast_##type, 0);		\
	if (is_constant_expression (a) && is_constant_expression (b))	\
          set_constant_expression (value);				\
	value = ast_binary_operation_hook (n, stack, a, b, value);	\
      }									\
      SET_VALUE (value, res);						\
      if ((value_flags (a) & unset) || (value_flags (b) & unset))	\
	value_set_flags (value, unset);					\
      return value;							\
    }									\
    else {								\
      int res;								\
      LOGICAL_OPERATION (res, va, vb);					\
      if (set) {							\
	Value * value = new_value (stack, n, (Ast *) &ast_int, 0);	\
	if (is_constant_expression (a) && is_constant_expression (b))	\
          set_constant_expression (value);				\
	ast_binary_operation_hook (n, stack, a, b, value);		\
	SET_VALUE (value, res);						\
	if ((value_flags (a) & unset) || (value_flags (b) & unset))	\
	  value_set_flags (value, unset);				\
	return value;							\
      }									\
    }									\
  } while (0)

#define INTEGER_OPERATION(type) do {					\
    type res;								\
    bool set;								\
    BITWISE_OPERATION (res, va, vb);					\
    if (set) {								\
      Value * value = new_value (stack, n, (Ast *) &ast_##type, 0);	\
      if (is_constant_expression (a) && is_constant_expression (b))	\
	set_constant_expression (value);				\
      value = ast_binary_operation_hook (n, stack, a, b, value);	\
      SET_VALUE (value, res);						\
      if ((value_flags (a) & unset) || (value_flags (b) & unset))	\
	value_set_flags (value, unset);					\
      return value;							\
    }									\
  } while (0)

  if (a->type->sym == sym_DOUBLE || b->type->sym == sym_DOUBLE) {
    double va, vb;
    GENERIC_OPERATION (double);
  }
  else if (a->type->sym == sym_FLOAT || b->type->sym == sym_FLOAT) {
    float va, vb;
    GENERIC_OPERATION (float);
  }
  else if (a->type->sym == sym_LONG || b->type->sym == sym_LONG) {
    long va, vb;
    GENERIC_OPERATION (long); 
    INTEGER_OPERATION (long);
  }
  else {
    int va, vb;
    GENERIC_OPERATION (int);
    INTEGER_OPERATION (int);
  }
  
  return not_implemented (NULL, n, stack);
}

#define has_value(identifier) (ast_terminal (identifier)->value == (identifier))

static
Value * constant_value (Ast * n, Stack * stack)
{
#if 0 // fixme: really necessary? apparently not
  if (ast_terminal (n)->value) // constants have a "static" value
    return ast_terminal (n)->value;
#endif
  switch (n->sym) {
  case sym_I_CONSTANT: {
    Value * v = new_value (stack, n, (Ast *) &ast_long, 0);
    set_constant_expression (v);
    if (ast_terminal (n)->start[0] == '\'') {
      // character constant
      if (ast_terminal (n)->start[1] == '\\') {
	switch (ast_terminal (n)->start[2]) {
	case '\'': value_data (v, long) = '\''; break;
	case '"': value_data (v, long) = '\"'; break;
	case '?': value_data (v, long) = '\?'; break;
	case '\\': value_data (v, long) = '\\'; break;
	case 'a': value_data (v, long) = '\a'; break;
	case 'b': value_data (v, long) = '\b'; break;
	case 'f': value_data (v, long) = '\f'; break;
	case 'n': value_data (v, long) = '\n'; break;
	case 'r': value_data (v, long) = '\r'; break;
	case 't': value_data (v, long) = '\t'; break;
	case 'v': value_data (v, long) = '\v'; break;
	default:
	  value_data (v, long) = atol (ast_terminal (n)->start + 2);
	}
      }
      else
	value_data (v, long) = ast_terminal (n)->start[1];
    }
    else
      value_data (v, long) = atol (ast_terminal (n)->start);
#if 0
    {
      AstTerminal * t = ast_terminal (n);
      if (t->file && !strcmp (t->file, BASILISK "/grid/events.h"))
	fprintf (stderr, "%s:%d: %s %p %p\n", t->file, t->line, t->start, n, v);
    }
#endif
    ast_terminal ((Ast *)v)->value = ast_terminal (n)->value = v;
    return v;
  }
  case sym_F_CONSTANT: {
    Value * v = new_value (stack, n, (Ast *) &ast_double, 0);
    set_constant_expression (v);
    value_data (v, double) = atof (ast_terminal (n)->start);
    ast_terminal ((Ast *)v)->value = ast_terminal (n)->value = v;
    return v;
  }
  case sym_STRING_LITERAL: case sym_FUNC_NAME: {
    Value * v = new_value (stack, n, (Ast *) &ast_char, 1);
    int size = strlen (ast_terminal (n)->start) - 1;
    Pointer p = {0};
    p.size = size*sizeof (char);
    p.p = allocate (((StackData *)stack_get_data (stack))->static_alloc, p.size);
    p.start = p.p;
    value_data (v, Pointer) = p;
    memcpy (value_data (v, char *), ast_terminal (n)->start + 1, (size - 1)*sizeof (char));
    ast_terminal ((Ast *)v)->value = ast_terminal (n)->value = v;
    return v;
  }
  default:
    assert (false);
  }
  return NULL;
}

static Ast * get_array_dimensions (Ast * direct_declarator, int symbol, Dimensions * d, int nd, Stack * stack)
{
  assert (d->dimension == NULL);
  d->dimension = allocate (stack_alloc (stack), (nd + 1)*sizeof (int));
  d->dimension[nd] = -1;
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
    Value * s;
    if (nelem) {
      s = run (nelem, stack);
      if (!s)
	return message (NULL, nelem, "undefined array size '%s'\n", warning_verbosity, stack);
      if (value_flags (s) & unset)
	return message (NULL, nelem, "array size '%s' not set\n", warning_verbosity, stack);
    }
    else if (ast_schema (direct_declarator->parent, symbol,
			 1, token_symbol (']')) ||
	     ast_schema (direct_declarator->parent, symbol,
			 2, token_symbol (']')))
      s = NULL;
    else
      // fixme: more complex size syntaxes not handled
      return not_implemented (NULL, direct_declarator, stack);
    if (s) {
      bool error = false;
      d->dimension[nd] = value_int (s, &error, stack);
      d->size *= d->dimension[nd];
      if (error)
	return NULL;
    }
    else {
      if (nd != 0)
	return message (NULL, direct_declarator,
			"only the first dimension of a multidimensional array can be undefined\n",
			1, stack);
      d->dimension[0] = 0;
      d->size = -1;
    }
    nd++;
    direct_declarator = direct_declarator->parent;
  }
  return direct_declarator;
}

static
Ast * direct_declarator_type (Ast * direct_declarator, Dimensions * d, Stack * stack)
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
    if (nd && !(direct_declarator = get_array_dimensions (direct_declarator, sym_direct_declarator, d, nd, stack)))
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

static
Ast * identifier_type (Ast * n, Dimensions * d, Stack * stack)
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

  return base_type (type, d, stack);
}

static
Ast * return_type (Ast * function_declaration, Dimensions * d, Stack * stack)
{
  if (!function_declaration)
    return NULL;
  
  Ast * declarator = ast_child (function_declaration, sym_declarator);
  while (declarator) {
    Ast * pointer = declarator;
    while ((pointer = ast_child (pointer, sym_pointer)))
      d->pointer++;
    if (ast_child (declarator, token_symbol ('[')))
      d->pointer++;
    if (ast_child (declarator, sym_direct_declarator))
      declarator = ast_child (declarator, sym_direct_declarator);
    else
      declarator = ast_child (declarator, sym_declarator);
  }
  
  Ast * specifiers = ast_child (function_declaration, sym_declaration_specifiers);
  assert (specifiers);
  Ast * type = ast_find (specifiers, sym_types);
  assert (type);
  type = type->child[0];
  if (type->sym == sym_struct_or_union_specifier &&
      !ast_child (type, sym_struct_declaration_list)) {
    if (ast_child (type, sym_YYerror))
      return NULL;
    type = ast_child (type, sym_generic_identifier)->child[0];
  }
  return base_type (type, d, stack);
}

static
Ast * type_name_type (Ast * type_name, Dimensions * d, Stack * stack)
{
  Ast * declarator = ast_child (type_name, sym_abstract_declarator);
  while (declarator) {
    Ast * pointer = declarator;
    while ((pointer = ast_child (pointer, sym_pointer)))
      d->pointer++;
    if (ast_child (declarator, token_symbol ('[')))
      d->pointer++;
    if (ast_child (declarator, sym_direct_abstract_declarator))
      declarator = ast_child (declarator, sym_direct_abstract_declarator);
    else
      declarator = ast_child (declarator, sym_abstract_declarator);
  }
  
  Ast * specifiers = ast_schema (type_name, sym_type_name,
				 0, sym_specifier_qualifier_list);
  assert (specifiers);
  Ast * type = ast_find (specifiers, sym_types);
  assert (type);
  type = type->child[0];
  
  if (type->sym == sym_struct_or_union_specifier &&
      !ast_child (type, sym_struct_declaration_list)) {
    if (ast_child (type, sym_YYerror))
      return NULL;
    type = ast_child (type, sym_generic_identifier)->child[0];
  }

  return base_type (type, d, stack);
}

static
Ast * struct_size (Ast * type,
		   const char * member, int index,
		   Dimensions * d,
		   int * size, Stack * stack)
{
  Ast * list = ast_child (type, sym_struct_declaration_list);
  while (!list) {
    type = base_type (ast_child (type, sym_generic_identifier)->child[0], d, stack);
    if (!type)
      return NULL;
    assert (type->sym == sym_struct_or_union_specifier);
    list = ast_child (type, sym_struct_declaration_list);
  }
  *size = 0;
  int i = 0;
  foreach_item_r (list, sym_struct_declaration, item) {
    Ast * list = ast_child (item, sym_struct_declarator_list);
    foreach_item_r (list, sym_struct_declarator, item) {
      if (ast_child (item, sym_constant_expression))
	return not_implemented (NULL, type, stack);
      Ast * declarator = ast_child (item, sym_declarator);
      Ast * direct_declarator = ast_child (declarator, sym_direct_declarator);
      while (direct_declarator->child[0]->sym == sym_direct_declarator ||
	     direct_declarator->child[0]->sym == token_symbol('('))
	direct_declarator = direct_declarator->child[0]->sym == sym_direct_declarator ?
	  direct_declarator->child[0] :
	  ast_child (direct_declarator->child[1], sym_direct_declarator);
      Dimensions dm = { .size = 1 };
      Ast * type = base_type (direct_declarator_type (direct_declarator, &dm, stack),
			      &dm, stack);
      Ast * identifier;
      if ((!member && i == index) ||
	  (member &&
	   ((identifier = ast_find (direct_declarator->child[0], sym_IDENTIFIER)) ||
	    (identifier = ast_find (direct_declarator->child[0], sym_TYPEDEF_NAME))) &&
	   !strcmp (member, ast_terminal (identifier)->start))) {
	*d = dm;
	return type;
      }
      *size += type_size (dm.pointer, type, stack)*dm.size;
      i++;
    }
  }
  return NULL;
}

static
Value * struct_member_value (Ast * n, Value * value, Ast * member, int index, Stack * stack)
{
  int offset;
  Dimensions d = { .size = 1 };
  Ast * type = struct_size (value->type, member ? ast_terminal (member)->start : NULL,
			    index, &d, &offset, stack);
  if (!type) {
    if (member)
      return message (NULL, member, "could not find type of '%s'\n", warning_verbosity, stack);
    return NULL;
  }

  if (d.dimension) {
    assert (d.dimension[0] >= 0);
    d.size *= type_size (d.pointer, type, stack);
    for (int * i = d.dimension; *i >= 0; i++)
      d.pointer++;
  }
  
  Value * v = new_value (stack, n, type, d.pointer);
  
  char * data;
  if (value->pointer) {
    if (value_flags (value) & unset)
      return message (NULL, n, "undefined structure pointer in '%s'\n", warning_verbosity, stack);
    data = value_data (value, char *) + offset;
    Pointer p = value_data (value, Pointer);
    if (data < p.start || data + v->size > p.start + p.size)
      return message (NULL, n, data < p.start ?
		      "structure pointer underflow in '%s'\n" :
		      "structure pointer overflow in '%s'\n",
		      1, stack);
    v->vscope = p.pscope;
  }
  else {
    data = ((char *)value->data.p) + offset;
    v->vscope = value->vscope;
  }

  if (d.dimension) {
    v->dimension = d.dimension;
    value_data (v, void *) = data;
    if (d.size > 0)
      value_data (v, Pointer).start = data,
	value_data (v, Pointer).size = d.size;
    else
      value_set_flags (v, unset);
    value_data (v, Pointer).pscope = v->vscope;
  }
  else {
    v->data.p = data;
    if (value->pointer) {
      v->data.start = value_data (value, Pointer).start;
      v->data.size = value_data (value, Pointer).size;
    }
    else {
      v->data.start = value->data.start;
      v->data.size = value->data.size;
    }
    assert (((char *)v->data.p) + v->size <= v->data.start + v->data.size);
  }
#if UNSET_ARRAY
  if (value->unset_array) {
    v->unset_array = value->unset_array;
    v->unset_member = member;
  }
#endif
  Ast * attributes;
  if (v->vscope == 0 && member && (value_flags (v) & unset) &&
      (attributes = ast_schema (ast_parent (value->type, sym_declaration), sym_declaration,
				1, sym_init_declarator_list,
				0, sym_init_declarator,
				0, sym_declarator,
				0, sym_direct_declarator,
				0, sym_generic_identifier,
				0, sym_IDENTIFIER)) &&
      !strcmp (ast_terminal (attributes)->start, "_Attributes") &&
      !strcmp (ast_terminal (member)->start, "freed"))
    value_unset_flags (v, unset);
  
  return v;
}

static
void print_expression (Ast * n, Value * value)
{
  AstTerminal * t = ast_left_terminal (n);
  if (t->file && t->line) {
    char * s = ast_str_append (n, NULL);
    fprintf (stderr, "%s:%d: %s: ", t->file, t->line, ast_crop_before (s));
    value_print (value, stderr);
    fputc ('\n', stderr);
    free (s);
  }
}

static
void print_function_call (Ast * n, FILE * fp, Stack * stack)
{
  char * s = ast_str_append (n->child[0], NULL);
  fprintf (fp, "%s", ast_crop_before (s));
  free (s);
  Ast * arguments = ast_schema (n, sym_function_call,
				2, sym_argument_expression_list);
  if (arguments) {
    fputs (" (", fp);
    bool first = true;
    foreach_item_r (arguments, sym_argument_expression_list_item, arg) {
      if (!first)
	fputs (", ", fp);
      first = false;
      value_print (run (arg, stack), fp);
    }
    fputc (')', fp);
  }
  else
    fputs ("()", fp);
}

static bool allocate_array (Stack * stack, Value * value, Ast * initializer)
{
  if (value->pointer && !value_data (value, void *)) {
    int size = 0;
    Ast * list = initializer->child[1];
    foreach_item (list, 2, i) size++;
    int pointer = value->pointer;
    if (value->dimension) {
      assert (value->dimension[0] == 0);
      value->dimension[0] = size;
      size = 1;
      for (int * i = value->dimension; *i >= 0; i++)
	pointer--, size *= *i;
      assert (size > 0);
      pointer++;
    }
    value_data (value, Pointer) = stack_allocate (stack, size*type_size (pointer - 1, value->type, stack));
    value_unset_flags (value, unset);
#if 0
    if (stack_verbosity (stack) > error_verbosity) {
      fprintf (stderr, "--- Array allocation ---\n");
      display_value (value);
      fprintf (stderr, "--- End Array allocation ---\n");
    }
#endif
    return true;
  }
  return false;
}

static
void initialize (Stack * stack, Ast * n, Value * value, Ast * initializer);

static
int initialize_struct_member (Stack * stack, Value * value, int index, Ast * initializer)
{
  int cindex = ast_child_index (initializer);
  if (cindex == 0 || cindex == 2) {
    Value * v = struct_member_value ((Ast *)value, value, NULL, index, stack);
    if (!v) {
      message (NULL, initializer, "too many elements in initializer?\n", error_verbosity, stack);
      return index + 1;
    }
    initialize (stack, initializer, v, initializer);
  }
  else { // designation
    Ast * designation = ast_child (initializer->parent, sym_designation);
    // fixme: this only handles simple designators i.e. '.identifier = ...'
    Ast * designator = ast_find (designation, sym_IDENTIFIER);
    if (!designator)
      designator = ast_find (designation, sym_TYPEDEF_NAME);
    Value * v = struct_member_value (designator, value, designator, -1, stack);
    if (!v) {
      message (NULL, designator, "unknown structure member '%s'\n", error_verbosity, stack);
      return index + 1;
    }
    initialize (stack, designator, v, initializer);
  }
  return index + 1;
}

static
void initialize (Stack * stack, Ast * n, Value * value, Ast * initializer)
{
  if (!value)
    return;
#if 0
  if (stack_verbosity (stack) > error_verbosity) {
    ast_print_tree (n, stderr, 0, 0, -1);
    display_value (value);
  }
#endif
  if (initializer && initializer->child[0]->sym == sym_assignment_expression)
    assign (n, value, run (initializer->child[0], stack), stack);
  else if (value->pointer) {
    int index = 0;
    if (initializer) {
      allocate_array (stack, value, initializer);
      Ast * list = initializer->child[1];
      foreach_item_r (list, sym_initializer, j) {
	int cindex = ast_child_index (j);
	if (cindex == 0 || cindex == 2) {
	  Value * v = array_member_value ((Ast *)value, value, index, false, stack);
	  if (!v) {
	    message (NULL, n, "too many elements in initializer?\n", error_verbosity, stack);
	    return;
	  }
	  initialize (stack, j, v, j);
	}
	index++;
      }
    }
  }
  else if (value->type->sym == sym_struct_or_union_specifier) {
    int index = 0;
    if (initializer) {
      Ast * list = initializer->child[1];
      foreach_item_r (list, sym_initializer, initializer)
	index = initialize_struct_member (stack, value, index, initializer);
    }

    /**
    We set to "zero", the remaining members of the structure. */

    Value * v;
    int verbosity = stack_verbosity (stack);
    ((StackData *)stack_get_data (stack))->verbosity = 0;
    while ((v = struct_member_value ((Ast *)value, value, NULL, index, stack))) {
      ((StackData *)stack_get_data (stack))->verbosity = verbosity;
      initialize (stack, n, v, NULL);
      ((StackData *)stack_get_data (stack))->verbosity = 0;
      index++;
    }
    ((StackData *)stack_get_data (stack))->verbosity = verbosity;
  }
  else if (!initializer) {
    Ast * parent = n->parent;
    n->parent = ((StackData *)stack_get_data (stack))->root;
    initializer = ast_new_constant (n, sym_I_CONSTANT, "0");
    n->parent = parent;
    assign (n, value, run (initializer, stack), stack);
    ast_destroy (initializer);
  }
  else
    message (NULL, n, "lists can only initialize structures or arrays\n", error_verbosity, stack);
}

static
Value * value_from_type (Ast * n, Ast * type, Dimensions * d, Ast * initializer, Stack * stack)
{

  /**
  Array allocation. */

  Value * v;
  if (d->dimension) {
    assert (d->dimension[0] >= 0);
    d->size *= type_size (d->pointer, type, stack);
    for (int * i = d->dimension; *i >= 0; i++)
      d->pointer++; 
    v = new_value (stack, n, type, d->pointer);
    v->dimension = d->dimension;
    if (d->size > 0)
      value_data (v, Pointer) = stack_allocate (stack, d->size);
    else
      value_set_flags (v, unset);
  }

  /**
  Standard allocation. */
  
  else {
    v = new_value (stack, n, type, d->pointer);
    if (type->sym == sym_function_definition && d->pointer == 1)
      value_data (v, Ast *) = n;
    else if (type == (Ast *) &ast_enum) {
      Ast * enumerator = ast_parent (n, sym_enumerator);
      Ast * constant = ast_child (enumerator, sym_constant_expression);
      if (constant)
	assign (n, v, run (constant, stack), stack);
      else {
	Ast * previous = ast_child (ast_child (ast_ancestor (enumerator, 1), sym_enumerator_list),
				    sym_enumerator);
	if (previous) {
	  Ast * identifier =
	    ast_identifier_declaration (stack, ast_terminal (ast_find (previous, sym_IDENTIFIER))->start);
	  assert (has_value (identifier));
	  value_data (v, int) = value_data ((Value *)identifier, int) + 1;
	}
      }
    }
    else if (!initializer)
      unset_value (v, stack);
  }

  /**
  Initialization */

  if (initializer)
    initialize (stack, initializer->parent, v, initializer);
  
  return v;
}

static
Value * identifier_value (Ast * n, Stack * stack)
{
  Dimensions d = { .size = 1 };
  Ast * type = identifier_type (n, &d, stack);
  if (!type)
    return message (NULL, n, "could not find type of '%s'\n", warning_verbosity, stack);

  if (type->sym == sym_IDENTIFIER) {
    type = ast_identifier_declaration (stack, ast_terminal (type)->start);
    if (!type || ast_ancestor (type, 2)->sym != sym_struct_or_union_specifier)
      return message (NULL, n, "could not find type of '%s'\n", warning_verbosity, stack);
    type = ast_ancestor (type, 2);
  }

  return value_from_type (n, type, &d, ast_child (ast_parent (n, sym_init_declarator), sym_initializer), stack);
}

static
Ast ** function_parameters (Ast * identifier)
{
  Ast * parameters = ast_find (ast_parent (identifier, sym_declarator), sym_direct_declarator,
			       2, sym_parameter_type_list,
			       0, sym_parameter_list);
  int n = 0;
  foreach_item (parameters, 2, param) n++;
  Ast ** params = malloc ((n + 1)*sizeof (Ast *));
  params[n] = NULL;
  foreach_item (parameters, 2, param)
    params[--n] = param;
  return params;
}

static
void * pointer_check (Ast * n, Value * v, long size, Stack * stack)
{
  Pointer p = value_data (v, Pointer);
  if ((char *) p.p < p.start || (char *) p.p + size > p.start + p.size)
    return message (NULL, n, (char *) p.p < p.start ?
		    "pointer underflow in '%s'\n" :
		    "pointer overflow in '%s'\n", error_verbosity, stack);
  return v;
}

typedef struct {
  size_t id, size;
} pmdata;

static void * pmalloc (pmdata * d, size_t size)
{
  assert (d != NULL);
  d->id = 1;
  d->size = size;
  return ((char *)d) + sizeof(pmdata);
}

static void * pfree (void * ptr, Ast * n, Stack * stack)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id != 1) {
    if (d->size == 0)
      return message (NULL, n, "possible double free\n", warning_verbosity, stack);
    return message (NULL, n, "corrupted pointer\n", warning_verbosity, stack);
  }
  d->id = 0;
  d->size = 0;
  return d;
}

static
Value * internal_functions_no_hook (Ast * call, Ast * identifier, Value ** parameters, Stack * stack, Value * value) {
  return value;
}

Value * (* ast_internal_functions_hook) (Ast * call, Ast * identifier, Value ** parameters,
					 Stack * stack, Value * value) =
  internal_functions_no_hook;

typedef struct {
  size_t size;
} FatPointer;

/**
Memory sizes above DYNAMIC_SIZE are allocated/freed using
`malloc/free` rather than using static allocation. */

#define DYNAMIC_SIZE (1024*1024)

static
void * static_malloc (size_t size, Stack * stack)
{
  if (!size)
    return NULL;
  FatPointer * p;
  if (size > DYNAMIC_SIZE)
    p = malloc (sizeof (FatPointer) + size);
  else {
    StackData * d = stack_get_data (stack);
    p = allocate (d->static_alloc, sizeof (FatPointer) + size);
  }
  p->size = size;
  return (void *) (((char *)p) + sizeof (FatPointer));
}

static
void static_free (void * ptr, Stack * stack)
{
  if (ptr) {
    FatPointer * p = (FatPointer *) (((char *)ptr) - sizeof (FatPointer));
    if (p->size > DYNAMIC_SIZE)
      free (p);
  }
}

static
void * static_calloc (size_t nmemb, size_t size, Stack * stack)
{
  if (!size)
    return NULL;
  void * p = static_malloc (size*nmemb, stack);
  memset (p, 0, size*nmemb);
  return p;
}

static
void * static_realloc (void * ptr, size_t size, Stack * stack)
{
  if (ptr && !size) {
    static_free (ptr, stack);
    return NULL;
  }
  void * p = static_malloc (size, stack);
  if (p && ptr) {
    FatPointer * f = (FatPointer *)(((char *)ptr) - sizeof (FatPointer));
    if (f->size < size) size = f->size;
    memcpy (p, ptr, size);
    static_free (ptr, stack);
  }
  return p;
}

static
bool continue_iterations (Ast * n, Stack * stack, int iter, int calls)
{
  if (iter <= 0) {
    message (NULL, n, "reached maximum number of iterations\n", warning_verbosity, stack);
    return false;
  }    
  return true;
#if 0  
  iter = maximum_iterations - iter;
  calls = calls - ((StackData *)stack_get_data (stack))->maxcalls;

  //  if (iter == maximum_iterations)
#if 0  
  fprintf (stderr, "%s:%d: %d #iterations %d\n",
	   ast_left_terminal (n)->file, ast_left_terminal (n)->line,
	   iter, 20000000/(calls/iter));
#endif
#if 1
  if (iter > 2 && 20000000/(calls/iter) < 1000) {
    if (!strncmp (ast_left_terminal (n)->file, "ast/", 4) && iter < 32)
      return true;
    message (NULL, n, "could reach maximum number of iterations\n", error_verbosity, stack);
    return false;
  }
#endif
  return true;
#endif
}

static
Value * internal_functions (Ast * call, Ast * identifier, Ast ** parameters, bool constant_arguments, Stack * stack)
{
  const char * name = ast_terminal (identifier)->start;
  if (!strcmp (name, "malloc")) {
    Value * params[] = { run (parameters[0], stack) };
    if (!params[0])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_void, 1);
    if (value_flags (params[0]) & unset)
      value_set_flags (value, unset);
    else {
      long size = value_data (params[0], long);
      if (size < 0) size = 0;
      value_data (value, void *) = memset (pmalloc (static_calloc (1, sizeof(pmdata) + size, stack), size), 0, size);
      value_data (value, Pointer).start = value_data (value, void *);
      value_data (value, Pointer).size = size;
      value_data (value, Pointer).pscope = 0; // global scope // value->vscope;
    }
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "calloc")) {
    Value * params[] = { run (parameters[0], stack), run (parameters[1], stack) };
    if (!params[0] || !params[1])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_void, 1);
    if ((value_flags (params[0]) & unset) || (value_flags (params[1]) & unset))
      value_set_flags (value, unset);
    else {
      long nmemb = value_data (params[0], long);
      long size = value_data (params[1], long);
      if (size < 0) size = 0;
      value_data (value, void *) = memset (pmalloc (static_malloc (sizeof(pmdata) + nmemb*size, stack), nmemb*size),
					   0, nmemb*size);
      value_data (value, Pointer).start = value_data (value, void *);
      value_data (value, Pointer).size = nmemb*size;    
      value_data (value, Pointer).pscope = 0; // global scope // value->vscope;
    }
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "realloc")) {
    Value * params[] = { run (parameters[0], stack), run (parameters[1], stack) };
    if (!params[0] || !params[1])
      return undefined (NULL, call, stack);
    if (!params[0]->pointer)
      return message (NULL, parameters[0], "'%s' is not a pointer\n", error_verbosity, stack);
    void * ptr = value_data (params[0], void *);
    if (ptr != value_data (params[0], Pointer).start)
      return message (NULL, parameters[0], "'%s' is not malloc'ed\n", error_verbosity, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_void, 1);
    if ((value_flags (params[0]) & unset) || (value_flags (params[1]) & unset))
      value_set_flags (value, unset);
    else {
      long size = value_data (params[1], long);
      if (size < 0) size = 0;
      value_data (value, void *) = pmalloc (static_realloc (pfree (ptr, call, stack), sizeof (pmdata) + size, stack), size);
      long oldsize = value_data (params[0], Pointer).size;
      if (size > oldsize)
	memset (value_data (value, char *) + oldsize, 0, size - oldsize);
      value_data (value, Pointer).start = value_data (value, void *);
      value_data (value, Pointer).size = size;
      value_data (value, Pointer).pscope = 0; // global scope // value->vscope;
    }
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "strdup")) {
    Value * params[] = { run (parameters[0], stack) };
    if (!params[0])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_char, 1);
    if (value_flags (params[0]) & unset)
      value_set_flags (value, unset);
    else {
      char * s = value_data (params[0], char *);
      if (!pointer_check (call, params[0], 1, stack))
	return NULL;
      value_data (value, char *) = pmalloc (static_malloc (sizeof (pmdata) + strlen(s) + 1, stack), strlen(s) + 1);
      strcpy (value_data (value, char *), s);
      value_data (value, Pointer).start = value_data (value, void *);
      value_data (value, Pointer).size = strlen (s) + 1;
      value_data (value, Pointer).pscope = 0; // global scope // value->vscope;
    }
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "free")) {
    Value * params[] = { run (parameters[0], stack) };
    if (!params[0])
      return undefined (NULL, call, stack);    
    if (!params[0]->pointer)
      return message (NULL, parameters[0], "'%s' is not a pointer\n", error_verbosity, stack);
    if (value_flags (params[0]) & unset)
      return NULL;
    void * ptr = value_data (params[0], void *);
    if (ptr != value_data (params[0], Pointer).start)
      return message (NULL, parameters[0], "'%s' is not malloc'ed\n", error_verbosity, stack);
    StackData * d = stack_get_data (stack);

    /**
    Remove non-local values whose addresses overlap with the memory
    being freed. */
    
    if (d->nonlocals) {
      Pointer p = value_data (params[0], Pointer);
      for (khiter_t k = kh_begin (d->nonlocals); k != kh_end (d->nonlocals); k++)
	if (kh_exist (d->nonlocals, k)) {
	  long l = kh_key (d->nonlocals, k) - (long) p.start;
	  if (l >= 0 && l < p.size) {
	    free (kh_value (d->nonlocals, k));
	    kh_del (INT64, d->nonlocals, k);
	  }
	}
    }
    
    static_free (pfree (ptr, call, stack), stack);
    return ast_internal_functions_hook (call, identifier, params, stack, NULL);
  }
  else if (!strcmp (name, "memset")) {
    Value * params[] = { run (parameters[0], stack), run (parameters[1], stack), run (parameters[2], stack) };
    if (!params[0] || !params[1] || !params[2])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_void, 1);
    if ((value_flags (params[0]) & unset) ||
	(value_flags (params[1]) & unset) ||
	(value_flags (params[2]) & unset))
      value_set_flags (value, unset);
    else {
      void * s = value_data (params[0], void *);      
      int c = value_data (params[1], int);
      long n = value_data (params[2], long);
      if (!pointer_check (call, params[0], n, stack))
	return NULL;
      memset (s, c, n);
      value_data (value, void *) = s;
    }
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "memcpy")) {
    Value * params[] = { run (parameters[0], stack), run (parameters[1], stack), run (parameters[2], stack) };
    if (!params[0] || !params[1] || !params[2])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_void, 1);
    if ((value_flags (params[0]) & unset) ||
	(value_flags (params[1]) & unset) ||
	(value_flags (params[2]) & unset))
      value_set_flags (value, unset);
    else {
      void * dest = value_data (params[0], void *);
      const void * src = value_data (params[1], void *);
      long n = value_data (params[2], long);
      if (!pointer_check (call, params[0], n, stack) ||
	  !pointer_check (call, params[1], n, stack))
	return NULL;
      memcpy (dest, src, n);
      value_data (value, void *) = dest;
    }
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "strlen")) {
    Value * params[] = { run (parameters[0], stack) };
    if (!params[0])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_long, 0);
    if (value_flags (params[0]) & unset)
      value_set_flags (value, unset);
    else {
      char * s = value_data (params[0], void *);
      if (!s)
	message (NULL, call, "strlen (null, stack)\n", warning_verbosity, stack);
      if (!pointer_check (call, params[0], 1, stack))
	return NULL;
      value_data (value, long) = s ? strlen (s) : 0;
    }
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "strcpy")) {
    Value * params[] = { run (parameters[0], stack), run (parameters[1], stack) };
    if (!params[0] || !params[1])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_char, 1);
    if ((value_flags (params[0]) & unset) ||
	(value_flags (params[1]) & unset))
      value_set_flags (value, unset);
    else {
      char * dest = value_data (params[0], char *);
      const char * src = value_data (params[1], char *);
      if (!pointer_check (call, params[1], 1, stack))
	return NULL;
      if (!pointer_check (call, params[0], 1, stack))
	return NULL;
      int len = strlen (src) + 1;
      if (!pointer_check (call, params[0], len, stack))
	return NULL;
      value_data (value, char *) = strcpy (dest, src);
      value_data (value, Pointer) = value_data (params[0], Pointer);
    }
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "strcat")) {
    Value * params[] = { run (parameters[0], stack), run (parameters[1], stack) };
    if (!params[0] || !params[1])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_char, 1);
    if ((value_flags (params[0]) & unset) ||
	(value_flags (params[1]) & unset))
      value_set_flags (value, unset);
    else {
      char * dest = value_data (params[0], char *);
      const char * src = value_data (params[1], char *);
      if (!pointer_check (call, params[1], 1, stack))
	return NULL;
      if (!pointer_check (call, params[0], 1, stack))
	return NULL;
      int len = strlen (dest) + strlen (src) + 1;
      if (!pointer_check (call, params[0], len, stack))
	return NULL;
      value_data (value, char *) = strcat (dest, src);
      value_data (value, Pointer) = value_data (params[0], Pointer);
    }
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "strcmp")) {
    Value * params[] = { run (parameters[0], stack), run (parameters[1], stack) };
    if (!params[0] || !params[1])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_int, 0);
    if ((value_flags (params[0]) & unset) ||
	(value_flags (params[1]) & unset))
      value_set_flags (value, unset);
    else {
      const char * s1 = value_data (params[0], void *);
      const char * s2 = value_data (params[1], void *);
      if (!s1 || !s2) {
	message (NULL, call, "strcmp (null, null, stack)\n", 2, stack);
	value_data (value, int) = 0;
      }
      else
	value_data (value, int) = strcmp (s1, s2);
    }
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "strncmp")) {
    Value * params[] = { run (parameters[0], stack), run (parameters[1], stack), run (parameters[2], stack) };
    if (!params[0] || !params[1] || !params[2])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_int, 0);
    if ((value_flags (params[0]) & unset) ||
	(value_flags (params[1]) & unset) ||
	(value_flags (params[2]) & unset))
      value_set_flags (value, unset);
    else {
      const char * s1 = value_data (params[0], void *);
      const char * s2 = value_data (params[1], void *);
      long size = value_data (params[2], long);
      if (!s1 || !s2) {
	message (NULL, call, "strncmp (null, null, stack)\n", 2, stack);
	value_data (value, int) = 0;
      }
      else
	value_data (value, int) = strncmp (s1, s2, size);
    }
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "atan2")) {
    Value * params[] = { run (parameters[0], stack), run (parameters[1], stack) };
    if (!params[0] || !params[1])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_double, 0);
    if (constant_arguments)
      set_constant_expression (value);
    if ((value_flags (params[0]) & unset) ||
	(value_flags (params[1]) & unset))
      value_set_flags (value, unset);
    else
      value_data (value, double) = atan2 (value_data (params[0], double), value_data (params[1], double));
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "pow")) {
    Value * params[] = { run (parameters[0], stack), run (parameters[1], stack) };
    if (!params[0] || !params[1])
      return undefined (NULL, call, stack);
    Value * value = new_value (stack, call, (Ast *)&ast_double, 0);
    if (constant_arguments)
      set_constant_expression (value);
    if ((value_flags (params[0]) & unset) ||
	(value_flags (params[1]) & unset))
      value_set_flags (value, unset);
    else
      value_data (value, double) = pow (value_data (params[0], double), value_data (params[1], double));
    return ast_internal_functions_hook (call, identifier, params, stack, value);
  }
  else if (!strcmp (name, "interpreter_verbosity"))
    ((StackData *)stack_get_data (stack))->verbosity = value_data (run (parameters[0], stack), int);
  else if (!strcmp (name, "interpreter_maximum_iterations"))
    maximum_iterations = value_data (run (parameters[0], stack), int);
  else if (!strcmp (name, "reset_field_value")) {
    Value * params[] = {
      run (parameters[0], stack), run (parameters[1], stack),
      run (parameters[2], stack), run (parameters[3], stack)
    };
    char * field = value_data (params[0], char *);
    memcpy (field, params[2]->data.p, params[2]->size);
    *((Flags *)(field + params[2]->size - sizeof (Flags))) |= unset;
    return ast_internal_functions_hook (call, identifier, params, stack, NULL);
  }
  else {
    struct {
      const char * name;
      double (* func) (double);
    } funcs[] = {
      {"fabs", fabs}, {"sqrt", sqrt}, {"exp", exp}, {"log", log}, {"log10", log10},
      {"sin", sin}, {"cos", cos}, {"tan", tan},
      {"asin", asin}, {"acos", acos}, {"atan", atan},
      {"sinh", sinh}, {"cosh", cosh}, {"tanh", tanh},
      {"asinh", sinh}, {"acosh", cosh}, {"atanh", tanh},
      {NULL}
    }, * i = funcs;
    for (; i->name; i++)
      if (!strcmp (name, i->name)) {
	Value * params[] = { run (parameters[0], stack) };
	if (!params[0])
	  return undefined (NULL, call, stack);
	Value * value = new_value (stack, call, (Ast *)&ast_double, 0);
	if (constant_arguments)
	  set_constant_expression (value);
	if (value_flags (params[0]) & unset)
	  value_set_flags (value, unset);
	else {
	  double x = value_data (params[0], double);
	  value_data (value, double) = i->func (x);
	}
	return ast_internal_functions_hook (call, identifier, params, stack, value);
      }
  }
  return NULL;
}

static Value * default_check (Stack * stack, Ast * n)
{
  Value * v = NULL;
  if (n->child)
    for (Ast ** c = n->child; *c; c++) {
      Value * i = run (*c, stack);
      if (i)
	v = i;
    }
  return v;
}

static
Stack * save_locals (Stack * stack)
{
#if 0
  fprintf (stderr, "========= save_locals ===========\n");
  ast_stack_print (stack, stderr);
  fprintf (stderr, "========= end save_locals ===========\n");
#endif
  Stack * locals = stack_new (sizeof (Ast *));
  Ast ** i;
  while ((i = stack_pop (stack)))
    if (*i) {
      stack_push (locals, i);
      if ((*i)->sym == sym_direct_declarator ||
	  (*i)->sym == sym_external_declaration)
	break;
    }
  return locals;
}

static
void restore_locals (Stack * stack, Stack * locals)
{
  void * prev = stack_set_push (stack, NULL);
  Ast ** i;
  while ((i = stack_pop (locals)))
    stack_push (stack, i);
  stack_set_push (stack, prev);
  stack_destroy (locals);
#if 0
  fprintf (stderr, "========= restore_locals ===========\n");
  ast_stack_print (stack, stderr);
  fprintf (stderr, "========= end restore_locals ===========\n");
#endif
}

static
void print_function_name (Ast * n)
{
  AstTerminal * t = ast_left_terminal (n);
  char * s = ast_str_append (n->child[0], NULL);
  fprintf (stderr, "%s:%d: %s", t->file, t->line, ast_crop_before (s));
  free (s);
}

static
void init_point_variables (Stack * stack)
{
  Ast * init = ast_parent (ast_identifier_declaration (stack, "_init_point_variables"),
			   sym_function_definition);
  assert (init);
  run (ast_child (init, sym_compound_statement), stack);
}
  
Value * (* ast_choose_hook) (Ast *, Stack *, Value *, Value *) = NULL;

Value * (* run) (Ast *, Stack *) = ast_run_node;

Value * ast_run_node (Ast * n, Stack * stack)
{
  if (!n)
    return NULL;
  if (((StackData *)stack_get_data (stack))->maxcalls <= 0)
    return message (NULL, n, "reached maximum number of interpreter calls\n", warning_verbosity, stack);
  ((StackData *)stack_get_data (stack))->maxcalls--;
  Value * value = NULL;
  Ast * scope = ast_push_declarations (n, stack);
  StackData d, * old = NULL;
  if (scope) {
    d = *((StackData *)stack_get_data (stack));
    d.alloc = new_allocator();
    d.scope++;
    old = stack_set_data (stack, &d);
  }
  
  switch (n->sym) {

  case sym_primary_expression:
    switch (n->child[0]->sym) {

    case sym_IDENTIFIER: {
      Ast * identifier = ast_identifier_declaration (stack, ast_terminal (n->child[0])->start);
      if (identifier) {
	if (has_value (identifier)) {
	  value = (Value *) identifier;
	  if (value->vscope == 0 && (value_flags (value) & unset) &&
	      !strcmp (ast_terminal (identifier)->start, "datasize"))
	    value_unset_flags (value, unset);
	}
	else
	  message (NULL, n->child[0], "no value for '%s'\n", warning_verbosity, stack);
      }
      else
	message (NULL, n->child[0], "'%s' undeclared\n", warning_verbosity, stack);
      break;
    }

    case sym_constant: case sym_string:
      value = constant_value (n->child[0]->child[0], stack);
      break;

    default:
      if (n->child[1]) // expression_error
	value = run (n->child[1], stack);
      else // generic_selection	
	value = default_check (stack, n);
      
    }
    break;

  case sym_postfix_expression:
    if (!n->child[1]) // primary_expression or function_call or array_access
      value = run (n->child[0], stack);
    else if (!n->child[2]) { // INC_OP and DEC_OP
      Value * a = run (n->child[0], stack);
      if (!a)
	return undefined (scope, n->child[0], stack);
      value = value_copy (a, stack_alloc (stack), true);
      unary_operation (n, n->child[1], a, stack);
    }
    else if (ast_schema (n, sym_postfix_expression, 1, token_symbol('.')) ||
	     ast_schema (n, sym_postfix_expression, 1, sym_PTR_OP)) {
      value = run (n->child[0], stack);
      if (value) {
	if (value->type->sym == sym_struct_or_union_specifier) {
	  if (value->pointer && ast_schema (n, sym_postfix_expression, 1, token_symbol('.')))
	    return message (scope, n->child[0], "'%s' is a structure pointer\n", error_verbosity, stack);
	  if (!value->pointer && ast_schema (n, sym_postfix_expression, 1, sym_PTR_OP))
	    return message (scope, n->child[0], "'%s' is a structure\n", error_verbosity, stack);
	  if (value->pointer && !value_data (value, void *))
	    return message (scope, n->child[0], "uninitialised structure pointer '%s'\n", warning_verbosity, stack);
	  Ast * member = ast_find (n->child[2], sym_generic_identifier)->child[0];
	  return struct_member_value (n, value, member, -1, stack);
	}
	else
	  return message (scope, n, "'%s' is not accessing a structure\n", error_verbosity, stack);
	return message (scope, n, "structure member not found in '%s'\n", error_verbosity, stack);
      }
#if 0
      else
	ast_print (n, stdout, 0);
#endif
    }
    else {
      Ast * type = ast_find (ast_schema (n, sym_postfix_expression,
					 1, sym_type_name,
					 0, sym_specifier_qualifier_list), sym_types);
      assert (type);
      type = base_type (type->child[0], NULL, stack);
      
      Dimensions d = { .size = 1 };
      Ast * declarator = ast_schema (n, sym_postfix_expression,
				     1, sym_type_name,
				     1, sym_abstract_declarator);
      // fixme: this is not very general, see also direct_declarator_type()
      for (Ast * c = ast_child (declarator, sym_pointer); c; c = ast_child (c, sym_pointer))
	d.pointer++;
      
      declarator = ast_child (declarator, sym_direct_abstract_declarator);
      
      if (declarator) {
	int nd = 1;
	Ast * array = declarator;
	while (array->child[1]->sym == token_symbol ('['))
	  array = array->child[0], nd++;
	if (array->child[0]->sym == token_symbol ('[') &&
	    !get_array_dimensions (array->child[0], sym_direct_abstract_declarator, &d, nd, stack))
	  return NULL;
      }
      
      value = value_from_type (n, type, &d, ast_child (n, sym_postfix_initializer), stack);
    }
    break;

  case sym_function_call:
    {
      Ast * identifier;
      if ((identifier = ast_schema (n, sym_function_call,
				    0, sym_postfix_expression,
				    0, sym_primary_expression,
				    0, sym_IDENTIFIER)) &&
	  !strcmp (ast_terminal (identifier)->start, "_register_fpointer")) {
	value = run (ast_find (n, sym_argument_expression_list_item), stack);
	break;
      }
    }
    value = run (n->child[0], stack);
    if (value) {
      Ast * identifier, * function_definition;
      bool unset_function_pointer = false;
      if (value->type->sym != sym_function_definition || value->pointer != 1) {
	if (((Ast *)value)->sym == sym_IDENTIFIER &&
	    (!strcmp (ast_terminal ((Ast *)value)->start, "val") ||
	     !strcmp (ast_terminal ((Ast *)value)->start, "depth") ||
	     !strcmp (ast_terminal ((Ast *)value)->start, "pid") ||
	     !strcmp (ast_terminal ((Ast *)value)->start, "tid") ||
	     !strcmp (ast_terminal ((Ast *)value)->start, "npe"))) {
	  identifier = (Ast *) value;
	  do {
	    identifier = ast_identifier_declaration_from_to (stack, ast_terminal (identifier)->start,
							     identifier, NULL);
	    assert (identifier);
	    function_definition = identifier;
	    while (function_definition &&
		   function_definition->sym != sym_declaration &&
		   function_definition->sym != sym_function_definition)
	      function_definition = function_definition->parent;
	  } while (!function_definition || function_definition->sym != sym_function_definition);
	}
	else
	  return message (scope, n->child[0], "called object '%s' is not a function\n", error_verbosity, stack);
      }
      else {
#if 0 // ignore unset function pointers
	if (value_flags (value) & unset)
	  return undefined (scope, n, stack);
#else
	if (value_flags (value) & unset)
	  unset_function_pointer = true;
#endif
	identifier = value_data (value, Ast *);
	if (!identifier)
	  return undefined (scope, n, stack);
	function_definition = identifier;
	while (function_definition &&
	       function_definition->sym != sym_declaration &&
	       function_definition->sym != sym_function_definition)
	  function_definition = function_definition->parent;
	if (!function_definition || function_definition->sym != sym_function_definition)
	  function_definition = ast_get_function_definition (stack, identifier, NULL);
      }
      if (!function_definition)
	return message (scope, identifier, "cannot find function definition for '%s'\n", warning_verbosity, stack);
      else {
	int * ncall = (int *)&ast_left_terminal (n)->value;
	if (*ncall >= maximum_recursion)
	  return message (scope, n, "maximum recursion reached for '%s'\n", warning_verbosity, stack);
	else {
	  Ast * parameters = ast_find (ast_schema (function_definition, sym_function_definition,
						   0, sym_function_declaration,
						   1, sym_declarator),
				       sym_direct_declarator,
				       2, sym_parameter_type_list,
				       0, sym_parameter_list);
	  Stack * locals;
	  bool constant_arguments = true;
	  if (!parameters) {
	    if (stack_verbosity (stack) >= function_verbosity) {
	      print_function_name (n);
	      fputs ("()\n", stderr);
	    }
	    locals = save_locals (stack);
	    stack_push (stack, &function_definition->parent);
	  }
	  else {
	    Ast * ellipsis = ast_schema (parameters->parent, sym_parameter_type_list,
					 2, sym_ELLIPSIS);
	    if (ellipsis)
	      return message (scope, n, "ellipsis in '%s' not implemented yet\n", warning_verbosity, stack);
	    int np = 0;
	    Ast * parameters1 = parameters; // fixme: because foreach_item modifies parameters
	    foreach_item (parameters1, 2, parameter)
	      if (ast_schema (parameter, sym_parameter_declaration,
			      1, sym_declarator))
		np++;
	    Ast * arguments = ast_schema (n, sym_function_call,
					  2, sym_argument_expression_list);
	    int narg = 0;
	    Ast * arguments1 = arguments; // fixme: because foreach_item modifies arguments
	    foreach_item (arguments1, 2, argument) narg++;
	    if ((ellipsis && narg < np) || (!ellipsis && narg != np))
	      return message (scope, n, "wrong number of arguments in '%s'\n", error_verbosity, stack);	      
	    if (!np) {
	      if (stack_verbosity (stack) >= function_verbosity) {
		print_function_name (n);
		fputs ("()\n", stderr);
	      }
	      locals = save_locals (stack);
	      stack_push (stack, &function_definition->parent);
	      ast_push_declarations (function_definition, stack);
	      Ast * function_declaration = ast_child (function_definition, sym_function_declaration);
	      run (function_declaration, stack);		
	    }
	    else {
	      Value * v[narg];
	      narg = 0;
	      Ast * arguments1 = arguments;
	      foreach_item_r (arguments, sym_argument_expression_list_item, argument) {
		v[narg] = run (argument, stack);
		if (ast_assign_hook)
		  v[narg] = ast_assign_hook (argument, v[narg], v[narg], stack);
		if (constant_arguments && !is_constant_expression (v[narg]))
		  constant_arguments = false;
		narg++;
	      }
	      if (stack_verbosity (stack) >= function_verbosity) {
		print_function_name (n);
		int i = 0;
		foreach_item_r (arguments1, sym_argument_expression_list_item, argument) {
		  fputs (i == 0 ? " (" : ", ", stderr);
		  value_print (v[i++], stderr);
		}
		fputs (")\n", stderr);
	      }
	      if (!strcmp (ast_terminal (identifier)->start, "display_value")) {
		assert (narg == 1);
		display_value (v[0]);
		return NULL;
	      }
	      locals = save_locals (stack);
	      StackData * d = stack_get_data (stack);
	      d->scope++;
	      stack_push (stack, &function_definition->parent);
	      ast_push_declarations (function_definition, stack);
	      if (ast_is_point_function (ast_schema (function_definition, sym_function_definition,
						     0, sym_function_declaration,
						     1, sym_declarator)) &&
		  !ast_is_stencil_function (function_definition))
		init_point_variables (stack);
	      Ast * function_declaration = ast_child (function_definition, sym_function_declaration);
	      run (function_declaration, stack);
	      foreach_item (parameters, 2, parameter)
		if (ast_schema (parameter, sym_parameter_declaration,
				1, sym_declarator) &&
		    !assign (parameter, run (parameter, stack), v[--narg], stack)) {
		  d->scope--;
		  ast_pop_scope (stack, function_definition->parent);
		  restore_locals (stack, locals);
		  return message (scope, parameter, "could not assign parameter '%s'\n", warning_verbosity, stack);
		}
	      d->scope--;
	      // fixme: deal with ellipsis....
	    }
	  }
	  if (!strcmp (ast_terminal (identifier)->file, "ast/interpreter/internal.h")) {
	    Ast ** parameters = function_parameters (identifier);
	    value = internal_functions (n, identifier, parameters, constant_arguments, stack);
	    free (parameters);
	  }
	  else {
	    (*ncall)++;
	    StackData * d = stack_get_data (stack);
	    int call = d->call;
	    d->call = d->scope;
	    run (ast_child (function_definition, sym_declaration_list), stack);
	    value = run (ast_child (function_definition, sym_compound_statement), stack);
	    if (value)
	      value_unset_flags (value, constant_expression);
	    ast_left_terminal (n)->value = NULL;
	    d->call = call;
	  }
	  if (value && ((Ast *)value)->sym == sym_jump_statement) {
	    value_set_parent (value, n);
	    if (unset_function_pointer)
	      unset_value (value, stack);
	  }
	  ast_pop_scope (stack, function_definition->parent);
	  restore_locals (stack, locals);
	  if (ast_terminal (identifier)->file &&
	      !strcmp (ast_terminal (identifier)->file, "ast/interpreter/overload.h") &&
	      !strcmp (ast_terminal (identifier)->start, "val")) {
	    value = array_member_value (n, value, 0, false, stack);
	    if (value)
	      unset_value (value, stack); // field values are always undefined
	  }
	}
      }
    }
    else {
      if (stack_verbosity (stack) >= warning_verbosity && is_new_message (n, stack)) {
	message (NULL, n, "undefined function call '", warning_verbosity, stack);
	print_function_call (n, stderr, stack);
	fputs ("'\n", stderr);
      }
      else {
	Ast * arguments = ast_schema (n, sym_function_call,
				      2, sym_argument_expression_list);
	if (arguments) {
	  foreach_item_r (arguments, sym_argument_expression_list_item, arg)
	    run (arg, stack);
	}
      }
    }
    break;

  case sym_array_access:
    if (n->child[2]) {
      Value * a = run (n->child[0], stack);
      if (!a) return undefined (scope, n->child[0], stack);
      if (!a->pointer && is_constant_expression (a)) // skips dimensions e.g. 3.1456 [3,-1]
	return a;
      Value * i = run (n->child[2], stack);
      if (!i) return undefined (scope, n->child[2], stack);
      bool error = false;
      int index = value_int (i, &error, stack);
      if (a->pointer && !error) {
#if 0	
	if (value_flags (i) & unset)
	  return message (scope, n->child[2], "undefined array index '%s'\n", warning_verbosity, stack);
#endif
	value = array_member_value (n, a, index, value_flags (i) & unset, stack);
      }
      else
	return message (scope, n, "inconsistent array access '%s'\n", error_verbosity, stack);
    }
    break;
    
  case sym_argument_expression_list_item:
    value = run (n->child[0], stack);
    break;
    
  case sym_unary_expression:
    if (!n->child[1]) // postfix_expression or new_field
      value = run (n->child[0], stack);
    else if (n->child[0]->sym == sym_SIZEOF) {
      if (n->child[1]->sym == token_symbol ('(')) {
	value = new_value (stack, n, (Ast *)&ast_int, 0);
	if (ast_schema (n, sym_unary_expression,
			2, sym_type_name,
			1, sym_abstract_declarator))
	  value_data (value, int) = type_size (1, NULL, stack);
	else {
	  Ast * type = ast_find (ast_schema (n, sym_unary_expression,
					     2, sym_type_name,
					     0, sym_specifier_qualifier_list), sym_types);
	  assert (type);
	  value_data (value, int) = type_size (0, type->child[0], stack);
	}
      }
      else { // SIZEOF unary_expression
	Value * v = run (n->child[1], stack);
	if (!v)
	  return undefined (scope, n->child[1], stack);
	value = new_value (stack, n, (Ast *)&ast_int, 0);
	value_data (value, int) = type_size (v->pointer, v->type, stack);
      }
    }
    else if (!n->child[2])
      value = unary_operation (n, n->child[0], run (n->child[1], stack), stack);
    else // ALIGNOF '(' type_name ')'
      return not_implemented (scope, n, stack);
    break;

  case sym_cast_expression:
    if (!n->child[1]) // unary_expression
      value = run (n->child[0], stack);
    else { // '(' type_name ')' cast_expression
      value = run (n->child[3], stack);
      if (value) {
	Dimensions d = { .size = 1 };
	Ast * type = type_name_type (n->child[1], &d, stack);
	if (value->pointer && !d.pointer)
	  return message (scope, n->child[1], "can only cast from pointers to pointers\n",
			  warning_verbosity, stack);
	Value * v = new_value (stack, n, type, d.pointer);
	assign (n, v, value, stack);
	value = v;
      }
    }
    break;

  case sym_relational_expression:
  case sym_equality_expression:
  case sym_and_expression:
  case sym_exclusive_or_expression:
  case sym_inclusive_or_expression:
  case sym_logical_and_expression:
  case sym_logical_or_expression:
  case sym_shift_expression:
  case sym_multiplicative_expression:
  case sym_additive_expression:
    if (!n->child[1])
      value = run (n->child[0], stack);
    else
      value = binary_operation (n, stack);
    break;
    
  case sym_conditional_expression:
    value = run (n->child[0], stack);
    if (n->child[1]) {
      if (!value)
	undefined (NULL, n->child[0], stack);
      else {
	bool error = false;
	bool cond = value_bool (value, &error, stack);
	if (error)
	  message (NULL, n->child[0], "could not evaluate condition '%s'\n", warning_verbosity, stack);
	else if (value_flags (value) & unset)
	  message (NULL, n->child[0], "undefined condition '%s'\n", warning_verbosity, stack);
	else	  
	  return run (n->child[cond ? 2 : 4], stack);
      }
      
      /**
      If the condition is undefined, we execute both branches. */

      value = run (n->child[2], stack);
      Value * v = run (n->child[4], stack);
      if (ast_choose_hook)
	value = ast_choose_hook (n, stack, value, v);
      if (value && ((Ast *)value)->sym == sym_IDENTIFIER)
	value = value_copy (value, stack_alloc (stack), true);
      value_set_flags (value, unset);
    }
    break;
    
  case sym_assignment_expression:
    if (!n->child[1]) // conditional_expression
      value = run (n->child[0], stack);
    else if (ast_schema (n, sym_assignment_expression,
			 1, sym_assignment_operator,
			 0, token_symbol('=')))
      value = assign (n, run (n->child[0], stack), run (n->child[2], stack), stack);
    else
      value = binary_operation (n, stack);
    if (stack_verbosity (stack) >= expression_verbosity && n->parent->sym == sym_direct_declarator)
      print_expression (n, value);
    break;

  case sym_parameter_declaration:
    value = run (ast_child (n, sym_declarator), stack);
    break;  

  case sym_direct_declarator:
    if (n->child[0]->sym == sym_direct_declarator)
      value = run (n->child[0], stack);
    else if (n->child[0]->sym == sym_generic_identifier) {
      Ast * identifier = ast_identifier_declaration (stack, ast_terminal (n->child[0]->child[0])->start);
      if (identifier && has_value (identifier))
	value = (Value *) identifier;
    }
    else
      value = run (ast_child (n, sym_declarator), stack);
    break;
    
  case sym_declarator:
    value = run (ast_child (n, sym_direct_declarator), stack);
    break;
    
  case sym_init_declarator:
    run (n->child[0], stack);
    break;

  case sym_expression:
    if (n->child[1])
      run (n->child[0], stack);
    value = run (ast_child (n, sym_assignment_expression), stack);
    if (stack_verbosity (stack) >= expression_verbosity && !n->child[1])
      print_expression (n, value);
    break;
    
  case sym_jump_statement:
    if (n->child[0]->sym == sym_RETURN) {
      if (n->child[1]->sym == sym_expression) {
	value = run (n->child[1], stack);
	if (value && ((Ast *)value)->sym == sym_IDENTIFIER)
	  value = value_copy (value, stack_alloc (stack), true);
	value_set_parent (value, n);
	if (ast_assign_hook)
	  value = ast_assign_hook (n, value, value, stack);
	/**
	Cast the value to the proper return type, if necessary. */
	if (value && value->type->sym != sym_struct_or_union_specifier) {
	  Dimensions d = { .size = 1 };
	  Ast * type = return_type (ast_schema (ast_parent (n, sym_function_definition), sym_function_definition,
						0, sym_function_declaration), &d, stack);
	  if (type && type->sym != value->type->sym) {
#if 0	    
	    ast_print_tree (type, stderr, 0, 0, -1);
	    ast_print_tree (value->type, stderr, 0, 0, -1);
#endif
	    Value * v = new_value (stack, n, type, d.pointer);
	    assign (n, v, value, stack);
	    value = v;
	  }
	}
      }
      else
	value = new_value (stack, n, (Ast *)&ast_void, 0);
#if 1
      StackData * d = stack_get_data (stack);
      if (d->call < d->conditional) // fixme: changed from d->call <= d->conditional
	value_set_flags (value, unset);
#endif
    }
    break;

  case sym_for_declaration_statement:
  case sym_iteration_statement: {
    int maxiter = maximum_iterations;
    int maxcalls = ((StackData *)stack_get_data (stack))->maxcalls;
    if (n->child[0]->sym == sym_WHILE) { // while ...
      Value * condition = run (n->child[2], stack);
      if (!condition)
	undefined (NULL, n->child[2], stack);
      else
	do {
	  bool error = false;
	  bool cond = value_bool (condition, &error, stack);
	  if (error) {
	    message (NULL, n->child[2], "could not evaluate condition '%s'\n", warning_verbosity, stack);
	    break;
	  }
	  if (value_flags (condition) & unset) {
	    message (NULL, n->child[2], "undefined condition '%s'\n", warning_verbosity, stack);
	    cond = false;
	  }
	  else if (!cond)
	    break;
	  value = run (n->child[4], stack);
	  if (!cond || (value && ((Ast *)value)->sym == sym_jump_statement))
	    break;
	  condition = run (n->child[2], stack);
	} while (condition && continue_iterations (n, stack, --maxiter, maxcalls));
    }
    else if (n->child[0]->sym == sym_DO) { // do ... while
      do {
	value = run (n->child[1], stack);
	if (value && ((Ast *)value)->sym == sym_jump_statement)
	  break;
	Value * condition = run (n->child[4], stack);
	if (!condition) {
	  undefined (NULL, n->child[4], stack);
	  break;
	}
	bool error = false;
	bool cond = value_bool (condition, &error, stack);
	if (error) {
	  message (NULL, n->child[4], "could not evaluate condition '%s'\n", warning_verbosity, stack);
	  break;
	}
	if (value_flags (condition) & unset) {
	  message (NULL, n->child[4], "undefined condition '%s'\n", warning_verbosity, stack);
	  break;
	}
	else if (!cond)
	  break;
      } while (continue_iterations (n, stack, --maxiter, maxcalls));
    }
    else if (n->child[0]->sym == sym_for_declaration_statement)
      value = run (n->child[0], stack);
    else { // for loop
      assert (n->child[0]->sym == sym_for_scope);
      run (n->child[2], stack);
      Value * condition = run (n->child[3], stack);
      if (!condition)
	undefined (NULL, n->child[3], stack);
      else
	do {
	  bool error = false;
	  bool cond = value_bool (condition, &error, stack);
	  if (error) {
	    message (NULL, n->child[3], "could not evaluate condition '%s'\n", warning_verbosity, stack);
	    break;
	  }
	  if (value_flags (condition) & unset) {
	    message (NULL, n->child[3], "undefined condition '%s'\n", warning_verbosity, stack);
	    cond = false;
	  }
	  else if (!cond)
	    break;
	  value = run (ast_child (n, sym_statement), stack);
	  if (value && ((Ast *)value)->sym == sym_jump_statement)
	    break;
	  Value * expr = run (ast_child (n, sym_expression), stack);
	  condition = run (n->child[3], stack);
	  if (!cond || !expr || !condition)
	    break;
	} while (continue_iterations (n, stack, --maxiter, maxcalls));
    }
    break;
  }

  case sym_selection_statement:
    if (n->child[0]->sym == sym_IF) {
      Value * condition = run (n->child[2], stack);
      if (condition && !(value_flags (condition) & unset)) {
	bool error = false;
	bool cond = value_bool (condition, &error, stack);
	if (!error) {
	  if (cond)
	    value = run (n->child[4], stack);
	  else if (n->child[5])
	    value = run (n->child[6], stack);
	}
      }
      else {
	
	/**
	If the condition is undefined, we execute both branches. */

	if (!condition)
	  undefined (NULL, n->child[2], stack);
	else
	  message (NULL, n->child[2], "undefined condition '%s'\n", warning_verbosity, stack);

	StackData * d = stack_get_data (stack);
	int conditional = d->conditional;
	d->conditional = d->scope;
	khash_t(INT64) * nonlocals = d->nonlocals;
	d->nonlocals = kh_init (INT64);

	value = run (n->child[4], stack); // execute if branch
	if (n->child[5]) { // else

	  /**
	  Non-local values could have been modified by the 'if'
	  branch, we restore their values before executing the 'else'
	  branch. */
	  
	  long l;
	  char * data;
	  kh_foreach (d->nonlocals, l, data, {
	      Value * v = (Value *) data;
	      if ((value_flags (v) & unset) ||
		  memcmp ((char *)l, data + sizeof (Value), v->size)) {
		
		/**
		We swap the current value and the value saved in the
		hash table i.e. we are back to the state before the
		'if' branch and the hash table now contains the value
		modified by the 'if' branch. */
		
		char buf[v->size];
		memcpy (buf, (char *)l, v->size);
		memcpy ((char *)l, data + sizeof (Value), v->size);
		memcpy (data + sizeof (Value), buf, v->size);
		v->data.p = NULL;
	      }
	    });
	  
	  value = run (n->child[6], stack); // execute else branch
	}

	/**
	We unset non-local values which have been modified by the 'if'
	or 'else' branches. */

	long l;
	char * data;
	kh_foreach (d->nonlocals, l, data, {
	    Value * v = (Value *) data;
	    if (!v->data.p || memcmp ((char *)l, data + sizeof (Value), v->size)) {
	      if (!v->data.p)
		v->data.p = (char *) l;
	      if (ast_choose_hook) {

		/**
		The value has been modified by the 'if' or 'else 'branch, we
		need to choose a value. */

		Value ifvalue = *v;
		ifvalue.data.p = data + sizeof (Value);
#if 0
		if (stack_verbosity (stack) > 1) {
		  fprintf (stderr, "if branch value\n");
		  display_value (&ifvalue);
		  fprintf (stderr, "else branch value\n");
		  display_value (v);
		}
#endif		  
		if (ast_choose_hook (n, stack, &ifvalue, v) != v)
		  memcpy (v->data.p, ifvalue.data.p, v->size);
#if 0
		if (stack_verbosity (stack) > 1) {
		  fprintf (stderr, "chosen value\n");
		  display_value (v);
		}
#endif
	      }

	      unset_value (v, stack);

	      /**
	      Characters cannot have flags because their size is
	      unity. This is only to simplify pointer operations and
	      could/should be lifted */
	      
	      if (!v->pointer && v->type->sym == sym_CHAR)
		message (NULL, ((Ast *)v), "characters are always defined/set\n", warning_verbosity, stack);
	    }
	    free (data);
	  });
	
	d->conditional = conditional;
	kh_destroy (INT64, d->nonlocals);
	d->nonlocals = nonlocals;
      }
    }
    else // fixme: SWITCH
      value = default_check (stack, n);
    break;

  case sym_statement:
    value = run (n->child[0], stack);
    if (value && ((Ast *)value)->sym != sym_jump_statement)
      value = NULL;
    break;

  case sym_block_item_list:
    value = run (n->child[0], stack);
    if (n->child[1]) {
      if (!value || ((Ast *)value)->sym != sym_jump_statement)
	value = run (n->child[1], stack);
      else if (value_flags (value) & unset) {
	Value * v = run (n->child[1], stack);
	if (ast_choose_hook && v && ((Ast *)v)->sym == sym_jump_statement) {
	  value = ast_choose_hook (n, stack, value, v);
	  if (!(value_flags (value) & unset)) {
	    if (value && ((Ast *)value)->sym == sym_IDENTIFIER)
	      value = value_copy (value, stack_alloc (stack), true);
	    value_set_flags (value, unset);
	  }
	}
      }
    }
    break;

  case sym_macro_statement: {
    if (ast_is_foreach_statement (n)) {
      if (ast_is_foreach_stencil (n)) {
	StackData * d = stack_get_data (stack);
	run (d->stencil, stack);
	default_check (stack, n);
	run (d->end_stencil, stack);
      }
      else {
	if (strcmp (ast_terminal (n->child[0])->start, "foreach_face_generic"))
	  init_point_variables (stack);
	default_check (stack, n);
      }
      break;
    }

    Ast * identifier = ast_schema (n, sym_macro_statement,
				   0, sym_MACRO);
    if (identifier && !strncmp (ast_terminal (identifier)->start, "is_face_", 8))
      init_point_variables (stack);
    value = default_check (stack, n);
    break;
  }
    
  default:
    value = default_check (stack, n);
  }
  
  if (scope) {
    ast_pop_scope (stack, scope);
    assert (old);
    old->warnings = d.warnings;
    old->maxcalls = d.maxcalls;
    old = stack_set_data (stack, old);
    if (value)
      value = value_copy (value, stack_alloc (stack), true);
    free_allocator (old->alloc);
  }

  return value;
}

void run_external_declarations (Ast * n, Stack * stack)
{
  if (!n)
    return;

  switch (n->sym) {

  case sym_root:
  case sym_external_declaration:
    run_external_declarations (n->child[0], stack);
    break;

  case sym_translation_unit:
    if (!n->child[1])
      run_external_declarations (n->child[0], stack);
    else if (!n->child[2]) {
      run_external_declarations (n->child[0], stack);
      run_external_declarations (n->child[1], stack);
    }
    break;

  case sym_function_definition:
    ast_pop_scope (stack, ast_push_declarations (n, stack));
    break;
    
  case sym_declaration:
    run (n, stack);
    break;
    
  }
}

static
void * push (Stack * stack, void * p)
{
  Ast * n = *((Ast **)p);
  if (n->sym == sym_IDENTIFIER) {
    Value * v = identifier_value (n, stack);
    if (v) {
#if 0
      ast_print (n, stdout, 0);
      display_value (v);
      fprintf (stderr, "%p\n", v);
#endif
      *((Ast **)p) = (Ast *) v;
    }
#if 0
    else
      ast_print_tree (ast_ancestor (n, 2), stderr, 0, 0, -1);
#endif
  }
  return p;
}

static AstRoot * add_external_declarations (AstRoot * root, const char * file)
{
  FILE * fp = fopen (file, "r");
  if (!fp) {
    perror (file);
    assert (fp);
  }
  Stack * stack = root->stack;
  void * prev = stack_set_push (stack, NULL);
  stack_push (stack, &root);
  AstRoot * d = ast_parse_file (fp, root);
  ast_pop_scope (stack, (Ast *) root);
  stack_set_push (stack, prev);
  fclose (fp);
  run_external_declarations ((Ast *) d, stack);
  return d;
}

void (* after_run) (Ast *, Stack *) = NULL;

int ast_run (AstRoot * root, Ast * n, int verbosity, int maxcalls, void * user_data)
{
  Stack * stack = root->stack;
  StackData d = {
    .alloc = root->alloc,
    .static_alloc = root->alloc,
    .verbosity = verbosity,
    .maxcalls = maxcalls,
    .messages = kh_init (INT64),
    .root = (Ast *) root,
    .stencil = ast_new_function_call ((Ast *) root, "_stencil"),
    .end_stencil = ast_new_function_call ((Ast *) root, "_end_stencil"),
    .data = user_data
  };
  void * data = stack_set_data (stack, &d);
  void * prev = stack_set_push (stack, push);
  stack_push (stack, &root);
  if (n) {
    AstRoot * interpreter = add_external_declarations (root, BASILISK "/ast/interpreter/declarations.h");
    AstRoot * internal = add_external_declarations (root, BASILISK "/ast/interpreter/internal.h");
    run_external_declarations ((Ast *) root, stack);
    AstRoot * overload = add_external_declarations (root, BASILISK "/ast/interpreter/overload.h");
    run (n, stack);
    if (after_run)
      after_run (n, stack);
    ast_destroy ((Ast *) interpreter);
    ast_destroy ((Ast *) internal);
    ast_destroy ((Ast *) overload);
  }
  else {
    run ((Ast *) root, stack);
    if (after_run)
      after_run ((Ast *) root, stack);
  }
  ast_pop_scope (stack, (Ast *) root);
  stack_set_push (stack, prev);
  stack_set_data (stack, data);
  kh_destroy (INT64, d.messages);
  ast_destroy (d.stencil);
  ast_destroy (d.end_stencil);
  return d.maxcalls;
}

#include "dimension.c"
