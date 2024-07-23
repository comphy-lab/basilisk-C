/* A Bison parser, made by GNU Bison 3.7.5.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30705

/* Bison version string.  */
#define YYBISON_VERSION "3.7.5"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 2

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* First part of user prologue.  */
#line 1 "parser.y"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
  
typedef struct _Allocator Allocator;  
  
struct _Allocator {
  void * m;
  long len, maxlen;
  Allocator * next;
};

static Allocator * new_allocator()
{
  Allocator * a = calloc (1, sizeof (Allocator));
  a->maxlen = 1 << 16;
  a->m = calloc (a->maxlen, 1);
  return a;
}

static void * allocate (Allocator * a, long size)
{
  Allocator * last = a;
  while (a && a->len + size >= a->maxlen)
    last = a, a = a->next;
  if (a == NULL)
    a = last->next = new_allocator();
  assert (a->len + size < a->maxlen);
  void * p = (void *)(((char *)a->m) + a->len);
  a->len += size;
  return p;
}

static void free_allocator (Allocator * a)
{
  free (a->m);
  if (a->next) free_allocator (a->next);
  free (a);
}
  
#include "parser.h"

static Node * new_value (Allocator * alloc, double val)
{
  Node * n = allocate (alloc, sizeof (Node));
  n->type = '1';
  n->d.value = val;
  return n;
}

static Node * new_expr (Allocator * alloc,
			char type, Node * e1, Node * e2, Node * e3)
{
  Node * n = allocate (alloc, sizeof (Node));
  n->type = type;
  n->e[0] = e1;
  n->e[1] = e2;
  n->e[2] = e3;
  return n;
}

static int yyparse (char ** input, Allocator * alloc, Node ** root);
static int yylex (Node ** lvalp, char ** input, Allocator * alloc);
static void yyerror (char ** input, Allocator * alloc,
		     Node ** root, char const *);

#line 143 "parser.tab.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif


/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    NUM = 258,                     /* NUM  */
    VAR = 259,                     /* VAR  */
    FNCT = 260,                    /* FNCT  */
    OR = 261,                      /* OR  */
    AND = 262,                     /* AND  */
    EQUAL = 263,                   /* EQUAL  */
    DIFF = 264,                    /* DIFF  */
    LE = 265,                      /* LE  */
    GE = 266,                      /* GE  */
    NEG = 267                      /* NEG  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef Node * YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yyparse (char ** input, Allocator * alloc, Node ** root);


/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_NUM = 3,                        /* NUM  */
  YYSYMBOL_VAR = 4,                        /* VAR  */
  YYSYMBOL_FNCT = 5,                       /* FNCT  */
  YYSYMBOL_6_ = 6,                         /* '?'  */
  YYSYMBOL_7_ = 7,                         /* ':'  */
  YYSYMBOL_OR = 8,                         /* OR  */
  YYSYMBOL_AND = 9,                        /* AND  */
  YYSYMBOL_EQUAL = 10,                     /* EQUAL  */
  YYSYMBOL_DIFF = 11,                      /* DIFF  */
  YYSYMBOL_12_ = 12,                       /* '<'  */
  YYSYMBOL_LE = 13,                        /* LE  */
  YYSYMBOL_14_ = 14,                       /* '>'  */
  YYSYMBOL_GE = 15,                        /* GE  */
  YYSYMBOL_16_ = 16,                       /* '-'  */
  YYSYMBOL_17_ = 17,                       /* '+'  */
  YYSYMBOL_18_ = 18,                       /* '*'  */
  YYSYMBOL_19_ = 19,                       /* '/'  */
  YYSYMBOL_NEG = 20,                       /* NEG  */
  YYSYMBOL_21_ = 21,                       /* '!'  */
  YYSYMBOL_22_ = 22,                       /* '^'  */
  YYSYMBOL_23_ = 23,                       /* '['  */
  YYSYMBOL_24_ = 24,                       /* ']'  */
  YYSYMBOL_25_ = 25,                       /* ','  */
  YYSYMBOL_26_ = 26,                       /* '('  */
  YYSYMBOL_27_ = 27,                       /* ')'  */
  YYSYMBOL_YYACCEPT = 28,                  /* $accept  */
  YYSYMBOL_input = 29,                     /* input  */
  YYSYMBOL_exp = 30                        /* exp  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int8 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

#if defined __GNUC__ && ! defined __ICC && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                            \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  15
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   257

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  28
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  3
/* YYNRULES -- Number of rules.  */
#define YYNRULES  27
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  58

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   267


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    21,     2,     2,     2,     2,     2,     2,
      26,    27,    18,    17,    25,    16,     2,    19,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     7,     2,
      12,     2,    14,     6,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    23,     2,    24,    22,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     8,     9,    10,    11,    13,    15,    20
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,    94,    94,    95,    99,   101,   103,   105,   109,   113,
     117,   121,   125,   129,   133,   137,   141,   143,   147,   151,
     155,   159,   163,   167,   171,   175,   179,   183
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "NUM", "VAR", "FNCT",
  "'?'", "':'", "OR", "AND", "EQUAL", "DIFF", "'<'", "LE", "'>'", "GE",
  "'-'", "'+'", "'*'", "'/'", "NEG", "'!'", "'^'", "'['", "']'", "','",
  "'('", "')'", "$accept", "input", "exp", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_int16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,    63,    58,   261,   262,
     263,   264,    60,   265,    62,   266,    45,    43,    42,    47,
     267,    33,    94,    91,    93,    44,    40,    41
};
#endif

#define YYPACT_NINF (-19)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-1)

#define yytable_value_is_error(Yyn) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
      25,   -19,   -19,   -14,   -18,    49,    49,    49,    10,   182,
      31,    49,     3,     3,    68,   -19,    49,    49,    49,    49,
      49,    49,    49,    49,    49,    49,    49,    49,    49,    49,
     -19,   108,    88,   -19,   167,   211,   224,   235,   235,    21,
      21,    21,    21,   -15,   -15,     3,     3,     3,   -19,    49,
     -19,    49,   128,   197,   -19,    49,   148,   -19
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       0,     3,     4,     5,     0,     0,     0,     0,     0,     2,
       0,     0,    15,    16,     0,     1,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       6,     0,     0,    27,     0,    24,    25,    22,    23,    19,
      20,    18,    21,    12,    11,    13,    14,    17,     7,     0,
      10,     0,     0,    26,     8,     0,     0,     9
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -19,   -19,    -5
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
       0,     8,     9
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int8 yytable[] =
{
      12,    13,    14,    27,    28,    31,    32,    29,    11,    10,
      15,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    29,     1,     0,     2,     3,
       4,     0,     0,     0,     2,     3,     4,    25,    26,    27,
      28,     5,     6,    29,    52,     0,    53,     5,     6,     0,
      56,     7,     2,     3,     4,    30,     0,     7,     0,     0,
       0,     0,     0,     0,     0,     5,     6,     0,     0,     0,
       0,     0,     0,     0,    16,     7,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,     0,     0,
      29,     0,     0,     0,    16,    33,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,     0,     0,
      29,     0,     0,     0,    16,    50,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,     0,     0,
      29,     0,    48,    49,    16,     0,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,     0,     0,
      29,     0,    54,    55,    16,     0,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,     0,     0,
      29,     0,    57,    16,    51,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,     0,    16,    29,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,     0,     0,    29,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,     0,     0,    29,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,     0,     0,    29,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,     0,     0,    29,    21,    22,    23,
      24,    25,    26,    27,    28,     0,     0,    29
};

static const yytype_int8 yycheck[] =
{
       5,     6,     7,    18,    19,    10,    11,    22,    26,    23,
       0,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    22,     1,    -1,     3,     4,
       5,    -1,    -1,    -1,     3,     4,     5,    16,    17,    18,
      19,    16,    17,    22,    49,    -1,    51,    16,    17,    -1,
      55,    26,     3,     4,     5,    24,    -1,    26,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    16,    17,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,     6,    26,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    -1,    -1,
      22,    -1,    -1,    -1,     6,    27,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    -1,    -1,
      22,    -1,    -1,    -1,     6,    27,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    -1,    -1,
      22,    -1,    24,    25,     6,    -1,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    -1,    -1,
      22,    -1,    24,    25,     6,    -1,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    -1,    -1,
      22,    -1,    24,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    -1,     6,    22,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    -1,    -1,    22,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    -1,    -1,    22,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    -1,    -1,    22,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    -1,    -1,    22,    12,    13,    14,
      15,    16,    17,    18,    19,    -1,    -1,    22
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     1,     3,     4,     5,    16,    17,    26,    29,    30,
      23,    26,    30,    30,    30,     0,     6,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    22,
      24,    30,    30,    27,    30,    30,    30,    30,    30,    30,
      30,    30,    30,    30,    30,    30,    30,    30,    24,    25,
      27,     7,    30,    30,    24,    25,    30,    24
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_int8 yyr1[] =
{
       0,    28,    29,    29,    30,    30,    30,    30,    30,    30,
      30,    30,    30,    30,    30,    30,    30,    30,    30,    30,
      30,    30,    30,    30,    30,    30,    30,    30
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     1,     1,     1,     1,     3,     4,     6,     8,
       4,     3,     3,     3,     3,     2,     2,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     5,     3
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (input, alloc, root, YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
# ifndef YY_LOCATION_PRINT
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif


# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value, input, alloc, root); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, char ** input, Allocator * alloc, Node ** root)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  YY_USE (input);
  YY_USE (alloc);
  YY_USE (root);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yykind < YYNTOKENS)
    YYPRINT (yyo, yytoknum[yykind], *yyvaluep);
# endif
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, char ** input, Allocator * alloc, Node ** root)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep, input, alloc, root);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule, char ** input, Allocator * alloc, Node ** root)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)], input, alloc, root);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, input, alloc, root); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep, char ** input, Allocator * alloc, Node ** root)
{
  YY_USE (yyvaluep);
  YY_USE (input);
  YY_USE (alloc);
  YY_USE (root);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}






/*----------.
| yyparse.  |
`----------*/

int
yyparse (char ** input, Allocator * alloc, Node ** root)
{
/* Lookahead token kind.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

    /* Number of syntax errors so far.  */
    int yynerrs = 0;

    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    goto yyexhaustedlab;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex (&yylval, input, alloc);
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 2: /* input: exp  */
#line 94 "parser.y"
         { *root = yyvsp[0];   }
#line 1284 "parser.tab.c"
    break;

  case 3: /* input: error  */
#line 95 "parser.y"
         { *root = NULL; YYABORT; }
#line 1290 "parser.tab.c"
    break;

  case 7: /* exp: VAR '[' exp ']'  */
#line 105 "parser.y"
                     {
  yyval = yyvsp[-3]; yyvsp[-3]->e[0] = yyvsp[-1];
}
#line 1298 "parser.tab.c"
    break;

  case 8: /* exp: VAR '[' exp ',' exp ']'  */
#line 109 "parser.y"
                          {
  yyval = yyvsp[-5]; yyvsp[-5]->e[0] = yyvsp[-3]; yyvsp[-5]->e[1] = yyvsp[-1];
}
#line 1306 "parser.tab.c"
    break;

  case 9: /* exp: VAR '[' exp ',' exp ',' exp ']'  */
#line 113 "parser.y"
                                  {
  yyval = yyvsp[-7]; yyvsp[-7]->e[0] = yyvsp[-5]; yyvsp[-7]->e[1] = yyvsp[-3]; yyvsp[-7]->e[2] = yyvsp[-1];
}
#line 1314 "parser.tab.c"
    break;

  case 10: /* exp: FNCT '(' exp ')'  */
#line 117 "parser.y"
                     {
  yyval = yyvsp[-3]; yyvsp[-3]->e[0] = yyvsp[-1];
}
#line 1322 "parser.tab.c"
    break;

  case 11: /* exp: exp '+' exp  */
#line 121 "parser.y"
                     {
  yyval = new_expr (alloc, '+', yyvsp[-2], yyvsp[0], NULL);
}
#line 1330 "parser.tab.c"
    break;

  case 12: /* exp: exp '-' exp  */
#line 125 "parser.y"
                     {
  yyval = new_expr (alloc, '-', yyvsp[-2], yyvsp[0], NULL);
}
#line 1338 "parser.tab.c"
    break;

  case 13: /* exp: exp '*' exp  */
#line 129 "parser.y"
                     {
  yyval = new_expr (alloc, '*', yyvsp[-2], yyvsp[0], NULL);
}
#line 1346 "parser.tab.c"
    break;

  case 14: /* exp: exp '/' exp  */
#line 133 "parser.y"
                     {
  yyval = new_expr (alloc, '/', yyvsp[-2], yyvsp[0], NULL);
}
#line 1354 "parser.tab.c"
    break;

  case 15: /* exp: '-' exp  */
#line 137 "parser.y"
                     {
  yyval = new_expr (alloc, 'm', yyvsp[0], NULL, NULL);
}
#line 1362 "parser.tab.c"
    break;

  case 17: /* exp: exp '^' exp  */
#line 143 "parser.y"
                     {
  yyval = new_expr (alloc, '^', yyvsp[-2], yyvsp[0], NULL);
}
#line 1370 "parser.tab.c"
    break;

  case 18: /* exp: exp '>' exp  */
#line 147 "parser.y"
                     {
  yyval = new_expr (alloc, '>', yyvsp[-2], yyvsp[0], NULL);
}
#line 1378 "parser.tab.c"
    break;

  case 19: /* exp: exp '<' exp  */
#line 151 "parser.y"
                     {
  yyval = new_expr (alloc, '<', yyvsp[-2], yyvsp[0], NULL);
}
#line 1386 "parser.tab.c"
    break;

  case 20: /* exp: exp LE exp  */
#line 155 "parser.y"
                    {
  yyval = new_expr (alloc, 'L', yyvsp[-2], yyvsp[0], NULL);
}
#line 1394 "parser.tab.c"
    break;

  case 21: /* exp: exp GE exp  */
#line 159 "parser.y"
                    {
  yyval = new_expr (alloc, 'G', yyvsp[-2], yyvsp[0], NULL);
}
#line 1402 "parser.tab.c"
    break;

  case 22: /* exp: exp EQUAL exp  */
#line 163 "parser.y"
                       {
  yyval = new_expr (alloc, '=', yyvsp[-2], yyvsp[0], NULL);
}
#line 1410 "parser.tab.c"
    break;

  case 23: /* exp: exp DIFF exp  */
#line 167 "parser.y"
                      {
  yyval = new_expr (alloc, 'i', yyvsp[-2], yyvsp[0], NULL);
}
#line 1418 "parser.tab.c"
    break;

  case 24: /* exp: exp OR exp  */
#line 171 "parser.y"
                    {
  yyval = new_expr (alloc, 'O', yyvsp[-2], yyvsp[0], NULL);
}
#line 1426 "parser.tab.c"
    break;

  case 25: /* exp: exp AND exp  */
#line 175 "parser.y"
                     {
  yyval = new_expr (alloc, 'A', yyvsp[-2], yyvsp[0], NULL);
}
#line 1434 "parser.tab.c"
    break;

  case 26: /* exp: exp '?' exp ':' exp  */
#line 179 "parser.y"
                      {
  yyval = new_expr (alloc, '?', yyvsp[-4], yyvsp[-2], yyvsp[0]);
}
#line 1442 "parser.tab.c"
    break;

  case 27: /* exp: '(' exp ')'  */
#line 183 "parser.y"
                     {
  yyval = yyvsp[-1];
}
#line 1450 "parser.tab.c"
    break;


#line 1454 "parser.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (input, alloc, root, YY_("syntax error"));
    }

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, input, alloc, root);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp, input, alloc, root);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;


#if !defined yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (input, alloc, root, YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturn;
#endif


/*-------------------------------------------------------.
| yyreturn -- parsing is finished, clean up and return.  |
`-------------------------------------------------------*/
yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, input, alloc, root);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp, input, alloc, root);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 187 "parser.y"


/* Called by yyparse on error.  */
static void
yyerror (char ** input, Allocator * alloc, Node ** root, char const *s)
{
  // printf ("%s\n", s);
}

struct init
{
  char const *fname;
  double (*fnct) (double);
};

static struct init arith_fncts[] =
{
  {"abs",   fabs},
  {"fabs",  fabs},
  {"exp",   exp},
  {"log",   log},
  {"ln",    log},
  {"log10", log10},
  {"sqrt",  sqrt},
  {"sin",   sin},
  {"cos",   cos},
  {"tan",   tan},
  {"asin",  asin},
  {"acos",  acos},
  {"atan",  atan},
  {"sinh",  sinh},
  {"cosh",  cosh},
  {"tanh",  tanh},
  {"asinh", asinh},
  {"acosh", acosh},
  {"atanh", atanh},
  {"erf",   erf},
  {"erfc",  erfc},
  {NULL, NULL}
};

static int
yylex (YYSTYPE *lvalp, char ** input, Allocator * alloc)
{
  int c;

  /* Ignore white space, get first nonwhite character.  */
  while ((c = *((*input)++)) == ' ' || c == '\t' || c == '\n' || c == '\t');

  if (c == '\0')
    return EOF;

  /* Char starts a number => parse the number.         */
  if (c == '.' || isdigit (c))
    {
      (*input)--;
      double val;
      int pos;
      if (sscanf (*input, "%lf%n", &val, &pos) == 1) {
	*lvalp = new_value (alloc, val);
	*input += pos;
      }
      else
	return EOF;
      return NUM;
    }

  /* Char starts an identifier => read the name.       */
  if (isalpha (c))
    {
      /* Initially make the buffer long enough
         for a 40-character symbol name.  */
      int length = 40;
      char * symbuf = (char *)malloc (length + 1);
      int i = 0;
      do
        {
          /* If buffer is full, make it bigger.        */
          if (i == length)
            {
              length *= 2;
              symbuf = (char *) realloc (symbuf, length + 1);
            }
          /* Add this character to the buffer.         */
          symbuf[i++] = c;
          /* Get another character.                    */
          c = *((*input)++);
        }
      while (isalnum (c) || c == '_' || c == '.');

      (*input)--;
      symbuf[i] = '\0';

      Node * n = allocate (alloc, sizeof (Node));
      
      struct init * p = arith_fncts;      
      int type = VAR;
      while (p->fname) {
	if (!strcmp (symbuf, p->fname)) {
	  n->type = 'f';
	  n->d.func = p->fnct;
	  type = FNCT;
	  break;
	}
	p++;
      }
      if (type == VAR) {
	n->type = 'v';
	n->d.id = allocate (alloc, strlen(symbuf) + 1);
	memcpy (n->d.id, symbuf, strlen(symbuf) + 1);
      }
      free (symbuf);
      *lvalp = n;
      return type;
    }

  if (c == '|' && **input == '|') {
    (*input)++;
    return OR;
  }

  if (c == '&' && **input == '&') {
    (*input)++;
    return AND;
  }

  if (c == '=' && **input == '=') {
    (*input)++;
    return EQUAL;
  }  

  if (c == '!' && **input == '=') {
    (*input)++;
    return DIFF;
  }

  if (c == '>' && **input == '=') {
    (*input)++;
    return GE;
  }

  if (c == '<' && **input == '=') {
    (*input)++;
    return LE;
  }
  
  /* Any other character is a token by itself.        */
  return c;
}

Node * copy_node (Node * n)
{
  Node * a = malloc (sizeof (Node));
  memcpy (a, n, sizeof (Node));
  if (a->type == 'v')
    a->d.id = strdup (n->d.id);
  int i;
  for (i = 0; i < 3; i++)
    if (a->e[i])
      a->e[i] = copy_node (a->e[i]);
  return a;
}

Node * parse_node (char * code)
{
  Node * root = NULL;
  Allocator * alloc = new_allocator();
  int status = yyparse (&code, alloc, &root);
  if (status == 0)
    root = copy_node (root);
  else
    root = NULL;
  free_allocator (alloc);
  return root;
}

void free_node (Node * n)
{
  int i;
  for (i = 0; i < 3; i++)
    if (n->e[i])
      free_node (n->e[i]);
  if (n->type == 'v')
    free (n->d.id);
  free (n);
}

void print_node (Node * n, FILE * fp)
{
  fprintf (fp, "n%p ", n);
  if (n->type == '1')
    fprintf (fp, "[label=\"%g\"];\n", n->d.value);
  else if (n->type == 'v')
    fprintf (fp, "[label=\"%s\"];\n", n->d.id);
  else if (n->type == 'f')
    fprintf (fp, "[label=\"%p\"];\n", n->d.func);
  else
    fprintf (fp, "[label=\"%c\", shape=box];\n", n->type);
  int i;
  for (i = 0; i < 3; i++)
    if (n->e[i])
      fprintf (fp, "n%p -> n%p;\n", n, n->e[i]);
  for (i = 0; i < 3; i++)
    if (n->e[i])
      print_node (n->e[i], fp);
}

void reset_node_type (Node * n, char type)
{
  n->type = type; free (n->d.id);
}

#if STANDALONE
int
main (int argc, char * argv[])
{
  assert (argc > 1);
  Node * root = parse_node (argv[1]);
  if (root) {
    printf ("digraph G {\n");
    print_node (root, stdout);
    printf ("}\n");
    free_node (root);
  }
  else
    fprintf (stderr, "syntax error\n");
}
#endif // STANDALONE
