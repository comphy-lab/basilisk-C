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
#line 54 "basilisk.yacc"

#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include "parser.h"

static Ast * ast_reduce (Allocator * alloc, int sym, Ast ** children, int n);
#define DEFAULT_ACTION(yyn)					\
  yyval = ast_reduce ((Allocator *)parse->alloc, yyr1[yyn], yyvsp, yyr2[yyn])
static int yyparse (AstRoot * parse, Ast * root);

#line 84 "basilisk.c"

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

#include "basilisk.h"
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_IDENTIFIER = 3,                 /* IDENTIFIER  */
  YYSYMBOL_I_CONSTANT = 4,                 /* I_CONSTANT  */
  YYSYMBOL_F_CONSTANT = 5,                 /* F_CONSTANT  */
  YYSYMBOL_STRING_LITERAL = 6,             /* STRING_LITERAL  */
  YYSYMBOL_FUNC_NAME = 7,                  /* FUNC_NAME  */
  YYSYMBOL_SIZEOF = 8,                     /* SIZEOF  */
  YYSYMBOL_PTR_OP = 9,                     /* PTR_OP  */
  YYSYMBOL_INC_OP = 10,                    /* INC_OP  */
  YYSYMBOL_DEC_OP = 11,                    /* DEC_OP  */
  YYSYMBOL_LEFT_OP = 12,                   /* LEFT_OP  */
  YYSYMBOL_RIGHT_OP = 13,                  /* RIGHT_OP  */
  YYSYMBOL_LE_OP = 14,                     /* LE_OP  */
  YYSYMBOL_GE_OP = 15,                     /* GE_OP  */
  YYSYMBOL_EQ_OP = 16,                     /* EQ_OP  */
  YYSYMBOL_NE_OP = 17,                     /* NE_OP  */
  YYSYMBOL_AND_OP = 18,                    /* AND_OP  */
  YYSYMBOL_OR_OP = 19,                     /* OR_OP  */
  YYSYMBOL_MUL_ASSIGN = 20,                /* MUL_ASSIGN  */
  YYSYMBOL_DIV_ASSIGN = 21,                /* DIV_ASSIGN  */
  YYSYMBOL_MOD_ASSIGN = 22,                /* MOD_ASSIGN  */
  YYSYMBOL_ADD_ASSIGN = 23,                /* ADD_ASSIGN  */
  YYSYMBOL_SUB_ASSIGN = 24,                /* SUB_ASSIGN  */
  YYSYMBOL_LEFT_ASSIGN = 25,               /* LEFT_ASSIGN  */
  YYSYMBOL_RIGHT_ASSIGN = 26,              /* RIGHT_ASSIGN  */
  YYSYMBOL_AND_ASSIGN = 27,                /* AND_ASSIGN  */
  YYSYMBOL_XOR_ASSIGN = 28,                /* XOR_ASSIGN  */
  YYSYMBOL_OR_ASSIGN = 29,                 /* OR_ASSIGN  */
  YYSYMBOL_TYPEDEF_NAME = 30,              /* TYPEDEF_NAME  */
  YYSYMBOL_ENUMERATION_CONSTANT = 31,      /* ENUMERATION_CONSTANT  */
  YYSYMBOL_TYPEDEF = 32,                   /* TYPEDEF  */
  YYSYMBOL_EXTERN = 33,                    /* EXTERN  */
  YYSYMBOL_STATIC = 34,                    /* STATIC  */
  YYSYMBOL_AUTO = 35,                      /* AUTO  */
  YYSYMBOL_REGISTER = 36,                  /* REGISTER  */
  YYSYMBOL_INLINE = 37,                    /* INLINE  */
  YYSYMBOL_CONST = 38,                     /* CONST  */
  YYSYMBOL_RESTRICT = 39,                  /* RESTRICT  */
  YYSYMBOL_VOLATILE = 40,                  /* VOLATILE  */
  YYSYMBOL_BOOL = 41,                      /* BOOL  */
  YYSYMBOL_CHAR = 42,                      /* CHAR  */
  YYSYMBOL_SHORT = 43,                     /* SHORT  */
  YYSYMBOL_INT = 44,                       /* INT  */
  YYSYMBOL_LONG = 45,                      /* LONG  */
  YYSYMBOL_SIGNED = 46,                    /* SIGNED  */
  YYSYMBOL_UNSIGNED = 47,                  /* UNSIGNED  */
  YYSYMBOL_FLOAT = 48,                     /* FLOAT  */
  YYSYMBOL_DOUBLE = 49,                    /* DOUBLE  */
  YYSYMBOL_VOID = 50,                      /* VOID  */
  YYSYMBOL_COMPLEX = 51,                   /* COMPLEX  */
  YYSYMBOL_IMAGINARY = 52,                 /* IMAGINARY  */
  YYSYMBOL_STRUCT = 53,                    /* STRUCT  */
  YYSYMBOL_UNION = 54,                     /* UNION  */
  YYSYMBOL_ENUM = 55,                      /* ENUM  */
  YYSYMBOL_ELLIPSIS = 56,                  /* ELLIPSIS  */
  YYSYMBOL_CASE = 57,                      /* CASE  */
  YYSYMBOL_DEFAULT = 58,                   /* DEFAULT  */
  YYSYMBOL_IF = 59,                        /* IF  */
  YYSYMBOL_ELSE = 60,                      /* ELSE  */
  YYSYMBOL_SWITCH = 61,                    /* SWITCH  */
  YYSYMBOL_WHILE = 62,                     /* WHILE  */
  YYSYMBOL_DO = 63,                        /* DO  */
  YYSYMBOL_FOR = 64,                       /* FOR  */
  YYSYMBOL_GOTO = 65,                      /* GOTO  */
  YYSYMBOL_CONTINUE = 66,                  /* CONTINUE  */
  YYSYMBOL_BREAK = 67,                     /* BREAK  */
  YYSYMBOL_RETURN = 68,                    /* RETURN  */
  YYSYMBOL_ALIGNAS = 69,                   /* ALIGNAS  */
  YYSYMBOL_ALIGNOF = 70,                   /* ALIGNOF  */
  YYSYMBOL_ATOMIC = 71,                    /* ATOMIC  */
  YYSYMBOL_GENERIC = 72,                   /* GENERIC  */
  YYSYMBOL_NORETURN = 73,                  /* NORETURN  */
  YYSYMBOL_STATIC_ASSERT = 74,             /* STATIC_ASSERT  */
  YYSYMBOL_THREAD_LOCAL = 75,              /* THREAD_LOCAL  */
  YYSYMBOL_MAYBECONST = 76,                /* MAYBECONST  */
  YYSYMBOL_NEW_FIELD = 77,                 /* NEW_FIELD  */
  YYSYMBOL_TRACE = 78,                     /* TRACE  */
  YYSYMBOL_FOREACH = 79,                   /* FOREACH  */
  YYSYMBOL_FOREACH_INNER = 80,             /* FOREACH_INNER  */
  YYSYMBOL_FOREACH_DIMENSION = 81,         /* FOREACH_DIMENSION  */
  YYSYMBOL_REDUCTION = 82,                 /* REDUCTION  */
  YYSYMBOL_83_ = 83,                       /* ';'  */
  YYSYMBOL_84_ = 84,                       /* '}'  */
  YYSYMBOL_85_ = 85,                       /* ')'  */
  YYSYMBOL_86_ = 86,                       /* '('  */
  YYSYMBOL_87_ = 87,                       /* ','  */
  YYSYMBOL_88_ = 88,                       /* ':'  */
  YYSYMBOL_89_ = 89,                       /* '.'  */
  YYSYMBOL_90_ = 90,                       /* '{'  */
  YYSYMBOL_91_ = 91,                       /* '['  */
  YYSYMBOL_92_ = 92,                       /* ']'  */
  YYSYMBOL_93_ = 93,                       /* '*'  */
  YYSYMBOL_94_ = 94,                       /* '&'  */
  YYSYMBOL_95_ = 95,                       /* '+'  */
  YYSYMBOL_96_ = 96,                       /* '-'  */
  YYSYMBOL_97_ = 97,                       /* '~'  */
  YYSYMBOL_98_ = 98,                       /* '!'  */
  YYSYMBOL_99_ = 99,                       /* '/'  */
  YYSYMBOL_100_ = 100,                     /* '%'  */
  YYSYMBOL_101_ = 101,                     /* '<'  */
  YYSYMBOL_102_ = 102,                     /* '>'  */
  YYSYMBOL_103_ = 103,                     /* '^'  */
  YYSYMBOL_104_ = 104,                     /* '|'  */
  YYSYMBOL_105_ = 105,                     /* '?'  */
  YYSYMBOL_106_ = 106,                     /* '='  */
  YYSYMBOL_YYACCEPT = 107,                 /* $accept  */
  YYSYMBOL_translation_unit = 108,         /* translation_unit  */
  YYSYMBOL_primary_expression = 109,       /* primary_expression  */
  YYSYMBOL_expression_error = 110,         /* expression_error  */
  YYSYMBOL_constant = 111,                 /* constant  */
  YYSYMBOL_enumeration_constant = 112,     /* enumeration_constant  */
  YYSYMBOL_string = 113,                   /* string  */
  YYSYMBOL_generic_selection = 114,        /* generic_selection  */
  YYSYMBOL_generic_assoc_list = 115,       /* generic_assoc_list  */
  YYSYMBOL_generic_association = 116,      /* generic_association  */
  YYSYMBOL_postfix_expression = 117,       /* postfix_expression  */
  YYSYMBOL_postfix_initializer = 118,      /* postfix_initializer  */
  YYSYMBOL_array_access = 119,             /* array_access  */
  YYSYMBOL_function_call = 120,            /* function_call  */
  YYSYMBOL_member_identifier = 121,        /* member_identifier  */
  YYSYMBOL_argument_expression_list = 122, /* argument_expression_list  */
  YYSYMBOL_argument_expression_list_item = 123, /* argument_expression_list_item  */
  YYSYMBOL_unary_expression = 124,         /* unary_expression  */
  YYSYMBOL_unary_operator = 125,           /* unary_operator  */
  YYSYMBOL_cast_expression = 126,          /* cast_expression  */
  YYSYMBOL_multiplicative_expression = 127, /* multiplicative_expression  */
  YYSYMBOL_additive_expression = 128,      /* additive_expression  */
  YYSYMBOL_shift_expression = 129,         /* shift_expression  */
  YYSYMBOL_relational_expression = 130,    /* relational_expression  */
  YYSYMBOL_equality_expression = 131,      /* equality_expression  */
  YYSYMBOL_and_expression = 132,           /* and_expression  */
  YYSYMBOL_exclusive_or_expression = 133,  /* exclusive_or_expression  */
  YYSYMBOL_inclusive_or_expression = 134,  /* inclusive_or_expression  */
  YYSYMBOL_logical_and_expression = 135,   /* logical_and_expression  */
  YYSYMBOL_logical_or_expression = 136,    /* logical_or_expression  */
  YYSYMBOL_conditional_expression = 137,   /* conditional_expression  */
  YYSYMBOL_assignment_expression = 138,    /* assignment_expression  */
  YYSYMBOL_assignment_operator = 139,      /* assignment_operator  */
  YYSYMBOL_expression = 140,               /* expression  */
  YYSYMBOL_constant_expression = 141,      /* constant_expression  */
  YYSYMBOL_declaration = 142,              /* declaration  */
  YYSYMBOL_declaration_specifiers = 143,   /* declaration_specifiers  */
  YYSYMBOL_init_declarator_list = 144,     /* init_declarator_list  */
  YYSYMBOL_init_declarator = 145,          /* init_declarator  */
  YYSYMBOL_storage_class_specifier = 146,  /* storage_class_specifier  */
  YYSYMBOL_type_specifier = 147,           /* type_specifier  */
  YYSYMBOL_types = 148,                    /* types  */
  YYSYMBOL_struct_or_union_specifier = 149, /* struct_or_union_specifier  */
  YYSYMBOL_struct_or_union = 150,          /* struct_or_union  */
  YYSYMBOL_struct_declaration_list = 151,  /* struct_declaration_list  */
  YYSYMBOL_struct_declaration = 152,       /* struct_declaration  */
  YYSYMBOL_specifier_qualifier_list = 153, /* specifier_qualifier_list  */
  YYSYMBOL_struct_declarator_list = 154,   /* struct_declarator_list  */
  YYSYMBOL_struct_declarator = 155,        /* struct_declarator  */
  YYSYMBOL_enum_specifier = 156,           /* enum_specifier  */
  YYSYMBOL_enumerator_list = 157,          /* enumerator_list  */
  YYSYMBOL_enumerator = 158,               /* enumerator  */
  YYSYMBOL_atomic_type_specifier = 159,    /* atomic_type_specifier  */
  YYSYMBOL_type_qualifier = 160,           /* type_qualifier  */
  YYSYMBOL_function_specifier = 161,       /* function_specifier  */
  YYSYMBOL_alignment_specifier = 162,      /* alignment_specifier  */
  YYSYMBOL_declarator = 163,               /* declarator  */
  YYSYMBOL_direct_declarator = 164,        /* direct_declarator  */
  YYSYMBOL_generic_identifier = 165,       /* generic_identifier  */
  YYSYMBOL_pointer = 166,                  /* pointer  */
  YYSYMBOL_type_qualifier_list = 167,      /* type_qualifier_list  */
  YYSYMBOL_parameter_type_list = 168,      /* parameter_type_list  */
  YYSYMBOL_parameter_list = 169,           /* parameter_list  */
  YYSYMBOL_parameter_declaration = 170,    /* parameter_declaration  */
  YYSYMBOL_identifier_list = 171,          /* identifier_list  */
  YYSYMBOL_type_name = 172,                /* type_name  */
  YYSYMBOL_abstract_declarator = 173,      /* abstract_declarator  */
  YYSYMBOL_direct_abstract_declarator = 174, /* direct_abstract_declarator  */
  YYSYMBOL_type_not_specified = 175,       /* type_not_specified  */
  YYSYMBOL_initializer = 176,              /* initializer  */
  YYSYMBOL_initializer_list = 177,         /* initializer_list  */
  YYSYMBOL_designation = 178,              /* designation  */
  YYSYMBOL_designator_list = 179,          /* designator_list  */
  YYSYMBOL_designator = 180,               /* designator  */
  YYSYMBOL_static_assert_declaration = 181, /* static_assert_declaration  */
  YYSYMBOL_statement = 182,                /* statement  */
  YYSYMBOL_labeled_statement = 183,        /* labeled_statement  */
  YYSYMBOL_compound_statement = 184,       /* compound_statement  */
  YYSYMBOL_185_1 = 185,                    /* @1  */
  YYSYMBOL_block_item_list = 186,          /* block_item_list  */
  YYSYMBOL_block_item = 187,               /* block_item  */
  YYSYMBOL_expression_statement = 188,     /* expression_statement  */
  YYSYMBOL_selection_statement = 189,      /* selection_statement  */
  YYSYMBOL_for_scope = 190,                /* for_scope  */
  YYSYMBOL_iteration_statement = 191,      /* iteration_statement  */
  YYSYMBOL_for_declaration_statement = 192, /* for_declaration_statement  */
  YYSYMBOL_jump_statement = 193,           /* jump_statement  */
  YYSYMBOL_external_declaration = 194,     /* external_declaration  */
  YYSYMBOL_function_declaration = 195,     /* function_declaration  */
  YYSYMBOL_function_definition = 196,      /* function_definition  */
  YYSYMBOL_declaration_list = 197,         /* declaration_list  */
  YYSYMBOL_basilisk_statements = 198,      /* basilisk_statements  */
  YYSYMBOL_macro_statement = 199,          /* macro_statement  */
  YYSYMBOL_foreach_statement = 200,        /* foreach_statement  */
  YYSYMBOL_foreach_parameters = 201,       /* foreach_parameters  */
  YYSYMBOL_foreach_parameter = 202,        /* foreach_parameter  */
  YYSYMBOL_reduction_list = 203,           /* reduction_list  */
  YYSYMBOL_reduction = 204,                /* reduction  */
  YYSYMBOL_reduction_operator = 205,       /* reduction_operator  */
  YYSYMBOL_reduction_array = 206,          /* reduction_array  */
  YYSYMBOL_foreach_inner_statement = 207,  /* foreach_inner_statement  */
  YYSYMBOL_foreach_dimension_statement = 208, /* foreach_dimension_statement  */
  YYSYMBOL_forin_declaration_statement = 209, /* forin_declaration_statement  */
  YYSYMBOL_forin_statement = 210,          /* forin_statement  */
  YYSYMBOL_forin_arguments = 211,          /* forin_arguments  */
  YYSYMBOL_event_definition = 212,         /* event_definition  */
  YYSYMBOL_event_parameters = 213,         /* event_parameters  */
  YYSYMBOL_event_parameter = 214,          /* event_parameter  */
  YYSYMBOL_boundary_definition = 215,      /* boundary_definition  */
  YYSYMBOL_external_foreach_dimension = 216, /* external_foreach_dimension  */
  YYSYMBOL_attribute = 217,                /* attribute  */
  YYSYMBOL_new_field = 218,                /* new_field  */
  YYSYMBOL_root = 219                      /* root  */
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
typedef yytype_int16 yy_state_t;

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
#define YYFINAL  197
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   3450

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  107
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  113
/* YYNRULES -- Number of rules.  */
#define YYNRULES  358
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  643

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   337


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
       2,     2,     2,    98,     2,     2,     2,   100,    94,     2,
      86,    85,    93,    95,    87,    96,    89,    99,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    88,    83,
     101,   106,   102,   105,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    91,     2,    92,   103,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    90,   104,    84,    97,     2,     2,     2,
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
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   105,   105,   106,   107,   108,   109,   113,   114,   115,
     116,   117,   121,   122,   126,   127,   128,   132,   136,   137,
     141,   145,   146,   150,   151,   155,   156,   157,   158,   159,
     160,   161,   162,   166,   167,   171,   172,   173,   177,   178,
     182,   186,   187,   188,   192,   193,   197,   198,   199,   200,
     201,   202,   203,   204,   208,   209,   210,   211,   212,   213,
     217,   218,   224,   225,   226,   227,   231,   232,   233,   237,
     238,   239,   243,   244,   245,   246,   247,   251,   252,   253,
     257,   258,   262,   263,   267,   268,   272,   273,   277,   278,
     282,   283,   287,   288,   289,   290,   291,   295,   296,   297,
     298,   299,   300,   301,   302,   303,   304,   305,   309,   310,
     314,   318,   319,   320,   324,   325,   326,   327,   328,   329,
     330,   331,   332,   333,   337,   338,   342,   343,   347,   348,
     349,   350,   351,   352,   353,   357,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   380,   381,   382,   383,   384,   388,   389,   393,
     394,   398,   399,   400,   404,   405,   406,   407,   411,   412,
     416,   417,   418,   422,   423,   424,   425,   426,   430,   431,
     435,   436,   440,   444,   445,   446,   447,   448,   452,   453,
     457,   458,   462,   463,   467,   468,   469,   470,   471,   472,
     473,   474,   475,   476,   477,   478,   479,   480,   481,   485,
     486,   490,   491,   492,   493,   497,   498,   502,   503,   507,
     508,   512,   513,   514,   515,   519,   520,   524,   525,   529,
     530,   531,   535,   536,   537,   538,   539,   540,   541,   542,
     543,   544,   545,   546,   547,   548,   549,   550,   551,   552,
     553,   554,   555,   559,   563,   564,   565,   569,   570,   571,
     572,   576,   580,   581,   585,   586,   590,   594,   595,   596,
     597,   598,   599,   600,   601,   605,   606,   607,   611,   613,
     613,   619,   620,   624,   625,   629,   630,   634,   635,   636,
     640,   644,   645,   646,   648,   650,   654,   656,   661,   662,
     663,   664,   665,   669,   670,   671,   672,   673,   674,   675,
     676,   680,   684,   686,   691,   692,   699,   700,   701,   702,
     703,   704,   708,   712,   713,   717,   718,   722,   723,   727,
     728,   732,   736,   737,   738,   742,   743,   747,   748,   752,
     753,   757,   762,   767,   768,   772,   776,   777,   778,   782,
     783,   784,   788,   792,   793,   797,   801,   802,   806
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
  "\"end of file\"", "error", "\"invalid token\"", "IDENTIFIER",
  "I_CONSTANT", "F_CONSTANT", "STRING_LITERAL", "FUNC_NAME", "SIZEOF",
  "PTR_OP", "INC_OP", "DEC_OP", "LEFT_OP", "RIGHT_OP", "LE_OP", "GE_OP",
  "EQ_OP", "NE_OP", "AND_OP", "OR_OP", "MUL_ASSIGN", "DIV_ASSIGN",
  "MOD_ASSIGN", "ADD_ASSIGN", "SUB_ASSIGN", "LEFT_ASSIGN", "RIGHT_ASSIGN",
  "AND_ASSIGN", "XOR_ASSIGN", "OR_ASSIGN", "TYPEDEF_NAME",
  "ENUMERATION_CONSTANT", "TYPEDEF", "EXTERN", "STATIC", "AUTO",
  "REGISTER", "INLINE", "CONST", "RESTRICT", "VOLATILE", "BOOL", "CHAR",
  "SHORT", "INT", "LONG", "SIGNED", "UNSIGNED", "FLOAT", "DOUBLE", "VOID",
  "COMPLEX", "IMAGINARY", "STRUCT", "UNION", "ENUM", "ELLIPSIS", "CASE",
  "DEFAULT", "IF", "ELSE", "SWITCH", "WHILE", "DO", "FOR", "GOTO",
  "CONTINUE", "BREAK", "RETURN", "ALIGNAS", "ALIGNOF", "ATOMIC", "GENERIC",
  "NORETURN", "STATIC_ASSERT", "THREAD_LOCAL", "MAYBECONST", "NEW_FIELD",
  "TRACE", "FOREACH", "FOREACH_INNER", "FOREACH_DIMENSION", "REDUCTION",
  "';'", "'}'", "')'", "'('", "','", "':'", "'.'", "'{'", "'['", "']'",
  "'*'", "'&'", "'+'", "'-'", "'~'", "'!'", "'/'", "'%'", "'<'", "'>'",
  "'^'", "'|'", "'?'", "'='", "$accept", "translation_unit",
  "primary_expression", "expression_error", "constant",
  "enumeration_constant", "string", "generic_selection",
  "generic_assoc_list", "generic_association", "postfix_expression",
  "postfix_initializer", "array_access", "function_call",
  "member_identifier", "argument_expression_list",
  "argument_expression_list_item", "unary_expression", "unary_operator",
  "cast_expression", "multiplicative_expression", "additive_expression",
  "shift_expression", "relational_expression", "equality_expression",
  "and_expression", "exclusive_or_expression", "inclusive_or_expression",
  "logical_and_expression", "logical_or_expression",
  "conditional_expression", "assignment_expression", "assignment_operator",
  "expression", "constant_expression", "declaration",
  "declaration_specifiers", "init_declarator_list", "init_declarator",
  "storage_class_specifier", "type_specifier", "types",
  "struct_or_union_specifier", "struct_or_union",
  "struct_declaration_list", "struct_declaration",
  "specifier_qualifier_list", "struct_declarator_list",
  "struct_declarator", "enum_specifier", "enumerator_list", "enumerator",
  "atomic_type_specifier", "type_qualifier", "function_specifier",
  "alignment_specifier", "declarator", "direct_declarator",
  "generic_identifier", "pointer", "type_qualifier_list",
  "parameter_type_list", "parameter_list", "parameter_declaration",
  "identifier_list", "type_name", "abstract_declarator",
  "direct_abstract_declarator", "type_not_specified", "initializer",
  "initializer_list", "designation", "designator_list", "designator",
  "static_assert_declaration", "statement", "labeled_statement",
  "compound_statement", "@1", "block_item_list", "block_item",
  "expression_statement", "selection_statement", "for_scope",
  "iteration_statement", "for_declaration_statement", "jump_statement",
  "external_declaration", "function_declaration", "function_definition",
  "declaration_list", "basilisk_statements", "macro_statement",
  "foreach_statement", "foreach_parameters", "foreach_parameter",
  "reduction_list", "reduction", "reduction_operator", "reduction_array",
  "foreach_inner_statement", "foreach_dimension_statement",
  "forin_declaration_statement", "forin_statement", "forin_arguments",
  "event_definition", "event_parameters", "event_parameter",
  "boundary_definition", "external_foreach_dimension", "attribute",
  "new_field", "root", YY_NULLPTR
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
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,    59,   125,    41,    40,    44,    58,    46,
     123,    91,    93,    42,    38,    43,    45,   126,    33,    47,
      37,    60,    62,    94,   124,    63,    61
};
#endif

#define YYPACT_NINF (-404)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-359)

#define yytable_value_is_error(Yyn) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    1067,   -33,    45,  -404,  -404,  -404,  -404,  2600,  2615,  2615,
     963,  -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,
    -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,
    -404,  -404,  -404,  -404,  -404,    68,    -8,    59,   102,   156,
    -404,   167,  -404,  -404,   112,  -404,   180,  1165,  -404,  -404,
    -404,  -404,  -404,  -404,   773,  -404,  -404,  -404,  -404,   314,
    -404,   -33,   669,  2711,  -404,   -50,    80,   235,    37,   201,
     194,   158,   202,   276,    36,  -404,   216,  -404,    61,  3254,
    3254,  -404,  -404,   147,  -404,  -404,  3254,  3254,  3254,   205,
    -404,  -404,  2952,  -404,  -404,  -404,  -404,  -404,  -404,  -404,
     312,   233,  -404,  -404,  1165,  -404,  -404,  1165,  -404,  -404,
    -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,
    -404,  2143,  -404,  -404,   335,   271,  1414,  3374,  3374,  2558,
    2711,   304,    14,  -404,   669,   264,  -404,   255,  3374,    78,
    3374,   286,   186,  -404,    26,  -404,  -404,  2018,    26,  2158,
    -404,  2143,  -404,  -404,  2711,  2711,  2711,  2711,  2711,  2711,
    2711,  2711,  2711,  2711,  2711,  2711,  2711,  2711,  2711,  2711,
    2711,  2711,  2558,  -404,  -404,    43,   342,    91,  -404,    -6,
     -58,  -404,   219,  -404,  -404,  -404,   687,   285,  -404,  -404,
    -404,  3296,   292,  -404,    61,  -404,  2952,  -404,  -404,   969,
     298,   307,   669,  2003,  -404,  -404,  -404,   283,   -47,  -404,
     335,  -404,   324,   326,   329,   341,   333,   340,   305,   343,
    3254,  -404,  2558,  -404,  2890,  1491,   -23,  -404,    -7,  -404,
    -404,  2571,  -404,  -404,  -404,  -404,  -404,  -404,  -404,   109,
    -404,  -404,  -404,  -404,   337,     9,  -404,  -404,  -404,  -404,
    -404,   -50,   -50,    80,    80,   235,   235,   235,   235,    37,
      37,   201,   194,   158,   202,   276,   195,  -404,   345,  -404,
    -404,  -404,   342,  -404,    43,  2187,  -404,  1535,  -404,   -58,
     347,  3060,  -404,   198,  -404,  2845,  3107,  2711,   327,  -404,
    -404,   349,  2711,   351,   359,   360,   361,  1222,  -404,    26,
     367,   371,  2264,   369,   374,   376,  -404,   123,  -404,   355,
    -404,  -404,  -404,   871,  -404,  -404,  -404,   400,  -404,  -404,
    -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,   397,   397,
      26,  2003,  2711,  -404,  -404,   132,  2187,   -53,  -404,  2711,
    -404,    71,   244,  -404,  -404,  -404,  -404,  3335,   451,  -404,
    3254,    43,  -404,  -404,  -404,   199,   404,   407,  -404,   406,
    1854,  -404,   405,   408,  1612,    -7,  -404,  -404,  1656,  -404,
    -404,  -404,  2143,  -404,  -404,  2711,  -404,  -404,  -404,  -404,
    -404,  -404,  2789,  1854,  -404,   410,   411,  1733,  -404,  -404,
    -404,  -404,  -404,  2711,   175,  -404,   418,   412,  3154,  -404,
     669,  -404,   -18,  -404,  -404,   420,  1222,   581,   581,  2558,
     669,   436,   416,  -404,  -404,  -404,   213,  2047,  2293,    81,
    -404,  1222,  -404,  -404,  1318,  -404,   246,   419,  -404,  1892,
    -404,  -404,  -404,  -404,  -404,  -404,  -404,   148,   424,   190,
    -404,   425,   429,  -404,  -404,  2725,   409,   181,  -404,  -404,
    3207,  -404,  -404,   427,  1854,  -404,  -404,  2558,  -404,   430,
    3011,  1854,  -404,   431,   432,  1777,  -404,  -404,   440,  -404,
    -404,   441,   252,   435,  1854,  -404,  -404,  2558,  -404,   443,
     444,  -404,  -404,  -404,   150,  2711,  -404,  -404,  2571,  2711,
    1222,  2711,  1222,  -404,   453,   454,   258,   442,  -404,  -404,
     458,  1222,  -404,   261,  -404,   463,  -404,  1222,   267,   455,
    1222,  -404,   157,  2322,    61,  2322,  -404,  1907,  -404,  -404,
    -404,  2187,  -404,  2558,  -404,  3335,  2558,   434,  -404,  2187,
    -404,  -404,  -404,   456,   457,  -404,  -404,   461,   459,  1854,
    -404,  -404,  2558,  -404,   460,  -404,  -404,  -404,   544,  -404,
     462,   464,  -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,
    -404,  -404,  -404,  1222,  1222,  1222,  2558,    42,  -404,  1222,
    2418,  -404,  -404,  1222,  1222,  -404,  2143,  2433,    10,  2462,
    -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,  -404,
    -404,  -404,   465,   467,  -404,  -404,  -404,  -404,   495,  -404,
    -404,   280,  -404,  -404,  -404,   472,  -404,  -404,  -404,  -404,
    -404,   255,   476,  1222,   308,  2143,  1222,   323,  -404,  -404,
    1222,   479,    26,  1222,  -404,  1222,   481,  -404,  1222,  -404,
    -404,   477,   482,  -404,  -404,  1222,  -404,   484,  -404,  -404,
    2558,   117,  -404
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_int16 yydefact[] =
{
       0,     0,     7,    14,    15,    18,    19,     0,     0,     0,
     151,    16,   128,   129,   130,   132,   133,   188,   183,   184,
     185,   145,   137,   138,   139,   140,   143,   144,   141,   142,
     136,   146,   147,   157,   158,     0,     0,     0,   186,     0,
     189,     0,   131,   187,   356,   134,     0,     0,    55,    54,
      56,    57,    58,    59,     0,    25,     8,     9,    11,    46,
      27,    26,    60,     0,    62,    66,    69,    72,    77,    80,
      82,    84,    86,    88,    90,    92,     0,   304,     0,   115,
     117,   135,   149,     0,   150,   148,   119,   121,   123,     0,
     113,     2,     0,   303,   305,   306,   307,   308,   309,    53,
       0,   279,   310,     7,     0,    26,    50,     0,    47,    48,
      98,    99,   100,   101,   102,   103,   104,   105,   106,   107,
      97,     0,   209,   210,     0,   177,     0,     0,     0,     0,
       0,     0,     0,    13,   151,     0,   108,    12,   165,   253,
     167,     0,     0,     3,     0,    30,    31,     0,     0,     0,
     322,     0,    60,    49,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   352,   253,     0,   214,     0,   124,   311,
     253,   194,     0,   151,   114,   116,     0,   154,   118,   120,
     122,     0,     0,   314,     0,   313,     0,     1,   278,     0,
       0,     0,     0,     0,    96,    95,    17,   181,     0,   178,
       0,   110,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    10,     0,   164,     0,     0,   230,   253,   231,   228,
     166,     0,     4,     5,     6,    29,    40,    38,    45,     0,
      41,    44,    28,    35,    55,     0,    94,    93,    63,    64,
      65,    67,    68,    70,    71,    75,    76,    73,    74,    78,
      79,    81,    83,    85,    87,    89,     0,   111,     0,   186,
     215,   213,   212,   253,     0,     0,   253,     0,   193,   253,
       0,     0,   159,     0,   163,     0,     0,     0,   127,   315,
     312,     0,     0,     0,     0,     0,     0,     0,   290,     0,
       0,     0,     0,     0,     0,     0,   285,     0,   283,     0,
     284,   267,   268,     0,   281,   269,   270,     0,   271,   295,
     272,   273,   316,   317,   318,   319,   320,   321,    51,     0,
       0,     0,     0,   256,   258,     0,     0,     0,   262,     0,
     173,     0,     0,   191,   190,    52,   182,     0,     0,   357,
       0,     0,   353,   109,   249,   224,     0,   218,   253,     0,
       0,   233,    55,     0,     0,   229,   227,   253,     0,    32,
      61,    39,    43,    36,    37,     0,   195,   216,   211,   112,
     125,   126,     0,     0,   196,    55,     0,     0,   192,   155,
     152,   160,   253,     0,     0,   168,   172,     0,     0,   355,
      60,   349,     0,   346,   274,     0,     0,     0,     0,     0,
     210,     0,     0,   299,   300,   301,     0,     0,     0,     0,
     286,     0,   280,   282,     0,   265,     0,     0,    33,     0,
     257,   261,   263,   180,   174,   179,   175,     0,     0,     0,
      21,     0,     0,   354,   311,     0,   221,   230,   223,   250,
       0,   219,   232,     0,     0,   234,   240,     0,   239,     0,
       0,     0,   241,    55,     0,     0,    42,    91,     0,   225,
     207,     0,     0,     0,     0,   197,   204,     0,   203,    55,
       0,   161,   170,   253,     0,     0,   156,   153,     0,     0,
       0,     0,     0,   277,     0,     0,     0,     0,   298,   302,
       0,     0,   327,     0,   325,   328,   329,     0,     0,     0,
       0,   275,     0,     0,     0,     0,   254,     0,   264,    34,
     260,     0,   176,     0,    20,     0,     0,     0,   151,     0,
     217,   253,   236,     0,     0,   238,   251,     0,     0,     0,
     242,   248,     0,   247,     0,   206,   205,   208,     0,   199,
       0,     0,   200,   202,   162,   169,   171,   351,   350,   348,
     345,   347,   276,     0,     0,     0,     0,     0,   323,     0,
       0,   330,   337,     0,     0,   339,     0,     0,   127,     0,
     255,   259,    24,    22,    23,   266,   222,   220,   235,   237,
     252,   244,     0,     0,   245,   226,   198,   201,   288,   289,
     291,     0,   334,   333,   332,     0,   324,   326,   338,   340,
     344,   343,     0,     0,     0,     0,     0,     0,   243,   246,
       0,     0,     0,     0,   296,     0,     0,   293,     0,   287,
     292,   335,     0,   342,   297,     0,   294,     0,   331,   341,
       0,     0,   336
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -404,  -404,  -404,   -86,  -404,  -404,  -404,  -404,  -404,    50,
     445,  -120,  -404,     0,   433,  -404,   206,   311,  -404,   -31,
     200,   257,   224,   259,   423,   415,   437,   426,   438,  -404,
    -119,    76,   -56,   -45,  -125,   -68,    34,  -404,   303,  -404,
     -21,  -404,  -404,  -404,  -176,  -264,   -96,  -404,    95,  -404,
     390,  -322,  -404,    55,  -404,  -404,   -62,  -174,     4,  -116,
    -211,  -357,  -404,   151,  -404,   357,  -129,  -223,  -133,  -266,
     273,  -399,  -404,   270,   -94,   137,  -404,    19,  -404,  -404,
     295,  -403,  -404,  -404,  -404,  -404,  -404,   556,  -404,  -208,
    -404,  -404,    35,  -404,  -404,    44,  -404,   114,  -404,  -404,
    -404,  -404,  -404,  -404,     2,  -404,  -404,   -90,  -404,  -404,
    -404,  -404,  -404
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
       0,    54,    55,   135,    56,   207,    57,    58,   439,   440,
      59,   369,    60,   105,   235,   239,   240,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,   136,   121,   307,   212,    77,   355,   177,   178,    79,
      80,    81,    82,    83,   281,   282,   139,   394,   395,    84,
     208,   209,    85,    86,    87,    88,   268,   180,   309,   182,
     272,   356,   357,   358,   472,   441,   359,   228,   229,   334,
     335,   336,   337,   338,    90,   310,   311,   312,   199,   313,
     314,   315,   316,   317,   318,   319,   320,    91,    92,    93,
     196,   321,   322,   323,   503,   504,   505,   506,   605,   632,
     324,   325,   326,   327,   612,    95,   402,   403,    96,    97,
      98,    99,   100
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      61,   204,   137,   365,    89,   217,   151,   211,   279,   381,
     227,   211,   352,   615,   364,   286,   179,   391,   219,   435,
     102,   515,   391,   226,   193,   471,   138,   238,   276,   122,
     521,   246,   153,   277,    78,    94,   330,   340,   332,   125,
     341,   267,   223,   154,   230,   122,   122,   278,  -209,   155,
     156,   161,   162,   431,    61,   171,   123,   101,    89,   137,
     271,   602,   137,   224,   122,   489,   387,   490,   225,   491,
     430,   122,   123,   123,   206,  -209,    76,  -127,   126,   367,
     150,  -127,   181,   138,   368,   509,   138,   187,    78,    94,
     283,   123,   284,   192,   366,   283,   222,   284,   123,   220,
     275,   374,   140,   537,   245,   138,   138,   138,   226,   398,
     577,   195,   579,   184,   185,   435,   275,   138,   521,   138,
     188,   189,   190,   248,   249,   250,   194,   266,   289,   175,
      76,   308,   288,  -209,   391,  -209,   176,   603,   163,   164,
     379,   172,   443,   382,   174,   127,   388,   175,   236,   454,
     122,   206,   236,   122,   176,   434,   378,   465,   124,   140,
     576,   102,   140,   520,   224,   138,   510,   405,   401,   225,
     138,   176,   474,   211,   273,   157,   158,   123,   274,   181,
     123,   140,   140,   140,   122,   283,   181,   284,   128,   283,
     283,   284,   284,   140,   371,   140,   372,   205,   181,    61,
     370,   122,   122,   131,   222,   216,   420,   427,   122,   642,
     222,   123,   288,   211,   433,   290,   428,   165,   166,   429,
     211,   396,   122,   241,   365,   451,   448,   247,   123,   123,
     194,   270,   522,   194,   460,   123,   175,   186,   393,   447,
     420,   140,   129,   176,   222,   308,   140,   159,   160,   123,
     539,   520,   238,   130,   351,   581,   467,   416,   483,   481,
     138,   168,   484,   586,   138,   138,   132,   445,   482,   232,
     233,   234,   225,   279,   211,   524,   101,   525,   181,   333,
     270,   392,   222,   375,   175,   445,   393,   181,   167,   444,
     225,   176,   176,   446,   170,   191,   499,    61,   353,   173,
     222,   363,   283,   412,   284,   175,   169,   103,     3,     4,
       5,     6,   197,    61,   144,   145,   146,   198,   106,   108,
     109,   494,   495,   144,   145,   146,   138,   377,   436,   447,
     516,   437,   270,   517,   425,    11,   140,   547,   206,   548,
     140,   140,   222,   565,   488,   222,   569,   194,   570,   221,
     554,   333,   573,   386,   222,   181,   513,   251,   252,   181,
     556,   210,   137,   137,   496,   621,   211,   222,   557,   558,
     401,   231,   401,   508,   152,   285,    39,   138,   287,   512,
      18,    19,    20,   328,   351,   255,   256,   257,   258,   339,
     107,   147,   329,   625,   148,   222,   149,   349,   587,   559,
     147,   561,   140,   148,   141,   149,    61,   333,   628,   343,
     222,   344,   333,   269,   345,   270,   253,   254,    43,   377,
     347,    61,   396,   270,   259,   260,   346,   348,   350,   373,
     376,   389,   404,   275,   411,   176,   453,   152,   270,   406,
     459,   152,   377,   421,   464,   407,   408,   409,   241,   181,
     413,   181,   578,   140,   414,   417,   610,   442,   514,   473,
     418,   200,   419,   480,   201,   152,   152,   152,   152,   152,
     152,   152,   152,   152,   152,   152,   152,   152,   152,   152,
     152,   152,   152,   213,   214,   215,   424,   203,   181,   449,
      61,   452,    61,   502,   450,   610,   486,   455,   497,   498,
     456,    61,   475,   476,   138,   333,   485,    61,   492,   377,
      61,   518,   523,   526,   527,   529,   270,   585,   181,   532,
     377,   601,   535,   540,   541,   545,   546,   549,   566,   377,
     533,   611,   614,   534,   617,   552,   553,   538,   563,   564,
     574,   544,   152,   493,   567,   500,   590,   595,   588,   589,
     550,   591,   594,   551,   596,   620,   597,   618,   511,   619,
     622,   623,   630,    61,    61,    61,   635,   638,   637,    61,
     611,   604,   640,    61,    61,   583,   218,   380,   466,   555,
     140,   242,   133,   262,   103,     3,     4,     5,     6,     7,
     261,     8,     9,   333,   377,   641,   264,   333,   400,   582,
     342,   531,   584,   152,   426,   333,   263,   432,   423,   265,
     143,   202,    11,    61,   607,   592,    61,   626,   593,   571,
      61,     0,     0,    61,     0,    61,   631,   560,    61,   562,
       0,     0,     0,     0,     0,    61,     0,     0,   568,     0,
       0,     0,     0,   152,   572,     0,   502,   575,     0,     0,
     152,    37,     0,    39,     0,     0,     0,     0,    44,     0,
       0,     0,     0,     0,     0,     0,     0,    47,     0,     0,
       0,     0,     0,     0,    48,    49,    50,    51,    52,    53,
       0,     0,     0,     0,     0,     0,   152,     0,   280,   110,
     111,   112,   113,   114,   115,   116,   117,   118,   119,     0,
     598,   599,   600,     0,   152,     0,   606,     0,     0,     0,
     608,   609,     0,     0,     0,     0,     0,   183,     0,     0,
       0,     0,     0,     0,     0,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,     0,     0,     0,     0,     0,     0,     0,
     624,     0,     0,   627,     0,     0,     0,   629,    38,     0,
     633,    41,   634,    43,     0,   636,     0,     0,     0,     0,
       0,     0,   639,  -358,   142,   120,     2,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   152,     0,     0,   152,
     400,     0,   400,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,     0,     0,    46,     0,     0,     0,     0,    47,
       0,     0,     0,     0,     0,     0,    48,    49,    50,    51,
      52,    53,   291,     0,     2,     3,     4,     5,     6,     7,
       0,     8,     9,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,     0,   292,   293,
     294,     0,   295,   296,   297,   298,   299,   300,   301,   302,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
     303,   304,   305,     0,   306,   422,     0,    47,     0,     0,
       0,   101,     0,     0,    48,    49,    50,    51,    52,    53,
     291,     0,     2,     3,     4,     5,     6,     7,     0,     8,
       9,     0,     0,   110,   111,   112,   113,   114,   115,   116,
     117,   118,   119,     0,     0,     0,     0,     0,     0,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,     0,   292,   293,   294,     0,
     295,   296,   297,   298,   299,   300,   301,   302,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,   303,   304,
     305,  -210,   306,  -210,     0,    47,     0,     0,     0,   101,
       0,     0,    48,    49,    50,    51,    52,    53,     1,   120,
       2,     3,     4,     5,     6,     7,     0,     8,     9,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,     0,     0,    46,     0,
       0,     0,     0,    47,     0,     0,     0,     0,     0,     0,
      48,    49,    50,    51,    52,    53,   133,     0,   103,     3,
       4,     5,     6,     7,     0,     8,     9,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   134,    11,     0,     0,     0,
       0,     0,     0,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,     0,     0,   291,     0,     2,     3,     4,     5,     6,
       7,     0,     8,     9,     0,    37,    38,    39,     0,     0,
       0,    43,    44,     0,     0,     0,     0,     0,     0,     0,
       0,    47,   410,    11,     0,     0,     0,     0,    48,    49,
      50,    51,    52,    53,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   292,
     293,   294,     0,   295,   296,   297,   298,   299,   300,   301,
     302,     0,    37,     0,    39,     0,     0,     0,     0,    44,
       0,   303,   304,   305,     0,   306,     0,     0,    47,     0,
       0,     0,   101,     0,     0,    48,    49,    50,    51,    52,
      53,   103,     3,     4,     5,     6,     7,     0,     8,     9,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   134,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,     0,     0,     0,
       0,   306,     0,     0,    47,     0,     0,     0,     0,     0,
       0,    48,    49,    50,    51,    52,    53,   103,     3,     4,
       5,     6,     7,     0,     8,     9,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   183,    11,     0,     0,     0,     0,
       0,     0,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    37,    38,    39,     0,     0,     0,
      43,    44,     0,     0,   103,     3,     4,     5,     6,     7,
      47,     8,     9,     0,     0,     0,     0,    48,    49,    50,
      51,    52,    53,     0,     0,     0,     0,     0,     0,     0,
       0,   202,    11,     0,     0,   360,     0,     0,     0,    18,
      19,    20,     0,     0,     0,     0,     0,     0,   103,     3,
       4,     5,     6,     7,     0,     8,     9,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    37,   269,    39,     0,   202,    11,    43,    44,   383,
       0,     0,     0,    18,    19,    20,     0,    47,     0,     0,
       0,     0,     0,   361,   362,    49,    50,    51,    52,    53,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    37,   269,    39,     0,     0,
       0,    43,    44,     0,     0,   103,     3,     4,     5,     6,
       7,    47,     8,     9,     0,     0,     0,   384,   385,    49,
      50,    51,    52,    53,     0,     0,     0,     0,     0,     0,
       0,     0,   202,    11,     0,     0,   457,     0,     0,     0,
      18,    19,    20,     0,     0,     0,     0,     0,     0,   103,
       3,     4,     5,     6,     7,     0,     8,     9,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    37,   269,    39,     0,   202,    11,    43,    44,
     461,     0,     0,     0,    18,    19,    20,     0,    47,     0,
       0,     0,     0,     0,   458,    48,    49,    50,    51,    52,
      53,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    37,   269,    39,     0,
       0,     0,    43,    44,     0,     0,   103,     3,     4,     5,
       6,     7,    47,     8,     9,     0,     0,     0,   462,   463,
      49,    50,    51,    52,    53,     0,     0,     0,     0,     0,
       0,     0,     0,   202,    11,     0,     0,   477,     0,     0,
       0,    18,    19,    20,     0,     0,     0,     0,     0,     0,
     103,     3,     4,     5,     6,     7,     0,     8,     9,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    37,   269,    39,     0,   202,    11,    43,
      44,   542,     0,     0,     0,    18,    19,    20,     0,    47,
       0,     0,     0,     0,     0,   478,   479,    49,    50,    51,
      52,    53,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    37,   269,    39,
       0,     0,     0,    43,    44,     0,     0,   103,     3,     4,
       5,     6,     7,    47,     8,     9,     0,     0,     0,   543,
      48,    49,    50,    51,    52,    53,     0,     0,     0,     0,
       0,     0,     0,     0,   202,    11,     0,     0,     0,     0,
       0,     0,    18,    19,    20,   103,     3,     4,     5,     6,
       7,     0,     8,     9,     0,     0,     0,     0,     0,     0,
     103,     3,     4,     5,     6,     7,     0,     8,     9,     0,
       0,     0,   202,    11,    37,   269,    39,     0,     0,     0,
      43,    44,     0,     0,     0,     0,     0,   202,    11,     0,
      47,     0,     0,     0,     0,     0,     0,    48,    49,    50,
      51,    52,    53,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    37,     0,    39,     0,     0,     0,     0,    44,
       0,     0,     0,     0,     0,     0,   519,    37,    47,    39,
       0,   330,   331,   332,    44,    48,    49,    50,    51,    52,
      53,   580,     0,    47,     0,     0,   330,   331,   332,     0,
      48,    49,    50,    51,    52,    53,   103,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,   103,     3,     4,     5,     6,     7,     0,     8,     9,
       0,     0,     0,   202,    11,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   202,    11,
     103,     3,     4,     5,     6,     7,     0,     8,     9,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    37,     0,    39,     0,   202,    11,     0,
      44,     0,     0,     0,     0,     0,     0,     0,    37,    47,
      39,     0,   330,   331,   332,    44,    48,    49,    50,    51,
      52,    53,     0,   237,    47,     0,     0,     0,   203,     0,
       0,    48,    49,    50,    51,    52,    53,    37,     0,    39,
       0,     0,     0,     0,    44,     0,     0,     0,     0,   500,
       0,     0,   501,    47,     0,     0,     0,     0,     0,     0,
      48,    49,    50,    51,    52,    53,   103,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,   103,     3,     4,     5,     6,     7,     0,     8,     9,
       0,     0,     0,   202,    11,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   202,    11,
     103,     3,     4,     5,     6,     7,     0,     8,     9,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    37,     0,    39,     0,   202,    11,     0,
      44,     0,     0,     0,     0,     0,     0,     0,    37,    47,
      39,     0,     0,   203,     0,    44,    48,    49,    50,    51,
      52,    53,     0,     0,    47,     0,     0,     0,     0,     0,
     243,   244,    49,    50,    51,    52,    53,    37,     0,    39,
       0,     0,     0,     0,    44,     0,     0,   103,     3,     4,
       5,     6,     7,    47,     8,     9,     0,   331,     0,     0,
      48,    49,    50,    51,    52,    53,     0,     0,     0,     0,
       0,     0,     0,     0,   202,    11,   103,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   202,    11,   103,     3,     4,     5,     6,
       7,     0,     8,     9,    37,     0,    39,     0,     0,     0,
       0,    44,     0,     0,     0,     0,     0,   415,     0,     0,
      47,     0,   202,    11,     0,     0,     0,    48,    49,    50,
      51,    52,    53,    37,     0,    39,     0,     0,     0,     0,
      44,     0,     0,     0,     0,     0,     0,     0,   507,    47,
       0,     0,     0,     0,     0,     0,    48,    49,    50,    51,
      52,    53,    37,     0,    39,     0,     0,     0,     0,    44,
       0,     0,     0,     0,     0,   306,     0,     0,    47,     0,
       0,     0,     0,     0,     0,    48,    49,    50,    51,    52,
      53,   103,     3,     4,     5,     6,     7,     0,     8,     9,
       0,     0,     0,     0,     0,     0,   103,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,   202,    11,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   202,    11,   103,     3,     4,     5,     6,
       7,     0,     8,     9,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    37,     0,
      39,     0,   202,    11,     0,    44,     0,     0,     0,     0,
     500,     0,     0,    37,    47,    39,     0,     0,     0,     0,
      44,    48,    49,    50,    51,    52,    53,     0,   613,    47,
       0,     0,     0,     0,     0,     0,    48,    49,    50,    51,
      52,    53,    37,     0,    39,     0,     0,     0,     0,    44,
       0,     0,     0,     0,     0,     0,     0,   616,    47,     0,
       0,     0,     0,     0,     0,    48,    49,    50,    51,    52,
      53,   103,     3,     4,     5,     6,     7,     0,     8,     9,
       0,     0,     0,     0,   103,     3,     4,     5,     6,     7,
       0,     8,     9,     0,     0,     0,     0,     0,   202,    11,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    11,   103,     3,     4,     5,     6,     7,     0,
       8,     9,     0,     0,     0,     0,     0,     0,   103,     3,
       4,     5,     6,     7,     0,     8,     9,     0,    37,     0,
      39,    11,     0,     0,     0,    44,     0,     0,     0,     0,
       0,    37,     0,    39,    47,     0,    11,     0,    44,     0,
       0,    48,    49,    50,    51,    52,    53,    47,     0,     0,
       0,   203,     0,     0,    48,    49,    50,    51,    52,    53,
      37,     0,    39,     0,     0,     0,     0,    44,     0,     0,
       0,     0,     0,     0,     0,    37,   104,    39,     0,     0,
       0,     0,    44,    48,    49,    50,    51,    52,    53,     0,
       0,   107,     0,     0,     0,     0,     0,     0,    48,    49,
      50,    51,    52,    53,   103,     3,     4,     5,     6,     7,
       0,     8,     9,     0,     0,     0,     0,     0,   122,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    11,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   528,     0,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    37,     0,    39,     0,     0,     0,     0,    44,     0,
     468,     0,   469,     0,    36,     0,    38,    47,    40,     0,
      42,    43,     0,    45,    48,    49,    50,    51,    52,    53,
     354,   445,     0,     0,     0,     0,   225,     0,   176,   183,
       0,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,     0,   397,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    36,     0,
      38,     0,    40,     0,    42,    43,     0,    45,     0,     0,
       0,     0,     0,     0,   470,   183,     0,     0,     0,     0,
       0,     0,     0,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    38,     0,     0,    41,
     183,    43,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    36,
       0,    38,     0,    40,     0,    42,    43,     0,    45,     0,
       0,     0,     0,     0,     0,   354,   224,     0,     0,     0,
       0,   225,   183,   176,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    36,     0,    38,     0,    40,    41,    42,    43,     0,
      45,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   183,   101,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      36,     0,    38,     0,    40,     0,    42,    43,     0,    45,
     183,     0,     0,     0,     0,     0,   536,     0,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    38,     0,     0,    41,     0,    43,   183,     0,     0,
       0,     0,     0,     0,   390,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    38,     0,
       0,    41,     0,    43,   183,     0,     0,     0,     0,     0,
       0,   399,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    38,     0,     0,    41,     0,
      43,     0,     0,     0,     0,     0,     0,   183,   487,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,   530,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    36,     0,    38,     0,
      40,     0,    42,    43,   183,    45,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    36,     0,    38,   183,    40,     0,    42,
      43,     0,    45,     0,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   183,     0,    38,     0,     0,
      41,     0,    43,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,     0,     0,   438,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   183,     0,    38,     0,     0,     0,
       0,    43,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    38,     0,     0,     0,     0,
      43
};

static const yytype_int16 yycheck[] =
{
       0,   121,    47,   226,     0,   130,    62,   126,   182,   275,
     139,   130,   220,     3,   225,   191,    78,   281,     4,   341,
       1,   424,   286,   139,    92,   382,    47,   147,    86,     3,
     429,   151,    63,    91,     0,     0,    89,    84,    91,    35,
      87,   174,   138,    93,   140,     3,     3,   180,     3,    99,
     100,    14,    15,   106,    54,    19,    30,    90,    54,   104,
     176,    19,   107,    86,     3,    83,   277,    85,    91,    87,
     336,     3,    30,    30,     3,    30,     0,    83,    86,    86,
      61,    87,    78,   104,    91,     4,   107,    83,    54,    54,
     186,    30,   186,    89,   227,   191,    87,   191,    30,    85,
     106,    92,    47,   460,   149,   126,   127,   128,   224,   285,
     513,    92,   515,    79,    80,   437,   106,   138,   517,   140,
      86,    87,    88,   154,   155,   156,    92,   172,   196,    86,
      54,   199,   194,    88,   398,    90,    93,    95,   101,   102,
     273,   105,   350,   276,    83,    86,   279,    86,   144,   360,
       3,     3,   148,     3,    93,    84,   272,   368,    90,   104,
       3,   142,   107,   429,    86,   186,    85,   292,   287,    91,
     191,    93,   383,   292,    83,    95,    96,    30,    87,   175,
      30,   126,   127,   128,     3,   281,   182,   281,    86,   285,
     286,   285,   286,   138,    85,   140,    87,   121,   194,   199,
     231,     3,     3,    91,    87,   129,    83,   332,     3,    92,
      87,    30,   274,   332,   339,   196,    84,    16,    17,    87,
     339,   283,     3,   147,   447,   358,   355,   151,    30,    30,
     196,   176,    84,   199,   367,    30,    86,    90,    88,   355,
      83,   186,    86,    93,    87,   313,   191,    12,    13,    30,
     461,   517,   372,    86,   220,   521,   375,   302,    83,   392,
     281,   103,    87,   529,   285,   286,    86,    86,   393,    83,
      84,    85,    91,   447,   393,    85,    90,    87,   274,   203,
     225,    83,    87,    88,    86,    86,    88,   283,    94,   351,
      91,    93,    93,   355,    18,    90,    83,   297,   222,    83,
      87,   225,   398,   299,   398,    86,   104,     3,     4,     5,
       6,     7,     0,   313,     9,    10,    11,    84,     7,     8,
       9,   407,   408,     9,    10,    11,   347,   272,    84,   445,
      84,    87,   277,    87,   330,    31,   281,    85,     3,    87,
     285,   286,    87,    85,   400,    87,    85,   313,    87,    85,
     483,   275,    85,   277,    87,   351,   424,   157,   158,   355,
     485,    90,   407,   408,   409,    85,   485,    87,   488,   488,
     489,    85,   491,   418,    63,    90,    72,   398,    86,   424,
      38,    39,    40,    85,   350,   161,   162,   163,   164,   106,
      86,    86,    85,    85,    89,    87,    91,    92,   531,   489,
      86,   491,   347,    89,    47,    91,   406,   331,    85,    85,
      87,    85,   336,    71,    85,   360,   159,   160,    76,   364,
      87,   421,   484,   368,   165,   166,    85,    87,    85,    92,
      85,    84,    83,   106,   297,    93,   360,   126,   383,    88,
     364,   130,   387,    88,   368,    86,    86,    86,   372,   445,
      83,   447,   514,   398,    83,    86,   576,     6,   424,   383,
      86,   104,    86,   387,   107,   154,   155,   156,   157,   158,
     159,   160,   161,   162,   163,   164,   165,   166,   167,   168,
     169,   170,   171,   126,   127,   128,    86,    90,   484,    85,
     490,    85,   492,   417,    87,   615,    84,    92,    62,    83,
      92,   501,    92,    92,   525,   429,    88,   507,    88,   454,
     510,    92,    88,    88,    85,   106,   461,    83,   514,    92,
     465,   566,    92,    92,    92,    85,    85,    92,    86,   474,
     454,   576,   577,   457,   579,    92,    92,   461,    85,    85,
      85,   465,   231,   406,    86,    82,    85,     3,    92,    92,
     474,    92,    92,   477,    92,    60,    92,    92,   421,    92,
      88,    85,    83,   563,   564,   565,    85,    85,    91,   569,
     615,   567,    88,   573,   574,   525,   131,   274,   372,   484,
     525,   148,     1,   168,     3,     4,     5,     6,     7,     8,
     167,    10,    11,   517,   539,   640,   170,   521,   287,   523,
     210,   450,   526,   292,   331,   529,   169,   337,   313,   171,
      54,    30,    31,   613,   570,   539,   616,   615,   542,   505,
     620,    -1,    -1,   623,    -1,   625,   622,   490,   628,   492,
      -1,    -1,    -1,    -1,    -1,   635,    -1,    -1,   501,    -1,
      -1,    -1,    -1,   332,   507,    -1,   570,   510,    -1,    -1,
     339,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    86,    -1,    -1,
      -1,    -1,    -1,    -1,    93,    94,    95,    96,    97,    98,
      -1,    -1,    -1,    -1,    -1,    -1,   375,    -1,     1,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    -1,
     563,   564,   565,    -1,   393,    -1,   569,    -1,    -1,    -1,
     573,   574,    -1,    -1,    -1,    -1,    -1,    30,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     613,    -1,    -1,   616,    -1,    -1,    -1,   620,    71,    -1,
     623,    74,   625,    76,    -1,   628,    -1,    -1,    -1,    -1,
      -1,    -1,   635,     0,     1,   106,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   485,    -1,    -1,   488,
     489,    -1,   491,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    55,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    69,    70,    71,    72,    73,    74,    75,    76,
      77,    78,    -1,    -1,    81,    -1,    -1,    -1,    -1,    86,
      -1,    -1,    -1,    -1,    -1,    -1,    93,    94,    95,    96,
      97,    98,     1,    -1,     3,     4,     5,     6,     7,     8,
      -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,    55,    -1,    57,    58,
      59,    -1,    61,    62,    63,    64,    65,    66,    67,    68,
      69,    70,    71,    72,    73,    74,    75,    76,    77,    78,
      79,    80,    81,    -1,    83,    84,    -1,    86,    -1,    -1,
      -1,    90,    -1,    -1,    93,    94,    95,    96,    97,    98,
       1,    -1,     3,     4,     5,     6,     7,     8,    -1,    10,
      11,    -1,    -1,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    -1,    -1,    -1,    -1,    -1,    -1,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      51,    52,    53,    54,    55,    -1,    57,    58,    59,    -1,
      61,    62,    63,    64,    65,    66,    67,    68,    69,    70,
      71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
      81,    88,    83,    90,    -1,    86,    -1,    -1,    -1,    90,
      -1,    -1,    93,    94,    95,    96,    97,    98,     1,   106,
       3,     4,     5,     6,     7,     8,    -1,    10,    11,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    69,    70,    71,    72,
      73,    74,    75,    76,    77,    78,    -1,    -1,    81,    -1,
      -1,    -1,    -1,    86,    -1,    -1,    -1,    -1,    -1,    -1,
      93,    94,    95,    96,    97,    98,     1,    -1,     3,     4,
       5,     6,     7,     8,    -1,    10,    11,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    30,    31,    -1,    -1,    -1,
      -1,    -1,    -1,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    -1,    -1,     1,    -1,     3,     4,     5,     6,     7,
       8,    -1,    10,    11,    -1,    70,    71,    72,    -1,    -1,
      -1,    76,    77,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    86,    30,    31,    -1,    -1,    -1,    -1,    93,    94,
      95,    96,    97,    98,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    57,
      58,    59,    -1,    61,    62,    63,    64,    65,    66,    67,
      68,    -1,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,
      -1,    79,    80,    81,    -1,    83,    -1,    -1,    86,    -1,
      -1,    -1,    90,    -1,    -1,    93,    94,    95,    96,    97,
      98,     3,     4,     5,     6,     7,     8,    -1,    10,    11,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    51,
      52,    53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    69,    70,    71,
      72,    73,    74,    75,    76,    77,    78,    -1,    -1,    -1,
      -1,    83,    -1,    -1,    86,    -1,    -1,    -1,    -1,    -1,
      -1,    93,    94,    95,    96,    97,    98,     3,     4,     5,
       6,     7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    30,    31,    -1,    -1,    -1,    -1,
      -1,    -1,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    70,    71,    72,    -1,    -1,    -1,
      76,    77,    -1,    -1,     3,     4,     5,     6,     7,     8,
      86,    10,    11,    -1,    -1,    -1,    -1,    93,    94,    95,
      96,    97,    98,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    30,    31,    -1,    -1,    34,    -1,    -1,    -1,    38,
      39,    40,    -1,    -1,    -1,    -1,    -1,    -1,     3,     4,
       5,     6,     7,     8,    -1,    10,    11,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    70,    71,    72,    -1,    30,    31,    76,    77,    34,
      -1,    -1,    -1,    38,    39,    40,    -1,    86,    -1,    -1,
      -1,    -1,    -1,    92,    93,    94,    95,    96,    97,    98,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    70,    71,    72,    -1,    -1,
      -1,    76,    77,    -1,    -1,     3,     4,     5,     6,     7,
       8,    86,    10,    11,    -1,    -1,    -1,    92,    93,    94,
      95,    96,    97,    98,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    30,    31,    -1,    -1,    34,    -1,    -1,    -1,
      38,    39,    40,    -1,    -1,    -1,    -1,    -1,    -1,     3,
       4,     5,     6,     7,     8,    -1,    10,    11,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    70,    71,    72,    -1,    30,    31,    76,    77,
      34,    -1,    -1,    -1,    38,    39,    40,    -1,    86,    -1,
      -1,    -1,    -1,    -1,    92,    93,    94,    95,    96,    97,
      98,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    70,    71,    72,    -1,
      -1,    -1,    76,    77,    -1,    -1,     3,     4,     5,     6,
       7,     8,    86,    10,    11,    -1,    -1,    -1,    92,    93,
      94,    95,    96,    97,    98,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    31,    -1,    -1,    34,    -1,    -1,
      -1,    38,    39,    40,    -1,    -1,    -1,    -1,    -1,    -1,
       3,     4,     5,     6,     7,     8,    -1,    10,    11,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    70,    71,    72,    -1,    30,    31,    76,
      77,    34,    -1,    -1,    -1,    38,    39,    40,    -1,    86,
      -1,    -1,    -1,    -1,    -1,    92,    93,    94,    95,    96,
      97,    98,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    70,    71,    72,
      -1,    -1,    -1,    76,    77,    -1,    -1,     3,     4,     5,
       6,     7,     8,    86,    10,    11,    -1,    -1,    -1,    92,
      93,    94,    95,    96,    97,    98,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    30,    31,    -1,    -1,    -1,    -1,
      -1,    -1,    38,    39,    40,     3,     4,     5,     6,     7,
       8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,
       3,     4,     5,     6,     7,     8,    -1,    10,    11,    -1,
      -1,    -1,    30,    31,    70,    71,    72,    -1,    -1,    -1,
      76,    77,    -1,    -1,    -1,    -1,    -1,    30,    31,    -1,
      86,    -1,    -1,    -1,    -1,    -1,    -1,    93,    94,    95,
      96,    97,    98,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,
      -1,    -1,    -1,    -1,    -1,    -1,    84,    70,    86,    72,
      -1,    89,    90,    91,    77,    93,    94,    95,    96,    97,
      98,    84,    -1,    86,    -1,    -1,    89,    90,    91,    -1,
      93,    94,    95,    96,    97,    98,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,     3,     4,     5,     6,     7,     8,    -1,    10,    11,
      -1,    -1,    -1,    30,    31,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,
       3,     4,     5,     6,     7,     8,    -1,    10,    11,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    70,    -1,    72,    -1,    30,    31,    -1,
      77,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    70,    86,
      72,    -1,    89,    90,    91,    77,    93,    94,    95,    96,
      97,    98,    -1,    85,    86,    -1,    -1,    -1,    90,    -1,
      -1,    93,    94,    95,    96,    97,    98,    70,    -1,    72,
      -1,    -1,    -1,    -1,    77,    -1,    -1,    -1,    -1,    82,
      -1,    -1,    85,    86,    -1,    -1,    -1,    -1,    -1,    -1,
      93,    94,    95,    96,    97,    98,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,     3,     4,     5,     6,     7,     8,    -1,    10,    11,
      -1,    -1,    -1,    30,    31,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,
       3,     4,     5,     6,     7,     8,    -1,    10,    11,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    70,    -1,    72,    -1,    30,    31,    -1,
      77,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    70,    86,
      72,    -1,    -1,    90,    -1,    77,    93,    94,    95,    96,
      97,    98,    -1,    -1,    86,    -1,    -1,    -1,    -1,    -1,
      92,    93,    94,    95,    96,    97,    98,    70,    -1,    72,
      -1,    -1,    -1,    -1,    77,    -1,    -1,     3,     4,     5,
       6,     7,     8,    86,    10,    11,    -1,    90,    -1,    -1,
      93,    94,    95,    96,    97,    98,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    30,    31,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    31,     3,     4,     5,     6,     7,
       8,    -1,    10,    11,    70,    -1,    72,    -1,    -1,    -1,
      -1,    77,    -1,    -1,    -1,    -1,    -1,    83,    -1,    -1,
      86,    -1,    30,    31,    -1,    -1,    -1,    93,    94,    95,
      96,    97,    98,    70,    -1,    72,    -1,    -1,    -1,    -1,
      77,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,    86,
      -1,    -1,    -1,    -1,    -1,    -1,    93,    94,    95,    96,
      97,    98,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,
      -1,    -1,    -1,    -1,    -1,    83,    -1,    -1,    86,    -1,
      -1,    -1,    -1,    -1,    -1,    93,    94,    95,    96,    97,
      98,     3,     4,     5,     6,     7,     8,    -1,    10,    11,
      -1,    -1,    -1,    -1,    -1,    -1,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    30,    31,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    31,     3,     4,     5,     6,     7,
       8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    70,    -1,
      72,    -1,    30,    31,    -1,    77,    -1,    -1,    -1,    -1,
      82,    -1,    -1,    70,    86,    72,    -1,    -1,    -1,    -1,
      77,    93,    94,    95,    96,    97,    98,    -1,    85,    86,
      -1,    -1,    -1,    -1,    -1,    -1,    93,    94,    95,    96,
      97,    98,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,    86,    -1,
      -1,    -1,    -1,    -1,    -1,    93,    94,    95,    96,    97,
      98,     3,     4,     5,     6,     7,     8,    -1,    10,    11,
      -1,    -1,    -1,    -1,     3,     4,     5,     6,     7,     8,
      -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    30,    31,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    31,     3,     4,     5,     6,     7,     8,    -1,
      10,    11,    -1,    -1,    -1,    -1,    -1,    -1,     3,     4,
       5,     6,     7,     8,    -1,    10,    11,    -1,    70,    -1,
      72,    31,    -1,    -1,    -1,    77,    -1,    -1,    -1,    -1,
      -1,    70,    -1,    72,    86,    -1,    31,    -1,    77,    -1,
      -1,    93,    94,    95,    96,    97,    98,    86,    -1,    -1,
      -1,    90,    -1,    -1,    93,    94,    95,    96,    97,    98,
      70,    -1,    72,    -1,    -1,    -1,    -1,    77,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    70,    86,    72,    -1,    -1,
      -1,    -1,    77,    93,    94,    95,    96,    97,    98,    -1,
      -1,    86,    -1,    -1,    -1,    -1,    -1,    -1,    93,    94,
      95,    96,    97,    98,     3,     4,     5,     6,     7,     8,
      -1,    10,    11,    -1,    -1,    -1,    -1,    -1,     3,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    31,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    30,    -1,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,    -1,
       1,    -1,     3,    -1,    69,    -1,    71,    86,    73,    -1,
      75,    76,    -1,    78,    93,    94,    95,    96,    97,    98,
      85,    86,    -1,    -1,    -1,    -1,    91,    -1,    93,    30,
      -1,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      51,    52,    53,    54,    55,    -1,     1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    69,    -1,
      71,    -1,    73,    -1,    75,    76,    -1,    78,    -1,    -1,
      -1,    -1,    -1,    -1,    85,    30,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    71,    -1,    -1,    74,
      30,    76,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    69,
      -1,    71,    -1,    73,    -1,    75,    76,    -1,    78,    -1,
      -1,    -1,    -1,    -1,    -1,    85,    86,    -1,    -1,    -1,
      -1,    91,    30,    93,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    69,    -1,    71,    -1,    73,    74,    75,    76,    -1,
      78,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    30,    90,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,    55,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      69,    -1,    71,    -1,    73,    -1,    75,    76,    -1,    78,
      30,    -1,    -1,    -1,    -1,    -1,    85,    -1,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    71,    -1,    -1,    74,    -1,    76,    30,    -1,    -1,
      -1,    -1,    -1,    -1,    84,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    71,    -1,
      -1,    74,    -1,    76,    30,    -1,    -1,    -1,    -1,    -1,
      -1,    84,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    71,    -1,    -1,    74,    -1,
      76,    -1,    -1,    -1,    -1,    -1,    -1,    30,    84,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    56,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    69,    -1,    71,    -1,
      73,    -1,    75,    76,    30,    78,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    69,    -1,    71,    30,    73,    -1,    75,
      76,    -1,    78,    -1,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    51,    52,    53,
      54,    55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    30,    -1,    71,    -1,    -1,
      74,    -1,    76,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    -1,    -1,    58,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    30,    -1,    71,    -1,    -1,    -1,
      -1,    76,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    71,    -1,    -1,    -1,    -1,
      76
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     5,     6,     7,     8,    10,    11,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    69,    70,    71,    72,
      73,    74,    75,    76,    77,    78,    81,    86,    93,    94,
      95,    96,    97,    98,   108,   109,   111,   113,   114,   117,
     119,   120,   124,   125,   126,   127,   128,   129,   130,   131,
     132,   133,   134,   135,   136,   137,   138,   142,   143,   146,
     147,   148,   149,   150,   156,   159,   160,   161,   162,   165,
     181,   194,   195,   196,   199,   212,   215,   216,   217,   218,
     219,    90,   184,     3,    86,   120,   124,    86,   124,   124,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
     106,   139,     3,    30,    90,   165,    86,    86,    86,    86,
      86,    91,    86,     1,    30,   110,   138,   140,   147,   153,
     160,   172,     1,   194,     9,    10,    11,    86,    89,    91,
     184,   139,   124,   126,    93,    99,   100,    95,    96,    12,
      13,    14,    15,   101,   102,    16,    17,    94,   103,   104,
      18,    19,   105,    83,    83,    86,    93,   144,   145,   163,
     164,   165,   166,    30,   143,   143,    90,   165,   143,   143,
     143,    90,   165,   142,   143,   184,   197,     0,    84,   185,
     172,   172,    30,    90,   118,   138,     3,   112,   157,   158,
      90,   137,   141,   172,   172,   172,   138,   141,   117,     4,
      85,    85,    87,   153,    86,    91,   166,   173,   174,   175,
     153,    85,    83,    84,    85,   121,   165,    85,   118,   122,
     123,   138,   121,    92,    93,   140,   118,   138,   126,   126,
     126,   127,   127,   128,   128,   129,   129,   129,   129,   130,
     130,   131,   132,   133,   134,   135,   140,   175,   163,    71,
     160,   166,   167,    83,    87,   106,    86,    91,   175,   164,
       1,   151,   152,   153,   181,    90,   151,    86,   163,   142,
     184,     1,    57,    58,    59,    61,    62,    63,    64,    65,
      66,    67,    68,    79,    80,    81,    83,   140,   142,   165,
     182,   183,   184,   186,   187,   188,   189,   190,   191,   192,
     193,   198,   199,   200,   207,   208,   209,   210,    85,    85,
      89,    90,    91,   138,   176,   177,   178,   179,   180,   106,
      84,    87,   157,    85,    85,    85,    85,    87,    87,    92,
      85,   143,   196,   138,    85,   143,   168,   169,   170,   173,
      34,    92,    93,   138,   167,   174,   175,    86,    91,   118,
     126,    85,    87,    92,    92,    88,    85,   160,   166,   175,
     145,   176,   175,    34,    92,    93,   138,   167,   175,    84,
      84,   152,    83,    88,   154,   155,   163,     1,   151,    84,
     124,   137,   213,   214,    83,   141,    88,    86,    86,    86,
      30,   182,   165,    83,    83,    83,   140,    86,    86,    86,
      83,    88,    84,   187,    86,   165,   177,   141,    84,    87,
     176,   106,   180,   141,    84,   158,    84,    87,    58,   115,
     116,   172,     6,   196,   163,    86,   163,   166,   173,    85,
      87,   175,    85,   138,   167,    92,    92,    34,    92,   138,
     175,    34,    92,    93,   138,   167,   123,   137,     1,     3,
      85,   168,   171,   138,   167,    92,    92,    34,    92,    93,
     138,   175,   141,    83,    87,    88,    84,    84,   139,    83,
      85,    87,    88,   182,   110,   110,   140,    62,    83,    83,
      82,    85,   138,   201,   202,   203,   204,    85,   140,     4,
      85,   182,   140,   142,   143,   188,    84,    87,    92,    84,
     176,   178,    84,    88,    85,    87,    88,    85,    30,   106,
      56,   170,    92,   138,   138,    92,    85,   168,   138,   167,
      92,    92,    34,    92,   138,    85,    85,    85,    87,    92,
     138,   138,    92,    92,   175,   155,   141,   118,   137,   214,
     182,   214,   182,    85,    85,    85,    86,    86,   182,    85,
      87,   204,   182,    85,    85,   182,     3,   188,   163,   188,
      84,   176,   138,   116,   138,    83,   176,   175,    92,    92,
      85,    92,   138,   138,    92,     3,    92,    92,   182,   182,
     182,   140,    19,    95,   165,   205,   182,   202,   182,   182,
     118,   140,   211,    85,   140,     3,    85,   140,    92,    92,
      60,    85,    88,    85,   182,    85,   211,   182,    85,   182,
      83,   165,   206,   182,   182,    85,   182,    91,    85,   182,
      88,   140,    92
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   107,   108,   108,   108,   108,   108,   109,   109,   109,
     109,   109,   110,   110,   111,   111,   111,   112,   113,   113,
     114,   115,   115,   116,   116,   117,   117,   117,   117,   117,
     117,   117,   117,   118,   118,   119,   119,   119,   120,   120,
     121,   122,   122,   122,   123,   123,   124,   124,   124,   124,
     124,   124,   124,   124,   125,   125,   125,   125,   125,   125,
     126,   126,   127,   127,   127,   127,   128,   128,   128,   129,
     129,   129,   130,   130,   130,   130,   130,   131,   131,   131,
     132,   132,   133,   133,   134,   134,   135,   135,   136,   136,
     137,   137,   138,   138,   138,   138,   138,   139,   139,   139,
     139,   139,   139,   139,   139,   139,   139,   139,   140,   140,
     141,   142,   142,   142,   143,   143,   143,   143,   143,   143,
     143,   143,   143,   143,   144,   144,   145,   145,   146,   146,
     146,   146,   146,   146,   146,   147,   148,   148,   148,   148,
     148,   148,   148,   148,   148,   148,   148,   148,   148,   148,
     148,   148,   149,   149,   149,   149,   149,   150,   150,   151,
     151,   152,   152,   152,   153,   153,   153,   153,   154,   154,
     155,   155,   155,   156,   156,   156,   156,   156,   157,   157,
     158,   158,   159,   160,   160,   160,   160,   160,   161,   161,
     162,   162,   163,   163,   164,   164,   164,   164,   164,   164,
     164,   164,   164,   164,   164,   164,   164,   164,   164,   165,
     165,   166,   166,   166,   166,   167,   167,   168,   168,   169,
     169,   170,   170,   170,   170,   171,   171,   172,   172,   173,
     173,   173,   174,   174,   174,   174,   174,   174,   174,   174,
     174,   174,   174,   174,   174,   174,   174,   174,   174,   174,
     174,   174,   174,   175,   176,   176,   176,   177,   177,   177,
     177,   178,   179,   179,   180,   180,   181,   182,   182,   182,
     182,   182,   182,   182,   182,   183,   183,   183,   184,   185,
     184,   186,   186,   187,   187,   188,   188,   189,   189,   189,
     190,   191,   191,   191,   191,   191,   192,   192,   193,   193,
     193,   193,   193,   194,   194,   194,   194,   194,   194,   194,
     194,   195,   196,   196,   197,   197,   198,   198,   198,   198,
     198,   198,   199,   200,   200,   201,   201,   202,   202,   203,
     203,   204,   205,   205,   205,   206,   206,   207,   207,   208,
     208,   209,   210,   211,   211,   212,   213,   213,   213,   214,
     214,   214,   215,   216,   216,   217,   218,   218,   219
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     1,     2,     3,     3,     3,     1,     1,     1,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       6,     1,     3,     3,     3,     1,     1,     1,     3,     3,
       2,     2,     4,     3,     4,     3,     4,     4,     3,     4,
       1,     1,     3,     2,     1,     1,     1,     2,     2,     2,
       2,     4,     4,     1,     1,     1,     1,     1,     1,     1,
       1,     4,     1,     3,     3,     3,     1,     3,     3,     1,
       3,     3,     1,     3,     3,     3,     3,     1,     3,     3,
       1,     3,     1,     3,     1,     3,     1,     3,     1,     3,
       1,     5,     1,     3,     3,     3,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     3,
       1,     3,     4,     1,     2,     1,     2,     1,     2,     1,
       2,     1,     2,     1,     1,     3,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     4,     5,     2,     4,     5,     1,     1,     1,
       2,     3,     4,     1,     2,     1,     2,     1,     1,     3,
       2,     3,     1,     4,     5,     5,     6,     2,     1,     3,
       3,     1,     4,     1,     1,     1,     1,     1,     1,     1,
       4,     4,     3,     2,     1,     3,     3,     4,     6,     5,
       5,     6,     5,     4,     4,     5,     5,     4,     5,     1,
       1,     3,     2,     2,     1,     1,     2,     3,     1,     2,
       4,     2,     4,     2,     1,     1,     3,     3,     2,     2,
       1,     1,     3,     2,     3,     5,     4,     5,     4,     3,
       3,     3,     4,     6,     5,     5,     6,     4,     4,     2,
       3,     4,     5,     0,     3,     4,     1,     2,     1,     4,
       3,     2,     1,     2,     3,     2,     7,     1,     1,     1,
       1,     1,     1,     1,     2,     3,     4,     3,     2,     0,
       4,     1,     2,     1,     1,     1,     2,     7,     5,     5,
       1,     5,     7,     6,     7,     1,     6,     7,     3,     2,
       2,     2,     3,     1,     1,     1,     1,     1,     1,     1,
       2,     2,     3,     2,     1,     2,     1,     1,     1,     1,
       1,     1,     2,     4,     5,     1,     3,     1,     1,     1,
       2,     6,     1,     1,     1,     1,     5,     4,     5,     4,
       5,     8,     7,     1,     1,     6,     1,     3,     3,     1,
       3,     3,     2,     4,     5,     4,     1,     4,     1
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
        yyerror (parse, root, YY_("syntax error: cannot back up")); \
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
                  Kind, Value, parse, root); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, AstRoot * parse, Ast * root)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  YY_USE (parse);
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
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, AstRoot * parse, Ast * root)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep, parse, root);
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
                 int yyrule, AstRoot * parse, Ast * root)
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
                       &yyvsp[(yyi + 1) - (yynrhs)], parse, root);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, parse, root); \
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
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep, AstRoot * parse, Ast * root)
{
  YY_USE (yyvaluep);
  YY_USE (parse);
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
yyparse (AstRoot * parse, Ast * root)
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
      yychar = yylex (&yylval, parse);
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


 YY_REDUCE_PRINT (yyn); DEFAULT_ACTION (yyn);
  switch (yyn)
    {
  case 4: /* translation_unit: translation_unit error ';'  */
#line 107 "basilisk.yacc"
                                      { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 2373 "basilisk.c"
    break;

  case 5: /* translation_unit: translation_unit error '}'  */
#line 108 "basilisk.yacc"
                                      { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 2379 "basilisk.c"
    break;

  case 6: /* translation_unit: translation_unit error ')'  */
#line 109 "basilisk.yacc"
                                      { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 2385 "basilisk.c"
    break;

  case 13: /* expression_error: error  */
#line 122 "basilisk.yacc"
                        { yyvsp[0]->sym = YYSYMBOL_YYerror; }
#line 2391 "basilisk.c"
    break;

  case 111: /* declaration: declaration_specifiers ';' type_not_specified  */
#line 318 "basilisk.yacc"
                                                                               { ast_push_declaration (parse->stack, yyval); }
#line 2397 "basilisk.c"
    break;

  case 112: /* declaration: declaration_specifiers init_declarator_list ';' type_not_specified  */
#line 319 "basilisk.yacc"
                                                                               { ast_push_declaration (parse->stack, yyval); }
#line 2403 "basilisk.c"
    break;

  case 135: /* type_specifier: types  */
#line 357 "basilisk.yacc"
                   { parse->type_already_specified = true; }
#line 2409 "basilisk.c"
    break;

  case 155: /* struct_or_union_specifier: struct_or_union '{' error '}'  */
#line 383 "basilisk.yacc"
                                                            { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 2415 "basilisk.c"
    break;

  case 156: /* struct_or_union_specifier: struct_or_union generic_identifier '{' error '}'  */
#line 384 "basilisk.yacc"
                                                            { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 2421 "basilisk.c"
    break;

  case 206: /* direct_declarator: direct_declarator '(' type_not_specified error ')'  */
#line 479 "basilisk.yacc"
                                                                                           { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 2427 "basilisk.c"
    break;

  case 253: /* type_not_specified: %empty  */
#line 559 "basilisk.yacc"
           { parse->type_already_specified = false; }
#line 2433 "basilisk.c"
    break;

  case 274: /* statement: error ';'  */
#line 601 "basilisk.yacc"
                     { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 2439 "basilisk.c"
    break;

  case 279: /* @1: %empty  */
#line 613 "basilisk.yacc"
                           { stack_push (parse->stack, &(yyvsp[0])); yyval->sym = YYSYMBOL_YYUNDEF; }
#line 2445 "basilisk.c"
    break;

  case 280: /* compound_statement: '{' @1 block_item_list '}'  */
#line 615 "basilisk.yacc"
                           { ast_pop_scope (parse->stack, yyvsp[-3]); }
#line 2451 "basilisk.c"
    break;

  case 290: /* for_scope: FOR  */
#line 640 "basilisk.yacc"
              { stack_push (parse->stack, &(yyval)); }
#line 2457 "basilisk.c"
    break;

  case 293: /* iteration_statement: for_scope '(' expression_statement expression_statement ')' statement  */
#line 647 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-5]); }
#line 2463 "basilisk.c"
    break;

  case 294: /* iteration_statement: for_scope '(' expression_statement expression_statement expression ')' statement  */
#line 649 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-6]); }
#line 2469 "basilisk.c"
    break;

  case 296: /* for_declaration_statement: for_scope '(' declaration expression_statement ')' statement  */
#line 655 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-5]); }
#line 2475 "basilisk.c"
    break;

  case 297: /* for_declaration_statement: for_scope '(' declaration expression_statement expression ')' statement  */
#line 657 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-6]); }
#line 2481 "basilisk.c"
    break;

  case 310: /* external_declaration: error compound_statement  */
#line 676 "basilisk.yacc"
                                                { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 2487 "basilisk.c"
    break;

  case 311: /* function_declaration: declaration_specifiers declarator  */
#line 680 "basilisk.yacc"
                                            { ast_push_function_definition (parse->stack, yyvsp[0]);  }
#line 2493 "basilisk.c"
    break;

  case 312: /* function_definition: function_declaration declaration_list compound_statement  */
#line 685 "basilisk.yacc"
                                                { ast_pop_scope (parse->stack, yyvsp[-2]->child[1]); }
#line 2499 "basilisk.c"
    break;

  case 313: /* function_definition: function_declaration compound_statement  */
#line 687 "basilisk.yacc"
                                                { ast_pop_scope (parse->stack, yyvsp[-1]->child[1]); }
#line 2505 "basilisk.c"
    break;

  case 341: /* forin_declaration_statement: for_scope '(' declaration_specifiers declarator IDENTIFIER forin_arguments ')' statement  */
#line 758 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-7]); }
#line 2511 "basilisk.c"
    break;

  case 342: /* forin_statement: for_scope '(' expression IDENTIFIER forin_arguments ')' statement  */
#line 763 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-6]); }
#line 2517 "basilisk.c"
    break;

  case 358: /* root: translation_unit  */
#line 806 "basilisk.yacc"
                           {
	  yyval = root;
	  yyval->sym = yyr1[yyn];
	  yyval->child = allocate ((Allocator *)parse->alloc, 2*sizeof(Ast *));
	  yyval->child[0] = yyvsp[0];
	  yyvsp[0]->parent = yyval;
	  yyval->child[1] = NULL;
        }
#line 2530 "basilisk.c"
    break;


#line 2534 "basilisk.c"

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
      yyerror (parse, root, YY_("syntax error"));
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
                      yytoken, &yylval, parse, root);
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
                  YY_ACCESSING_SYMBOL (yystate), yyvsp, parse, root);
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
  yyerror (parse, root, YY_("memory exhausted"));
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
                  yytoken, &yylval, parse, root);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp, parse, root);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 816 "basilisk.yacc"


/**
# Parsing functions */

/* Called by yyparse on error.  */
void
yyerror (AstRoot * parse, Ast * root, char const *s)
{
#if 0
  fprintf (stderr, "%d: %s near '", *line, s);
  char * s1 = *input - 1;
  while (!strchr("}{;\n", *s1)) s1--;
  s1++;
  while (strchr(" \t", *s1)) s1++;
  for (; s1 < *input; s1++)
    fputc (*s1, stderr);
  fputs ("'\n", stderr);
#endif
}

static char * copy_range (const char * start, const char * end, long offset)
{
  char * c = NULL;
  int len = end - start;
  if (len > 0) {
    char * s = c = malloc (len + 1);
    for (const char * i = start; i < end; i++, s++)
      *s = *(i + offset);
    *s = '\0';
  }
  return c;
}

static const char * copy_strings (const char * i, Ast * n, long offset)
{
  AstTerminal * t = ast_terminal (n);
  if (t) {
    t->before = copy_range (i, t->start, offset);
    if (t->start > i)
      i = t->start;

    t->start = copy_range (i, t->after + 1, offset);
    if (t->after + 1 > i)
      i = t->after + 1;
    t->after = NULL;
  }
  else
    for (Ast ** c = n->child; *c; c++)
      i = copy_strings (i, *c, offset);
  return i;
}

static void remove_child (Ast * c)
{
  if (!c->parent)
    return;
  Ast ** i = c->parent->child;
  for (; *i && *i != c; i++);
  assert (*i == c);
  for (; *i; i++)
    *i = *(i + 1);
  c->parent = NULL;
}

static Ast * ast_reduce (Allocator * alloc, int sym, Ast ** children, int n)
{
  Ast * ast = allocate (alloc, sizeof(Ast));
  memset (ast, 0, sizeof(Ast));
  ast->sym = sym;
  int ndef = 0;
  for (int i = 0; i < n; i++) {
    Ast * c = children[i + 1 - n];
    if (c->sym != YYSYMBOL_YYUNDEF && c->sym != YYSYMBOL_type_not_specified)
      ndef++;
  }
  ast->child = allocate (alloc, (ndef + 1)*sizeof(Ast *));
  ndef = 0;
  for (int i = 0; i < n; i++) {
    Ast * c = children[i + 1 - n];
    if (c->sym != YYSYMBOL_YYUNDEF && c->sym != YYSYMBOL_type_not_specified) {
      if (c->parent)
	remove_child (c);
      c->parent = ast;
      ast->child[ndef++] = c;
    }
    else
      assert (!c->parent);
  }
  ast->child[ndef] = NULL;
  return ast;
}

static Stack * stack_internalize (Stack * stack)
{
  Ast ** n;
  for (int i = 0; (n = stack_index (stack, i)); i++)
    if ((*n)->sym == YYSYMBOL_IDENTIFIER) {
      AstTerminal * t = ast_terminal (*n);
      char * after = t->start + strlen (t->start) - 1;
      if (t->after != NULL && t->after != after) {
	fprintf (stderr, "%s:%d: %s after: %s\n",
		 t->file, t->line, t->start, t->after);
	abort();
      }
      t->after = after;
    }
  return stack;
}

static void stack_externalize (Stack * stack)
{
  Ast ** n;
  for (int i = 0; (n = stack_index (stack, i)); i++)
    if ((*n)->sym == YYSYMBOL_IDENTIFIER) {
      AstTerminal * t = ast_terminal(*n);
      if (t->after != NULL) {
	if (t->after[1] != '\0') {
	  
	  /**
	  This is a declaration which has not been through
	  copy_strings() i.e. which is not connected to the root, due
	  to a syntax error which lead to the corresponding branch being
	  discarded. We set the symbol to UNDEF. */

	  t->start = t->before = NULL;
	  (*n)->sym = YYSYMBOL_YYUNDEF;
	}
	t->after = NULL;
      }
    }
}

AstRoot * ast_parse (const char * code, AstRoot * parent)
{
  AstRoot parse;
  parse.alloc = parent ? parent->alloc : new_allocator();
  parse.stack = parent ? parent->stack : stack_new (sizeof (Ast *));
  parse.file = allocate (parse.alloc, strlen ("<basilisk>") + 1);
  strcpy ((char *) parse.file, "<basilisk>");
  stack_internalize (parse.stack);
  parse.type_already_specified = false;
  extern void lexer_setup (char * buffer, size_t len);
  size_t len = strlen (code) + 1;
  char * buffer = malloc (len + 1);
  memcpy (buffer, code, len);
  buffer[len] = '\0';
  lexer_setup (buffer, len + 1);
  //  yydebug = 1;
  AstRoot * root = allocate ((Allocator *)parse.alloc, sizeof(AstRoot));
  memset (root, 0, sizeof(AstRoot));
  stack_push (parse.stack, &root);
  yyparse (&parse, (Ast *) root);
  if (((Ast *)root)->child) {
    const char * i = copy_strings (buffer, (Ast *) root, code - buffer);
    const char * end = i; while (*end != '\0') end++;
    root->after = copy_range (i, end, code - buffer);
    root->alloc = parent ? NULL : parse.alloc;
    root->stack = parent ? NULL : parse.stack;
    stack_externalize (parse.stack);
  }
  else {
    root = NULL;
    if (parent)
      stack_externalize (parse.stack);
    else {
      free_allocator (parse.alloc);
      stack_destroy (parse.stack);
    }
  }
  free (buffer);
  yylex_destroy();
  return root;
}

int token_symbol (int token)
{
  return YYTRANSLATE (token);
}

const char * symbol_name (int sym)
{
  return yytname[sym];
}
