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
#line 59 "basilisk.yacc"

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
  YYSYMBOL_FOREACH_DIMENSION = 79,         /* FOREACH_DIMENSION  */
  YYSYMBOL_REDUCTION = 80,                 /* REDUCTION  */
  YYSYMBOL_MACRO = 81,                     /* MACRO  */
  YYSYMBOL_ELLIPSIS_MACRO = 82,            /* ELLIPSIS_MACRO  */
  YYSYMBOL_MACRODEF = 83,                  /* MACRODEF  */
  YYSYMBOL_foreach_statement = 84,         /* foreach_statement  */
  YYSYMBOL_85_ = 85,                       /* ';'  */
  YYSYMBOL_86_ = 86,                       /* '}'  */
  YYSYMBOL_87_ = 87,                       /* ')'  */
  YYSYMBOL_88_ = 88,                       /* '('  */
  YYSYMBOL_89_ = 89,                       /* ','  */
  YYSYMBOL_90_ = 90,                       /* ':'  */
  YYSYMBOL_91_ = 91,                       /* '.'  */
  YYSYMBOL_92_ = 92,                       /* '{'  */
  YYSYMBOL_93_ = 93,                       /* '['  */
  YYSYMBOL_94_ = 94,                       /* ']'  */
  YYSYMBOL_95_ = 95,                       /* '*'  */
  YYSYMBOL_96_ = 96,                       /* '&'  */
  YYSYMBOL_97_ = 97,                       /* '+'  */
  YYSYMBOL_98_ = 98,                       /* '-'  */
  YYSYMBOL_99_ = 99,                       /* '~'  */
  YYSYMBOL_100_ = 100,                     /* '!'  */
  YYSYMBOL_101_ = 101,                     /* '/'  */
  YYSYMBOL_102_ = 102,                     /* '%'  */
  YYSYMBOL_103_ = 103,                     /* '<'  */
  YYSYMBOL_104_ = 104,                     /* '>'  */
  YYSYMBOL_105_ = 105,                     /* '^'  */
  YYSYMBOL_106_ = 106,                     /* '|'  */
  YYSYMBOL_107_ = 107,                     /* '?'  */
  YYSYMBOL_108_ = 108,                     /* '='  */
  YYSYMBOL_YYACCEPT = 109,                 /* $accept  */
  YYSYMBOL_translation_unit = 110,         /* translation_unit  */
  YYSYMBOL_primary_expression = 111,       /* primary_expression  */
  YYSYMBOL_expression_error = 112,         /* expression_error  */
  YYSYMBOL_constant = 113,                 /* constant  */
  YYSYMBOL_enumeration_constant = 114,     /* enumeration_constant  */
  YYSYMBOL_string = 115,                   /* string  */
  YYSYMBOL_generic_selection = 116,        /* generic_selection  */
  YYSYMBOL_generic_assoc_list = 117,       /* generic_assoc_list  */
  YYSYMBOL_generic_association = 118,      /* generic_association  */
  YYSYMBOL_postfix_expression = 119,       /* postfix_expression  */
  YYSYMBOL_postfix_initializer = 120,      /* postfix_initializer  */
  YYSYMBOL_array_access = 121,             /* array_access  */
  YYSYMBOL_function_call = 122,            /* function_call  */
  YYSYMBOL_member_identifier = 123,        /* member_identifier  */
  YYSYMBOL_argument_expression_list = 124, /* argument_expression_list  */
  YYSYMBOL_argument_expression_list_item = 125, /* argument_expression_list_item  */
  YYSYMBOL_unary_expression = 126,         /* unary_expression  */
  YYSYMBOL_unary_operator = 127,           /* unary_operator  */
  YYSYMBOL_cast_expression = 128,          /* cast_expression  */
  YYSYMBOL_multiplicative_expression = 129, /* multiplicative_expression  */
  YYSYMBOL_additive_expression = 130,      /* additive_expression  */
  YYSYMBOL_shift_expression = 131,         /* shift_expression  */
  YYSYMBOL_relational_expression = 132,    /* relational_expression  */
  YYSYMBOL_equality_expression = 133,      /* equality_expression  */
  YYSYMBOL_and_expression = 134,           /* and_expression  */
  YYSYMBOL_exclusive_or_expression = 135,  /* exclusive_or_expression  */
  YYSYMBOL_inclusive_or_expression = 136,  /* inclusive_or_expression  */
  YYSYMBOL_logical_and_expression = 137,   /* logical_and_expression  */
  YYSYMBOL_logical_or_expression = 138,    /* logical_or_expression  */
  YYSYMBOL_conditional_expression = 139,   /* conditional_expression  */
  YYSYMBOL_assignment_expression = 140,    /* assignment_expression  */
  YYSYMBOL_assignment_operator = 141,      /* assignment_operator  */
  YYSYMBOL_expression = 142,               /* expression  */
  YYSYMBOL_constant_expression = 143,      /* constant_expression  */
  YYSYMBOL_declaration = 144,              /* declaration  */
  YYSYMBOL_declaration_specifiers = 145,   /* declaration_specifiers  */
  YYSYMBOL_init_declarator_list = 146,     /* init_declarator_list  */
  YYSYMBOL_init_declarator = 147,          /* init_declarator  */
  YYSYMBOL_storage_class_specifier = 148,  /* storage_class_specifier  */
  YYSYMBOL_type_specifier = 149,           /* type_specifier  */
  YYSYMBOL_types = 150,                    /* types  */
  YYSYMBOL_struct_or_union_specifier = 151, /* struct_or_union_specifier  */
  YYSYMBOL_struct_or_union = 152,          /* struct_or_union  */
  YYSYMBOL_struct_declaration_list = 153,  /* struct_declaration_list  */
  YYSYMBOL_struct_declaration = 154,       /* struct_declaration  */
  YYSYMBOL_specifier_qualifier_list = 155, /* specifier_qualifier_list  */
  YYSYMBOL_struct_declarator_list = 156,   /* struct_declarator_list  */
  YYSYMBOL_struct_declarator = 157,        /* struct_declarator  */
  YYSYMBOL_enum_specifier = 158,           /* enum_specifier  */
  YYSYMBOL_enumerator_list = 159,          /* enumerator_list  */
  YYSYMBOL_enumerator = 160,               /* enumerator  */
  YYSYMBOL_atomic_type_specifier = 161,    /* atomic_type_specifier  */
  YYSYMBOL_type_qualifier = 162,           /* type_qualifier  */
  YYSYMBOL_function_specifier = 163,       /* function_specifier  */
  YYSYMBOL_alignment_specifier = 164,      /* alignment_specifier  */
  YYSYMBOL_declarator = 165,               /* declarator  */
  YYSYMBOL_direct_declarator = 166,        /* direct_declarator  */
  YYSYMBOL_generic_identifier = 167,       /* generic_identifier  */
  YYSYMBOL_pointer = 168,                  /* pointer  */
  YYSYMBOL_type_qualifier_list = 169,      /* type_qualifier_list  */
  YYSYMBOL_parameter_type_list = 170,      /* parameter_type_list  */
  YYSYMBOL_parameter_list = 171,           /* parameter_list  */
  YYSYMBOL_parameter_declaration = 172,    /* parameter_declaration  */
  YYSYMBOL_identifier_list = 173,          /* identifier_list  */
  YYSYMBOL_type_name = 174,                /* type_name  */
  YYSYMBOL_abstract_declarator = 175,      /* abstract_declarator  */
  YYSYMBOL_direct_abstract_declarator = 176, /* direct_abstract_declarator  */
  YYSYMBOL_type_not_specified = 177,       /* type_not_specified  */
  YYSYMBOL_initializer = 178,              /* initializer  */
  YYSYMBOL_initializer_list = 179,         /* initializer_list  */
  YYSYMBOL_designation = 180,              /* designation  */
  YYSYMBOL_designator_list = 181,          /* designator_list  */
  YYSYMBOL_designator = 182,               /* designator  */
  YYSYMBOL_static_assert_declaration = 183, /* static_assert_declaration  */
  YYSYMBOL_statement = 184,                /* statement  */
  YYSYMBOL_labeled_statement = 185,        /* labeled_statement  */
  YYSYMBOL_compound_statement = 186,       /* compound_statement  */
  YYSYMBOL_187_1 = 187,                    /* @1  */
  YYSYMBOL_block_item_list = 188,          /* block_item_list  */
  YYSYMBOL_block_item = 189,               /* block_item  */
  YYSYMBOL_expression_statement = 190,     /* expression_statement  */
  YYSYMBOL_selection_statement = 191,      /* selection_statement  */
  YYSYMBOL_for_scope = 192,                /* for_scope  */
  YYSYMBOL_iteration_statement = 193,      /* iteration_statement  */
  YYSYMBOL_for_declaration_statement = 194, /* for_declaration_statement  */
  YYSYMBOL_jump_statement = 195,           /* jump_statement  */
  YYSYMBOL_external_declaration = 196,     /* external_declaration  */
  YYSYMBOL_function_declaration = 197,     /* function_declaration  */
  YYSYMBOL_function_definition = 198,      /* function_definition  */
  YYSYMBOL_declaration_list = 199,         /* declaration_list  */
  YYSYMBOL_basilisk_statements = 200,      /* basilisk_statements  */
  YYSYMBOL_macro_statement = 201,          /* macro_statement  */
  YYSYMBOL_reduction_list = 202,           /* reduction_list  */
  YYSYMBOL_reduction = 203,                /* reduction  */
  YYSYMBOL_reduction_operator = 204,       /* reduction_operator  */
  YYSYMBOL_reduction_array = 205,          /* reduction_array  */
  YYSYMBOL_foreach_dimension_statement = 206, /* foreach_dimension_statement  */
  YYSYMBOL_forin_declaration_statement = 207, /* forin_declaration_statement  */
  YYSYMBOL_forin_statement = 208,          /* forin_statement  */
  YYSYMBOL_forin_arguments = 209,          /* forin_arguments  */
  YYSYMBOL_event_definition = 210,         /* event_definition  */
  YYSYMBOL_event_parameters = 211,         /* event_parameters  */
  YYSYMBOL_event_parameter = 212,          /* event_parameter  */
  YYSYMBOL_boundary_definition = 213,      /* boundary_definition  */
  YYSYMBOL_external_foreach_dimension = 214, /* external_foreach_dimension  */
  YYSYMBOL_attribute = 215,                /* attribute  */
  YYSYMBOL_new_field = 216,                /* new_field  */
  YYSYMBOL_root = 217                      /* root  */
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
#define YYFINAL  203
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   3359

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  109
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  109
/* YYNRULES -- Number of rules.  */
#define YYNRULES  355
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  636

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   339


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
       2,     2,     2,   100,     2,     2,     2,   102,    96,     2,
      88,    87,    95,    97,    89,    98,    91,   101,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    90,    85,
     103,   108,   104,   107,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    93,     2,    94,   105,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    92,   106,    86,    99,     2,     2,     2,
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
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   116,   116,   117,   118,   119,   120,   124,   125,   126,
     127,   128,   129,   130,   134,   135,   139,   140,   141,   145,
     149,   150,   154,   158,   159,   163,   164,   168,   169,   170,
     171,   172,   173,   174,   175,   179,   180,   184,   185,   186,
     190,   191,   195,   199,   200,   201,   205,   206,   210,   211,
     212,   213,   214,   215,   216,   217,   221,   222,   223,   224,
     225,   226,   230,   231,   237,   238,   239,   240,   244,   245,
     246,   250,   251,   252,   256,   257,   258,   259,   260,   264,
     265,   266,   270,   271,   275,   276,   280,   281,   285,   286,
     290,   291,   295,   296,   300,   301,   302,   303,   304,   308,
     309,   310,   311,   312,   313,   314,   315,   316,   317,   318,
     322,   323,   327,   331,   332,   333,   337,   338,   339,   340,
     341,   342,   343,   344,   345,   346,   350,   351,   355,   356,
     360,   361,   362,   363,   364,   365,   366,   367,   371,   375,
     376,   377,   378,   379,   380,   381,   382,   383,   384,   385,
     386,   387,   388,   389,   390,   394,   395,   396,   397,   398,
     402,   403,   407,   408,   412,   413,   414,   418,   419,   420,
     421,   425,   426,   430,   431,   432,   436,   437,   438,   439,
     440,   444,   445,   449,   450,   454,   458,   459,   460,   461,
     462,   466,   467,   471,   472,   476,   477,   481,   482,   483,
     484,   485,   486,   487,   488,   489,   490,   491,   492,   493,
     494,   495,   499,   500,   501,   505,   506,   507,   508,   512,
     513,   517,   518,   522,   523,   527,   528,   529,   530,   531,
     532,   536,   537,   541,   542,   546,   547,   548,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   576,
     580,   581,   582,   586,   587,   588,   589,   593,   597,   598,
     602,   603,   607,   611,   612,   613,   614,   615,   616,   617,
     618,   622,   623,   624,   628,   630,   630,   636,   637,   641,
     642,   646,   647,   651,   652,   653,   657,   661,   662,   663,
     665,   667,   671,   673,   678,   679,   680,   681,   682,   686,
     687,   688,   689,   690,   691,   692,   696,   700,   702,   707,
     708,   715,   716,   717,   718,   719,   723,   724,   728,   729,
     733,   737,   738,   739,   743,   744,   748,   749,   753,   758,
     763,   764,   768,   772,   773,   774,   778,   779,   780,   784,
     788,   789,   793,   797,   798,   802
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
  "TRACE", "FOREACH_DIMENSION", "REDUCTION", "MACRO", "ELLIPSIS_MACRO",
  "MACRODEF", "foreach_statement", "';'", "'}'", "')'", "'('", "','",
  "':'", "'.'", "'{'", "'['", "']'", "'*'", "'&'", "'+'", "'-'", "'~'",
  "'!'", "'/'", "'%'", "'<'", "'>'", "'^'", "'|'", "'?'", "'='", "$accept",
  "translation_unit", "primary_expression", "expression_error", "constant",
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
  "reduction_list", "reduction", "reduction_operator", "reduction_array",
  "foreach_dimension_statement", "forin_declaration_statement",
  "forin_statement", "forin_arguments", "event_definition",
  "event_parameters", "event_parameter", "boundary_definition",
  "external_foreach_dimension", "attribute", "new_field", "root", YY_NULLPTR
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
     335,   336,   337,   338,   339,    59,   125,    41,    40,    44,
      58,    46,   123,    91,    93,    42,    38,    43,    45,   126,
      33,    47,    37,    60,    62,    94,   124,    63,    61
};
#endif

#define YYPACT_NINF (-410)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-356)

#define yytable_value_is_error(Yyn) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    1013,   -67,   163,  -410,  -410,  -410,  -410,  2500,  2547,  2547,
     552,  -410,  -410,  -410,  -410,  -410,  -410,  -410,  -410,  -410,
    -410,  -410,  -410,  -410,  -410,  -410,  -410,  -410,  -410,  -410,
    -410,  -410,  -410,  -410,  -410,    53,   -28,    -7,     4,    22,
    -410,    28,  -410,  -410,     2,  -410,    54,    80,   175,  -410,
    1113,  -410,  -410,  -410,  -410,  -410,  -410,   709,  -410,  -410,
    -410,  -410,   363,  -410,  -410,    17,  2579,  -410,   -29,    79,
     295,     8,   194,    90,    -2,   154,   178,    -1,  -410,   138,
    -410,   218,  3197,  3197,  -410,  -410,   198,  -410,  -410,  3197,
    3197,  3197,   271,  -410,  -410,  2881,  -410,   188,  -410,  -410,
    -410,  -410,  -410,  -410,   278,   203,  -410,  -410,  -410,  1113,
    -410,  1113,  -410,  -410,  -410,  -410,  -410,  -410,  -410,  -410,
    -410,  -410,  -410,  -410,  -410,  2111,  -410,  -410,  -410,   312,
     243,  1390,   586,   586,  2421,  2579,   421,    50,    49,  -410,
      17,   254,  -410,   264,   586,   171,   586,   304,   117,  -410,
     262,  -410,  -410,  1963,   262,  2142,  2111,  -410,  -410,  2579,
    2579,  2579,  2579,  2579,  2579,  2579,  2579,  2579,  2579,  2579,
    2579,  2579,  2579,  2579,  2579,  2579,  2579,  2421,  -410,  -410,
     228,   474,   -15,  -410,   -32,   -38,  -410,   321,  -410,  -410,
    -410,   647,   302,  -410,  -410,  -410,  3244,   313,  -410,   218,
    -410,  2881,  -410,  -410,  -410,   913,   317,   323,    17,  1932,
    -410,  -410,  -410,   306,   -54,  -410,   312,  -410,   335,   345,
     346,   349,   356,   361,   396,   366,  3197,  -410,  -410,  -410,
     369,  -410,  2421,  -410,  2817,  1443,    -8,  -410,   102,  -410,
    -410,  2468,  -410,  -410,  -410,  -410,  -410,  -410,  -410,   135,
    -410,  -410,  -410,  -410,   367,   182,  -410,  -410,  -410,  -410,
    -410,   -29,   -29,    79,    79,   295,   295,   295,   295,     8,
       8,   194,    90,    -2,   154,   178,   281,  -410,   379,  -410,
    -410,  -410,   474,  -410,   228,  2190,  -410,  1495,  -410,   -38,
     381,  2992,  -410,   162,  -410,  2772,  3041,  2579,   360,  -410,
    -410,   385,  2579,   383,   386,   388,   389,  1168,  -410,   262,
     393,   395,  2221,   394,   207,  -410,  -410,   232,  -410,   391,
    -410,  -410,  -410,   813,  -410,  -410,  -410,   397,  -410,  -410,
    -410,  -410,  -410,  -410,  -410,  -410,   402,   402,   262,  1932,
    2579,  -410,  -410,    15,  2190,   179,  -410,  2579,  -410,    62,
     275,  -410,  -410,  -410,  -410,  3283,   465,  -410,  3197,   228,
    -410,   262,  -410,   378,  -410,    48,   401,   409,  -410,   412,
    1755,  -410,   406,   411,  1547,   102,  -410,  -410,  1599,  -410,
    -410,  -410,  2111,  -410,  -410,  2579,  -410,  -410,  -410,  -410,
    -410,  -410,  2717,  1755,  -410,   413,   414,  1651,  -410,  -410,
    -410,  -410,  -410,  2579,   233,  -410,   416,   417,  3090,  -410,
      17,  -410,   215,  -410,  -410,   420,  1168,  1213,  1213,  2421,
      17,   454,   432,  -410,  -410,  -410,   242,    71,  2011,  -410,
    1168,  -410,  -410,  1292,  -410,   301,   424,  -410,  1803,  -410,
    -410,  -410,  -410,  -410,  -410,  -410,    74,   429,   239,  -410,
     431,   435,  -410,  -410,   433,   436,  2063,  2650,   419,   266,
    -410,  -410,  3145,  -410,  -410,   441,  1755,  -410,  -410,  2421,
    -410,   442,  2942,  1755,  -410,   444,   445,  1703,  -410,  -410,
     455,  -410,  -410,   459,   253,   453,  1755,  -410,  -410,  2421,
    -410,   457,   461,  -410,  -410,  -410,   224,  2579,  -410,  -410,
    2468,  2579,  1168,  2579,  1168,  -410,   462,   469,   308,   460,
    -410,  -410,   471,  1168,  1168,   311,  -410,    68,  2269,   218,
    2269,  -410,  1834,  -410,  -410,  -410,  2190,  -410,  2421,  -410,
    3283,  2421,   475,   439,  -410,  -410,  -410,  -410,  2190,  -410,
    -410,  -410,   468,   472,  -410,  -410,   476,   473,  1755,  -410,
    -410,  2421,  -410,   488,  -410,  -410,  -410,   556,  -410,   489,
     490,  -410,  -410,  -410,  -410,  -410,  -410,  -410,  -410,  -410,
    -410,  -410,  1168,  1168,  1168,  2421,  1168,  -410,  -410,  1168,
    2111,  2321,    12,  2369,  -410,  -410,  -410,  -410,  -410,  -410,
    2421,  -410,  -410,  -410,  -410,  -410,  -410,   491,   493,  -410,
    -410,  -410,  -410,   510,  -410,  -410,   324,  -410,  -410,  -410,
     264,   501,  1168,   328,  2111,  1168,   331,   186,  -410,  -410,
    1168,   522,  1168,  -410,  1168,   521,  -410,  1168,  -410,  -410,
    -410,  -410,  -410,  1168,  -410,  -410
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_int16 yydefact[] =
{
       0,     0,     7,    16,    17,    20,    21,     0,     0,     0,
     154,    18,   130,   131,   132,   134,   135,   191,   186,   187,
     188,   148,   140,   141,   142,   143,   146,   147,   144,   145,
     139,   149,   150,   160,   161,     0,     0,     0,   189,     0,
     192,     0,   133,   190,   353,   136,     0,     0,     8,   137,
       0,    57,    56,    58,    59,    60,    61,     0,    27,     9,
      10,    12,    48,    29,    28,    62,     0,    64,    68,    71,
      74,    79,    82,    84,    86,    88,    90,    92,    94,     0,
     310,     0,   117,   119,   138,   152,     0,   153,   151,   121,
     123,   125,     0,   115,     2,     0,   309,    13,   328,   311,
     312,   313,   314,    55,     0,   285,   315,     7,     8,     0,
      52,     0,    49,    50,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   109,    99,     0,   212,   213,   214,     0,
     180,     0,     0,     0,     0,     0,     0,     0,     0,    15,
     154,     0,   110,    14,   168,   259,   170,     0,     0,     3,
       0,    32,    33,     0,     0,     0,     0,    62,    51,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   349,   259,
       0,   218,     0,   126,   316,   259,   197,     0,   154,   116,
     118,     0,   157,   120,   122,   124,     0,     0,   319,     0,
     318,     0,   329,     1,   284,     0,     0,     0,     0,     0,
      98,    97,    19,   184,     0,   181,     0,   112,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   333,   332,   331,
       0,    11,     0,   167,     0,     0,   236,   259,   237,   234,
     169,     0,     4,     5,     6,    31,    42,    40,    47,     0,
      43,    46,    30,    37,    57,     0,    96,    95,    65,    66,
      67,    69,    70,    72,    73,    77,    78,    75,    76,    80,
      81,    83,    85,    87,    89,    91,     0,   113,     0,   189,
     219,   217,   216,   259,     0,     0,   259,     0,   196,   259,
       0,     0,   162,     0,   166,     0,     0,     0,   129,   320,
     317,     0,     0,     0,     0,     0,     0,     0,   296,     0,
       0,     0,     0,     0,     8,   321,   291,     0,   289,     0,
     290,   273,   274,     0,   287,   275,   276,     0,   277,   301,
     278,   279,   322,   323,   324,   325,    53,     0,     0,     0,
       0,   262,   264,     0,     0,     0,   268,     0,   176,     0,
       0,   194,   193,    54,   185,     0,     0,   354,     0,     0,
     350,     0,   111,     0,   255,   228,     0,   222,   259,     0,
       0,   239,    57,     0,     0,   235,   233,   259,     0,    34,
      63,    41,    45,    38,    39,     0,   198,   220,   215,   114,
     127,   128,     0,     0,   199,    57,     0,     0,   195,   158,
     155,   163,   259,     0,     0,   171,   175,     0,     0,   352,
      62,   346,     0,   343,   280,     0,     0,     0,     0,     0,
     213,     0,     0,   305,   306,   307,     0,     0,     0,   292,
       0,   286,   288,     0,   271,     0,     0,    35,     0,   263,
     267,   269,   183,   177,   182,   178,     0,     0,     0,    23,
       0,     0,   351,   316,   334,     0,     0,     0,   225,   236,
     227,   256,     0,   223,   238,     0,     0,   240,   246,     0,
     245,     0,     0,     0,   247,    57,     0,     0,    44,    93,
       0,   231,   210,     0,     0,     0,     0,   200,   207,     0,
     206,    57,     0,   164,   173,   259,     0,     0,   159,   156,
       0,     0,     0,     0,     0,   283,     0,     0,     0,     0,
     304,   308,     0,     0,     0,     0,   281,     0,     0,     0,
       0,   260,     0,   270,    36,   266,     0,   179,     0,    22,
       0,     0,     0,     0,   330,   230,   229,   154,     0,   221,
     259,   242,     0,     0,   244,   257,     0,     0,     0,   248,
     254,     0,   253,     0,   209,   208,   211,     0,   202,     0,
       0,   203,   205,   165,   172,   174,   348,   347,   345,   342,
     344,   282,     0,     0,     0,     0,     0,   336,   326,     0,
       0,     0,   129,     0,   261,   265,    26,    24,    25,   272,
       0,   226,   224,   241,   243,   258,   250,     0,     0,   251,
     232,   201,   204,   294,   295,   297,     0,   337,   327,   341,
     340,     0,     0,     0,     0,     0,     0,     0,   249,   252,
       0,     0,     0,   302,     0,     0,   299,     0,   335,   293,
     298,   339,   303,     0,   300,   338
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -410,  -410,  -410,  -134,  -410,  -410,  -410,  -410,  -410,    84,
     481,  -123,  -410,  -410,   456,   181,   237,   430,  -410,   -52,
     223,   265,   189,   277,   448,   449,   447,   470,   467,  -410,
    -122,     6,   -62,   -50,  -105,   -84,     7,  -410,   339,  -410,
      41,  -410,  -410,  -410,  -169,  -275,   -97,  -410,   150,  -410,
     434,  -318,  -410,    38,  -410,  -410,   -77,  -182,     1,  -119,
    -187,  -308,  -410,   185,  -410,   332,  -135,  -219,  -151,  -277,
     310,  -409,  -410,   307,  -172,    39,  -410,    19,  -410,  -410,
     330,  -364,  -410,  -410,  -410,  -410,  -410,   598,  -410,  -214,
    -410,  -410,  -410,  -410,   559,  -410,  -410,  -410,  -410,  -410,
      44,  -410,  -410,   -80,  -410,  -410,  -410,  -410,  -410
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
       0,    57,    58,   141,    59,   213,    60,    61,   448,   449,
      62,   248,    63,    64,   245,   249,   250,    65,    66,    67,
      68,    69,    70,    71,    72,    73,    74,    75,    76,    77,
      78,   142,   125,   317,   218,    80,   365,   182,   183,    82,
      83,    84,    85,    86,   291,   292,   145,   404,   405,    87,
     214,   215,    88,    89,    90,    91,   278,   185,   319,   187,
     282,   366,   367,   368,   484,   450,   369,   238,   239,   342,
     343,   344,   345,   346,    93,   320,   321,   322,   205,   323,
     324,   325,   326,   327,   328,   329,   330,    94,    95,    96,
     201,   331,   332,    97,    98,   230,   455,   333,   334,   335,
     611,    99,   412,   413,   100,   101,   102,   103,   104
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     143,    92,   210,   156,   184,   289,    79,    81,   391,   217,
     237,   198,   360,   217,   158,   614,   401,   375,   176,   294,
     106,   401,   166,   167,   294,   105,   236,   296,   277,   526,
     223,   444,   348,   256,   288,   349,   130,   114,   115,   116,
     117,   118,   119,   120,   121,   122,   123,   233,   374,   240,
     286,   126,   126,  -129,   225,   287,   126,  -129,    92,   143,
     131,   143,   281,    79,    81,   212,   159,   439,   227,   520,
     283,   580,   160,   161,   284,   512,   285,   212,   127,   127,
     234,   132,   186,   127,   483,   235,   376,   192,   146,   189,
     190,   144,   133,   197,   293,   136,   193,   194,   195,   293,
     397,   437,   199,   173,   438,   255,   177,   258,   259,   260,
     134,   168,   169,   526,   200,   236,   135,   299,   379,   294,
     285,   318,   298,   294,   294,   124,   408,   276,   444,   128,
     128,   211,   389,   401,   128,   392,   457,   226,   398,   229,
     222,   235,   137,   181,   452,   129,   228,   146,   443,   146,
     144,   246,   144,   429,   581,   246,   583,   232,   513,   251,
     527,   525,   257,   388,   546,   126,  -212,   106,   138,   146,
     146,   146,   144,   144,   144,   411,   162,   163,  -214,   536,
     217,   186,   146,   466,   146,   144,   172,   144,   186,   380,
     377,   477,   127,  -212,   293,   378,   175,   415,   293,   293,
     186,   126,   242,   243,   244,  -214,   486,   298,   199,   105,
     170,   171,   199,   379,   379,   341,   406,   463,   217,   280,
     300,   126,   381,   178,   382,   217,   472,   126,   127,   146,
     460,   126,   144,   359,   146,   436,   294,   144,   362,   318,
     375,   373,   442,   128,  -212,   525,   459,   402,   127,   585,
     180,   493,   403,  -212,   127,  -212,  -214,   181,   127,   234,
     174,   591,   426,   479,   235,   126,   181,  -214,    47,   126,
     338,   232,   340,   280,   126,   232,   384,   289,   203,   128,
     628,   217,   453,   506,   507,   186,   548,   440,   458,   204,
     191,   341,   127,   396,   186,   428,   127,  -214,   494,   128,
     501,   127,   502,   179,   503,   128,   180,   164,   165,   128,
     422,   293,   180,   181,   403,   212,   180,   429,   495,   181,
     387,   232,   496,   181,   126,   280,   529,   511,   530,   146,
     199,   232,   144,   146,   146,   216,   144,   144,   459,   434,
     556,   231,   557,   128,   563,   341,   421,   128,   500,   518,
     341,   127,   128,   232,   457,   265,   266,   267,   268,   235,
     186,   445,   454,   196,   446,   359,   186,   143,   143,   508,
     232,   385,   150,   151,   152,   217,   465,   566,   567,   411,
     471,   411,   147,   517,   476,   261,   262,   521,   251,   592,
     522,   241,   565,   146,   295,   574,   144,   232,   579,   485,
     382,   297,   128,   492,   336,   150,   151,   152,   280,   180,
     337,   621,   387,   232,   347,   624,   280,   232,   627,   406,
     232,   568,   351,   570,   107,     3,     4,     5,     6,   263,
     264,   280,   352,   353,   251,   387,   354,   110,   112,   113,
     519,   206,   582,   207,   341,   355,   146,   269,   270,   144,
     356,   153,    11,   358,   154,   505,   155,   609,   186,   361,
     186,   383,   341,   219,   220,   221,   386,   399,   285,   516,
     414,   451,   542,   416,   417,   543,   418,   419,   423,   547,
     424,   430,   427,   553,   153,   433,   456,   154,   461,   155,
     357,   609,   559,    39,   209,   560,   157,   186,   462,   464,
     467,    47,   108,   498,   387,   468,   497,   487,   488,   111,
     504,   280,    18,    19,    20,   387,   509,   510,   523,   528,
     186,   531,   532,   534,   387,   606,   533,   538,   341,   590,
     610,   613,   341,   616,   586,   541,   544,   588,   549,   550,
     617,   569,   554,   571,   341,   279,   555,   558,   575,   572,
      43,   561,   577,   578,   597,   562,   573,   598,   576,   600,
     589,   157,   593,   595,   610,   157,   594,   596,   146,   181,
     620,   144,   114,   115,   116,   117,   118,   119,   120,   121,
     122,   123,   599,   601,   602,   618,   387,   619,   622,   157,
     157,   157,   157,   157,   157,   157,   157,   157,   157,   157,
     157,   157,   157,   157,   157,   157,   157,   630,   633,   515,
     252,   603,   604,   605,   587,   607,   188,   224,   608,   478,
     271,   273,   272,   390,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,  -213,   275,  -213,   274,   564,   540,   290,   435,
     350,   623,   441,   432,   626,   149,   202,    38,   625,   629,
     124,   631,    43,   632,     0,     0,   634,     0,     0,     0,
       0,   157,   635,     0,     0,     0,     0,   188,     0,     0,
       0,     0,     0,     0,     0,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,     0,     0,     0,     0,     0,     0,  -355,
     148,     0,     2,     3,     4,     5,     6,     7,    38,     8,
       9,    41,     0,    43,     0,     0,     0,   410,     0,     0,
       0,     0,   157,     0,     0,     0,     0,     0,     0,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,     0,     0,     0,     0,     0,
     157,     0,     0,     0,     0,     0,     0,   157,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,     0,    49,     0,     0,     0,     0,    50,     0,     0,
       0,     0,     0,     0,    51,    52,    53,    54,    55,    56,
       0,     0,     0,     0,   301,   157,     2,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,     0,     0,   157,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,     0,
     302,   303,   304,     0,   305,   306,   307,   308,   309,   310,
     311,   312,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,   313,    47,   314,   315,    49,     0,   316,   431,
       0,    50,     0,     0,     0,   105,     0,     0,    51,    52,
      53,    54,    55,    56,   301,     0,     2,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,   157,     0,     0,
     157,   410,     0,   410,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,     0,
     302,   303,   304,     0,   305,   306,   307,   308,   309,   310,
     311,   312,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,   313,    47,   314,   315,    49,     0,   316,     0,
       0,    50,     0,     0,     0,   105,     0,     0,    51,    52,
      53,    54,    55,    56,     1,     0,     2,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,     0,    49,     0,     0,     0,
       0,    50,     0,     0,     0,     0,     0,     0,    51,    52,
      53,    54,    55,    56,   139,     0,   107,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   140,    11,     0,     0,     0,     0,     0,
       0,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,   301,
       0,     2,     3,     4,     5,     6,     7,     0,     8,     9,
       0,     0,     0,    37,    38,    39,     0,     0,     0,    43,
      44,     0,     0,    47,   108,     0,     0,     0,   420,    11,
       0,    50,     0,     0,     0,     0,     0,     0,    51,    52,
      53,    54,    55,    56,   139,     0,   107,     3,     4,     5,
       6,     7,     0,     8,     9,   302,   303,   304,     0,   305,
     306,   307,   308,   309,   310,   311,   312,     0,    37,     0,
      39,     0,     0,   208,    11,    44,     0,   313,    47,   314,
     315,     0,     0,   316,     0,     0,    50,     0,     0,     0,
     105,     0,     0,    51,    52,    53,    54,    55,    56,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    37,     0,    39,     0,     0,     0,     0,
      44,     0,     0,    47,   108,   107,     3,     4,     5,     6,
       7,    50,     8,     9,     0,     0,     0,     0,    51,    52,
      53,    54,    55,    56,     0,     0,     0,     0,     0,     0,
       0,     0,   140,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,     0,    47,   108,     0,    49,     0,   316,     0,     0,
      50,     0,     0,     0,     0,     0,     0,    51,    52,    53,
      54,    55,    56,   107,     3,     4,     5,     6,     7,     0,
       8,     9,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     188,    11,     0,     0,     0,     0,     0,     0,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,   107,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
      37,    38,    39,     0,     0,     0,    43,    44,     0,     0,
      47,   108,     0,   208,    11,     0,     0,   370,    50,     0,
       0,    18,    19,    20,     0,    51,    52,    53,    54,    55,
      56,     0,     0,     0,     0,     0,     0,     0,   107,     3,
       4,     5,     6,     7,     0,     8,     9,     0,     0,     0,
       0,     0,     0,    37,   279,    39,     0,     0,     0,    43,
      44,     0,     0,    47,   108,   208,    11,     0,     0,   393,
       0,    50,     0,    18,    19,    20,     0,   371,   372,    52,
      53,    54,    55,    56,     0,     0,     0,     0,     0,     0,
     107,     3,     4,     5,     6,     7,     0,     8,     9,     0,
       0,     0,     0,     0,     0,    37,   279,    39,     0,     0,
       0,    43,    44,     0,     0,    47,   108,   208,    11,     0,
       0,   469,     0,    50,     0,    18,    19,    20,     0,   394,
     395,    52,    53,    54,    55,    56,     0,     0,     0,     0,
       0,     0,   107,     3,     4,     5,     6,     7,     0,     8,
       9,     0,     0,     0,     0,     0,     0,    37,   279,    39,
       0,     0,     0,    43,    44,     0,     0,    47,   108,   208,
      11,     0,     0,   473,     0,    50,     0,    18,    19,    20,
       0,   470,    51,    52,    53,    54,    55,    56,     0,     0,
       0,     0,     0,     0,   107,     3,     4,     5,     6,     7,
       0,     8,     9,     0,     0,     0,     0,     0,     0,    37,
     279,    39,     0,     0,     0,    43,    44,     0,     0,    47,
     108,   208,    11,     0,     0,   489,     0,    50,     0,    18,
      19,    20,     0,   474,   475,    52,    53,    54,    55,    56,
       0,     0,     0,     0,     0,     0,   107,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,    37,   279,    39,     0,     0,     0,    43,    44,     0,
       0,    47,   108,   208,    11,     0,     0,   551,     0,    50,
       0,    18,    19,    20,     0,   490,   491,    52,    53,    54,
      55,    56,     0,     0,     0,     0,     0,     0,   107,     3,
       4,     5,     6,     7,     0,     8,     9,     0,     0,     0,
       0,     0,     0,    37,   279,    39,     0,     0,     0,    43,
      44,     0,     0,    47,   108,   208,    11,     0,     0,     0,
       0,    50,     0,    18,    19,    20,     0,   552,    51,    52,
      53,    54,    55,    56,     0,     0,   107,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    37,   279,    39,     0,     0,
       0,    43,    44,   208,    11,    47,   108,   107,     3,     4,
       5,     6,     7,    50,     8,     9,     0,     0,     0,     0,
      51,    52,    53,    54,    55,    56,     0,     0,     0,     0,
       0,     0,     0,     0,   208,    11,     0,     0,     0,     0,
       0,     0,     0,    37,     0,    39,     0,     0,     0,     0,
      44,     0,     0,    47,   108,     0,     0,     0,     0,   524,
       0,    50,     0,     0,   338,   339,   340,     0,    51,    52,
      53,    54,    55,    56,    37,     0,    39,     0,     0,     0,
       0,    44,     0,     0,    47,   108,     0,     0,     0,     0,
     584,     0,    50,     0,     0,   338,   339,   340,     0,    51,
      52,    53,    54,    55,    56,   107,     3,     4,     5,     6,
       7,     0,     8,     9,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   208,    11,     0,     0,   107,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   208,    11,     0,     0,     0,     0,     0,
       0,     0,    37,     0,    39,     0,     0,     0,     0,    44,
       0,     0,    47,   108,   107,     3,     4,     5,     6,     7,
      50,     8,     9,   338,   339,   340,     0,    51,    52,    53,
      54,    55,    56,    37,     0,    39,     0,     0,     0,     0,
      44,   208,    11,    47,   108,     0,     0,     0,     0,     0,
     247,    50,     0,     0,     0,   209,     0,     0,    51,    52,
      53,    54,    55,    56,     0,     0,   107,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,    37,     0,    39,     0,     0,     0,     0,    44,     0,
       0,    47,   108,   208,    11,     0,     0,     0,   514,    50,
       0,     0,     0,   209,     0,     0,    51,    52,    53,    54,
      55,    56,     0,     0,   107,     3,     4,     5,     6,     7,
       0,     8,     9,     0,     0,     0,     0,     0,     0,     0,
     535,     0,     0,    37,     0,    39,     0,     0,     0,     0,
      44,   208,    11,    47,   108,   107,     3,     4,     5,     6,
       7,    50,     8,     9,     0,   339,     0,     0,    51,    52,
      53,    54,    55,    56,     0,     0,     0,     0,     0,     0,
       0,     0,   208,    11,     0,     0,     0,     0,     0,     0,
       0,    37,     0,    39,     0,     0,     0,     0,    44,     0,
       0,    47,   108,   107,     3,     4,     5,     6,     7,    50,
       8,     9,     0,   209,     0,     0,    51,    52,    53,    54,
      55,    56,    37,     0,    39,     0,     0,     0,     0,    44,
     208,    11,    47,   108,   107,     3,     4,     5,     6,     7,
      50,     8,     9,     0,     0,     0,   253,   254,    52,    53,
      54,    55,    56,     0,     0,     0,     0,     0,     0,     0,
       0,   208,    11,     0,     0,     0,     0,     0,     0,     0,
      37,     0,    39,     0,     0,     0,     0,    44,     0,     0,
      47,   108,   107,     3,     4,     5,     6,     7,    50,     8,
       9,     0,   339,     0,     0,    51,    52,    53,    54,    55,
      56,    37,     0,    39,     0,     0,     0,     0,    44,   208,
      11,    47,   108,     0,     0,     0,   425,     0,     0,    50,
       0,     0,     0,     0,     0,     0,    51,    52,    53,    54,
      55,    56,     0,     0,   107,     3,     4,     5,     6,     7,
       0,     8,     9,     0,     0,     0,     0,     0,     0,    37,
       0,    39,     0,     0,     0,     0,    44,     0,     0,    47,
     108,   208,    11,     0,   316,     0,     0,    50,     0,     0,
       0,     0,     0,     0,    51,    52,    53,    54,    55,    56,
       0,     0,   107,     3,     4,     5,     6,     7,     0,     8,
       9,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    37,     0,    39,     0,     0,     0,     0,    44,   208,
      11,    47,   108,     0,     0,     0,     0,     0,   612,    50,
       0,     0,     0,     0,     0,     0,    51,    52,    53,    54,
      55,    56,     0,     0,   107,     3,     4,     5,     6,     7,
       0,     8,     9,     0,     0,     0,     0,     0,     0,    37,
       0,    39,     0,     0,     0,     0,    44,     0,     0,    47,
     108,   208,    11,     0,     0,     0,   615,    50,     0,     0,
       0,     0,     0,     0,    51,    52,    53,    54,    55,    56,
       0,   107,     3,     4,     5,     6,     7,     0,     8,     9,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    37,     0,    39,     0,     0,     0,     0,    44,    11,
       0,    47,   108,   107,     3,     4,     5,     6,     7,    50,
       8,     9,     0,     0,     0,     0,    51,    52,    53,    54,
      55,    56,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    11,     0,     0,     0,     0,     0,     0,    37,     0,
      39,     0,     0,     0,     0,    44,     0,     0,    47,   108,
     107,     3,     4,     5,     6,     7,    50,     8,     9,     0,
     209,     0,     0,    51,    52,    53,    54,    55,    56,     0,
      37,     0,    39,     0,     0,     0,     0,    44,    11,     0,
      47,   108,   107,     3,     4,     5,     6,     7,   109,     8,
       9,     0,     0,     0,     0,    51,    52,    53,    54,    55,
      56,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      11,     0,     0,     0,     0,     0,     0,    37,     0,    39,
       0,     0,     0,     0,    44,     0,     0,    47,   108,     0,
       0,     0,     0,     0,     0,   111,     0,     0,     0,     0,
       0,     0,    51,    52,    53,    54,    55,    56,     0,    37,
       0,    39,     0,   126,     0,     0,    44,     0,     0,    47,
     108,     0,     0,     0,     0,     0,     0,    50,     0,     0,
       0,     0,     0,     0,    51,    52,    53,    54,    55,    56,
     537,     0,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   363,   480,    36,
     481,    38,     0,    40,     0,    42,    43,     0,    45,     0,
       0,   128,     0,    49,     0,     0,     0,   364,   457,     0,
       0,     0,     0,   235,     0,   181,     0,   188,     0,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,   407,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   363,     0,    36,     0,    38,     0,
      40,     0,    42,    43,     0,    45,     0,     0,     0,     0,
      49,     0,   188,     0,   482,     0,     0,     0,     0,     0,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    38,     0,     0,    41,   188,    43,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   363,     0,    36,     0,    38,     0,
      40,     0,    42,    43,     0,    45,     0,     0,     0,     0,
      49,     0,     0,     0,   364,   234,     0,     0,     0,     0,
     235,   188,   181,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      36,     0,    38,     0,    40,    41,    42,    43,     0,    45,
       0,     0,     0,     0,    49,     0,     0,     0,     0,     0,
       0,     0,   188,   105,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   363,
       0,    36,     0,    38,     0,    40,     0,    42,    43,     0,
      45,     0,   188,     0,     0,    49,     0,     0,     0,   545,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    38,     0,     0,    41,     0,    43,     0,
       0,   188,     0,     0,     0,     0,     0,     0,   400,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    38,     0,     0,    41,     0,    43,     0,     0,
     188,     0,     0,     0,     0,     0,     0,   409,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    38,     0,     0,    41,     0,    43,     0,     0,     0,
       0,     0,     0,     0,     0,   188,   499,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,   539,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   363,     0,    36,     0,    38,     0,    40,     0,
      42,    43,     0,    45,     0,     0,     0,   188,    49,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    36,     0,    38,     0,
      40,     0,    42,    43,   188,    45,     0,     0,     0,     0,
      49,     0,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   188,     0,    38,     0,     0,    41,     0,
      43,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,     0,
       0,   447,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    38,     0,     0,     0,     0,    43
};

static const yytype_int16 yycheck[] =
{
      50,     0,   125,    65,    81,   187,     0,     0,   285,   131,
     145,    95,   226,   135,    66,     3,   291,   236,    19,   191,
       1,   296,    14,    15,   196,    92,   145,   196,   179,   438,
     135,   349,    86,   156,   185,    89,    35,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,   144,   235,   146,
      88,     3,     3,    85,     4,    93,     3,    89,    57,   109,
      88,   111,   181,    57,    57,     3,    95,   344,    19,   433,
      85,     3,   101,   102,    89,     4,   108,     3,    30,    30,
      88,    88,    81,    30,   392,    93,   237,    86,    50,    82,
      83,    50,    88,    92,   191,    93,    89,    90,    91,   196,
     287,    86,    95,   105,    89,   155,   107,   159,   160,   161,
      88,   103,   104,   522,    95,   234,    88,   201,   241,   291,
     108,   205,   199,   295,   296,   108,   295,   177,   446,    81,
      81,   125,   283,   408,    81,   286,    88,    87,   289,   138,
     134,    93,    88,    95,   358,    92,    97,   109,    86,   111,
     109,   150,   111,    85,   518,   154,   520,    89,    87,   153,
      86,   438,   156,   282,   472,     3,     3,   148,    88,   131,
     132,   133,   131,   132,   133,   297,    97,    98,     3,   456,
     302,   180,   144,   370,   146,   144,    96,   146,   187,   241,
      88,   378,    30,    30,   291,    93,    18,   302,   295,   296,
     199,     3,    85,    86,    87,    30,   393,   284,   201,    92,
      16,    17,   205,   336,   337,   209,   293,   368,   340,   181,
     201,     3,    87,    85,    89,   347,   377,     3,    30,   191,
     365,     3,   191,   226,   196,   340,   408,   196,   232,   323,
     459,   235,   347,    81,    81,   522,   365,    85,    30,   526,
      88,   402,    90,    90,    30,    92,    81,    95,    30,    88,
     106,   538,   312,   385,    93,     3,    95,    92,    80,     3,
      91,    89,    93,   235,     3,    89,    94,   459,     0,    81,
      94,   403,   359,   417,   418,   284,   473,   108,   365,    86,
      92,   285,    30,   287,   293,    88,    30,    90,   403,    81,
      85,    30,    87,    85,    89,    81,    88,    12,    13,    81,
     309,   408,    88,    95,    90,     3,    88,    85,    85,    95,
     282,    89,    89,    95,     3,   287,    87,    85,    89,   291,
     323,    89,   291,   295,   296,    92,   295,   296,   457,   338,
      87,    87,    89,    81,   495,   339,   307,    81,   410,   433,
     344,    30,    81,    89,    88,   166,   167,   168,   169,    93,
     359,    86,   361,    92,    89,   358,   365,   417,   418,   419,
      89,    90,     9,    10,    11,   497,   370,   500,   500,   501,
     374,   503,    50,   433,   378,   162,   163,    86,   382,   540,
      89,    87,   497,   355,    92,    87,   355,    89,    87,   393,
      89,    88,    81,   397,    87,     9,    10,    11,   370,    88,
      87,    87,   374,    89,   108,    87,   378,    89,    87,   496,
      89,   501,    87,   503,     3,     4,     5,     6,     7,   164,
     165,   393,    87,    87,   428,   397,    87,     7,     8,     9,
     433,   109,   519,   111,   438,    89,   408,   170,   171,   408,
      89,    88,    31,    87,    91,   416,    93,   580,   457,    90,
     459,    94,   456,   131,   132,   133,    87,    86,   108,   430,
      85,     6,   466,    90,    88,   469,    88,    88,    85,   473,
      85,    90,    88,   477,    88,    88,   108,    91,    87,    93,
      94,   614,   486,    72,    92,   489,    66,   496,    89,    87,
      94,    80,    81,    86,   466,    94,    90,    94,    94,    88,
      90,   473,    38,    39,    40,   477,    62,    85,    94,    90,
     519,    90,    87,    87,   486,   575,    93,   108,   522,    90,
     580,   581,   526,   583,   528,    94,    94,   531,    94,    94,
     590,   502,    87,   504,   538,    71,    87,    94,    88,    87,
      76,    94,   513,   514,   548,    94,    87,   551,    87,     3,
      85,   131,    94,    87,   614,   135,    94,    94,   530,    95,
      60,   530,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    94,    94,    94,    94,   548,    94,    87,   159,
     160,   161,   162,   163,   164,   165,   166,   167,   168,   169,
     170,   171,   172,   173,   174,   175,   176,    85,    87,   428,
     154,   572,   573,   574,   530,   576,    30,   136,   579,   382,
     172,   174,   173,   284,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    51,    52,    53,
      54,    55,    90,   176,    92,   175,   496,   462,     1,   339,
     216,   612,   345,   323,   615,    57,    97,    71,   614,   620,
     108,   622,    76,   624,    -1,    -1,   627,    -1,    -1,    -1,
      -1,   241,   633,    -1,    -1,    -1,    -1,    30,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,     0,
       1,    -1,     3,     4,     5,     6,     7,     8,    71,    10,
      11,    74,    -1,    76,    -1,    -1,    -1,   297,    -1,    -1,
      -1,    -1,   302,    -1,    -1,    -1,    -1,    -1,    -1,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      51,    52,    53,    54,    55,    -1,    -1,    -1,    -1,    -1,
     340,    -1,    -1,    -1,    -1,    -1,    -1,   347,    69,    70,
      71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
      81,    -1,    83,    -1,    -1,    -1,    -1,    88,    -1,    -1,
      -1,    -1,    -1,    -1,    95,    96,    97,    98,    99,   100,
      -1,    -1,    -1,    -1,     1,   385,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   403,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    55,    -1,
      57,    58,    59,    -1,    61,    62,    63,    64,    65,    66,
      67,    68,    69,    70,    71,    72,    73,    74,    75,    76,
      77,    78,    79,    80,    81,    82,    83,    -1,    85,    86,
      -1,    88,    -1,    -1,    -1,    92,    -1,    -1,    95,    96,
      97,    98,    99,   100,     1,    -1,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,   497,    -1,    -1,
     500,   501,    -1,   503,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    55,    -1,
      57,    58,    59,    -1,    61,    62,    63,    64,    65,    66,
      67,    68,    69,    70,    71,    72,    73,    74,    75,    76,
      77,    78,    79,    80,    81,    82,    83,    -1,    85,    -1,
      -1,    88,    -1,    -1,    -1,    92,    -1,    -1,    95,    96,
      97,    98,    99,   100,     1,    -1,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    31,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    55,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    69,    70,    71,    72,    73,    74,    75,    76,
      77,    78,    79,    80,    81,    -1,    83,    -1,    -1,    -1,
      -1,    88,    -1,    -1,    -1,    -1,    -1,    -1,    95,    96,
      97,    98,    99,   100,     1,    -1,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    31,    -1,    -1,    -1,    -1,    -1,
      -1,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    55,     1,
      -1,     3,     4,     5,     6,     7,     8,    -1,    10,    11,
      -1,    -1,    -1,    70,    71,    72,    -1,    -1,    -1,    76,
      77,    -1,    -1,    80,    81,    -1,    -1,    -1,    30,    31,
      -1,    88,    -1,    -1,    -1,    -1,    -1,    -1,    95,    96,
      97,    98,    99,   100,     1,    -1,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    57,    58,    59,    -1,    61,
      62,    63,    64,    65,    66,    67,    68,    -1,    70,    -1,
      72,    -1,    -1,    30,    31,    77,    -1,    79,    80,    81,
      82,    -1,    -1,    85,    -1,    -1,    88,    -1,    -1,    -1,
      92,    -1,    -1,    95,    96,    97,    98,    99,   100,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    70,    -1,    72,    -1,    -1,    -1,    -1,
      77,    -1,    -1,    80,    81,     3,     4,     5,     6,     7,
       8,    88,    10,    11,    -1,    -1,    -1,    -1,    95,    96,
      97,    98,    99,   100,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    69,    70,    71,    72,    73,    74,    75,    76,    77,
      78,    -1,    80,    81,    -1,    83,    -1,    85,    -1,    -1,
      88,    -1,    -1,    -1,    -1,    -1,    -1,    95,    96,    97,
      98,    99,   100,     3,     4,     5,     6,     7,     8,    -1,
      10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      30,    31,    -1,    -1,    -1,    -1,    -1,    -1,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      70,    71,    72,    -1,    -1,    -1,    76,    77,    -1,    -1,
      80,    81,    -1,    30,    31,    -1,    -1,    34,    88,    -1,
      -1,    38,    39,    40,    -1,    95,    96,    97,    98,    99,
     100,    -1,    -1,    -1,    -1,    -1,    -1,    -1,     3,     4,
       5,     6,     7,     8,    -1,    10,    11,    -1,    -1,    -1,
      -1,    -1,    -1,    70,    71,    72,    -1,    -1,    -1,    76,
      77,    -1,    -1,    80,    81,    30,    31,    -1,    -1,    34,
      -1,    88,    -1,    38,    39,    40,    -1,    94,    95,    96,
      97,    98,    99,   100,    -1,    -1,    -1,    -1,    -1,    -1,
       3,     4,     5,     6,     7,     8,    -1,    10,    11,    -1,
      -1,    -1,    -1,    -1,    -1,    70,    71,    72,    -1,    -1,
      -1,    76,    77,    -1,    -1,    80,    81,    30,    31,    -1,
      -1,    34,    -1,    88,    -1,    38,    39,    40,    -1,    94,
      95,    96,    97,    98,    99,   100,    -1,    -1,    -1,    -1,
      -1,    -1,     3,     4,     5,     6,     7,     8,    -1,    10,
      11,    -1,    -1,    -1,    -1,    -1,    -1,    70,    71,    72,
      -1,    -1,    -1,    76,    77,    -1,    -1,    80,    81,    30,
      31,    -1,    -1,    34,    -1,    88,    -1,    38,    39,    40,
      -1,    94,    95,    96,    97,    98,    99,   100,    -1,    -1,
      -1,    -1,    -1,    -1,     3,     4,     5,     6,     7,     8,
      -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    70,
      71,    72,    -1,    -1,    -1,    76,    77,    -1,    -1,    80,
      81,    30,    31,    -1,    -1,    34,    -1,    88,    -1,    38,
      39,    40,    -1,    94,    95,    96,    97,    98,    99,   100,
      -1,    -1,    -1,    -1,    -1,    -1,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,    70,    71,    72,    -1,    -1,    -1,    76,    77,    -1,
      -1,    80,    81,    30,    31,    -1,    -1,    34,    -1,    88,
      -1,    38,    39,    40,    -1,    94,    95,    96,    97,    98,
      99,   100,    -1,    -1,    -1,    -1,    -1,    -1,     3,     4,
       5,     6,     7,     8,    -1,    10,    11,    -1,    -1,    -1,
      -1,    -1,    -1,    70,    71,    72,    -1,    -1,    -1,    76,
      77,    -1,    -1,    80,    81,    30,    31,    -1,    -1,    -1,
      -1,    88,    -1,    38,    39,    40,    -1,    94,    95,    96,
      97,    98,    99,   100,    -1,    -1,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    70,    71,    72,    -1,    -1,
      -1,    76,    77,    30,    31,    80,    81,     3,     4,     5,
       6,     7,     8,    88,    10,    11,    -1,    -1,    -1,    -1,
      95,    96,    97,    98,    99,   100,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    30,    31,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    70,    -1,    72,    -1,    -1,    -1,    -1,
      77,    -1,    -1,    80,    81,    -1,    -1,    -1,    -1,    86,
      -1,    88,    -1,    -1,    91,    92,    93,    -1,    95,    96,
      97,    98,    99,   100,    70,    -1,    72,    -1,    -1,    -1,
      -1,    77,    -1,    -1,    80,    81,    -1,    -1,    -1,    -1,
      86,    -1,    88,    -1,    -1,    91,    92,    93,    -1,    95,
      96,    97,    98,    99,   100,     3,     4,     5,     6,     7,
       8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    30,    31,    -1,    -1,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    31,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,
      -1,    -1,    80,    81,     3,     4,     5,     6,     7,     8,
      88,    10,    11,    91,    92,    93,    -1,    95,    96,    97,
      98,    99,   100,    70,    -1,    72,    -1,    -1,    -1,    -1,
      77,    30,    31,    80,    81,    -1,    -1,    -1,    -1,    -1,
      87,    88,    -1,    -1,    -1,    92,    -1,    -1,    95,    96,
      97,    98,    99,   100,    -1,    -1,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,    -1,
      -1,    80,    81,    30,    31,    -1,    -1,    -1,    87,    88,
      -1,    -1,    -1,    92,    -1,    -1,    95,    96,    97,    98,
      99,   100,    -1,    -1,     3,     4,     5,     6,     7,     8,
      -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      67,    -1,    -1,    70,    -1,    72,    -1,    -1,    -1,    -1,
      77,    30,    31,    80,    81,     3,     4,     5,     6,     7,
       8,    88,    10,    11,    -1,    92,    -1,    -1,    95,    96,
      97,    98,    99,   100,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    30,    31,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,    -1,
      -1,    80,    81,     3,     4,     5,     6,     7,     8,    88,
      10,    11,    -1,    92,    -1,    -1,    95,    96,    97,    98,
      99,   100,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,
      30,    31,    80,    81,     3,     4,     5,     6,     7,     8,
      88,    10,    11,    -1,    -1,    -1,    94,    95,    96,    97,
      98,    99,   100,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    30,    31,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      70,    -1,    72,    -1,    -1,    -1,    -1,    77,    -1,    -1,
      80,    81,     3,     4,     5,     6,     7,     8,    88,    10,
      11,    -1,    92,    -1,    -1,    95,    96,    97,    98,    99,
     100,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,    30,
      31,    80,    81,    -1,    -1,    -1,    85,    -1,    -1,    88,
      -1,    -1,    -1,    -1,    -1,    -1,    95,    96,    97,    98,
      99,   100,    -1,    -1,     3,     4,     5,     6,     7,     8,
      -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    70,
      -1,    72,    -1,    -1,    -1,    -1,    77,    -1,    -1,    80,
      81,    30,    31,    -1,    85,    -1,    -1,    88,    -1,    -1,
      -1,    -1,    -1,    -1,    95,    96,    97,    98,    99,   100,
      -1,    -1,     3,     4,     5,     6,     7,     8,    -1,    10,
      11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,    30,
      31,    80,    81,    -1,    -1,    -1,    -1,    -1,    87,    88,
      -1,    -1,    -1,    -1,    -1,    -1,    95,    96,    97,    98,
      99,   100,    -1,    -1,     3,     4,     5,     6,     7,     8,
      -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    70,
      -1,    72,    -1,    -1,    -1,    -1,    77,    -1,    -1,    80,
      81,    30,    31,    -1,    -1,    -1,    87,    88,    -1,    -1,
      -1,    -1,    -1,    -1,    95,    96,    97,    98,    99,   100,
      -1,     3,     4,     5,     6,     7,     8,    -1,    10,    11,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,    31,
      -1,    80,    81,     3,     4,     5,     6,     7,     8,    88,
      10,    11,    -1,    -1,    -1,    -1,    95,    96,    97,    98,
      99,   100,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    31,    -1,    -1,    -1,    -1,    -1,    -1,    70,    -1,
      72,    -1,    -1,    -1,    -1,    77,    -1,    -1,    80,    81,
       3,     4,     5,     6,     7,     8,    88,    10,    11,    -1,
      92,    -1,    -1,    95,    96,    97,    98,    99,   100,    -1,
      70,    -1,    72,    -1,    -1,    -1,    -1,    77,    31,    -1,
      80,    81,     3,     4,     5,     6,     7,     8,    88,    10,
      11,    -1,    -1,    -1,    -1,    95,    96,    97,    98,    99,
     100,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      31,    -1,    -1,    -1,    -1,    -1,    -1,    70,    -1,    72,
      -1,    -1,    -1,    -1,    77,    -1,    -1,    80,    81,    -1,
      -1,    -1,    -1,    -1,    -1,    88,    -1,    -1,    -1,    -1,
      -1,    -1,    95,    96,    97,    98,    99,   100,    -1,    70,
      -1,    72,    -1,     3,    -1,    -1,    77,    -1,    -1,    80,
      81,    -1,    -1,    -1,    -1,    -1,    -1,    88,    -1,    -1,
      -1,    -1,    -1,    -1,    95,    96,    97,    98,    99,   100,
      30,    -1,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    67,     1,    69,
       3,    71,    -1,    73,    -1,    75,    76,    -1,    78,    -1,
      -1,    81,    -1,    83,    -1,    -1,    -1,    87,    88,    -1,
      -1,    -1,    -1,    93,    -1,    95,    -1,    30,    -1,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,     1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    67,    -1,    69,    -1,    71,    -1,
      73,    -1,    75,    76,    -1,    78,    -1,    -1,    -1,    -1,
      83,    -1,    30,    -1,    87,    -1,    -1,    -1,    -1,    -1,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    71,    -1,    -1,    74,    30,    76,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    67,    -1,    69,    -1,    71,    -1,
      73,    -1,    75,    76,    -1,    78,    -1,    -1,    -1,    -1,
      83,    -1,    -1,    -1,    87,    88,    -1,    -1,    -1,    -1,
      93,    30,    95,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,    55,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      69,    -1,    71,    -1,    73,    74,    75,    76,    -1,    78,
      -1,    -1,    -1,    -1,    83,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    30,    92,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    67,
      -1,    69,    -1,    71,    -1,    73,    -1,    75,    76,    -1,
      78,    -1,    30,    -1,    -1,    83,    -1,    -1,    -1,    87,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    71,    -1,    -1,    74,    -1,    76,    -1,
      -1,    30,    -1,    -1,    -1,    -1,    -1,    -1,    86,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,    55,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    71,    -1,    -1,    74,    -1,    76,    -1,    -1,
      30,    -1,    -1,    -1,    -1,    -1,    -1,    86,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    71,    -1,    -1,    74,    -1,    76,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    30,    86,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    67,    -1,    69,    -1,    71,    -1,    73,    -1,
      75,    76,    -1,    78,    -1,    -1,    -1,    30,    83,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    69,    -1,    71,    -1,
      73,    -1,    75,    76,    30,    78,    -1,    -1,    -1,    -1,
      83,    -1,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    -1,    71,    -1,    -1,    74,    -1,
      76,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    55,    -1,
      -1,    58,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    71,    -1,    -1,    -1,    -1,    76
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     5,     6,     7,     8,    10,    11,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    69,    70,    71,    72,
      73,    74,    75,    76,    77,    78,    79,    80,    81,    83,
      88,    95,    96,    97,    98,    99,   100,   110,   111,   113,
     115,   116,   119,   121,   122,   126,   127,   128,   129,   130,
     131,   132,   133,   134,   135,   136,   137,   138,   139,   140,
     144,   145,   148,   149,   150,   151,   152,   158,   161,   162,
     163,   164,   167,   183,   196,   197,   198,   202,   203,   210,
     213,   214,   215,   216,   217,    92,   186,     3,    81,    88,
     126,    88,   126,   126,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,   108,   141,     3,    30,    81,    92,
     167,    88,    88,    88,    88,    88,    93,    88,    88,     1,
      30,   112,   140,   142,   149,   155,   162,   174,     1,   196,
       9,    10,    11,    88,    91,    93,   141,   126,   128,    95,
     101,   102,    97,    98,    12,    13,    14,    15,   103,   104,
      16,    17,    96,   105,   106,    18,    19,   107,    85,    85,
      88,    95,   146,   147,   165,   166,   167,   168,    30,   145,
     145,    92,   167,   145,   145,   145,    92,   167,   144,   145,
     186,   199,   203,     0,    86,   187,   174,   174,    30,    92,
     120,   140,     3,   114,   159,   160,    92,   139,   143,   174,
     174,   174,   140,   143,   119,     4,    87,    19,    97,   167,
     204,    87,    89,   155,    88,    93,   168,   175,   176,   177,
     155,    87,    85,    86,    87,   123,   167,    87,   120,   124,
     125,   140,   123,    94,    95,   142,   120,   140,   128,   128,
     128,   129,   129,   130,   130,   131,   131,   131,   131,   132,
     132,   133,   134,   135,   136,   137,   142,   177,   165,    71,
     162,   168,   169,    85,    89,   108,    88,    93,   177,   166,
       1,   153,   154,   155,   183,    92,   153,    88,   165,   144,
     186,     1,    57,    58,    59,    61,    62,    63,    64,    65,
      66,    67,    68,    79,    81,    82,    85,   142,   144,   167,
     184,   185,   186,   188,   189,   190,   191,   192,   193,   194,
     195,   200,   201,   206,   207,   208,    87,    87,    91,    92,
      93,   140,   178,   179,   180,   181,   182,   108,    86,    89,
     159,    87,    87,    87,    87,    89,    89,    94,    87,   145,
     198,    90,   140,    67,    87,   145,   170,   171,   172,   175,
      34,    94,    95,   140,   169,   176,   177,    88,    93,   120,
     128,    87,    89,    94,    94,    90,    87,   162,   168,   177,
     147,   178,   177,    34,    94,    95,   140,   169,   177,    86,
      86,   154,    85,    90,   156,   157,   165,     1,   153,    86,
     126,   139,   211,   212,    85,   143,    90,    88,    88,    88,
      30,   184,   167,    85,    85,    85,   142,    88,    88,    85,
      90,    86,   189,    88,   167,   179,   143,    86,    89,   178,
     108,   182,   143,    86,   160,    86,    89,    58,   117,   118,
     174,     6,   198,   165,   167,   205,   108,    88,   165,   168,
     175,    87,    89,   177,    87,   140,   169,    94,    94,    34,
      94,   140,   177,    34,    94,    95,   140,   169,   125,   139,
       1,     3,    87,   170,   173,   140,   169,    94,    94,    34,
      94,    95,   140,   177,   143,    85,    89,    90,    86,    86,
     141,    85,    87,    89,    90,   184,   112,   112,   142,    62,
      85,    85,     4,    87,    87,   124,   184,   142,   144,   145,
     190,    86,    89,    94,    86,   178,   180,    86,    90,    87,
      89,    90,    87,    93,    87,    67,   178,    30,   108,    56,
     172,    94,   140,   140,    94,    87,   170,   140,   169,    94,
      94,    34,    94,   140,    87,    87,    87,    89,    94,   140,
     140,    94,    94,   177,   157,   143,   120,   139,   212,   184,
     212,   184,    87,    87,    87,    88,    87,   184,   184,    87,
       3,   190,   165,   190,    86,   178,   140,   118,   140,    85,
      90,   178,   177,    94,    94,    87,    94,   140,   140,    94,
       3,    94,    94,   184,   184,   184,   142,   184,   184,   120,
     142,   209,    87,   142,     3,    87,   142,   142,    94,    94,
      60,    87,    87,   184,    87,   209,   184,    87,    94,   184,
      85,   184,   184,    87,   184,   184
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   109,   110,   110,   110,   110,   110,   111,   111,   111,
     111,   111,   111,   111,   112,   112,   113,   113,   113,   114,
     115,   115,   116,   117,   117,   118,   118,   119,   119,   119,
     119,   119,   119,   119,   119,   120,   120,   121,   121,   121,
     122,   122,   123,   124,   124,   124,   125,   125,   126,   126,
     126,   126,   126,   126,   126,   126,   127,   127,   127,   127,
     127,   127,   128,   128,   129,   129,   129,   129,   130,   130,
     130,   131,   131,   131,   132,   132,   132,   132,   132,   133,
     133,   133,   134,   134,   135,   135,   136,   136,   137,   137,
     138,   138,   139,   139,   140,   140,   140,   140,   140,   141,
     141,   141,   141,   141,   141,   141,   141,   141,   141,   141,
     142,   142,   143,   144,   144,   144,   145,   145,   145,   145,
     145,   145,   145,   145,   145,   145,   146,   146,   147,   147,
     148,   148,   148,   148,   148,   148,   148,   148,   149,   150,
     150,   150,   150,   150,   150,   150,   150,   150,   150,   150,
     150,   150,   150,   150,   150,   151,   151,   151,   151,   151,
     152,   152,   153,   153,   154,   154,   154,   155,   155,   155,
     155,   156,   156,   157,   157,   157,   158,   158,   158,   158,
     158,   159,   159,   160,   160,   161,   162,   162,   162,   162,
     162,   163,   163,   164,   164,   165,   165,   166,   166,   166,
     166,   166,   166,   166,   166,   166,   166,   166,   166,   166,
     166,   166,   167,   167,   167,   168,   168,   168,   168,   169,
     169,   170,   170,   171,   171,   172,   172,   172,   172,   172,
     172,   173,   173,   174,   174,   175,   175,   175,   176,   176,
     176,   176,   176,   176,   176,   176,   176,   176,   176,   176,
     176,   176,   176,   176,   176,   176,   176,   176,   176,   177,
     178,   178,   178,   179,   179,   179,   179,   180,   181,   181,
     182,   182,   183,   184,   184,   184,   184,   184,   184,   184,
     184,   185,   185,   185,   186,   187,   186,   188,   188,   189,
     189,   190,   190,   191,   191,   191,   192,   193,   193,   193,
     193,   193,   194,   194,   195,   195,   195,   195,   195,   196,
     196,   196,   196,   196,   196,   196,   197,   198,   198,   199,
     199,   200,   200,   200,   200,   200,   201,   201,   202,   202,
     203,   204,   204,   204,   205,   205,   206,   206,   207,   208,
     209,   209,   210,   211,   211,   211,   212,   212,   212,   213,
     214,   214,   215,   216,   216,   217
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     1,     2,     3,     3,     3,     1,     1,     1,
       1,     3,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     6,     1,     3,     3,     3,     1,     1,     1,
       3,     3,     2,     2,     4,     3,     4,     3,     4,     4,
       3,     4,     1,     1,     3,     2,     1,     1,     1,     2,
       2,     2,     2,     4,     4,     1,     1,     1,     1,     1,
       1,     1,     1,     4,     1,     3,     3,     3,     1,     3,
       3,     1,     3,     3,     1,     3,     3,     3,     3,     1,
       3,     3,     1,     3,     1,     3,     1,     3,     1,     3,
       1,     3,     1,     5,     1,     3,     3,     3,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     3,     1,     3,     4,     1,     2,     1,     2,     1,
       2,     1,     2,     1,     2,     1,     1,     3,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     4,     5,     2,     4,     5,
       1,     1,     1,     2,     3,     4,     1,     2,     1,     2,
       1,     1,     3,     2,     3,     1,     4,     5,     5,     6,
       2,     1,     3,     3,     1,     4,     1,     1,     1,     1,
       1,     1,     1,     4,     4,     3,     2,     1,     3,     3,
       4,     6,     5,     5,     6,     5,     4,     4,     5,     5,
       4,     5,     1,     1,     1,     3,     2,     2,     1,     1,
       2,     3,     1,     2,     4,     2,     4,     2,     1,     3,
       3,     1,     3,     3,     2,     2,     1,     1,     3,     2,
       3,     5,     4,     5,     4,     3,     3,     3,     4,     6,
       5,     5,     6,     4,     4,     2,     3,     4,     5,     0,
       3,     4,     1,     2,     1,     4,     3,     2,     1,     2,
       3,     2,     7,     1,     1,     1,     1,     1,     1,     1,
       2,     3,     4,     3,     2,     0,     4,     1,     2,     1,
       1,     1,     2,     7,     5,     5,     1,     5,     7,     6,
       7,     1,     6,     7,     3,     2,     2,     2,     3,     1,
       1,     1,     1,     1,     1,     2,     2,     3,     2,     1,
       2,     1,     1,     1,     1,     1,     4,     5,     1,     2,
       6,     1,     1,     1,     1,     5,     4,     5,     8,     7,
       1,     1,     6,     1,     3,     3,     1,     3,     3,     2,
       4,     5,     4,     1,     4,     1
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
  switch (yykind)
    {
    case YYSYMBOL_IDENTIFIER: /* IDENTIFIER  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 1961 "basilisk.c"
        break;

    case YYSYMBOL_I_CONSTANT: /* I_CONSTANT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 1970 "basilisk.c"
        break;

    case YYSYMBOL_F_CONSTANT: /* F_CONSTANT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 1979 "basilisk.c"
        break;

    case YYSYMBOL_STRING_LITERAL: /* STRING_LITERAL  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 1988 "basilisk.c"
        break;

    case YYSYMBOL_FUNC_NAME: /* FUNC_NAME  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 1997 "basilisk.c"
        break;

    case YYSYMBOL_SIZEOF: /* SIZEOF  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2006 "basilisk.c"
        break;

    case YYSYMBOL_PTR_OP: /* PTR_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2015 "basilisk.c"
        break;

    case YYSYMBOL_INC_OP: /* INC_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2024 "basilisk.c"
        break;

    case YYSYMBOL_DEC_OP: /* DEC_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2033 "basilisk.c"
        break;

    case YYSYMBOL_LEFT_OP: /* LEFT_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2042 "basilisk.c"
        break;

    case YYSYMBOL_RIGHT_OP: /* RIGHT_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2051 "basilisk.c"
        break;

    case YYSYMBOL_LE_OP: /* LE_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2060 "basilisk.c"
        break;

    case YYSYMBOL_GE_OP: /* GE_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2069 "basilisk.c"
        break;

    case YYSYMBOL_EQ_OP: /* EQ_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2078 "basilisk.c"
        break;

    case YYSYMBOL_NE_OP: /* NE_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2087 "basilisk.c"
        break;

    case YYSYMBOL_AND_OP: /* AND_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2096 "basilisk.c"
        break;

    case YYSYMBOL_OR_OP: /* OR_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2105 "basilisk.c"
        break;

    case YYSYMBOL_MUL_ASSIGN: /* MUL_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2114 "basilisk.c"
        break;

    case YYSYMBOL_DIV_ASSIGN: /* DIV_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2123 "basilisk.c"
        break;

    case YYSYMBOL_MOD_ASSIGN: /* MOD_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2132 "basilisk.c"
        break;

    case YYSYMBOL_ADD_ASSIGN: /* ADD_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2141 "basilisk.c"
        break;

    case YYSYMBOL_SUB_ASSIGN: /* SUB_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2150 "basilisk.c"
        break;

    case YYSYMBOL_LEFT_ASSIGN: /* LEFT_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2159 "basilisk.c"
        break;

    case YYSYMBOL_RIGHT_ASSIGN: /* RIGHT_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2168 "basilisk.c"
        break;

    case YYSYMBOL_AND_ASSIGN: /* AND_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2177 "basilisk.c"
        break;

    case YYSYMBOL_XOR_ASSIGN: /* XOR_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2186 "basilisk.c"
        break;

    case YYSYMBOL_OR_ASSIGN: /* OR_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2195 "basilisk.c"
        break;

    case YYSYMBOL_TYPEDEF_NAME: /* TYPEDEF_NAME  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2204 "basilisk.c"
        break;

    case YYSYMBOL_ENUMERATION_CONSTANT: /* ENUMERATION_CONSTANT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2213 "basilisk.c"
        break;

    case YYSYMBOL_TYPEDEF: /* TYPEDEF  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2222 "basilisk.c"
        break;

    case YYSYMBOL_EXTERN: /* EXTERN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2231 "basilisk.c"
        break;

    case YYSYMBOL_STATIC: /* STATIC  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2240 "basilisk.c"
        break;

    case YYSYMBOL_AUTO: /* AUTO  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2249 "basilisk.c"
        break;

    case YYSYMBOL_REGISTER: /* REGISTER  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2258 "basilisk.c"
        break;

    case YYSYMBOL_INLINE: /* INLINE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2267 "basilisk.c"
        break;

    case YYSYMBOL_CONST: /* CONST  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2276 "basilisk.c"
        break;

    case YYSYMBOL_RESTRICT: /* RESTRICT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2285 "basilisk.c"
        break;

    case YYSYMBOL_VOLATILE: /* VOLATILE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2294 "basilisk.c"
        break;

    case YYSYMBOL_BOOL: /* BOOL  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2303 "basilisk.c"
        break;

    case YYSYMBOL_CHAR: /* CHAR  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2312 "basilisk.c"
        break;

    case YYSYMBOL_SHORT: /* SHORT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2321 "basilisk.c"
        break;

    case YYSYMBOL_INT: /* INT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2330 "basilisk.c"
        break;

    case YYSYMBOL_LONG: /* LONG  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2339 "basilisk.c"
        break;

    case YYSYMBOL_SIGNED: /* SIGNED  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2348 "basilisk.c"
        break;

    case YYSYMBOL_UNSIGNED: /* UNSIGNED  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2357 "basilisk.c"
        break;

    case YYSYMBOL_FLOAT: /* FLOAT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2366 "basilisk.c"
        break;

    case YYSYMBOL_DOUBLE: /* DOUBLE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2375 "basilisk.c"
        break;

    case YYSYMBOL_VOID: /* VOID  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2384 "basilisk.c"
        break;

    case YYSYMBOL_COMPLEX: /* COMPLEX  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2393 "basilisk.c"
        break;

    case YYSYMBOL_IMAGINARY: /* IMAGINARY  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2402 "basilisk.c"
        break;

    case YYSYMBOL_STRUCT: /* STRUCT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2411 "basilisk.c"
        break;

    case YYSYMBOL_UNION: /* UNION  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2420 "basilisk.c"
        break;

    case YYSYMBOL_ENUM: /* ENUM  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2429 "basilisk.c"
        break;

    case YYSYMBOL_ELLIPSIS: /* ELLIPSIS  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2438 "basilisk.c"
        break;

    case YYSYMBOL_CASE: /* CASE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2447 "basilisk.c"
        break;

    case YYSYMBOL_DEFAULT: /* DEFAULT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2456 "basilisk.c"
        break;

    case YYSYMBOL_IF: /* IF  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2465 "basilisk.c"
        break;

    case YYSYMBOL_ELSE: /* ELSE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2474 "basilisk.c"
        break;

    case YYSYMBOL_SWITCH: /* SWITCH  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2483 "basilisk.c"
        break;

    case YYSYMBOL_WHILE: /* WHILE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2492 "basilisk.c"
        break;

    case YYSYMBOL_DO: /* DO  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2501 "basilisk.c"
        break;

    case YYSYMBOL_FOR: /* FOR  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2510 "basilisk.c"
        break;

    case YYSYMBOL_GOTO: /* GOTO  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2519 "basilisk.c"
        break;

    case YYSYMBOL_CONTINUE: /* CONTINUE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2528 "basilisk.c"
        break;

    case YYSYMBOL_BREAK: /* BREAK  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2537 "basilisk.c"
        break;

    case YYSYMBOL_RETURN: /* RETURN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2546 "basilisk.c"
        break;

    case YYSYMBOL_ALIGNAS: /* ALIGNAS  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2555 "basilisk.c"
        break;

    case YYSYMBOL_ALIGNOF: /* ALIGNOF  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2564 "basilisk.c"
        break;

    case YYSYMBOL_ATOMIC: /* ATOMIC  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2573 "basilisk.c"
        break;

    case YYSYMBOL_GENERIC: /* GENERIC  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2582 "basilisk.c"
        break;

    case YYSYMBOL_NORETURN: /* NORETURN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2591 "basilisk.c"
        break;

    case YYSYMBOL_STATIC_ASSERT: /* STATIC_ASSERT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2600 "basilisk.c"
        break;

    case YYSYMBOL_THREAD_LOCAL: /* THREAD_LOCAL  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2609 "basilisk.c"
        break;

    case YYSYMBOL_MAYBECONST: /* MAYBECONST  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2618 "basilisk.c"
        break;

    case YYSYMBOL_NEW_FIELD: /* NEW_FIELD  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2627 "basilisk.c"
        break;

    case YYSYMBOL_TRACE: /* TRACE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2636 "basilisk.c"
        break;

    case YYSYMBOL_FOREACH_DIMENSION: /* FOREACH_DIMENSION  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2645 "basilisk.c"
        break;

    case YYSYMBOL_REDUCTION: /* REDUCTION  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2654 "basilisk.c"
        break;

    case YYSYMBOL_MACRO: /* MACRO  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2663 "basilisk.c"
        break;

    case YYSYMBOL_ELLIPSIS_MACRO: /* ELLIPSIS_MACRO  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2672 "basilisk.c"
        break;

    case YYSYMBOL_MACRODEF: /* MACRODEF  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2681 "basilisk.c"
        break;

    case YYSYMBOL_foreach_statement: /* foreach_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2690 "basilisk.c"
        break;

    case YYSYMBOL_85_: /* ';'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2699 "basilisk.c"
        break;

    case YYSYMBOL_86_: /* '}'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2708 "basilisk.c"
        break;

    case YYSYMBOL_87_: /* ')'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2717 "basilisk.c"
        break;

    case YYSYMBOL_88_: /* '('  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2726 "basilisk.c"
        break;

    case YYSYMBOL_89_: /* ','  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2735 "basilisk.c"
        break;

    case YYSYMBOL_90_: /* ':'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2744 "basilisk.c"
        break;

    case YYSYMBOL_91_: /* '.'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2753 "basilisk.c"
        break;

    case YYSYMBOL_92_: /* '{'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2762 "basilisk.c"
        break;

    case YYSYMBOL_93_: /* '['  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2771 "basilisk.c"
        break;

    case YYSYMBOL_94_: /* ']'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2780 "basilisk.c"
        break;

    case YYSYMBOL_95_: /* '*'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2789 "basilisk.c"
        break;

    case YYSYMBOL_96_: /* '&'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2798 "basilisk.c"
        break;

    case YYSYMBOL_97_: /* '+'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2807 "basilisk.c"
        break;

    case YYSYMBOL_98_: /* '-'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2816 "basilisk.c"
        break;

    case YYSYMBOL_99_: /* '~'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2825 "basilisk.c"
        break;

    case YYSYMBOL_100_: /* '!'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2834 "basilisk.c"
        break;

    case YYSYMBOL_101_: /* '/'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2843 "basilisk.c"
        break;

    case YYSYMBOL_102_: /* '%'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2852 "basilisk.c"
        break;

    case YYSYMBOL_103_: /* '<'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2861 "basilisk.c"
        break;

    case YYSYMBOL_104_: /* '>'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2870 "basilisk.c"
        break;

    case YYSYMBOL_105_: /* '^'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2879 "basilisk.c"
        break;

    case YYSYMBOL_106_: /* '|'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2888 "basilisk.c"
        break;

    case YYSYMBOL_107_: /* '?'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2897 "basilisk.c"
        break;

    case YYSYMBOL_108_: /* '='  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2906 "basilisk.c"
        break;

    case YYSYMBOL_translation_unit: /* translation_unit  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2915 "basilisk.c"
        break;

    case YYSYMBOL_primary_expression: /* primary_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2924 "basilisk.c"
        break;

    case YYSYMBOL_expression_error: /* expression_error  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2933 "basilisk.c"
        break;

    case YYSYMBOL_constant: /* constant  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2942 "basilisk.c"
        break;

    case YYSYMBOL_enumeration_constant: /* enumeration_constant  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2951 "basilisk.c"
        break;

    case YYSYMBOL_string: /* string  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2960 "basilisk.c"
        break;

    case YYSYMBOL_generic_selection: /* generic_selection  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2969 "basilisk.c"
        break;

    case YYSYMBOL_generic_assoc_list: /* generic_assoc_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2978 "basilisk.c"
        break;

    case YYSYMBOL_generic_association: /* generic_association  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2987 "basilisk.c"
        break;

    case YYSYMBOL_postfix_expression: /* postfix_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2996 "basilisk.c"
        break;

    case YYSYMBOL_postfix_initializer: /* postfix_initializer  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3005 "basilisk.c"
        break;

    case YYSYMBOL_array_access: /* array_access  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3014 "basilisk.c"
        break;

    case YYSYMBOL_function_call: /* function_call  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3023 "basilisk.c"
        break;

    case YYSYMBOL_member_identifier: /* member_identifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3032 "basilisk.c"
        break;

    case YYSYMBOL_argument_expression_list: /* argument_expression_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3041 "basilisk.c"
        break;

    case YYSYMBOL_argument_expression_list_item: /* argument_expression_list_item  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3050 "basilisk.c"
        break;

    case YYSYMBOL_unary_expression: /* unary_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3059 "basilisk.c"
        break;

    case YYSYMBOL_unary_operator: /* unary_operator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3068 "basilisk.c"
        break;

    case YYSYMBOL_cast_expression: /* cast_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3077 "basilisk.c"
        break;

    case YYSYMBOL_multiplicative_expression: /* multiplicative_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3086 "basilisk.c"
        break;

    case YYSYMBOL_additive_expression: /* additive_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3095 "basilisk.c"
        break;

    case YYSYMBOL_shift_expression: /* shift_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3104 "basilisk.c"
        break;

    case YYSYMBOL_relational_expression: /* relational_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3113 "basilisk.c"
        break;

    case YYSYMBOL_equality_expression: /* equality_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3122 "basilisk.c"
        break;

    case YYSYMBOL_and_expression: /* and_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3131 "basilisk.c"
        break;

    case YYSYMBOL_exclusive_or_expression: /* exclusive_or_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3140 "basilisk.c"
        break;

    case YYSYMBOL_inclusive_or_expression: /* inclusive_or_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3149 "basilisk.c"
        break;

    case YYSYMBOL_logical_and_expression: /* logical_and_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3158 "basilisk.c"
        break;

    case YYSYMBOL_logical_or_expression: /* logical_or_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3167 "basilisk.c"
        break;

    case YYSYMBOL_conditional_expression: /* conditional_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3176 "basilisk.c"
        break;

    case YYSYMBOL_assignment_expression: /* assignment_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3185 "basilisk.c"
        break;

    case YYSYMBOL_assignment_operator: /* assignment_operator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3194 "basilisk.c"
        break;

    case YYSYMBOL_expression: /* expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3203 "basilisk.c"
        break;

    case YYSYMBOL_constant_expression: /* constant_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3212 "basilisk.c"
        break;

    case YYSYMBOL_declaration: /* declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3221 "basilisk.c"
        break;

    case YYSYMBOL_declaration_specifiers: /* declaration_specifiers  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3230 "basilisk.c"
        break;

    case YYSYMBOL_init_declarator_list: /* init_declarator_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3239 "basilisk.c"
        break;

    case YYSYMBOL_init_declarator: /* init_declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3248 "basilisk.c"
        break;

    case YYSYMBOL_storage_class_specifier: /* storage_class_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3257 "basilisk.c"
        break;

    case YYSYMBOL_type_specifier: /* type_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3266 "basilisk.c"
        break;

    case YYSYMBOL_types: /* types  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3275 "basilisk.c"
        break;

    case YYSYMBOL_struct_or_union_specifier: /* struct_or_union_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3284 "basilisk.c"
        break;

    case YYSYMBOL_struct_or_union: /* struct_or_union  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3293 "basilisk.c"
        break;

    case YYSYMBOL_struct_declaration_list: /* struct_declaration_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3302 "basilisk.c"
        break;

    case YYSYMBOL_struct_declaration: /* struct_declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3311 "basilisk.c"
        break;

    case YYSYMBOL_specifier_qualifier_list: /* specifier_qualifier_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3320 "basilisk.c"
        break;

    case YYSYMBOL_struct_declarator_list: /* struct_declarator_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3329 "basilisk.c"
        break;

    case YYSYMBOL_struct_declarator: /* struct_declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3338 "basilisk.c"
        break;

    case YYSYMBOL_enum_specifier: /* enum_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3347 "basilisk.c"
        break;

    case YYSYMBOL_enumerator_list: /* enumerator_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3356 "basilisk.c"
        break;

    case YYSYMBOL_enumerator: /* enumerator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3365 "basilisk.c"
        break;

    case YYSYMBOL_atomic_type_specifier: /* atomic_type_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3374 "basilisk.c"
        break;

    case YYSYMBOL_type_qualifier: /* type_qualifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3383 "basilisk.c"
        break;

    case YYSYMBOL_function_specifier: /* function_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3392 "basilisk.c"
        break;

    case YYSYMBOL_alignment_specifier: /* alignment_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3401 "basilisk.c"
        break;

    case YYSYMBOL_declarator: /* declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3410 "basilisk.c"
        break;

    case YYSYMBOL_direct_declarator: /* direct_declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3419 "basilisk.c"
        break;

    case YYSYMBOL_generic_identifier: /* generic_identifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3428 "basilisk.c"
        break;

    case YYSYMBOL_pointer: /* pointer  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3437 "basilisk.c"
        break;

    case YYSYMBOL_type_qualifier_list: /* type_qualifier_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3446 "basilisk.c"
        break;

    case YYSYMBOL_parameter_type_list: /* parameter_type_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3455 "basilisk.c"
        break;

    case YYSYMBOL_parameter_list: /* parameter_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3464 "basilisk.c"
        break;

    case YYSYMBOL_parameter_declaration: /* parameter_declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3473 "basilisk.c"
        break;

    case YYSYMBOL_identifier_list: /* identifier_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3482 "basilisk.c"
        break;

    case YYSYMBOL_type_name: /* type_name  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3491 "basilisk.c"
        break;

    case YYSYMBOL_abstract_declarator: /* abstract_declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3500 "basilisk.c"
        break;

    case YYSYMBOL_direct_abstract_declarator: /* direct_abstract_declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3509 "basilisk.c"
        break;

    case YYSYMBOL_type_not_specified: /* type_not_specified  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3518 "basilisk.c"
        break;

    case YYSYMBOL_initializer: /* initializer  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3527 "basilisk.c"
        break;

    case YYSYMBOL_initializer_list: /* initializer_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3536 "basilisk.c"
        break;

    case YYSYMBOL_designation: /* designation  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3545 "basilisk.c"
        break;

    case YYSYMBOL_designator_list: /* designator_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3554 "basilisk.c"
        break;

    case YYSYMBOL_designator: /* designator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3563 "basilisk.c"
        break;

    case YYSYMBOL_static_assert_declaration: /* static_assert_declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3572 "basilisk.c"
        break;

    case YYSYMBOL_statement: /* statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3581 "basilisk.c"
        break;

    case YYSYMBOL_labeled_statement: /* labeled_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3590 "basilisk.c"
        break;

    case YYSYMBOL_compound_statement: /* compound_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3599 "basilisk.c"
        break;

    case YYSYMBOL_187_1: /* @1  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3608 "basilisk.c"
        break;

    case YYSYMBOL_block_item_list: /* block_item_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3617 "basilisk.c"
        break;

    case YYSYMBOL_block_item: /* block_item  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3626 "basilisk.c"
        break;

    case YYSYMBOL_expression_statement: /* expression_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3635 "basilisk.c"
        break;

    case YYSYMBOL_selection_statement: /* selection_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3644 "basilisk.c"
        break;

    case YYSYMBOL_for_scope: /* for_scope  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3653 "basilisk.c"
        break;

    case YYSYMBOL_iteration_statement: /* iteration_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3662 "basilisk.c"
        break;

    case YYSYMBOL_for_declaration_statement: /* for_declaration_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3671 "basilisk.c"
        break;

    case YYSYMBOL_jump_statement: /* jump_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3680 "basilisk.c"
        break;

    case YYSYMBOL_external_declaration: /* external_declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3689 "basilisk.c"
        break;

    case YYSYMBOL_function_declaration: /* function_declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3698 "basilisk.c"
        break;

    case YYSYMBOL_function_definition: /* function_definition  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3707 "basilisk.c"
        break;

    case YYSYMBOL_declaration_list: /* declaration_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3716 "basilisk.c"
        break;

    case YYSYMBOL_basilisk_statements: /* basilisk_statements  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3725 "basilisk.c"
        break;

    case YYSYMBOL_macro_statement: /* macro_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3734 "basilisk.c"
        break;

    case YYSYMBOL_reduction_list: /* reduction_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3743 "basilisk.c"
        break;

    case YYSYMBOL_reduction: /* reduction  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3752 "basilisk.c"
        break;

    case YYSYMBOL_reduction_operator: /* reduction_operator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3761 "basilisk.c"
        break;

    case YYSYMBOL_reduction_array: /* reduction_array  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3770 "basilisk.c"
        break;

    case YYSYMBOL_foreach_dimension_statement: /* foreach_dimension_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3779 "basilisk.c"
        break;

    case YYSYMBOL_forin_declaration_statement: /* forin_declaration_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3788 "basilisk.c"
        break;

    case YYSYMBOL_forin_statement: /* forin_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3797 "basilisk.c"
        break;

    case YYSYMBOL_forin_arguments: /* forin_arguments  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3806 "basilisk.c"
        break;

    case YYSYMBOL_event_definition: /* event_definition  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3815 "basilisk.c"
        break;

    case YYSYMBOL_event_parameters: /* event_parameters  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3824 "basilisk.c"
        break;

    case YYSYMBOL_event_parameter: /* event_parameter  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3833 "basilisk.c"
        break;

    case YYSYMBOL_boundary_definition: /* boundary_definition  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3842 "basilisk.c"
        break;

    case YYSYMBOL_external_foreach_dimension: /* external_foreach_dimension  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3851 "basilisk.c"
        break;

    case YYSYMBOL_attribute: /* attribute  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3860 "basilisk.c"
        break;

    case YYSYMBOL_new_field: /* new_field  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3869 "basilisk.c"
        break;

    case YYSYMBOL_root: /* root  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3878 "basilisk.c"
        break;

      default:
        break;
    }
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
#line 118 "basilisk.yacc"
                                      { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4274 "basilisk.c"
    break;

  case 5: /* translation_unit: translation_unit error '}'  */
#line 119 "basilisk.yacc"
                                      { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4280 "basilisk.c"
    break;

  case 6: /* translation_unit: translation_unit error ')'  */
#line 120 "basilisk.yacc"
                                      { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4286 "basilisk.c"
    break;

  case 15: /* expression_error: error  */
#line 135 "basilisk.yacc"
                        { yyvsp[0]->sym = YYSYMBOL_YYerror; }
#line 4292 "basilisk.c"
    break;

  case 113: /* declaration: declaration_specifiers ';' type_not_specified  */
#line 331 "basilisk.yacc"
                                                                               { ast_push_declaration (parse->stack, yyval); }
#line 4298 "basilisk.c"
    break;

  case 114: /* declaration: declaration_specifiers init_declarator_list ';' type_not_specified  */
#line 332 "basilisk.yacc"
                                                                               { ast_push_declaration (parse->stack, yyval); }
#line 4304 "basilisk.c"
    break;

  case 138: /* type_specifier: types  */
#line 371 "basilisk.yacc"
                   { parse->type_already_specified = true; }
#line 4310 "basilisk.c"
    break;

  case 158: /* struct_or_union_specifier: struct_or_union '{' error '}'  */
#line 397 "basilisk.yacc"
                                                            { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4316 "basilisk.c"
    break;

  case 159: /* struct_or_union_specifier: struct_or_union generic_identifier '{' error '}'  */
#line 398 "basilisk.yacc"
                                                            { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4322 "basilisk.c"
    break;

  case 209: /* direct_declarator: direct_declarator '(' type_not_specified error ')'  */
#line 493 "basilisk.yacc"
                                                                                           { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4328 "basilisk.c"
    break;

  case 259: /* type_not_specified: %empty  */
#line 576 "basilisk.yacc"
           { parse->type_already_specified = false; }
#line 4334 "basilisk.c"
    break;

  case 280: /* statement: error ';'  */
#line 618 "basilisk.yacc"
                     { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4340 "basilisk.c"
    break;

  case 285: /* @1: %empty  */
#line 630 "basilisk.yacc"
                           { stack_push (parse->stack, &(yyvsp[0])); yyval->sym = YYSYMBOL_YYUNDEF; }
#line 4346 "basilisk.c"
    break;

  case 286: /* compound_statement: '{' @1 block_item_list '}'  */
#line 632 "basilisk.yacc"
                           { ast_pop_scope (parse->stack, yyvsp[-3]); }
#line 4352 "basilisk.c"
    break;

  case 296: /* for_scope: FOR  */
#line 657 "basilisk.yacc"
              { stack_push (parse->stack, &(yyval)); }
#line 4358 "basilisk.c"
    break;

  case 299: /* iteration_statement: for_scope '(' expression_statement expression_statement ')' statement  */
#line 664 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-5]); }
#line 4364 "basilisk.c"
    break;

  case 300: /* iteration_statement: for_scope '(' expression_statement expression_statement expression ')' statement  */
#line 666 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-6]); }
#line 4370 "basilisk.c"
    break;

  case 302: /* for_declaration_statement: for_scope '(' declaration expression_statement ')' statement  */
#line 672 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-5]); }
#line 4376 "basilisk.c"
    break;

  case 303: /* for_declaration_statement: for_scope '(' declaration expression_statement expression ')' statement  */
#line 674 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-6]); }
#line 4382 "basilisk.c"
    break;

  case 315: /* external_declaration: error compound_statement  */
#line 692 "basilisk.yacc"
                                                { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4388 "basilisk.c"
    break;

  case 316: /* function_declaration: declaration_specifiers declarator  */
#line 696 "basilisk.yacc"
                                            { ast_push_function_definition (parse->stack, yyvsp[0]);  }
#line 4394 "basilisk.c"
    break;

  case 317: /* function_definition: function_declaration declaration_list compound_statement  */
#line 701 "basilisk.yacc"
                                                { ast_pop_scope (parse->stack, yyvsp[-2]->child[1]); }
#line 4400 "basilisk.c"
    break;

  case 318: /* function_definition: function_declaration compound_statement  */
#line 703 "basilisk.yacc"
                                                { ast_pop_scope (parse->stack, yyvsp[-1]->child[1]); }
#line 4406 "basilisk.c"
    break;

  case 338: /* forin_declaration_statement: for_scope '(' declaration_specifiers declarator IDENTIFIER forin_arguments ')' statement  */
#line 754 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-7]); }
#line 4412 "basilisk.c"
    break;

  case 339: /* forin_statement: for_scope '(' expression IDENTIFIER forin_arguments ')' statement  */
#line 759 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-6]); }
#line 4418 "basilisk.c"
    break;

  case 355: /* root: translation_unit  */
#line 802 "basilisk.yacc"
                           {
	  yyval = root;
	  yyval->sym = yyr1[yyn];
	  yyval->child = allocate ((Allocator *)parse->alloc, 2*sizeof(Ast *));
	  yyval->child[0] = yyvsp[0];
	  yyvsp[0]->parent = yyval;
	  yyval->child[1] = NULL;
        }
#line 4431 "basilisk.c"
    break;


#line 4435 "basilisk.c"

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

#line 812 "basilisk.yacc"


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
