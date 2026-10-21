/* A Bison parser, made by GNU Bison 3.8.2.  */

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
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

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
#define YYBISON 30802

/* Bison version string.  */
#define YYBISON_VERSION "3.8.2"

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
  YYSYMBOL_ELLIPSIS_EXPRESSION = 83,       /* ELLIPSIS_EXPRESSION  */
  YYSYMBOL_MACRODEF = 84,                  /* MACRODEF  */
  YYSYMBOL_foreach_statement = 85,         /* foreach_statement  */
  YYSYMBOL_86_ = 86,                       /* ';'  */
  YYSYMBOL_87_ = 87,                       /* '}'  */
  YYSYMBOL_88_ = 88,                       /* ')'  */
  YYSYMBOL_89_ = 89,                       /* '('  */
  YYSYMBOL_90_ = 90,                       /* ','  */
  YYSYMBOL_91_ = 91,                       /* ':'  */
  YYSYMBOL_92_ = 92,                       /* '.'  */
  YYSYMBOL_93_ = 93,                       /* '{'  */
  YYSYMBOL_94_ = 94,                       /* '['  */
  YYSYMBOL_95_ = 95,                       /* ']'  */
  YYSYMBOL_96_ = 96,                       /* '*'  */
  YYSYMBOL_97_ = 97,                       /* '&'  */
  YYSYMBOL_98_ = 98,                       /* '+'  */
  YYSYMBOL_99_ = 99,                       /* '-'  */
  YYSYMBOL_100_ = 100,                     /* '~'  */
  YYSYMBOL_101_ = 101,                     /* '!'  */
  YYSYMBOL_102_ = 102,                     /* '/'  */
  YYSYMBOL_103_ = 103,                     /* '%'  */
  YYSYMBOL_104_ = 104,                     /* '<'  */
  YYSYMBOL_105_ = 105,                     /* '>'  */
  YYSYMBOL_106_ = 106,                     /* '^'  */
  YYSYMBOL_107_ = 107,                     /* '|'  */
  YYSYMBOL_108_ = 108,                     /* '?'  */
  YYSYMBOL_109_ = 109,                     /* '='  */
  YYSYMBOL_YYACCEPT = 110,                 /* $accept  */
  YYSYMBOL_translation_unit = 111,         /* translation_unit  */
  YYSYMBOL_primary_expression = 112,       /* primary_expression  */
  YYSYMBOL_expression_error = 113,         /* expression_error  */
  YYSYMBOL_constant = 114,                 /* constant  */
  YYSYMBOL_enumeration_constant = 115,     /* enumeration_constant  */
  YYSYMBOL_string = 116,                   /* string  */
  YYSYMBOL_generic_selection = 117,        /* generic_selection  */
  YYSYMBOL_generic_assoc_list = 118,       /* generic_assoc_list  */
  YYSYMBOL_generic_association = 119,      /* generic_association  */
  YYSYMBOL_postfix_expression = 120,       /* postfix_expression  */
  YYSYMBOL_postfix_initializer = 121,      /* postfix_initializer  */
  YYSYMBOL_array_access = 122,             /* array_access  */
  YYSYMBOL_function_call = 123,            /* function_call  */
  YYSYMBOL_member_identifier = 124,        /* member_identifier  */
  YYSYMBOL_argument_expression_list = 125, /* argument_expression_list  */
  YYSYMBOL_argument_expression_list_item = 126, /* argument_expression_list_item  */
  YYSYMBOL_unary_expression = 127,         /* unary_expression  */
  YYSYMBOL_unary_operator = 128,           /* unary_operator  */
  YYSYMBOL_cast_expression = 129,          /* cast_expression  */
  YYSYMBOL_multiplicative_expression = 130, /* multiplicative_expression  */
  YYSYMBOL_additive_expression = 131,      /* additive_expression  */
  YYSYMBOL_shift_expression = 132,         /* shift_expression  */
  YYSYMBOL_relational_expression = 133,    /* relational_expression  */
  YYSYMBOL_equality_expression = 134,      /* equality_expression  */
  YYSYMBOL_and_expression = 135,           /* and_expression  */
  YYSYMBOL_exclusive_or_expression = 136,  /* exclusive_or_expression  */
  YYSYMBOL_inclusive_or_expression = 137,  /* inclusive_or_expression  */
  YYSYMBOL_logical_and_expression = 138,   /* logical_and_expression  */
  YYSYMBOL_logical_or_expression = 139,    /* logical_or_expression  */
  YYSYMBOL_conditional_expression = 140,   /* conditional_expression  */
  YYSYMBOL_assignment_expression = 141,    /* assignment_expression  */
  YYSYMBOL_assignment_operator = 142,      /* assignment_operator  */
  YYSYMBOL_expression = 143,               /* expression  */
  YYSYMBOL_constant_expression = 144,      /* constant_expression  */
  YYSYMBOL_declaration = 145,              /* declaration  */
  YYSYMBOL_declaration_specifiers = 146,   /* declaration_specifiers  */
  YYSYMBOL_init_declarator_list = 147,     /* init_declarator_list  */
  YYSYMBOL_init_declarator = 148,          /* init_declarator  */
  YYSYMBOL_storage_class_specifier = 149,  /* storage_class_specifier  */
  YYSYMBOL_type_specifier = 150,           /* type_specifier  */
  YYSYMBOL_types = 151,                    /* types  */
  YYSYMBOL_struct_or_union_specifier = 152, /* struct_or_union_specifier  */
  YYSYMBOL_struct_or_union = 153,          /* struct_or_union  */
  YYSYMBOL_struct_declaration_list = 154,  /* struct_declaration_list  */
  YYSYMBOL_struct_declaration = 155,       /* struct_declaration  */
  YYSYMBOL_specifier_qualifier_list = 156, /* specifier_qualifier_list  */
  YYSYMBOL_struct_declarator_list = 157,   /* struct_declarator_list  */
  YYSYMBOL_struct_declarator = 158,        /* struct_declarator  */
  YYSYMBOL_enum_specifier = 159,           /* enum_specifier  */
  YYSYMBOL_enumerator_list = 160,          /* enumerator_list  */
  YYSYMBOL_enumerator = 161,               /* enumerator  */
  YYSYMBOL_atomic_type_specifier = 162,    /* atomic_type_specifier  */
  YYSYMBOL_type_qualifier = 163,           /* type_qualifier  */
  YYSYMBOL_function_specifier = 164,       /* function_specifier  */
  YYSYMBOL_alignment_specifier = 165,      /* alignment_specifier  */
  YYSYMBOL_declarator = 166,               /* declarator  */
  YYSYMBOL_direct_declarator = 167,        /* direct_declarator  */
  YYSYMBOL_generic_identifier = 168,       /* generic_identifier  */
  YYSYMBOL_pointer = 169,                  /* pointer  */
  YYSYMBOL_type_qualifier_list = 170,      /* type_qualifier_list  */
  YYSYMBOL_parameter_type_list = 171,      /* parameter_type_list  */
  YYSYMBOL_parameter_list = 172,           /* parameter_list  */
  YYSYMBOL_parameter_declaration = 173,    /* parameter_declaration  */
  YYSYMBOL_identifier_list = 174,          /* identifier_list  */
  YYSYMBOL_type_name = 175,                /* type_name  */
  YYSYMBOL_abstract_declarator = 176,      /* abstract_declarator  */
  YYSYMBOL_direct_abstract_declarator = 177, /* direct_abstract_declarator  */
  YYSYMBOL_type_not_specified = 178,       /* type_not_specified  */
  YYSYMBOL_initializer = 179,              /* initializer  */
  YYSYMBOL_initializer_list = 180,         /* initializer_list  */
  YYSYMBOL_designation = 181,              /* designation  */
  YYSYMBOL_designator_list = 182,          /* designator_list  */
  YYSYMBOL_designator = 183,               /* designator  */
  YYSYMBOL_static_assert_declaration = 184, /* static_assert_declaration  */
  YYSYMBOL_statement = 185,                /* statement  */
  YYSYMBOL_labeled_statement = 186,        /* labeled_statement  */
  YYSYMBOL_compound_statement = 187,       /* compound_statement  */
  YYSYMBOL_188_1 = 188,                    /* @1  */
  YYSYMBOL_block_item_list = 189,          /* block_item_list  */
  YYSYMBOL_block_item = 190,               /* block_item  */
  YYSYMBOL_expression_statement = 191,     /* expression_statement  */
  YYSYMBOL_selection_statement = 192,      /* selection_statement  */
  YYSYMBOL_for_scope = 193,                /* for_scope  */
  YYSYMBOL_iteration_statement = 194,      /* iteration_statement  */
  YYSYMBOL_for_declaration_statement = 195, /* for_declaration_statement  */
  YYSYMBOL_jump_statement = 196,           /* jump_statement  */
  YYSYMBOL_external_declaration = 197,     /* external_declaration  */
  YYSYMBOL_function_declaration = 198,     /* function_declaration  */
  YYSYMBOL_function_definition = 199,      /* function_definition  */
  YYSYMBOL_declaration_list = 200,         /* declaration_list  */
  YYSYMBOL_basilisk_statements = 201,      /* basilisk_statements  */
  YYSYMBOL_macro_statement = 202,          /* macro_statement  */
  YYSYMBOL_reduction_list = 203,           /* reduction_list  */
  YYSYMBOL_reduction = 204,                /* reduction  */
  YYSYMBOL_reduction_operator = 205,       /* reduction_operator  */
  YYSYMBOL_reduction_array = 206,          /* reduction_array  */
  YYSYMBOL_foreach_dimension_statement = 207, /* foreach_dimension_statement  */
  YYSYMBOL_forin_declaration_statement = 208, /* forin_declaration_statement  */
  YYSYMBOL_forin_statement = 209,          /* forin_statement  */
  YYSYMBOL_forin_arguments = 210,          /* forin_arguments  */
  YYSYMBOL_event_definition = 211,         /* event_definition  */
  YYSYMBOL_event_parameters = 212,         /* event_parameters  */
  YYSYMBOL_event_parameter = 213,          /* event_parameter  */
  YYSYMBOL_boundary_definition = 214,      /* boundary_definition  */
  YYSYMBOL_external_foreach_dimension = 215, /* external_foreach_dimension  */
  YYSYMBOL_attribute = 216,                /* attribute  */
  YYSYMBOL_new_field = 217,                /* new_field  */
  YYSYMBOL_root = 218                      /* root  */
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

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
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
#define YYFINAL  204
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   3461

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  110
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  109
/* YYNRULES -- Number of rules.  */
#define YYNRULES  356
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  637

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   340


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
       2,     2,     2,   101,     2,     2,     2,   103,    97,     2,
      89,    88,    96,    98,    90,    99,    92,   102,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    91,    86,
     104,   109,   105,   108,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    94,     2,    95,   106,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    93,   107,    87,   100,     2,     2,     2,
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
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   116,   116,   117,   118,   119,   120,   124,   125,   126,
     127,   128,   129,   130,   131,   135,   136,   140,   141,   142,
     146,   150,   151,   155,   159,   160,   164,   165,   169,   170,
     171,   172,   173,   174,   175,   176,   180,   181,   185,   186,
     187,   191,   192,   196,   200,   201,   202,   206,   207,   211,
     212,   213,   214,   215,   216,   217,   218,   222,   223,   224,
     225,   226,   227,   231,   232,   238,   239,   240,   241,   245,
     246,   247,   251,   252,   253,   257,   258,   259,   260,   261,
     265,   266,   267,   271,   272,   276,   277,   281,   282,   286,
     287,   291,   292,   296,   297,   301,   302,   303,   304,   305,
     309,   310,   311,   312,   313,   314,   315,   316,   317,   318,
     319,   323,   324,   328,   332,   333,   334,   338,   339,   340,
     341,   342,   343,   344,   345,   346,   347,   351,   352,   356,
     357,   361,   362,   363,   364,   365,   366,   367,   368,   372,
     376,   377,   378,   379,   380,   381,   382,   383,   384,   385,
     386,   387,   388,   389,   390,   391,   395,   396,   397,   398,
     399,   403,   404,   408,   409,   413,   414,   415,   419,   420,
     421,   422,   426,   427,   431,   432,   433,   437,   438,   439,
     440,   441,   445,   446,   450,   451,   455,   459,   460,   461,
     462,   463,   467,   468,   472,   473,   477,   478,   482,   483,
     484,   485,   486,   487,   488,   489,   490,   491,   492,   493,
     494,   495,   496,   500,   501,   502,   506,   507,   508,   509,
     513,   514,   518,   519,   523,   524,   528,   529,   530,   531,
     532,   533,   537,   538,   542,   543,   547,   548,   549,   553,
     554,   555,   556,   557,   558,   559,   560,   561,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     577,   581,   582,   583,   587,   588,   589,   590,   594,   598,
     599,   603,   604,   608,   612,   613,   614,   615,   616,   617,
     618,   619,   623,   624,   625,   629,   631,   631,   637,   638,
     642,   643,   647,   648,   652,   653,   654,   658,   662,   663,
     664,   666,   668,   672,   674,   679,   680,   681,   682,   683,
     687,   688,   689,   690,   691,   692,   693,   697,   701,   703,
     708,   709,   716,   717,   718,   719,   720,   724,   725,   729,
     730,   734,   738,   739,   740,   744,   745,   749,   750,   754,
     759,   764,   765,   769,   773,   774,   775,   779,   780,   781,
     785,   789,   790,   794,   798,   799,   803
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
  "ELLIPSIS_EXPRESSION", "MACRODEF", "foreach_statement", "';'", "'}'",
  "')'", "'('", "','", "':'", "'.'", "'{'", "'['", "']'", "'*'", "'&'",
  "'+'", "'-'", "'~'", "'!'", "'/'", "'%'", "'<'", "'>'", "'^'", "'|'",
  "'?'", "'='", "$accept", "translation_unit", "primary_expression",
  "expression_error", "constant", "enumeration_constant", "string",
  "generic_selection", "generic_assoc_list", "generic_association",
  "postfix_expression", "postfix_initializer", "array_access",
  "function_call", "member_identifier", "argument_expression_list",
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

#define YYPACT_NINF (-375)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-357)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     981,   -15,    76,  -375,  -375,  -375,  -375,  2537,  2591,  2591,
     493,  -375,  -375,  -375,  -375,  -375,  -375,  -375,  -375,  -375,
    -375,  -375,  -375,  -375,  -375,  -375,  -375,  -375,  -375,  -375,
    -375,  -375,  -375,  -375,  -375,    51,   -55,    -5,   120,   139,
    -375,   151,  -375,  -375,   153,  -375,   161,   176,   314,  -375,
    -375,  1082,  -375,  -375,  -375,  -375,  -375,  -375,   678,  -375,
    -375,  -375,  -375,   556,  -375,  -375,   567,  2636,  -375,   -59,
     249,   343,    34,   346,   200,   201,   219,   317,    17,  -375,
     253,  -375,   168,  3260,  3260,  -375,  -375,   328,  -375,  -375,
    3260,  3260,  3260,   334,  -375,  -375,  2938,  -375,   279,  -375,
    -375,  -375,  -375,  -375,  -375,   366,   256,  -375,  -375,  -375,
    1082,  -375,  1082,  -375,  -375,  -375,  -375,  -375,  -375,  -375,
    -375,  -375,  -375,  -375,  -375,  -375,  2096,  -375,  -375,  -375,
     368,   280,  1436,  3385,  3385,  2438,  2636,   575,    14,    52,
    -375,   567,   301,  -375,   303,  3385,   106,  3385,   310,   352,
    -375,    62,  -375,  -375,  1943,    62,  2141,  2096,  -375,  -375,
    2636,  2636,  2636,  2636,  2636,  2636,  2636,  2636,  2636,  2636,
    2636,  2636,  2636,  2636,  2636,  2636,  2636,  2636,  2438,  -375,
    -375,   330,    -9,   -40,  -375,   -58,   -67,  -375,   177,  -375,
    -375,  -375,  1246,   313,  -375,  -375,  -375,  3307,   302,  -375,
     168,  -375,  2938,  -375,  -375,  -375,   880,   324,   326,   567,
    1181,  -375,  -375,  -375,   307,    92,  -375,   368,  -375,   336,
     347,   348,   354,   353,   357,   458,   361,  3260,  -375,  -375,
    -375,   327,  -375,  2438,  -375,  2873,  1490,   -54,  -375,   -42,
    -375,  -375,  2492,  -375,  -375,  -375,  -375,  -375,  -375,  -375,
     104,  -375,  -375,  -375,  -375,   356,   122,  -375,  -375,  -375,
    -375,  -375,   -59,   -59,   249,   249,   343,   343,   343,   343,
      34,    34,   346,   200,   201,   219,   317,   295,  -375,   365,
    -375,  -375,  -375,    -9,  -375,   330,  2195,  -375,  1538,  -375,
     -67,   367,  3051,  -375,    77,  -375,  2828,  3101,  2636,   350,
    -375,  -375,   370,  2636,   369,   372,   375,   376,  1137,  -375,
      62,   380,   384,  2240,   382,   198,  -375,  -375,   -33,  -375,
     383,  -375,  -375,  -375,   779,  -375,  -375,  -375,   391,  -375,
    -375,  -375,  -375,  -375,  -375,  -375,  -375,   388,   388,    62,
    1181,  2636,  -375,  -375,   116,  2195,   -71,  -375,  2636,  -375,
      67,   129,  -375,  -375,  -375,  -375,  3346,   467,  -375,  3260,
     330,  -375,    62,  -375,   373,  -375,   298,   397,   398,  -375,
     401,  1796,  -375,   395,   406,  1592,   -42,  -375,  -375,  1640,
    -375,  -375,  -375,  2096,  -375,  -375,  2636,  -375,  -375,  -375,
    -375,  -375,  -375,  2773,  1796,  -375,   408,   409,  1694,  -375,
    -375,  -375,  -375,  -375,  2636,   -26,  -375,   414,   422,  3151,
    -375,   567,  -375,   -14,  -375,  -375,   419,  1137,  1238,  1238,
    2438,   567,   461,   438,  -375,  -375,  -375,   136,    90,  1997,
    -375,  1137,  -375,  -375,  1337,  -375,   147,   431,  -375,  1844,
    -375,  -375,  -375,  -375,  -375,  -375,  -375,   101,   444,   206,
    -375,   445,   441,  -375,  -375,   447,   450,  2042,  2708,   433,
      53,  -375,  -375,  3207,  -375,  -375,   448,  1796,  -375,  -375,
    2438,  -375,   451,  3000,  1796,  -375,   453,   460,  1742,  -375,
    -375,   457,  -375,  -375,   468,   210,   463,  1796,  -375,  -375,
    2438,  -375,   464,   473,  -375,  -375,  -375,   233,  2636,  -375,
    -375,  2492,  2636,  1137,  2636,  1137,  -375,   475,   481,   214,
     482,  -375,  -375,   484,  1137,  1137,   223,  -375,    66,  2294,
     168,  2294,  -375,  1898,  -375,  -375,  -375,  2195,  -375,  2438,
    -375,  3346,  2438,   489,   486,  -375,  -375,  -375,  -375,  2195,
    -375,  -375,  -375,   488,   490,  -375,  -375,   509,   503,  1796,
    -375,  -375,  2438,  -375,   504,  -375,  -375,  -375,   597,  -375,
     506,   508,  -375,  -375,  -375,  -375,  -375,  -375,  -375,  -375,
    -375,  -375,  -375,  1137,  1137,  1137,  2438,  1137,  -375,  -375,
    1137,  2096,  2339,    21,  2393,  -375,  -375,  -375,  -375,  -375,
    -375,  2438,  -375,  -375,  -375,  -375,  -375,  -375,   512,   515,
    -375,  -375,  -375,  -375,   544,  -375,  -375,   228,  -375,  -375,
    -375,   303,   524,  1137,   235,  2096,  1137,   242,   166,  -375,
    -375,  1137,   527,  1137,  -375,  1137,   526,  -375,  1137,  -375,
    -375,  -375,  -375,  -375,  1137,  -375,  -375
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_int16 yydefact[] =
{
       0,     0,     7,    17,    18,    21,    22,     0,     0,     0,
     155,    19,   131,   132,   133,   135,   136,   192,   187,   188,
     189,   149,   141,   142,   143,   144,   147,   148,   145,   146,
     140,   150,   151,   161,   162,     0,     0,     0,   190,     0,
     193,     0,   134,   191,   354,   137,     0,     0,     8,     9,
     138,     0,    58,    57,    59,    60,    61,    62,     0,    28,
      10,    11,    13,    49,    30,    29,    63,     0,    65,    69,
      72,    75,    80,    83,    85,    87,    89,    91,    93,    95,
       0,   311,     0,   118,   120,   139,   153,     0,   154,   152,
     122,   124,   126,     0,   116,     2,     0,   310,    14,   329,
     312,   313,   314,   315,    56,     0,   286,   316,     7,     8,
       0,    53,     0,    50,    51,   101,   102,   103,   104,   105,
     106,   107,   108,   109,   110,   100,     0,   213,   214,   215,
       0,   181,     0,     0,     0,     0,     0,     0,     0,     0,
      16,   155,     0,   111,    15,   169,   260,   171,     0,     0,
       3,     0,    33,    34,     0,     0,     0,     0,    63,    52,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   350,
     260,     0,   219,     0,   127,   317,   260,   198,     0,   155,
     117,   119,     0,   158,   121,   123,   125,     0,     0,   320,
       0,   319,     0,   330,     1,   285,     0,     0,     0,     0,
       0,    99,    98,    20,   185,     0,   182,     0,   113,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   334,   333,
     332,     0,    12,     0,   168,     0,     0,   237,   260,   238,
     235,   170,     0,     4,     5,     6,    32,    43,    41,    48,
       0,    44,    47,    31,    38,    58,     0,    97,    96,    66,
      67,    68,    70,    71,    73,    74,    78,    79,    76,    77,
      81,    82,    84,    86,    88,    90,    92,     0,   114,     0,
     190,   220,   218,   217,   260,     0,     0,   260,     0,   197,
     260,     0,     0,   163,     0,   167,     0,     0,     0,   130,
     321,   318,     0,     0,     0,     0,     0,     0,     0,   297,
       0,     0,     0,     0,     0,     8,   322,   292,     0,   290,
       0,   291,   274,   275,     0,   288,   276,   277,     0,   278,
     302,   279,   280,   323,   324,   325,   326,    54,     0,     0,
       0,     0,   263,   265,     0,     0,     0,   269,     0,   177,
       0,     0,   195,   194,    55,   186,     0,     0,   355,     0,
       0,   351,     0,   112,     0,   256,   229,     0,   223,   260,
       0,     0,   240,    58,     0,     0,   236,   234,   260,     0,
      35,    64,    42,    46,    39,    40,     0,   199,   221,   216,
     115,   128,   129,     0,     0,   200,    58,     0,     0,   196,
     159,   156,   164,   260,     0,     0,   172,   176,     0,     0,
     353,    63,   347,     0,   344,   281,     0,     0,     0,     0,
       0,   214,     0,     0,   306,   307,   308,     0,     0,     0,
     293,     0,   287,   289,     0,   272,     0,     0,    36,     0,
     264,   268,   270,   184,   178,   183,   179,     0,     0,     0,
      24,     0,     0,   352,   317,   335,     0,     0,     0,   226,
     237,   228,   257,     0,   224,   239,     0,     0,   241,   247,
       0,   246,     0,     0,     0,   248,    58,     0,     0,    45,
      94,     0,   232,   211,     0,     0,     0,     0,   201,   208,
       0,   207,    58,     0,   165,   174,   260,     0,     0,   160,
     157,     0,     0,     0,     0,     0,   284,     0,     0,     0,
       0,   305,   309,     0,     0,     0,     0,   282,     0,     0,
       0,     0,   261,     0,   271,    37,   267,     0,   180,     0,
      23,     0,     0,     0,     0,   331,   231,   230,   155,     0,
     222,   260,   243,     0,     0,   245,   258,     0,     0,     0,
     249,   255,     0,   254,     0,   210,   209,   212,     0,   203,
       0,     0,   204,   206,   166,   173,   175,   349,   348,   346,
     343,   345,   283,     0,     0,     0,     0,     0,   337,   327,
       0,     0,     0,   130,     0,   262,   266,    27,    25,    26,
     273,     0,   227,   225,   242,   244,   259,   251,     0,     0,
     252,   233,   202,   205,   295,   296,   298,     0,   338,   328,
     342,   341,     0,     0,     0,     0,     0,     0,     0,   250,
     253,     0,     0,     0,   303,     0,     0,   300,     0,   336,
     294,   299,   340,   304,     0,   301,   339
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -375,  -375,  -375,   -22,  -375,  -375,  -375,  -375,  -375,    84,
     479,  -124,  -375,  -375,   462,   189,   239,   107,  -375,   -41,
     237,   263,   329,   260,   452,   454,   449,   455,   456,  -375,
      72,     5,   -65,   -51,  -116,   -93,    93,  -375,   342,  -375,
     -44,  -375,  -375,  -375,  -185,  -281,  -106,  -375,   132,  -375,
     413,  -325,  -375,   -37,  -375,  -375,   -74,  -175,    42,  -137,
    -230,  -374,  -375,   169,  -375,   374,  -131,  -227,   -63,  -149,
     294,  -278,  -375,   289,    13,    46,  -375,    16,  -375,  -375,
     312,  -349,  -375,  -375,  -375,  -375,  -375,   579,  -375,  -223,
    -375,  -375,  -375,  -375,   540,  -375,  -375,  -375,  -375,  -375,
      24,  -375,  -375,  -151,  -375,  -375,  -375,  -375,  -375
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
       0,    58,    59,   142,    60,   214,    61,    62,   449,   450,
      63,   249,    64,    65,   246,   250,   251,    66,    67,    68,
      69,    70,    71,    72,    73,    74,    75,    76,    77,    78,
      79,   143,   126,   318,   219,    81,   366,   183,   184,    83,
      84,    85,    86,    87,   292,   293,   146,   405,   406,    88,
     215,   216,    89,    90,    91,    92,   279,   186,   320,   188,
     283,   367,   368,   369,   485,   451,   370,   239,   240,   343,
     344,   345,   346,   347,    94,   321,   322,   323,   206,   324,
     325,   326,   327,   328,   329,   330,   331,    95,    96,    97,
     202,   332,   333,    98,    99,   231,   456,   334,   335,   336,
     612,   100,   413,   414,   101,   102,   103,   104,   105
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     144,   157,   211,   199,   361,    80,   375,   145,   185,   237,
     376,   402,   297,   290,   147,   238,   402,   107,   226,   484,
     224,   339,   287,   341,   615,   445,   159,   288,  -130,    18,
      19,    20,  -130,   257,   132,   235,   177,   160,   441,   234,
     236,   241,    93,   161,   162,   282,   284,   378,   167,   168,
     285,   286,   379,   430,   127,   127,   127,   233,   398,   144,
     496,   144,   280,    80,   497,   127,   145,    43,   145,   581,
     213,   228,   502,   147,   503,   147,   504,   131,   106,  -213,
     127,   128,   128,   128,   133,   521,   294,   182,   145,   145,
     145,   294,   128,    82,   513,   147,   147,   147,   237,   547,
      93,   145,   227,   145,   213,   256,  -213,   128,   147,   300,
     147,   409,   201,   319,   111,   113,   114,   278,   380,   259,
     260,   261,   445,   289,   187,   178,   299,   277,   402,   193,
     286,   212,   129,   129,   129,   198,   453,   392,   169,   170,
     223,   467,   458,   129,   130,   281,   389,   236,   145,   478,
     229,    82,   430,   145,   444,   147,   233,  -213,   129,   252,
     147,   527,   258,   403,   487,   107,   181,  -213,   404,  -213,
     582,   127,   584,   182,   158,   377,   190,   191,   514,   349,
     127,   230,   350,   194,   195,   196,   294,   416,   528,   200,
     294,   294,   382,   247,   383,   235,   440,   247,   128,   281,
     236,   381,   182,   438,   218,   295,   439,   128,   218,   134,
     295,   299,   233,   380,   380,   342,   446,   385,   301,   447,
     407,   390,   512,   187,   393,   437,   233,   399,   135,   460,
     187,   319,   443,   376,   522,   461,   127,   523,   363,   158,
     136,   374,   187,   158,   549,   527,   388,   137,   145,   129,
     138,   281,   145,   145,   180,   147,   233,   181,   129,   147,
     147,   629,   427,   128,   182,   139,   181,   158,   158,   158,
     158,   158,   158,   158,   158,   158,   158,   158,   158,   158,
     158,   158,   158,   158,   158,   290,   454,   429,   495,  -215,
     526,   342,   459,   397,   530,   200,   531,   173,   557,   200,
     558,   127,   575,   294,   233,   295,   464,   174,   537,   295,
     295,   580,   145,   383,   129,   473,   622,  -215,   233,   147,
     360,   460,   181,   625,   404,   233,   175,   187,   128,   182,
     628,   127,   233,   127,   281,   176,   187,   127,   388,   179,
     494,   519,   281,   205,  -215,   342,   501,   163,   164,   158,
     342,   569,   423,   571,   422,   165,   166,   281,   128,    47,
     128,   388,   171,   172,   128,   145,   204,   144,   144,   509,
     412,   213,   147,   217,   526,   218,   466,   567,   586,   129,
     472,   435,   566,   518,   477,   233,   386,   458,   252,   232,
     592,   298,   236,   233,   182,  -215,   507,   508,   242,   486,
     262,   263,   187,   493,   455,   411,   296,  -215,   187,   129,
     158,   129,   337,   218,   338,   129,   348,   200,   362,   181,
     218,   192,   295,   407,   352,   148,   182,   197,   264,   265,
     388,   270,   271,   564,   252,   353,   354,   281,   243,   244,
     245,   388,   355,   356,   342,   106,   583,   357,   158,   359,
     388,   384,   360,   387,   400,   158,   415,   610,   480,   286,
     417,   418,   342,   506,   419,   420,   424,   151,   152,   153,
     425,   428,   543,   452,   431,   544,   218,   517,   593,   548,
     434,   210,   457,   554,   207,   462,   208,   145,   463,   465,
     468,   610,   560,   158,   147,   561,   266,   267,   268,   269,
     187,   469,   187,   488,   489,   498,   220,   221,   222,   499,
     505,   158,   388,   115,   116,   117,   118,   119,   120,   121,
     122,   123,   124,   510,   511,   607,   524,   520,   342,   533,
     611,   614,   342,   617,   587,   529,   532,   589,   535,   187,
     618,   534,   539,   542,   342,   555,   545,   154,   550,   570,
     155,   572,   156,   358,   598,   551,   556,   599,   559,   562,
     578,   579,   187,   573,   611,   151,   152,   153,   563,   574,
     218,   576,   577,   568,   412,   590,   412,   591,   108,     3,
       4,     5,     6,   594,  -214,   595,  -214,   115,   116,   117,
     118,   119,   120,   121,   122,   123,   124,   596,   597,   600,
     601,   602,   125,   603,   621,   158,    11,   619,   158,   411,
     620,   411,   623,   631,   634,   588,   225,   253,   516,   604,
     605,   606,   479,   608,   274,   272,   609,   391,   273,   565,
     351,   275,   541,   276,   436,   442,   433,   150,   203,   626,
       0,     0,     0,     0,     0,   154,     0,    39,   155,     0,
     156,     0,     0,     0,     0,    47,   109,     0,    49,   624,
       0,     0,   627,     0,   112,     0,     0,   630,     0,   632,
       0,   633,     0,     0,   635,     0,   125,     0,  -356,   149,
     636,     2,     3,     4,     5,     6,     7,     0,     8,     9,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
       0,    49,    50,     0,     0,     0,     0,    51,     0,     0,
       0,     0,     0,     0,    52,    53,    54,    55,    56,    57,
     302,     0,     2,     3,     4,     5,     6,     7,     0,     8,
       9,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,     0,   303,   304,   305,     0,
     306,   307,   308,   309,   310,   311,   312,   313,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,   314,    47,
     315,   316,    49,    50,     0,   317,   432,     0,    51,     0,
       0,     0,   106,     0,     0,    52,    53,    54,    55,    56,
      57,   302,     0,     2,     3,     4,     5,     6,     7,     0,
       8,     9,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,     0,   303,   304,   305,
       0,   306,   307,   308,   309,   310,   311,   312,   313,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,   314,
      47,   315,   316,    49,    50,     0,   317,     0,     0,    51,
       0,     0,     0,   106,     0,     0,    52,    53,    54,    55,
      56,    57,     1,     0,     2,     3,     4,     5,     6,     7,
       0,     8,     9,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,     0,    49,    50,     0,     0,     0,     0,
      51,     0,     0,     0,     0,     0,     0,    52,    53,    54,
      55,    56,    57,   140,     0,   108,     3,     4,     5,     6,
       7,     0,     8,     9,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   141,    11,     0,     0,     0,     0,     0,     0,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,   302,     0,
       2,     3,     4,     5,     6,     7,     0,     8,     9,     0,
       0,     0,    37,    38,    39,     0,     0,     0,    43,    44,
       0,     0,    47,   109,     0,    49,     0,   421,    11,     0,
       0,    51,     0,     0,     0,     0,     0,     0,    52,    53,
      54,    55,    56,    57,   108,     3,     4,     5,     6,     7,
       0,     8,     9,     0,   303,   304,   305,     0,   306,   307,
     308,   309,   310,   311,   312,   313,     0,    37,     0,    39,
       0,   209,    11,     0,    44,     0,   314,    47,   315,   316,
      49,     0,     0,   317,     0,     0,    51,     0,     0,     0,
     106,     0,     0,    52,    53,    54,    55,    56,    57,   140,
       0,   108,     3,     4,     5,     6,     7,   291,     8,     9,
       0,    37,     0,    39,     0,     0,     0,     0,    44,     0,
       0,    47,   109,     0,    49,     0,     0,     0,   209,    11,
      51,     0,     0,   339,   340,   341,   189,    52,    53,    54,
      55,    56,    57,     0,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,     0,     0,     0,     0,     0,     0,    37,     0,
      39,     0,     0,     0,     0,    44,     0,    38,    47,   109,
      41,    49,    43,     0,     0,     0,     0,    51,     0,     0,
       0,     0,     0,     0,    52,    53,    54,    55,    56,    57,
     108,     3,     4,     5,     6,     7,     0,     8,     9,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   141,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,     0,    47,   109,     0,
      49,    50,     0,   317,     0,     0,    51,     0,     0,     0,
       0,     0,     0,    52,    53,    54,    55,    56,    57,   108,
       3,     4,     5,     6,     7,     0,     8,     9,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   189,    11,     0,     0,
       0,     0,     0,     0,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,     0,   108,     3,     4,     5,     6,     7,     0,
       8,     9,     0,     0,     0,     0,    37,    38,    39,     0,
       0,     0,    43,    44,     0,     0,    47,   109,     0,    49,
     209,    11,     0,     0,   371,    51,     0,     0,    18,    19,
      20,     0,    52,    53,    54,    55,    56,    57,     0,     0,
       0,   108,     3,     4,     5,     6,     7,     0,     8,     9,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      37,   280,    39,     0,     0,     0,    43,    44,   209,    11,
      47,   109,   394,    49,     0,     0,    18,    19,    20,    51,
       0,     0,     0,     0,     0,   372,   373,    53,    54,    55,
      56,    57,     0,     0,     0,   108,     3,     4,     5,     6,
       7,     0,     8,     9,     0,     0,     0,     0,    37,   280,
      39,     0,     0,     0,    43,    44,     0,     0,    47,   109,
       0,    49,   209,    11,     0,     0,   470,    51,     0,     0,
      18,    19,    20,   395,   396,    53,    54,    55,    56,    57,
       0,     0,     0,   108,     3,     4,     5,     6,     7,     0,
       8,     9,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    37,   280,    39,     0,     0,     0,    43,    44,
     209,    11,    47,   109,   474,    49,     0,     0,    18,    19,
      20,    51,     0,     0,     0,     0,     0,   471,    52,    53,
      54,    55,    56,    57,     0,     0,     0,   108,     3,     4,
       5,     6,     7,     0,     8,     9,     0,     0,     0,     0,
      37,   280,    39,     0,     0,     0,    43,    44,     0,     0,
      47,   109,     0,    49,   209,    11,     0,     0,   490,    51,
       0,     0,    18,    19,    20,   475,   476,    53,    54,    55,
      56,    57,     0,     0,     0,   108,     3,     4,     5,     6,
       7,     0,     8,     9,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    37,   280,    39,     0,     0,     0,
      43,    44,   209,    11,    47,   109,   552,    49,     0,     0,
      18,    19,    20,    51,     0,     0,     0,     0,     0,   491,
     492,    53,    54,    55,    56,    57,     0,     0,     0,   108,
       3,     4,     5,     6,     7,     0,     8,     9,     0,     0,
       0,     0,    37,   280,    39,     0,     0,     0,    43,    44,
       0,     0,    47,   109,     0,    49,   209,    11,     0,     0,
       0,    51,     0,     0,    18,    19,    20,   553,    52,    53,
      54,    55,    56,    57,     0,     0,     0,   108,     3,     4,
       5,     6,     7,     0,     8,     9,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    37,   280,    39,     0,
       0,     0,    43,    44,   209,    11,    47,   109,     0,    49,
       0,     0,     0,     0,     0,    51,     0,     0,     0,     0,
       0,     0,    52,    53,    54,    55,    56,    57,     0,     0,
       0,   108,     3,     4,     5,     6,     7,     0,     8,     9,
       0,     0,     0,     0,    37,     0,    39,     0,     0,     0,
       0,    44,     0,     0,    47,   109,     0,    49,   209,    11,
       0,   525,     0,    51,     0,     0,   339,   340,   341,     0,
      52,    53,    54,    55,    56,    57,   108,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    37,     0,
      39,     0,     0,   209,    11,    44,     0,     0,    47,   109,
       0,    49,     0,     0,     0,   585,     0,    51,     0,     0,
     339,   340,   341,     0,    52,    53,    54,    55,    56,    57,
     108,     3,     4,     5,     6,     7,     0,     8,     9,     0,
       0,     0,     0,    37,     0,    39,     0,     0,     0,     0,
      44,     0,     0,    47,   109,     0,    49,   209,    11,     0,
       0,   248,    51,     0,     0,     0,   210,     0,     0,    52,
      53,    54,    55,    56,    57,   108,     3,     4,     5,     6,
       7,     0,     8,     9,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    37,     0,    39,
       0,     0,   209,    11,    44,     0,     0,    47,   109,     0,
      49,     0,     0,     0,     0,   515,    51,     0,     0,     0,
     210,     0,     0,    52,    53,    54,    55,    56,    57,   108,
       3,     4,     5,     6,     7,     0,     8,     9,     0,   536,
       0,     0,    37,     0,    39,     0,     0,     0,     0,    44,
       0,     0,    47,   109,     0,    49,   209,    11,     0,     0,
       0,    51,     0,     0,     0,   340,     0,     0,    52,    53,
      54,    55,    56,    57,   108,     3,     4,     5,     6,     7,
       0,     8,     9,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    37,     0,    39,     0,
       0,   209,    11,    44,     0,     0,    47,   109,     0,    49,
       0,     0,     0,     0,     0,    51,     0,     0,     0,   210,
       0,     0,    52,    53,    54,    55,    56,    57,   108,     3,
       4,     5,     6,     7,     0,     8,     9,     0,     0,     0,
       0,    37,     0,    39,     0,     0,     0,     0,    44,     0,
       0,    47,   109,     0,    49,   209,    11,     0,     0,     0,
      51,     0,     0,     0,     0,     0,   254,   255,    53,    54,
      55,    56,    57,   108,     3,     4,     5,     6,     7,     0,
       8,     9,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    37,     0,    39,     0,     0,
     209,    11,    44,     0,     0,    47,   109,     0,    49,     0,
       0,     0,     0,     0,    51,     0,     0,     0,   340,     0,
       0,    52,    53,    54,    55,    56,    57,   108,     3,     4,
       5,     6,     7,     0,     8,     9,     0,     0,     0,     0,
      37,     0,    39,     0,     0,     0,     0,    44,     0,     0,
      47,   109,     0,    49,   209,    11,   426,     0,     0,    51,
       0,     0,     0,     0,     0,     0,    52,    53,    54,    55,
      56,    57,   108,     3,     4,     5,     6,     7,     0,     8,
       9,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    37,     0,    39,     0,     0,   209,
      11,    44,     0,     0,    47,   109,     0,    49,     0,     0,
     317,     0,     0,    51,     0,     0,     0,     0,     0,     0,
      52,    53,    54,    55,    56,    57,   108,     3,     4,     5,
       6,     7,     0,     8,     9,     0,     0,     0,     0,    37,
       0,    39,     0,     0,     0,     0,    44,     0,     0,    47,
     109,     0,    49,   209,    11,     0,     0,   613,    51,     0,
       0,     0,     0,     0,     0,    52,    53,    54,    55,    56,
      57,   108,     3,     4,     5,     6,     7,     0,     8,     9,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    37,     0,    39,     0,     0,   209,    11,
      44,     0,     0,    47,   109,     0,    49,     0,     0,     0,
       0,   616,    51,     0,     0,     0,     0,     0,     0,    52,
      53,    54,    55,    56,    57,   108,     3,     4,     5,     6,
       7,     0,     8,     9,     0,     0,     0,     0,    37,     0,
      39,     0,     0,     0,     0,    44,     0,     0,    47,   109,
       0,    49,     0,    11,     0,     0,     0,    51,     0,     0,
       0,     0,     0,     0,    52,    53,    54,    55,    56,    57,
     108,     3,     4,     5,     6,     7,     0,     8,     9,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    37,     0,    39,     0,     0,     0,    11,    44,
       0,     0,    47,   109,     0,    49,     0,     0,     0,     0,
       0,    51,     0,     0,     0,   210,     0,     0,    52,    53,
      54,    55,    56,    57,   108,     3,     4,     5,     6,     7,
       0,     8,     9,     0,     0,     0,     0,    37,     0,    39,
       0,     0,     0,     0,    44,     0,     0,    47,   109,     0,
      49,     0,    11,     0,     0,     0,   110,     0,     0,     0,
       0,     0,     0,    52,    53,    54,    55,    56,    57,   108,
       3,     4,     5,     6,     7,     0,     8,     9,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    37,     0,    39,     0,     0,     0,    11,    44,     0,
       0,    47,   109,     0,    49,     0,     0,     0,     0,     0,
     112,     0,     0,     0,     0,     0,     0,    52,    53,    54,
      55,    56,    57,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    37,     0,    39,     0,
       0,   127,     0,    44,     0,     0,    47,   109,     0,    49,
       0,     0,     0,     0,     0,    51,     0,     0,     0,     0,
       0,     0,    52,    53,    54,    55,    56,    57,   538,     0,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   481,   364,   482,    36,     0,    38,
       0,    40,     0,    42,    43,     0,    45,     0,     0,   129,
       0,     0,    50,     0,     0,     0,   365,   458,     0,     0,
       0,     0,   236,   189,   182,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,   408,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     364,     0,    36,     0,    38,     0,    40,     0,    42,    43,
       0,    45,     0,     0,     0,     0,     0,    50,   189,     0,
       0,   483,     0,     0,     0,     0,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    38,
       0,     0,    41,   189,    43,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     364,     0,    36,     0,    38,     0,    40,     0,    42,    43,
       0,    45,     0,     0,     0,     0,     0,    50,     0,     0,
       0,   365,   235,     0,     0,     0,     0,   236,   189,   182,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    36,     0,    38,
       0,    40,    41,    42,    43,     0,    45,     0,     0,     0,
       0,     0,    50,     0,     0,     0,     0,     0,     0,     0,
     189,   106,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   364,     0,    36,
       0,    38,     0,    40,     0,    42,    43,     0,    45,     0,
       0,   189,     0,     0,    50,     0,     0,     0,   546,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    38,     0,     0,    41,     0,    43,     0,     0,
       0,   189,     0,     0,     0,     0,     0,     0,   401,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    38,     0,     0,    41,     0,    43,     0,     0,
       0,   189,     0,     0,     0,     0,     0,     0,   410,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    38,     0,     0,    41,     0,    43,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   189,   500,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,   540,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   364,     0,    36,     0,    38,     0,
      40,     0,    42,    43,     0,    45,     0,     0,     0,     0,
     189,    50,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    36,
       0,    38,     0,    40,     0,    42,    43,   189,    45,     0,
       0,     0,     0,     0,    50,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   189,     0,    38,     0,
       0,    41,     0,    43,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,     0,     0,   448,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   189,     0,    38,     0,     0,
       0,     0,    43,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    38,     0,     0,     0,
       0,    43
};

static const yytype_int16 yycheck[] =
{
      51,    66,   126,    96,   227,     0,   236,    51,    82,   146,
     237,   292,   197,   188,    51,   146,   297,     1,     4,   393,
     136,    92,    89,    94,     3,   350,    67,    94,    86,    38,
      39,    40,    90,   157,    89,    89,    19,    96,   109,   145,
      94,   147,     0,   102,   103,   182,    86,    89,    14,    15,
      90,   109,    94,    86,     3,     3,     3,    90,   288,   110,
      86,   112,    71,    58,    90,     3,   110,    76,   112,     3,
       3,    19,    86,   110,    88,   112,    90,    35,    93,     3,
       3,    30,    30,    30,    89,   434,   192,    96,   132,   133,
     134,   197,    30,     0,     4,   132,   133,   134,   235,   473,
      58,   145,    88,   147,     3,   156,    30,    30,   145,   202,
     147,   296,    96,   206,     7,     8,     9,   180,   242,   160,
     161,   162,   447,   186,    82,   108,   200,   178,   409,    87,
     109,   126,    81,    81,    81,    93,   359,   286,   104,   105,
     135,   371,    89,    81,    93,   182,   283,    94,   192,   379,
      98,    58,    86,   197,    87,   192,    90,    81,    81,   154,
     197,   439,   157,    86,   394,   149,    89,    91,    91,    93,
     519,     3,   521,    96,    67,   238,    83,    84,    88,    87,
       3,   139,    90,    90,    91,    92,   292,   303,    87,    96,
     296,   297,    88,   151,    90,    89,   345,   155,    30,   236,
      94,   242,    96,    87,   132,   192,    90,    30,   136,    89,
     197,   285,    90,   337,   338,   210,    87,    95,   202,    90,
     294,   284,    86,   181,   287,   341,    90,   290,    89,   366,
     188,   324,   348,   460,    87,   366,     3,    90,   233,   132,
      89,   236,   200,   136,   474,   523,   283,    94,   292,    81,
      89,   288,   296,   297,    86,   292,    90,    89,    81,   296,
     297,    95,   313,    30,    96,    89,    89,   160,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,   171,   172,
     173,   174,   175,   176,   177,   460,   360,    89,   404,    91,
     439,   286,   366,   288,    88,   202,    90,    97,    88,   206,
      90,     3,    88,   409,    90,   292,   369,   106,   457,   296,
     297,    88,   356,    90,    81,   378,    88,     3,    90,   356,
     227,   458,    89,    88,    91,    90,   107,   285,    30,    96,
      88,     3,    90,     3,   371,    18,   294,     3,   375,    86,
     403,   434,   379,    87,    30,   340,   411,    98,    99,   242,
     345,   502,   310,   504,   308,    12,    13,   394,    30,    80,
      30,   398,    16,    17,    30,   409,     0,   418,   419,   420,
     298,     3,   409,    93,   523,   303,   371,   501,   527,    81,
     375,   339,   498,   434,   379,    90,    91,    89,   383,    88,
     539,    89,    94,    90,    96,    81,   418,   419,    88,   394,
     163,   164,   360,   398,   362,   298,    93,    93,   366,    81,
     303,    81,    88,   341,    88,    81,   109,   324,    91,    89,
     348,    93,   409,   497,    88,    51,    96,    93,   165,   166,
     467,   171,   172,   496,   429,    88,    88,   474,    86,    87,
      88,   478,    88,    90,   439,    93,   520,    90,   341,    88,
     487,    95,   359,    88,    87,   348,    86,   581,   386,   109,
      91,    89,   457,   417,    89,    89,    86,     9,    10,    11,
      86,    89,   467,     6,    91,   470,   404,   431,   541,   474,
      89,    93,   109,   478,   110,    88,   112,   531,    90,    88,
      95,   615,   487,   386,   531,   490,   167,   168,   169,   170,
     458,    95,   460,    95,    95,    91,   132,   133,   134,    87,
      91,   404,   549,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    62,    86,   576,    95,   434,   523,    88,
     581,   582,   527,   584,   529,    91,    91,   532,    88,   497,
     591,    94,   109,    95,   539,    88,    95,    89,    95,   503,
      92,   505,    94,    95,   549,    95,    88,   552,    95,    95,
     514,   515,   520,    88,   615,     9,    10,    11,    95,    88,
     498,    89,    88,   501,   502,    86,   504,    91,     3,     4,
       5,     6,     7,    95,    91,    95,    93,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    88,    95,    95,
       3,    95,   109,    95,    60,   498,    31,    95,   501,   502,
      95,   504,    88,    86,    88,   531,   137,   155,   429,   573,
     574,   575,   383,   577,   175,   173,   580,   285,   174,   497,
     217,   176,   463,   177,   340,   346,   324,    58,    98,   615,
      -1,    -1,    -1,    -1,    -1,    89,    -1,    72,    92,    -1,
      94,    -1,    -1,    -1,    -1,    80,    81,    -1,    83,   613,
      -1,    -1,   616,    -1,    89,    -1,    -1,   621,    -1,   623,
      -1,   625,    -1,    -1,   628,    -1,   109,    -1,     0,     1,
     634,     3,     4,     5,     6,     7,     8,    -1,    10,    11,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    51,
      52,    53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    69,    70,    71,
      72,    73,    74,    75,    76,    77,    78,    79,    80,    81,
      -1,    83,    84,    -1,    -1,    -1,    -1,    89,    -1,    -1,
      -1,    -1,    -1,    -1,    96,    97,    98,    99,   100,   101,
       1,    -1,     3,     4,     5,     6,     7,     8,    -1,    10,
      11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      51,    52,    53,    54,    55,    -1,    57,    58,    59,    -1,
      61,    62,    63,    64,    65,    66,    67,    68,    69,    70,
      71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
      81,    82,    83,    84,    -1,    86,    87,    -1,    89,    -1,
      -1,    -1,    93,    -1,    -1,    96,    97,    98,    99,   100,
     101,     1,    -1,     3,     4,     5,     6,     7,     8,    -1,
      10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    -1,    57,    58,    59,
      -1,    61,    62,    63,    64,    65,    66,    67,    68,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    78,    79,
      80,    81,    82,    83,    84,    -1,    86,    -1,    -1,    89,
      -1,    -1,    -1,    93,    -1,    -1,    96,    97,    98,    99,
     100,   101,     1,    -1,     3,     4,     5,     6,     7,     8,
      -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,    55,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      69,    70,    71,    72,    73,    74,    75,    76,    77,    78,
      79,    80,    81,    -1,    83,    84,    -1,    -1,    -1,    -1,
      89,    -1,    -1,    -1,    -1,    -1,    -1,    96,    97,    98,
      99,   100,   101,     1,    -1,     3,     4,     5,     6,     7,
       8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    30,    31,    -1,    -1,    -1,    -1,    -1,    -1,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,     1,    -1,
       3,     4,     5,     6,     7,     8,    -1,    10,    11,    -1,
      -1,    -1,    70,    71,    72,    -1,    -1,    -1,    76,    77,
      -1,    -1,    80,    81,    -1,    83,    -1,    30,    31,    -1,
      -1,    89,    -1,    -1,    -1,    -1,    -1,    -1,    96,    97,
      98,    99,   100,   101,     3,     4,     5,     6,     7,     8,
      -1,    10,    11,    -1,    57,    58,    59,    -1,    61,    62,
      63,    64,    65,    66,    67,    68,    -1,    70,    -1,    72,
      -1,    30,    31,    -1,    77,    -1,    79,    80,    81,    82,
      83,    -1,    -1,    86,    -1,    -1,    89,    -1,    -1,    -1,
      93,    -1,    -1,    96,    97,    98,    99,   100,   101,     1,
      -1,     3,     4,     5,     6,     7,     8,     1,    10,    11,
      -1,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,    -1,
      -1,    80,    81,    -1,    83,    -1,    -1,    -1,    30,    31,
      89,    -1,    -1,    92,    93,    94,    30,    96,    97,    98,
      99,   100,   101,    -1,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    51,    52,    53,
      54,    55,    -1,    -1,    -1,    -1,    -1,    -1,    70,    -1,
      72,    -1,    -1,    -1,    -1,    77,    -1,    71,    80,    81,
      74,    83,    76,    -1,    -1,    -1,    -1,    89,    -1,    -1,
      -1,    -1,    -1,    -1,    96,    97,    98,    99,   100,   101,
       3,     4,     5,     6,     7,     8,    -1,    10,    11,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    69,    70,    71,    72,
      73,    74,    75,    76,    77,    78,    -1,    80,    81,    -1,
      83,    84,    -1,    86,    -1,    -1,    89,    -1,    -1,    -1,
      -1,    -1,    -1,    96,    97,    98,    99,   100,   101,     3,
       4,     5,     6,     7,     8,    -1,    10,    11,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    30,    31,    -1,    -1,
      -1,    -1,    -1,    -1,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    51,    52,    53,
      54,    55,    -1,     3,     4,     5,     6,     7,     8,    -1,
      10,    11,    -1,    -1,    -1,    -1,    70,    71,    72,    -1,
      -1,    -1,    76,    77,    -1,    -1,    80,    81,    -1,    83,
      30,    31,    -1,    -1,    34,    89,    -1,    -1,    38,    39,
      40,    -1,    96,    97,    98,    99,   100,   101,    -1,    -1,
      -1,     3,     4,     5,     6,     7,     8,    -1,    10,    11,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      70,    71,    72,    -1,    -1,    -1,    76,    77,    30,    31,
      80,    81,    34,    83,    -1,    -1,    38,    39,    40,    89,
      -1,    -1,    -1,    -1,    -1,    95,    96,    97,    98,    99,
     100,   101,    -1,    -1,    -1,     3,     4,     5,     6,     7,
       8,    -1,    10,    11,    -1,    -1,    -1,    -1,    70,    71,
      72,    -1,    -1,    -1,    76,    77,    -1,    -1,    80,    81,
      -1,    83,    30,    31,    -1,    -1,    34,    89,    -1,    -1,
      38,    39,    40,    95,    96,    97,    98,    99,   100,   101,
      -1,    -1,    -1,     3,     4,     5,     6,     7,     8,    -1,
      10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    70,    71,    72,    -1,    -1,    -1,    76,    77,
      30,    31,    80,    81,    34,    83,    -1,    -1,    38,    39,
      40,    89,    -1,    -1,    -1,    -1,    -1,    95,    96,    97,
      98,    99,   100,   101,    -1,    -1,    -1,     3,     4,     5,
       6,     7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,
      70,    71,    72,    -1,    -1,    -1,    76,    77,    -1,    -1,
      80,    81,    -1,    83,    30,    31,    -1,    -1,    34,    89,
      -1,    -1,    38,    39,    40,    95,    96,    97,    98,    99,
     100,   101,    -1,    -1,    -1,     3,     4,     5,     6,     7,
       8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    70,    71,    72,    -1,    -1,    -1,
      76,    77,    30,    31,    80,    81,    34,    83,    -1,    -1,
      38,    39,    40,    89,    -1,    -1,    -1,    -1,    -1,    95,
      96,    97,    98,    99,   100,   101,    -1,    -1,    -1,     3,
       4,     5,     6,     7,     8,    -1,    10,    11,    -1,    -1,
      -1,    -1,    70,    71,    72,    -1,    -1,    -1,    76,    77,
      -1,    -1,    80,    81,    -1,    83,    30,    31,    -1,    -1,
      -1,    89,    -1,    -1,    38,    39,    40,    95,    96,    97,
      98,    99,   100,   101,    -1,    -1,    -1,     3,     4,     5,
       6,     7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    70,    71,    72,    -1,
      -1,    -1,    76,    77,    30,    31,    80,    81,    -1,    83,
      -1,    -1,    -1,    -1,    -1,    89,    -1,    -1,    -1,    -1,
      -1,    -1,    96,    97,    98,    99,   100,   101,    -1,    -1,
      -1,     3,     4,     5,     6,     7,     8,    -1,    10,    11,
      -1,    -1,    -1,    -1,    70,    -1,    72,    -1,    -1,    -1,
      -1,    77,    -1,    -1,    80,    81,    -1,    83,    30,    31,
      -1,    87,    -1,    89,    -1,    -1,    92,    93,    94,    -1,
      96,    97,    98,    99,   100,   101,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    70,    -1,
      72,    -1,    -1,    30,    31,    77,    -1,    -1,    80,    81,
      -1,    83,    -1,    -1,    -1,    87,    -1,    89,    -1,    -1,
      92,    93,    94,    -1,    96,    97,    98,    99,   100,   101,
       3,     4,     5,     6,     7,     8,    -1,    10,    11,    -1,
      -1,    -1,    -1,    70,    -1,    72,    -1,    -1,    -1,    -1,
      77,    -1,    -1,    80,    81,    -1,    83,    30,    31,    -1,
      -1,    88,    89,    -1,    -1,    -1,    93,    -1,    -1,    96,
      97,    98,    99,   100,   101,     3,     4,     5,     6,     7,
       8,    -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    70,    -1,    72,
      -1,    -1,    30,    31,    77,    -1,    -1,    80,    81,    -1,
      83,    -1,    -1,    -1,    -1,    88,    89,    -1,    -1,    -1,
      93,    -1,    -1,    96,    97,    98,    99,   100,   101,     3,
       4,     5,     6,     7,     8,    -1,    10,    11,    -1,    67,
      -1,    -1,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,
      -1,    -1,    80,    81,    -1,    83,    30,    31,    -1,    -1,
      -1,    89,    -1,    -1,    -1,    93,    -1,    -1,    96,    97,
      98,    99,   100,   101,     3,     4,     5,     6,     7,     8,
      -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    70,    -1,    72,    -1,
      -1,    30,    31,    77,    -1,    -1,    80,    81,    -1,    83,
      -1,    -1,    -1,    -1,    -1,    89,    -1,    -1,    -1,    93,
      -1,    -1,    96,    97,    98,    99,   100,   101,     3,     4,
       5,     6,     7,     8,    -1,    10,    11,    -1,    -1,    -1,
      -1,    70,    -1,    72,    -1,    -1,    -1,    -1,    77,    -1,
      -1,    80,    81,    -1,    83,    30,    31,    -1,    -1,    -1,
      89,    -1,    -1,    -1,    -1,    -1,    95,    96,    97,    98,
      99,   100,   101,     3,     4,     5,     6,     7,     8,    -1,
      10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    70,    -1,    72,    -1,    -1,
      30,    31,    77,    -1,    -1,    80,    81,    -1,    83,    -1,
      -1,    -1,    -1,    -1,    89,    -1,    -1,    -1,    93,    -1,
      -1,    96,    97,    98,    99,   100,   101,     3,     4,     5,
       6,     7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,
      70,    -1,    72,    -1,    -1,    -1,    -1,    77,    -1,    -1,
      80,    81,    -1,    83,    30,    31,    86,    -1,    -1,    89,
      -1,    -1,    -1,    -1,    -1,    -1,    96,    97,    98,    99,
     100,   101,     3,     4,     5,     6,     7,     8,    -1,    10,
      11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    70,    -1,    72,    -1,    -1,    30,
      31,    77,    -1,    -1,    80,    81,    -1,    83,    -1,    -1,
      86,    -1,    -1,    89,    -1,    -1,    -1,    -1,    -1,    -1,
      96,    97,    98,    99,   100,   101,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    -1,    -1,    -1,    -1,    70,
      -1,    72,    -1,    -1,    -1,    -1,    77,    -1,    -1,    80,
      81,    -1,    83,    30,    31,    -1,    -1,    88,    89,    -1,
      -1,    -1,    -1,    -1,    -1,    96,    97,    98,    99,   100,
     101,     3,     4,     5,     6,     7,     8,    -1,    10,    11,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    70,    -1,    72,    -1,    -1,    30,    31,
      77,    -1,    -1,    80,    81,    -1,    83,    -1,    -1,    -1,
      -1,    88,    89,    -1,    -1,    -1,    -1,    -1,    -1,    96,
      97,    98,    99,   100,   101,     3,     4,     5,     6,     7,
       8,    -1,    10,    11,    -1,    -1,    -1,    -1,    70,    -1,
      72,    -1,    -1,    -1,    -1,    77,    -1,    -1,    80,    81,
      -1,    83,    -1,    31,    -1,    -1,    -1,    89,    -1,    -1,
      -1,    -1,    -1,    -1,    96,    97,    98,    99,   100,   101,
       3,     4,     5,     6,     7,     8,    -1,    10,    11,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    70,    -1,    72,    -1,    -1,    -1,    31,    77,
      -1,    -1,    80,    81,    -1,    83,    -1,    -1,    -1,    -1,
      -1,    89,    -1,    -1,    -1,    93,    -1,    -1,    96,    97,
      98,    99,   100,   101,     3,     4,     5,     6,     7,     8,
      -1,    10,    11,    -1,    -1,    -1,    -1,    70,    -1,    72,
      -1,    -1,    -1,    -1,    77,    -1,    -1,    80,    81,    -1,
      83,    -1,    31,    -1,    -1,    -1,    89,    -1,    -1,    -1,
      -1,    -1,    -1,    96,    97,    98,    99,   100,   101,     3,
       4,     5,     6,     7,     8,    -1,    10,    11,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    70,    -1,    72,    -1,    -1,    -1,    31,    77,    -1,
      -1,    80,    81,    -1,    83,    -1,    -1,    -1,    -1,    -1,
      89,    -1,    -1,    -1,    -1,    -1,    -1,    96,    97,    98,
      99,   100,   101,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    70,    -1,    72,    -1,
      -1,     3,    -1,    77,    -1,    -1,    80,    81,    -1,    83,
      -1,    -1,    -1,    -1,    -1,    89,    -1,    -1,    -1,    -1,
      -1,    -1,    96,    97,    98,    99,   100,   101,    30,    -1,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    51,
      52,    53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,     1,    67,     3,    69,    -1,    71,
      -1,    73,    -1,    75,    76,    -1,    78,    -1,    -1,    81,
      -1,    -1,    84,    -1,    -1,    -1,    88,    89,    -1,    -1,
      -1,    -1,    94,    30,    96,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    55,     1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      67,    -1,    69,    -1,    71,    -1,    73,    -1,    75,    76,
      -1,    78,    -1,    -1,    -1,    -1,    -1,    84,    30,    -1,
      -1,    88,    -1,    -1,    -1,    -1,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    51,
      52,    53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    71,
      -1,    -1,    74,    30,    76,    32,    33,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,    52,    53,    54,    55,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      67,    -1,    69,    -1,    71,    -1,    73,    -1,    75,    76,
      -1,    78,    -1,    -1,    -1,    -1,    -1,    84,    -1,    -1,
      -1,    88,    89,    -1,    -1,    -1,    -1,    94,    30,    96,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    51,
      52,    53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    69,    -1,    71,
      -1,    73,    74,    75,    76,    -1,    78,    -1,    -1,    -1,
      -1,    -1,    84,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      30,    93,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    67,    -1,    69,
      -1,    71,    -1,    73,    -1,    75,    76,    -1,    78,    -1,
      -1,    30,    -1,    -1,    84,    -1,    -1,    -1,    88,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,    55,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    71,    -1,    -1,    74,    -1,    76,    -1,    -1,
      -1,    30,    -1,    -1,    -1,    -1,    -1,    -1,    87,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,    55,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    71,    -1,    -1,    74,    -1,    76,    -1,    -1,
      -1,    30,    -1,    -1,    -1,    -1,    -1,    -1,    87,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    54,    55,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    71,    -1,    -1,    74,    -1,    76,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,    87,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    56,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    67,    -1,    69,    -1,    71,    -1,
      73,    -1,    75,    76,    -1,    78,    -1,    -1,    -1,    -1,
      30,    84,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    69,
      -1,    71,    -1,    73,    -1,    75,    76,    30,    78,    -1,
      -1,    -1,    -1,    -1,    84,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,    71,    -1,
      -1,    74,    -1,    76,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    50,    51,    52,    53,
      54,    55,    -1,    -1,    58,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    30,    -1,    71,    -1,    -1,
      -1,    -1,    76,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    71,    -1,    -1,    -1,
      -1,    76
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     5,     6,     7,     8,    10,    11,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,    52,    53,    54,    55,    69,    70,    71,    72,
      73,    74,    75,    76,    77,    78,    79,    80,    81,    83,
      84,    89,    96,    97,    98,    99,   100,   101,   111,   112,
     114,   116,   117,   120,   122,   123,   127,   128,   129,   130,
     131,   132,   133,   134,   135,   136,   137,   138,   139,   140,
     141,   145,   146,   149,   150,   151,   152,   153,   159,   162,
     163,   164,   165,   168,   184,   197,   198,   199,   203,   204,
     211,   214,   215,   216,   217,   218,    93,   187,     3,    81,
      89,   127,    89,   127,   127,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,   109,   142,     3,    30,    81,
      93,   168,    89,    89,    89,    89,    89,    94,    89,    89,
       1,    30,   113,   141,   143,   150,   156,   163,   175,     1,
     197,     9,    10,    11,    89,    92,    94,   142,   127,   129,
      96,   102,   103,    98,    99,    12,    13,    14,    15,   104,
     105,    16,    17,    97,   106,   107,    18,    19,   108,    86,
      86,    89,    96,   147,   148,   166,   167,   168,   169,    30,
     146,   146,    93,   168,   146,   146,   146,    93,   168,   145,
     146,   187,   200,   204,     0,    87,   188,   175,   175,    30,
      93,   121,   141,     3,   115,   160,   161,    93,   140,   144,
     175,   175,   175,   141,   144,   120,     4,    88,    19,    98,
     168,   205,    88,    90,   156,    89,    94,   169,   176,   177,
     178,   156,    88,    86,    87,    88,   124,   168,    88,   121,
     125,   126,   141,   124,    95,    96,   143,   121,   141,   129,
     129,   129,   130,   130,   131,   131,   132,   132,   132,   132,
     133,   133,   134,   135,   136,   137,   138,   143,   178,   166,
      71,   163,   169,   170,    86,    90,   109,    89,    94,   178,
     167,     1,   154,   155,   156,   184,    93,   154,    89,   166,
     145,   187,     1,    57,    58,    59,    61,    62,    63,    64,
      65,    66,    67,    68,    79,    81,    82,    86,   143,   145,
     168,   185,   186,   187,   189,   190,   191,   192,   193,   194,
     195,   196,   201,   202,   207,   208,   209,    88,    88,    92,
      93,    94,   141,   179,   180,   181,   182,   183,   109,    87,
      90,   160,    88,    88,    88,    88,    90,    90,    95,    88,
     146,   199,    91,   141,    67,    88,   146,   171,   172,   173,
     176,    34,    95,    96,   141,   170,   177,   178,    89,    94,
     121,   129,    88,    90,    95,    95,    91,    88,   163,   169,
     178,   148,   179,   178,    34,    95,    96,   141,   170,   178,
      87,    87,   155,    86,    91,   157,   158,   166,     1,   154,
      87,   127,   140,   212,   213,    86,   144,    91,    89,    89,
      89,    30,   185,   168,    86,    86,    86,   143,    89,    89,
      86,    91,    87,   190,    89,   168,   180,   144,    87,    90,
     179,   109,   183,   144,    87,   161,    87,    90,    58,   118,
     119,   175,     6,   199,   166,   168,   206,   109,    89,   166,
     169,   176,    88,    90,   178,    88,   141,   170,    95,    95,
      34,    95,   141,   178,    34,    95,    96,   141,   170,   126,
     140,     1,     3,    88,   171,   174,   141,   170,    95,    95,
      34,    95,    96,   141,   178,   144,    86,    90,    91,    87,
      87,   142,    86,    88,    90,    91,   185,   113,   113,   143,
      62,    86,    86,     4,    88,    88,   125,   185,   143,   145,
     146,   191,    87,    90,    95,    87,   179,   181,    87,    91,
      88,    90,    91,    88,    94,    88,    67,   179,    30,   109,
      56,   173,    95,   141,   141,    95,    88,   171,   141,   170,
      95,    95,    34,    95,   141,    88,    88,    88,    90,    95,
     141,   141,    95,    95,   178,   158,   144,   121,   140,   213,
     185,   213,   185,    88,    88,    88,    89,    88,   185,   185,
      88,     3,   191,   166,   191,    87,   179,   141,   119,   141,
      86,    91,   179,   178,    95,    95,    88,    95,   141,   141,
      95,     3,    95,    95,   185,   185,   185,   143,   185,   185,
     121,   143,   210,    88,   143,     3,    88,   143,   143,    95,
      95,    60,    88,    88,   185,    88,   210,   185,    88,    95,
     185,    86,   185,   185,    88,   185,   185
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_uint8 yyr1[] =
{
       0,   110,   111,   111,   111,   111,   111,   112,   112,   112,
     112,   112,   112,   112,   112,   113,   113,   114,   114,   114,
     115,   116,   116,   117,   118,   118,   119,   119,   120,   120,
     120,   120,   120,   120,   120,   120,   121,   121,   122,   122,
     122,   123,   123,   124,   125,   125,   125,   126,   126,   127,
     127,   127,   127,   127,   127,   127,   127,   128,   128,   128,
     128,   128,   128,   129,   129,   130,   130,   130,   130,   131,
     131,   131,   132,   132,   132,   133,   133,   133,   133,   133,
     134,   134,   134,   135,   135,   136,   136,   137,   137,   138,
     138,   139,   139,   140,   140,   141,   141,   141,   141,   141,
     142,   142,   142,   142,   142,   142,   142,   142,   142,   142,
     142,   143,   143,   144,   145,   145,   145,   146,   146,   146,
     146,   146,   146,   146,   146,   146,   146,   147,   147,   148,
     148,   149,   149,   149,   149,   149,   149,   149,   149,   150,
     151,   151,   151,   151,   151,   151,   151,   151,   151,   151,
     151,   151,   151,   151,   151,   151,   152,   152,   152,   152,
     152,   153,   153,   154,   154,   155,   155,   155,   156,   156,
     156,   156,   157,   157,   158,   158,   158,   159,   159,   159,
     159,   159,   160,   160,   161,   161,   162,   163,   163,   163,
     163,   163,   164,   164,   165,   165,   166,   166,   167,   167,
     167,   167,   167,   167,   167,   167,   167,   167,   167,   167,
     167,   167,   167,   168,   168,   168,   169,   169,   169,   169,
     170,   170,   171,   171,   172,   172,   173,   173,   173,   173,
     173,   173,   174,   174,   175,   175,   176,   176,   176,   177,
     177,   177,   177,   177,   177,   177,   177,   177,   177,   177,
     177,   177,   177,   177,   177,   177,   177,   177,   177,   177,
     178,   179,   179,   179,   180,   180,   180,   180,   181,   182,
     182,   183,   183,   184,   185,   185,   185,   185,   185,   185,
     185,   185,   186,   186,   186,   187,   188,   187,   189,   189,
     190,   190,   191,   191,   192,   192,   192,   193,   194,   194,
     194,   194,   194,   195,   195,   196,   196,   196,   196,   196,
     197,   197,   197,   197,   197,   197,   197,   198,   199,   199,
     200,   200,   201,   201,   201,   201,   201,   202,   202,   203,
     203,   204,   205,   205,   205,   206,   206,   207,   207,   208,
     209,   210,   210,   211,   212,   212,   212,   213,   213,   213,
     214,   215,   215,   216,   217,   217,   218
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     1,     2,     3,     3,     3,     1,     1,     1,
       1,     1,     3,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     6,     1,     3,     3,     3,     1,     1,
       1,     3,     3,     2,     2,     4,     3,     4,     3,     4,
       4,     3,     4,     1,     1,     3,     2,     1,     1,     1,
       2,     2,     2,     2,     4,     4,     1,     1,     1,     1,
       1,     1,     1,     1,     4,     1,     3,     3,     3,     1,
       3,     3,     1,     3,     3,     1,     3,     3,     3,     3,
       1,     3,     3,     1,     3,     1,     3,     1,     3,     1,
       3,     1,     3,     1,     5,     1,     3,     3,     3,     3,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     3,     1,     3,     4,     1,     2,     1,     2,
       1,     2,     1,     2,     1,     2,     1,     1,     3,     3,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     4,     5,     2,     4,
       5,     1,     1,     1,     2,     3,     4,     1,     2,     1,
       2,     1,     1,     3,     2,     3,     1,     4,     5,     5,
       6,     2,     1,     3,     3,     1,     4,     1,     1,     1,
       1,     1,     1,     1,     4,     4,     3,     2,     1,     3,
       3,     4,     6,     5,     5,     6,     5,     4,     4,     5,
       5,     4,     5,     1,     1,     1,     3,     2,     2,     1,
       1,     2,     3,     1,     2,     4,     2,     4,     2,     1,
       3,     3,     1,     3,     3,     2,     2,     1,     1,     3,
       2,     3,     5,     4,     5,     4,     3,     3,     3,     4,
       6,     5,     5,     6,     4,     4,     2,     3,     4,     5,
       0,     3,     4,     1,     2,     1,     4,     3,     2,     1,
       2,     3,     2,     7,     1,     1,     1,     1,     1,     1,
       1,     2,     3,     4,     3,     2,     0,     4,     1,     2,
       1,     1,     1,     2,     7,     5,     5,     1,     5,     7,
       6,     7,     1,     6,     7,     3,     2,     2,     2,     3,
       1,     1,     1,     1,     1,     1,     2,     2,     3,     2,
       1,     2,     1,     1,     1,     1,     1,     4,     5,     1,
       2,     6,     1,     1,     1,     1,     5,     4,     5,     8,
       7,     1,     1,     6,     1,     3,     3,     1,     3,     3,
       2,     4,     5,     4,     1,     4,     1
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


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
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yykind)
    {
    case YYSYMBOL_IDENTIFIER: /* IDENTIFIER  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 1966 "basilisk.c"
        break;

    case YYSYMBOL_I_CONSTANT: /* I_CONSTANT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 1975 "basilisk.c"
        break;

    case YYSYMBOL_F_CONSTANT: /* F_CONSTANT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 1984 "basilisk.c"
        break;

    case YYSYMBOL_STRING_LITERAL: /* STRING_LITERAL  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 1993 "basilisk.c"
        break;

    case YYSYMBOL_FUNC_NAME: /* FUNC_NAME  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2002 "basilisk.c"
        break;

    case YYSYMBOL_SIZEOF: /* SIZEOF  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2011 "basilisk.c"
        break;

    case YYSYMBOL_PTR_OP: /* PTR_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2020 "basilisk.c"
        break;

    case YYSYMBOL_INC_OP: /* INC_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2029 "basilisk.c"
        break;

    case YYSYMBOL_DEC_OP: /* DEC_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2038 "basilisk.c"
        break;

    case YYSYMBOL_LEFT_OP: /* LEFT_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2047 "basilisk.c"
        break;

    case YYSYMBOL_RIGHT_OP: /* RIGHT_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2056 "basilisk.c"
        break;

    case YYSYMBOL_LE_OP: /* LE_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2065 "basilisk.c"
        break;

    case YYSYMBOL_GE_OP: /* GE_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2074 "basilisk.c"
        break;

    case YYSYMBOL_EQ_OP: /* EQ_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2083 "basilisk.c"
        break;

    case YYSYMBOL_NE_OP: /* NE_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2092 "basilisk.c"
        break;

    case YYSYMBOL_AND_OP: /* AND_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2101 "basilisk.c"
        break;

    case YYSYMBOL_OR_OP: /* OR_OP  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2110 "basilisk.c"
        break;

    case YYSYMBOL_MUL_ASSIGN: /* MUL_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2119 "basilisk.c"
        break;

    case YYSYMBOL_DIV_ASSIGN: /* DIV_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2128 "basilisk.c"
        break;

    case YYSYMBOL_MOD_ASSIGN: /* MOD_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2137 "basilisk.c"
        break;

    case YYSYMBOL_ADD_ASSIGN: /* ADD_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2146 "basilisk.c"
        break;

    case YYSYMBOL_SUB_ASSIGN: /* SUB_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2155 "basilisk.c"
        break;

    case YYSYMBOL_LEFT_ASSIGN: /* LEFT_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2164 "basilisk.c"
        break;

    case YYSYMBOL_RIGHT_ASSIGN: /* RIGHT_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2173 "basilisk.c"
        break;

    case YYSYMBOL_AND_ASSIGN: /* AND_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2182 "basilisk.c"
        break;

    case YYSYMBOL_XOR_ASSIGN: /* XOR_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2191 "basilisk.c"
        break;

    case YYSYMBOL_OR_ASSIGN: /* OR_ASSIGN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2200 "basilisk.c"
        break;

    case YYSYMBOL_TYPEDEF_NAME: /* TYPEDEF_NAME  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2209 "basilisk.c"
        break;

    case YYSYMBOL_ENUMERATION_CONSTANT: /* ENUMERATION_CONSTANT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2218 "basilisk.c"
        break;

    case YYSYMBOL_TYPEDEF: /* TYPEDEF  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2227 "basilisk.c"
        break;

    case YYSYMBOL_EXTERN: /* EXTERN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2236 "basilisk.c"
        break;

    case YYSYMBOL_STATIC: /* STATIC  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2245 "basilisk.c"
        break;

    case YYSYMBOL_AUTO: /* AUTO  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2254 "basilisk.c"
        break;

    case YYSYMBOL_REGISTER: /* REGISTER  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2263 "basilisk.c"
        break;

    case YYSYMBOL_INLINE: /* INLINE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2272 "basilisk.c"
        break;

    case YYSYMBOL_CONST: /* CONST  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2281 "basilisk.c"
        break;

    case YYSYMBOL_RESTRICT: /* RESTRICT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2290 "basilisk.c"
        break;

    case YYSYMBOL_VOLATILE: /* VOLATILE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2299 "basilisk.c"
        break;

    case YYSYMBOL_BOOL: /* BOOL  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2308 "basilisk.c"
        break;

    case YYSYMBOL_CHAR: /* CHAR  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2317 "basilisk.c"
        break;

    case YYSYMBOL_SHORT: /* SHORT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2326 "basilisk.c"
        break;

    case YYSYMBOL_INT: /* INT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2335 "basilisk.c"
        break;

    case YYSYMBOL_LONG: /* LONG  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2344 "basilisk.c"
        break;

    case YYSYMBOL_SIGNED: /* SIGNED  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2353 "basilisk.c"
        break;

    case YYSYMBOL_UNSIGNED: /* UNSIGNED  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2362 "basilisk.c"
        break;

    case YYSYMBOL_FLOAT: /* FLOAT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2371 "basilisk.c"
        break;

    case YYSYMBOL_DOUBLE: /* DOUBLE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2380 "basilisk.c"
        break;

    case YYSYMBOL_VOID: /* VOID  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2389 "basilisk.c"
        break;

    case YYSYMBOL_COMPLEX: /* COMPLEX  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2398 "basilisk.c"
        break;

    case YYSYMBOL_IMAGINARY: /* IMAGINARY  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2407 "basilisk.c"
        break;

    case YYSYMBOL_STRUCT: /* STRUCT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2416 "basilisk.c"
        break;

    case YYSYMBOL_UNION: /* UNION  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2425 "basilisk.c"
        break;

    case YYSYMBOL_ENUM: /* ENUM  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2434 "basilisk.c"
        break;

    case YYSYMBOL_ELLIPSIS: /* ELLIPSIS  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2443 "basilisk.c"
        break;

    case YYSYMBOL_CASE: /* CASE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2452 "basilisk.c"
        break;

    case YYSYMBOL_DEFAULT: /* DEFAULT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2461 "basilisk.c"
        break;

    case YYSYMBOL_IF: /* IF  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2470 "basilisk.c"
        break;

    case YYSYMBOL_ELSE: /* ELSE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2479 "basilisk.c"
        break;

    case YYSYMBOL_SWITCH: /* SWITCH  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2488 "basilisk.c"
        break;

    case YYSYMBOL_WHILE: /* WHILE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2497 "basilisk.c"
        break;

    case YYSYMBOL_DO: /* DO  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2506 "basilisk.c"
        break;

    case YYSYMBOL_FOR: /* FOR  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2515 "basilisk.c"
        break;

    case YYSYMBOL_GOTO: /* GOTO  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2524 "basilisk.c"
        break;

    case YYSYMBOL_CONTINUE: /* CONTINUE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2533 "basilisk.c"
        break;

    case YYSYMBOL_BREAK: /* BREAK  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2542 "basilisk.c"
        break;

    case YYSYMBOL_RETURN: /* RETURN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2551 "basilisk.c"
        break;

    case YYSYMBOL_ALIGNAS: /* ALIGNAS  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2560 "basilisk.c"
        break;

    case YYSYMBOL_ALIGNOF: /* ALIGNOF  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2569 "basilisk.c"
        break;

    case YYSYMBOL_ATOMIC: /* ATOMIC  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2578 "basilisk.c"
        break;

    case YYSYMBOL_GENERIC: /* GENERIC  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2587 "basilisk.c"
        break;

    case YYSYMBOL_NORETURN: /* NORETURN  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2596 "basilisk.c"
        break;

    case YYSYMBOL_STATIC_ASSERT: /* STATIC_ASSERT  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2605 "basilisk.c"
        break;

    case YYSYMBOL_THREAD_LOCAL: /* THREAD_LOCAL  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2614 "basilisk.c"
        break;

    case YYSYMBOL_MAYBECONST: /* MAYBECONST  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2623 "basilisk.c"
        break;

    case YYSYMBOL_NEW_FIELD: /* NEW_FIELD  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2632 "basilisk.c"
        break;

    case YYSYMBOL_TRACE: /* TRACE  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2641 "basilisk.c"
        break;

    case YYSYMBOL_FOREACH_DIMENSION: /* FOREACH_DIMENSION  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2650 "basilisk.c"
        break;

    case YYSYMBOL_REDUCTION: /* REDUCTION  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2659 "basilisk.c"
        break;

    case YYSYMBOL_MACRO: /* MACRO  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2668 "basilisk.c"
        break;

    case YYSYMBOL_ELLIPSIS_MACRO: /* ELLIPSIS_MACRO  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2677 "basilisk.c"
        break;

    case YYSYMBOL_ELLIPSIS_EXPRESSION: /* ELLIPSIS_EXPRESSION  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2686 "basilisk.c"
        break;

    case YYSYMBOL_MACRODEF: /* MACRODEF  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2695 "basilisk.c"
        break;

    case YYSYMBOL_foreach_statement: /* foreach_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2704 "basilisk.c"
        break;

    case YYSYMBOL_86_: /* ';'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2713 "basilisk.c"
        break;

    case YYSYMBOL_87_: /* '}'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2722 "basilisk.c"
        break;

    case YYSYMBOL_88_: /* ')'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2731 "basilisk.c"
        break;

    case YYSYMBOL_89_: /* '('  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2740 "basilisk.c"
        break;

    case YYSYMBOL_90_: /* ','  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2749 "basilisk.c"
        break;

    case YYSYMBOL_91_: /* ':'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2758 "basilisk.c"
        break;

    case YYSYMBOL_92_: /* '.'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2767 "basilisk.c"
        break;

    case YYSYMBOL_93_: /* '{'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2776 "basilisk.c"
        break;

    case YYSYMBOL_94_: /* '['  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2785 "basilisk.c"
        break;

    case YYSYMBOL_95_: /* ']'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2794 "basilisk.c"
        break;

    case YYSYMBOL_96_: /* '*'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2803 "basilisk.c"
        break;

    case YYSYMBOL_97_: /* '&'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2812 "basilisk.c"
        break;

    case YYSYMBOL_98_: /* '+'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2821 "basilisk.c"
        break;

    case YYSYMBOL_99_: /* '-'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2830 "basilisk.c"
        break;

    case YYSYMBOL_100_: /* '~'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2839 "basilisk.c"
        break;

    case YYSYMBOL_101_: /* '!'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2848 "basilisk.c"
        break;

    case YYSYMBOL_102_: /* '/'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2857 "basilisk.c"
        break;

    case YYSYMBOL_103_: /* '%'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2866 "basilisk.c"
        break;

    case YYSYMBOL_104_: /* '<'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2875 "basilisk.c"
        break;

    case YYSYMBOL_105_: /* '>'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2884 "basilisk.c"
        break;

    case YYSYMBOL_106_: /* '^'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2893 "basilisk.c"
        break;

    case YYSYMBOL_107_: /* '|'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2902 "basilisk.c"
        break;

    case YYSYMBOL_108_: /* '?'  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2911 "basilisk.c"
        break;

    case YYSYMBOL_109_: /* '='  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2920 "basilisk.c"
        break;

    case YYSYMBOL_translation_unit: /* translation_unit  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2929 "basilisk.c"
        break;

    case YYSYMBOL_primary_expression: /* primary_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2938 "basilisk.c"
        break;

    case YYSYMBOL_expression_error: /* expression_error  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2947 "basilisk.c"
        break;

    case YYSYMBOL_constant: /* constant  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2956 "basilisk.c"
        break;

    case YYSYMBOL_enumeration_constant: /* enumeration_constant  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2965 "basilisk.c"
        break;

    case YYSYMBOL_string: /* string  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2974 "basilisk.c"
        break;

    case YYSYMBOL_generic_selection: /* generic_selection  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2983 "basilisk.c"
        break;

    case YYSYMBOL_generic_assoc_list: /* generic_assoc_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 2992 "basilisk.c"
        break;

    case YYSYMBOL_generic_association: /* generic_association  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3001 "basilisk.c"
        break;

    case YYSYMBOL_postfix_expression: /* postfix_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3010 "basilisk.c"
        break;

    case YYSYMBOL_postfix_initializer: /* postfix_initializer  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3019 "basilisk.c"
        break;

    case YYSYMBOL_array_access: /* array_access  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3028 "basilisk.c"
        break;

    case YYSYMBOL_function_call: /* function_call  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3037 "basilisk.c"
        break;

    case YYSYMBOL_member_identifier: /* member_identifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3046 "basilisk.c"
        break;

    case YYSYMBOL_argument_expression_list: /* argument_expression_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3055 "basilisk.c"
        break;

    case YYSYMBOL_argument_expression_list_item: /* argument_expression_list_item  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3064 "basilisk.c"
        break;

    case YYSYMBOL_unary_expression: /* unary_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3073 "basilisk.c"
        break;

    case YYSYMBOL_unary_operator: /* unary_operator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3082 "basilisk.c"
        break;

    case YYSYMBOL_cast_expression: /* cast_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3091 "basilisk.c"
        break;

    case YYSYMBOL_multiplicative_expression: /* multiplicative_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3100 "basilisk.c"
        break;

    case YYSYMBOL_additive_expression: /* additive_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3109 "basilisk.c"
        break;

    case YYSYMBOL_shift_expression: /* shift_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3118 "basilisk.c"
        break;

    case YYSYMBOL_relational_expression: /* relational_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3127 "basilisk.c"
        break;

    case YYSYMBOL_equality_expression: /* equality_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3136 "basilisk.c"
        break;

    case YYSYMBOL_and_expression: /* and_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3145 "basilisk.c"
        break;

    case YYSYMBOL_exclusive_or_expression: /* exclusive_or_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3154 "basilisk.c"
        break;

    case YYSYMBOL_inclusive_or_expression: /* inclusive_or_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3163 "basilisk.c"
        break;

    case YYSYMBOL_logical_and_expression: /* logical_and_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3172 "basilisk.c"
        break;

    case YYSYMBOL_logical_or_expression: /* logical_or_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3181 "basilisk.c"
        break;

    case YYSYMBOL_conditional_expression: /* conditional_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3190 "basilisk.c"
        break;

    case YYSYMBOL_assignment_expression: /* assignment_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3199 "basilisk.c"
        break;

    case YYSYMBOL_assignment_operator: /* assignment_operator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3208 "basilisk.c"
        break;

    case YYSYMBOL_expression: /* expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3217 "basilisk.c"
        break;

    case YYSYMBOL_constant_expression: /* constant_expression  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3226 "basilisk.c"
        break;

    case YYSYMBOL_declaration: /* declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3235 "basilisk.c"
        break;

    case YYSYMBOL_declaration_specifiers: /* declaration_specifiers  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3244 "basilisk.c"
        break;

    case YYSYMBOL_init_declarator_list: /* init_declarator_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3253 "basilisk.c"
        break;

    case YYSYMBOL_init_declarator: /* init_declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3262 "basilisk.c"
        break;

    case YYSYMBOL_storage_class_specifier: /* storage_class_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3271 "basilisk.c"
        break;

    case YYSYMBOL_type_specifier: /* type_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3280 "basilisk.c"
        break;

    case YYSYMBOL_types: /* types  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3289 "basilisk.c"
        break;

    case YYSYMBOL_struct_or_union_specifier: /* struct_or_union_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3298 "basilisk.c"
        break;

    case YYSYMBOL_struct_or_union: /* struct_or_union  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3307 "basilisk.c"
        break;

    case YYSYMBOL_struct_declaration_list: /* struct_declaration_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3316 "basilisk.c"
        break;

    case YYSYMBOL_struct_declaration: /* struct_declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3325 "basilisk.c"
        break;

    case YYSYMBOL_specifier_qualifier_list: /* specifier_qualifier_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3334 "basilisk.c"
        break;

    case YYSYMBOL_struct_declarator_list: /* struct_declarator_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3343 "basilisk.c"
        break;

    case YYSYMBOL_struct_declarator: /* struct_declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3352 "basilisk.c"
        break;

    case YYSYMBOL_enum_specifier: /* enum_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3361 "basilisk.c"
        break;

    case YYSYMBOL_enumerator_list: /* enumerator_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3370 "basilisk.c"
        break;

    case YYSYMBOL_enumerator: /* enumerator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3379 "basilisk.c"
        break;

    case YYSYMBOL_atomic_type_specifier: /* atomic_type_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3388 "basilisk.c"
        break;

    case YYSYMBOL_type_qualifier: /* type_qualifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3397 "basilisk.c"
        break;

    case YYSYMBOL_function_specifier: /* function_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3406 "basilisk.c"
        break;

    case YYSYMBOL_alignment_specifier: /* alignment_specifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3415 "basilisk.c"
        break;

    case YYSYMBOL_declarator: /* declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3424 "basilisk.c"
        break;

    case YYSYMBOL_direct_declarator: /* direct_declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3433 "basilisk.c"
        break;

    case YYSYMBOL_generic_identifier: /* generic_identifier  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3442 "basilisk.c"
        break;

    case YYSYMBOL_pointer: /* pointer  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3451 "basilisk.c"
        break;

    case YYSYMBOL_type_qualifier_list: /* type_qualifier_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3460 "basilisk.c"
        break;

    case YYSYMBOL_parameter_type_list: /* parameter_type_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3469 "basilisk.c"
        break;

    case YYSYMBOL_parameter_list: /* parameter_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3478 "basilisk.c"
        break;

    case YYSYMBOL_parameter_declaration: /* parameter_declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3487 "basilisk.c"
        break;

    case YYSYMBOL_identifier_list: /* identifier_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3496 "basilisk.c"
        break;

    case YYSYMBOL_type_name: /* type_name  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3505 "basilisk.c"
        break;

    case YYSYMBOL_abstract_declarator: /* abstract_declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3514 "basilisk.c"
        break;

    case YYSYMBOL_direct_abstract_declarator: /* direct_abstract_declarator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3523 "basilisk.c"
        break;

    case YYSYMBOL_type_not_specified: /* type_not_specified  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3532 "basilisk.c"
        break;

    case YYSYMBOL_initializer: /* initializer  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3541 "basilisk.c"
        break;

    case YYSYMBOL_initializer_list: /* initializer_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3550 "basilisk.c"
        break;

    case YYSYMBOL_designation: /* designation  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3559 "basilisk.c"
        break;

    case YYSYMBOL_designator_list: /* designator_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3568 "basilisk.c"
        break;

    case YYSYMBOL_designator: /* designator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3577 "basilisk.c"
        break;

    case YYSYMBOL_static_assert_declaration: /* static_assert_declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3586 "basilisk.c"
        break;

    case YYSYMBOL_statement: /* statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3595 "basilisk.c"
        break;

    case YYSYMBOL_labeled_statement: /* labeled_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3604 "basilisk.c"
        break;

    case YYSYMBOL_compound_statement: /* compound_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3613 "basilisk.c"
        break;

    case YYSYMBOL_188_1: /* @1  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3622 "basilisk.c"
        break;

    case YYSYMBOL_block_item_list: /* block_item_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3631 "basilisk.c"
        break;

    case YYSYMBOL_block_item: /* block_item  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3640 "basilisk.c"
        break;

    case YYSYMBOL_expression_statement: /* expression_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3649 "basilisk.c"
        break;

    case YYSYMBOL_selection_statement: /* selection_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3658 "basilisk.c"
        break;

    case YYSYMBOL_for_scope: /* for_scope  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3667 "basilisk.c"
        break;

    case YYSYMBOL_iteration_statement: /* iteration_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3676 "basilisk.c"
        break;

    case YYSYMBOL_for_declaration_statement: /* for_declaration_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3685 "basilisk.c"
        break;

    case YYSYMBOL_jump_statement: /* jump_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3694 "basilisk.c"
        break;

    case YYSYMBOL_external_declaration: /* external_declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3703 "basilisk.c"
        break;

    case YYSYMBOL_function_declaration: /* function_declaration  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3712 "basilisk.c"
        break;

    case YYSYMBOL_function_definition: /* function_definition  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3721 "basilisk.c"
        break;

    case YYSYMBOL_declaration_list: /* declaration_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3730 "basilisk.c"
        break;

    case YYSYMBOL_basilisk_statements: /* basilisk_statements  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3739 "basilisk.c"
        break;

    case YYSYMBOL_macro_statement: /* macro_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3748 "basilisk.c"
        break;

    case YYSYMBOL_reduction_list: /* reduction_list  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3757 "basilisk.c"
        break;

    case YYSYMBOL_reduction: /* reduction  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3766 "basilisk.c"
        break;

    case YYSYMBOL_reduction_operator: /* reduction_operator  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3775 "basilisk.c"
        break;

    case YYSYMBOL_reduction_array: /* reduction_array  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3784 "basilisk.c"
        break;

    case YYSYMBOL_foreach_dimension_statement: /* foreach_dimension_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3793 "basilisk.c"
        break;

    case YYSYMBOL_forin_declaration_statement: /* forin_declaration_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3802 "basilisk.c"
        break;

    case YYSYMBOL_forin_statement: /* forin_statement  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3811 "basilisk.c"
        break;

    case YYSYMBOL_forin_arguments: /* forin_arguments  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3820 "basilisk.c"
        break;

    case YYSYMBOL_event_definition: /* event_definition  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3829 "basilisk.c"
        break;

    case YYSYMBOL_event_parameters: /* event_parameters  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3838 "basilisk.c"
        break;

    case YYSYMBOL_event_parameter: /* event_parameter  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3847 "basilisk.c"
        break;

    case YYSYMBOL_boundary_definition: /* boundary_definition  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3856 "basilisk.c"
        break;

    case YYSYMBOL_external_foreach_dimension: /* external_foreach_dimension  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3865 "basilisk.c"
        break;

    case YYSYMBOL_attribute: /* attribute  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3874 "basilisk.c"
        break;

    case YYSYMBOL_new_field: /* new_field  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3883 "basilisk.c"
        break;

    case YYSYMBOL_root: /* root  */
#line 54 "basilisk.yacc"
         {
  fputc ('\n', stderr);
  ast_print_tree ((*yyvaluep), stderr, "  ", 0, -1);
}
#line 3892 "basilisk.c"
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
    YYNOMEM;
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
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
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
#line 4290 "basilisk.c"
    break;

  case 5: /* translation_unit: translation_unit error '}'  */
#line 119 "basilisk.yacc"
                                      { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4296 "basilisk.c"
    break;

  case 6: /* translation_unit: translation_unit error ')'  */
#line 120 "basilisk.yacc"
                                      { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4302 "basilisk.c"
    break;

  case 16: /* expression_error: error  */
#line 136 "basilisk.yacc"
                        { yyvsp[0]->sym = YYSYMBOL_YYerror; }
#line 4308 "basilisk.c"
    break;

  case 114: /* declaration: declaration_specifiers ';' type_not_specified  */
#line 332 "basilisk.yacc"
                                                                               { ast_push_declaration (parse->stack, yyval); }
#line 4314 "basilisk.c"
    break;

  case 115: /* declaration: declaration_specifiers init_declarator_list ';' type_not_specified  */
#line 333 "basilisk.yacc"
                                                                               { ast_push_declaration (parse->stack, yyval); }
#line 4320 "basilisk.c"
    break;

  case 139: /* type_specifier: types  */
#line 372 "basilisk.yacc"
                   { parse->type_already_specified = true; }
#line 4326 "basilisk.c"
    break;

  case 159: /* struct_or_union_specifier: struct_or_union '{' error '}'  */
#line 398 "basilisk.yacc"
                                                            { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4332 "basilisk.c"
    break;

  case 160: /* struct_or_union_specifier: struct_or_union generic_identifier '{' error '}'  */
#line 399 "basilisk.yacc"
                                                            { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4338 "basilisk.c"
    break;

  case 210: /* direct_declarator: direct_declarator '(' type_not_specified error ')'  */
#line 494 "basilisk.yacc"
                                                                                           { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4344 "basilisk.c"
    break;

  case 260: /* type_not_specified: %empty  */
#line 577 "basilisk.yacc"
           { parse->type_already_specified = false; }
#line 4350 "basilisk.c"
    break;

  case 281: /* statement: error ';'  */
#line 619 "basilisk.yacc"
                     { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4356 "basilisk.c"
    break;

  case 286: /* @1: %empty  */
#line 631 "basilisk.yacc"
                           { stack_push (parse->stack, &(yyvsp[0])); yyval->sym = YYSYMBOL_YYUNDEF; }
#line 4362 "basilisk.c"
    break;

  case 287: /* compound_statement: '{' @1 block_item_list '}'  */
#line 633 "basilisk.yacc"
                           { ast_pop_scope (parse->stack, yyvsp[-3]); }
#line 4368 "basilisk.c"
    break;

  case 297: /* for_scope: FOR  */
#line 658 "basilisk.yacc"
              { stack_push (parse->stack, &(yyval)); }
#line 4374 "basilisk.c"
    break;

  case 300: /* iteration_statement: for_scope '(' expression_statement expression_statement ')' statement  */
#line 665 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-5]); }
#line 4380 "basilisk.c"
    break;

  case 301: /* iteration_statement: for_scope '(' expression_statement expression_statement expression ')' statement  */
#line 667 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-6]); }
#line 4386 "basilisk.c"
    break;

  case 303: /* for_declaration_statement: for_scope '(' declaration expression_statement ')' statement  */
#line 673 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-5]); }
#line 4392 "basilisk.c"
    break;

  case 304: /* for_declaration_statement: for_scope '(' declaration expression_statement expression ')' statement  */
#line 675 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-6]); }
#line 4398 "basilisk.c"
    break;

  case 316: /* external_declaration: error compound_statement  */
#line 693 "basilisk.yacc"
                                                { yyvsp[-1]->sym = YYSYMBOL_YYerror; }
#line 4404 "basilisk.c"
    break;

  case 317: /* function_declaration: declaration_specifiers declarator  */
#line 697 "basilisk.yacc"
                                            { ast_push_function_definition (parse->stack, yyvsp[0]);  }
#line 4410 "basilisk.c"
    break;

  case 318: /* function_definition: function_declaration declaration_list compound_statement  */
#line 702 "basilisk.yacc"
                                                { ast_pop_scope (parse->stack, yyvsp[-2]->child[1]); }
#line 4416 "basilisk.c"
    break;

  case 319: /* function_definition: function_declaration compound_statement  */
#line 704 "basilisk.yacc"
                                                { ast_pop_scope (parse->stack, yyvsp[-1]->child[1]); }
#line 4422 "basilisk.c"
    break;

  case 339: /* forin_declaration_statement: for_scope '(' declaration_specifiers declarator IDENTIFIER forin_arguments ')' statement  */
#line 755 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-7]); }
#line 4428 "basilisk.c"
    break;

  case 340: /* forin_statement: for_scope '(' expression IDENTIFIER forin_arguments ')' statement  */
#line 760 "basilisk.yacc"
                    { ast_pop_scope (parse->stack, yyvsp[-6]); }
#line 4434 "basilisk.c"
    break;

  case 356: /* root: translation_unit  */
#line 803 "basilisk.yacc"
                           {
	  yyval = root;
	  yyval->sym = yyr1[yyn];
	  yyval->child = allocate ((Allocator *)parse->alloc, 2*sizeof(Ast *));
	  yyval->child[0] = yyvsp[0];
	  yyvsp[0]->parent = yyval;
	  yyval->child[1] = NULL;
        }
#line 4447 "basilisk.c"
    break;


#line 4451 "basilisk.c"

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
  ++yynerrs;

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
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (parse, root, YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
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

#line 813 "basilisk.yacc"


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
