# 2 "ast/defaults.h"

/**
# Default typedef/enum declarations 

To avoid having to parse large (and non-standard) system header
files. */

typedef int extrae_type_t, extrae_value_t;
typedef void QFILE;
struct timeval { long tv_sec, tv_usec; };

/**
For Python/SWIG */

typedef void PyObject;

/**
## From MPI */

typedef void MPI_Datatype, MPI_Request, MPI_Comm, MPI_Op, MPI_Aint;
typedef int MPI_Status;
typedef long long MPI_Offset;
typedef struct MPIR_Info *MPI_Info;

/**
## From OpenGL */

typedef int GLint, GLenum;
typedef float GLfloat;
typedef char GLubyte;

/**
## From standard C libraries */

typedef int bool;
typedef long ssize_t, size_t, clock_t, ptrdiff_t;
typedef long int64_t, int32_t, uint32_t, uint16_t, uint64_t;
typedef void va_list, FILE;
typedef unsigned char uint8_t;
typedef char int8_t;
typedef short int16_t;
typedef unsigned short uint16_t;
typedef int int32_t;
typedef unsigned int uint32_t;
typedef long int64_t;
typedef unsigned long uint64_t;
typedef double time_t;

/**
## Tricks for AST

The following are declarations semantically equivalent to their real
implementations (which are often C preprocessor macros). */

enum AstBoolean { false, true };

/**
# Stencils

Need to know about these implicitly declared variables/macros. */

int BGHOSTS, o_stencil;
double HUGE;
void * NULL;

/**
# Functions supported by GLSL */

bool is_constant();
void dimensional();
void NOT_UNUSED();
void neighborp();

double val_diagonal();

double abs (double x);
double acos (double x);
double acosh (double x);
double asin (double x);
double asinh (double x);
double atan (double x);
double atanh (double x);
double ceil (double x);
double cos (double x);
double cosh (double x);
double exp (double x);
double fabs(double x);
double log (double x);
double log2 (double x);
double fmax (double a, double b);
double fmin (double a, double b);
double mix (double a, double b, double c);
double mod (double a, double b);
double modf (double a, double * b);
double pow (double a, double b);
double round (double x);
double sin (double x);
double sinh (double x);
double sqrt (double x);
double tan (double x);
double tanh (double x);
double trunc (double x);
