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

/**
## Tricks for AST

The following are declarations semantically equivalent to their real
implementations (which are often C preprocessor macros). */

enum AstBoolean { false, true };

/**
# Stencils

Need to know about these implicitly declared variables/macros. */

void point;
int BGHOSTS, o_stencil;
double HUGE;
void * NULL;

void _Variables() {
  double x, y, z;
  double Delta;
  int level;
}
