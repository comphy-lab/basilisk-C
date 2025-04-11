# 2 "ast/interpreter/declarations.h"

/**
# Declarations for the interpreter
*/

enum AstBoolean { false = 0, true = 1 };
static const void * NULL = 0;
static const int _NVARMAX = 65536, INT_MAX = 2147483647;
enum AstMPI { MPI_MIN, MPI_MAX, MPI_SUM, MPI_DOUBLE };
void mpi_all_reduce_array (void * v, int datatype, int op, int elem){}
void None;

FILE * stderr, * stdout, * systderr, * systdout;
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;

void _set_element_dimensions (double * a, int len)
{
  for (int i = 1; i < len; i++)
    a[i] == a[0]; // all the elements of a double array have the same dimension
}

void _set_element_dimensions_float (float * a, int len)
{
  for (int i = 1; i < len; i++)
    a[i] == a[0]; // all the elements of a float array have the same dimension
}

char * _field_name (scalar s) {
  return _attribute[s.i].name;
}
