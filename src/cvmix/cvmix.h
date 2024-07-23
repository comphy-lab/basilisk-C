/**
# Interface with the Community Ocean and Vertical Mixing (CVMix) project

Community Ocean Vertical Mixing
([CVMix](https://github.com/CVMix/CVMix-description/raw/master/cvmix.pdf
)) is a software package that aims to provide transparent, robust,
flexible, well documented, shared Fortran source code for use in
parameterizing vertical mixing processes in numerical ocean models.

This header file provides a C-language interface for Basilisk (and
other projects, since the interface does not depend on Basilisk
itself).

## Installation

[CVMix](https://github.com/CVMix) needs to be
[installed](http://cvmix.github.io/#building) first and its source
code must be accessible, since the C interface is generated
automatically from the Fortran sources.

Once this done the C interface can be built using this [Makefile]()
and the correct location for the sources and/or Fortran compiler. For
example:

~~~bash
cd $BASILISK/src/cvmix/
CVMIX=$HOME/local/cvmix F90=gfortran FCFLAGS="-Wall -O2" make
~~~

## Interoperability with Fortran 90

Aside from standard considerations on compatibility of different types
between C and Fortran (as documented
[here](https://northstar-www.dartmouth.edu/doc/solaris-forte/manuals/fortran/prog_guide/11_cfort.html)
for example), the main difficulty is interoperability of ["assumed
size
arrays"](https://thinkingeek.com/2017/01/14/gfortran-array-descriptor/).

The approach taken here is a bit of a hack and should be replaced with
a more portable way of doing things [when it becomes
available](https://gcc.gnu.org/onlinedocs/gfortran/Further-Interoperability-of-Fortran-with-C.html#Further-Interoperability-of-Fortran-with-C).

For the moment, it relies on the assumption that the "shaped array
descriptor" works in a similar manner to that of
[gfortran](https://gcc.gnu.org/onlinedocs/gfortran). These assumptions are:

1. The first element of the array descriptor is a pointer to the data
stored in the array.

2. The size occupied by the array descriptor is smaller than or equal
to the size of the data structure below (i.e. `sizeof(cvmix_1d)`).

3. The only types used by CVMix are `integer`, `cvmix_r8`,
`character`, `logical` and derived types composed of these.

4. Fortran derived types never use more memory than the corresponding
C structure.

If any of these assumptions is violated, the most likely result will
be a low-level memory fault (i.e. segmentation fault, stack smashing
etc.). */

#pragma autolink -L$BASILISK/cvmix -lcvmixc -lgfortran

/**
We define a data structure describing the array descriptor of gfortran
as documented
[here](https://github.com/gcc-mirror/gcc/blob/master/gcc/fortran/trans-types.c#L1259). Note
that the goal is only to define a `cvmix_1d` structure which occupies
the right amount of memory so that conditions 2 and 4 above are
verified. The actual members are never used, with the exception of
`gfc_array_descriptor.a` which corresponds to assumption 1 above. */

typedef int indexing;

struct dtype_type
{
  size_t elem_len;
  int version;
  signed char rank;
  signed char type;
  signed short attribute;
};

struct descriptor_dimension
{
  indexing stride;
  indexing lbound;
  indexing ubound;
};

struct gfc_array_descriptor
{
  double * a;
  indexing offset;
  struct dtype_type dtype;
  struct descriptor_dimension dim[1];
  int padding[2]; // this seems necessary!!!
};

typedef double cvmix_r8;
typedef int logical;
typedef int integer;
typedef struct gfc_array_descriptor cvmix_1d;
typedef cvmix_1d cvmix_nd; // for n-rank arrays (not used)

typedef struct {
  double * a;
  indexing offset;
  struct dtype_type dtype;
  struct descriptor_dimension dim[2];
  int padding[2]; // this seems necessary!!!
} cvmix_2d;

cvmix_r8 cvmix_zero = 0., cvmix_one = 1.;
  
#define strlencheck(s) (s != NULL ? strlen(s) : 0)

#include "kinds_and_types.h"
#include "put_get.h"

extern void cvmix_redirect_stdout_ (void);
extern void allocate1d_ (const int * len, cvmix_1d * mem);
extern void deallocate1d_ (cvmix_1d * mem);
extern void allocate2d_ (const int * len1, const int * len2, cvmix_2d * mem);
extern void deallocate2d_ (cvmix_2d * mem);
extern void sizeofall_ (void);

extern void cvmix_deallocate_ (cvmix_data_type * CVmix_vars);

cvmix_1d cvmix_allocate_1d (int len)
{
  cvmix_1d a;
  allocate1d_ (&len, &a);
  return a;
}

void cvmix_deallocate_1d (cvmix_1d a)
{
  deallocate1d_ (&a);
}

void cvmix_1d_print (cvmix_1d * a)
{
  fprintf (stderr, "a: %p offset: %d dtype: %ld %d %d %d %d dim: %d %d %d\n",
	   (void *) a->a, a->offset, a->dtype.elem_len,
	   a->dtype.version, a->dtype.rank, a->dtype.type,
	   a->dtype.attribute,
	   a->dim[0].stride, a->dim[0].lbound, a->dim[0].ubound);  
}

cvmix_2d cvmix_allocate_2d (int len1, int len2)
{
  cvmix_2d a;
  allocate2d_ (&len1, &len2, &a);
  return a;
}

void cvmix_deallocate_2d (cvmix_2d a)
{
  deallocate2d_ (&a);
}

#define cvmix_deallocate(a) cvmix_deallocate_(a) 
#define cvmix_2d(c,len1,i,j) (c).a[(i) + (j)*(len1)]
