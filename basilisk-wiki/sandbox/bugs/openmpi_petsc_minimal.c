/**
Note: this possible issue was discovered by Robert Arslanbekov.

The [Petsc library](https://www.petsc.org) is built with a version of MPI specified when it is configured/installed. A straightforward installation may lead to OpenMPI+Petsc. If you then try to link Petsc to a different implementation of MPI, it will prevent you. 

If you try to use OpenMPI+Petsc in a Basilisk program, the qcc compiler stops early with syntax errors. You can modify the header mpi.h and then Basilisk+OpenMPI+Petsc will run. 

A minimal example is:
*/

#include "grid/quadtree.h"
#include "run.h"

#include "petscksp.h"

int main()
{
  size (1.[0]);
  init_grid(N);
  run();
}

event init(i = 0) {
  PetscInitialize(NULL, NULL, NULL, 0);
  printf("[init] i = %d\n", i);
}

/**

Compile with e.g.

`qcc -std=c99 -Wall -O3 -D_XOPEN_SOURCE=700 ${MPI_INCLUDE} ${PETSC_INCLUDE} ${target} -L{PETSC_LIB} -L{MPI_LIB} -lpetsc -lmpi -lm -otest`

and run with e.g.

`LD_LIBRARY_PATH=${MPI_LIB}:${PETSC_LIB} ./test`

where you might have 

`MPI_INCLUDE="-I/usr/lib/x86_64-linux-gnu/openmpi/include"`

`MPI_LIB="/usr/lib/x86_64-linux-gnu/openmpi/lib"`

`PETSC_INCLUDE="-I/home/user/petsc/include -I/home/user/petsc/arch-linux-c-debug/include"`

`PETSC_LIB="/home/user/petsc/arch-linux-c-debug/lib"`

This should produce compiler errors that can be spot-checked by editing the file `mpi.h` in the MPI_INCLUDE directory, by responding to errors produced by string-printing messages of the form:

`/usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h:2881:20: warning: ‘__error_’ attribute directive ignored [-Wattributes]`

`/usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h:2881:126: error: expected ‘,’ or ‘;’ before ‘)’ token`

This was tested for OpenMPI versions 4.1.6 and 5.0.3. For version 5.0.3, the following changes to `mpi.h` were sufficient:

* replace everywhere `__mpi_interface_deprecated__(<string>)` with `__mpi_interface_deprecated__()`
* remove the `#if` block at line 2805: `#if (!OMPI_OMIT_MPI1_COMPAT_DECLS || OMPI_BUILDING)...` which is only activated if Open MPI was configured with `--enable-mpi1-compatibility`, and was deactivated in this case.

The bug/issue is whether or not the qcc compiler can be updated to avoid having to modify headers in the OpenMPI library, in order to use Basilisk+OpenMPI+Petsc. 



**Removed block.** The following is the block at line 2805 of `mpi.h` in OpenMPI version 5.0.3:

*/

#if (!OMPI_OMIT_MPI1_COMPAT_DECLS || OMPI_BUILDING)
/*
 * Removed typedefs.  These typedefs are only available if Open MPI
 * was configured with --enable-mpi1-compatibility.
 *
 * These typedefs were formally removed from the MPI specification
 * and should no longer be used in MPI applications.
 *
 * Even though MPI_Handler_function is removed, we do not use the
 * attributes marking it as such, because otherwise the compiler
 * will warn for all the functions that are declared using them
 * (e.g., MPI_Errhandler_create).
 */
typedef void (MPI_Handler_function)(MPI_Comm *, int *, ...);
/* MPI_Handler_function was removed in MPI-3.0; use MPI_Comm_use_errhandler_function instead. */

/*
 * Removed prototypes.  These prototypes are only available if Open
 * MPI was configured with --enable-mpi1-compatibility.
 *
 * These functions were formally removed from the MPI specification
 * and should no longer be used in MPI applications.
 */
OMPI_DECLSPEC  int MPI_Address(void *location, MPI_Aint *address)
                   __mpi_interface_removed__(MPI_Address, MPI_Get_address);
OMPI_DECLSPEC  int PMPI_Address(void *location, MPI_Aint *address)
                   __mpi_interface_removed__(PMPI_Address, PMPI_Get_address);
OMPI_DECLSPEC  int MPI_Errhandler_create(MPI_Handler_function *function,
                                         MPI_Errhandler *errhandler)
                   __mpi_interface_removed__(MPI_Errhandler_create, MPI_Comm_create_errhandler);
OMPI_DECLSPEC  int PMPI_Errhandler_create(MPI_Handler_function *function,
                                          MPI_Errhandler *errhandler)
                   __mpi_interface_removed__(PMPI_Errhandler_create, PMPI_Comm_create_errhandler);
OMPI_DECLSPEC  int MPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler)
                   __mpi_interface_removed__(MPI_Errhandler_get, MPI_Comm_get_errhandler);
OMPI_DECLSPEC  int PMPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler)
                   __mpi_interface_removed__(PMPI_Errhandler_get, PMPI_Comm_get_errhandler);
OMPI_DECLSPEC  int MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler)
                   __mpi_interface_removed__(MPI_Errhandler_set, MPI_Comm_set_errhandler);
OMPI_DECLSPEC  int PMPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler)
                   __mpi_interface_removed__(PMPI_Errhandler_set, PMPI_Comm_set_errhandler);
OMPI_DECLSPEC  int MPI_Type_extent(MPI_Datatype type, MPI_Aint *extent)
                   __mpi_interface_removed__(MPI_Type_extent, MPI_Type_get_extent);
OMPI_DECLSPEC  int PMPI_Type_extent(MPI_Datatype type, MPI_Aint *extent)
                   __mpi_interface_removed__(PMPI_Type_extent, PMPI_Type_get_extent);
OMPI_DECLSPEC  int MPI_Type_hindexed(int count, int array_of_blocklengths[],
                                     MPI_Aint array_of_displacements[],
                                     MPI_Datatype oldtype, MPI_Datatype *newtype)
                   __mpi_interface_removed__(MPI_Type_hindexed, MPI_Type_create_hindexed);
OMPI_DECLSPEC  int PMPI_Type_hindexed(int count, int array_of_blocklengths[],
                                      MPI_Aint array_of_displacements[],
                                      MPI_Datatype oldtype, MPI_Datatype *newtype)
                   __mpi_interface_removed__(PMPI_Type_hindexed, PMPI_Type_create_hindexed);
OMPI_DECLSPEC  int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride,
                                    MPI_Datatype oldtype, MPI_Datatype *newtype)
                   __mpi_interface_removed__(MPI_Type_hvector, MPI_Type_create_hvector);
OMPI_DECLSPEC  int PMPI_Type_hvector(int count, int blocklength, MPI_Aint stride,
                                     MPI_Datatype oldtype, MPI_Datatype *newtype)
                   __mpi_interface_removed__(PMPI_Type_hvector, PMPI_Type_create_hvector);
OMPI_DECLSPEC  int MPI_Type_lb(MPI_Datatype type, MPI_Aint *lb)
                   __mpi_interface_removed__(MPI_Type_lb, MPI_Type_get_extent);
OMPI_DECLSPEC  int PMPI_Type_lb(MPI_Datatype type, MPI_Aint *lb)
                   __mpi_interface_removed__(PMPI_Type_lb, PMPI_Type_get_extent);
OMPI_DECLSPEC  int MPI_Type_struct(int count, int array_of_blocklengths[],
                                   MPI_Aint array_of_displacements[],
                                   MPI_Datatype array_of_types[],
                                   MPI_Datatype *newtype)
                   __mpi_interface_removed__(MPI_Type_struct, MPI_Type_create_struct);
OMPI_DECLSPEC  int PMPI_Type_struct(int count, int array_of_blocklengths[],
                                    MPI_Aint array_of_displacements[],
                                    MPI_Datatype array_of_types[],
                                    MPI_Datatype *newtype)
                   __mpi_interface_removed__(PMPI_Type_struct, PMPI_Type_create_struct);
OMPI_DECLSPEC  int MPI_Type_ub(MPI_Datatype mtype, MPI_Aint *ub)
                   __mpi_interface_removed__(MPI_Type_ub, MPI_Type_get_extent);
OMPI_DECLSPEC  int PMPI_Type_ub(MPI_Datatype mtype, MPI_Aint *ub)
                   __mpi_interface_removed__(PMPI_Type_ub, PMPI_Type_get_extent);
#endif /* !OMPI_OMIT_MPI1_COMPAT_DECLS */


