/**
# Dumps are mangled when using 3D multigrid over a certain level

When using 3D multigrid over a certain level, dumps are mangled
and you get the following error:

```
restore(): error: the number of processes don't match: 0 != 1024
```

despite having used the right amount of processors.

*/

#include "navier-stokes/multigrid3D.h"
#define LEVEL 8 // You need to change this to 10 or 11 to see the bug

int main(){

  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << (LEVEL);
  init_grid(N);

  scalar p[];
  foreach(){
    p[] = x;
  }

  dump("restore.dump");
  restore("restore.dump");
}

/** 

The issue can be traced output.h. The offset is a `long` but the multiplication 
`index[] * cell_size` is not cast to `long` before multiplication, which can
lead to an integer overflow when `index[]` and `cell_size` are large. A race 
condition when writing the header and the data can also lead to mangled dumps.

~~~diff
--- a/output.h
+++ b/output.h
@@ -1154,10 +1154,13 @@
   strcpy (name, file);
   if (!unbuffered)
     strcat (name, "~");
-  FILE * fh = fopen (name, "w");
-  if (fh == NULL) {
-    perror (name);
-    exit (1);    
+  FILE * fh = NULL;
+  if (pid() == 0) {
+    fh = fopen (name, "w");
+    if (fh == NULL) {
+      perror (name);
+      exit (1);    
+    }
   }
 
   scalar * dlist = dump_list (list, zero);
@@ -1172,8 +1175,20 @@
   MPI_Barrier (MPI_COMM_WORLD);
 #endif
 
-  if (pid() == 0)
+  if (pid() == 0) {
     dump_header (fh, &header, slist);
+    fflush (fh);
+  }
+  
+  MPI_Barrier (MPI_COMM_WORLD);
+
+  if (pid() != 0) {
+    fh = fopen (name, "r+");
+    if (fh == NULL) {
+      perror (name);
+      exit (1);
+    }
+  }
   
   scalar index = {-1};
   
@@ -1190,7 +1205,7 @@
   foreach_cell() {
     // fixme: this won't work when combining MPI and mask()
     if (is_local(cell)) {
-      long offset = sizeofheader + index[]*cell_size;
+      long offset = sizeofheader + (long)index[] * (long)cell_size;
       if (pos != offset) {
 	fseek (fh, offset, SEEK_SET);
 	pos = offset;
@@ -1323,8 +1338,8 @@
 
 #if MULTIGRID_MPI
   long cell_size = sizeof(unsigned) + header.len*sizeof(double);
-  long offset = pid()*((1 << dimension*(header.depth + 1)) - 1)/
-    ((1 << dimension) - 1)*cell_size;
+  long offset = pid()*((1L << dimension*(header.depth + 1)) - 1)/
+    ((1L << dimension) - 1)*cell_size;
   if (fseek (fp, offset, SEEK_CUR) < 0) {
     perror ("restore(): error while seeking");
     exit (1);
--- a/grid/multigrid-mpi.h
+++ b/grid/multigrid-mpi.h
@@ -232,9 +232,9 @@
 {
   long i;
   if (leaves)
-    i = pid()*(1 << dimension*depth());
+    i = (long)pid() * (1L << (dimension * depth()));
   else
-    i = pid()*((1 << dimension*(depth() + 1)) - 1)/((1 << dimension) - 1);
+    i = (long)pid() * ((1L << (dimension * (depth() + 1))) - 1) / ((1L << dimension) - 1);
   foreach_cell() {
     if (!leaves || is_leaf(cell))
~~~

*/