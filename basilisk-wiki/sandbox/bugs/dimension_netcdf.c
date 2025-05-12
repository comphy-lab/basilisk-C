/**
 The code Below does not compile with 

qcc dimension_netcdf.c -lm -lnetcdf

(hangs forever). But it does with 

qcc -disable-dimensions dimension_netcdf.c -lm -lnetcdf

<div class="message">
This is because qcc does not like (when emulating the code) the
redundant names for 'struct NC_Dispatch' and 'typedef NC_Dispatch'. I
did not even realize that this was possible in C... 

The "proper" way to write this would be something like:

~~~c
typedef struct _NC_Dispatch NC_Dispatch;
~~~

but since this is legal (confusing) C, Basilisk should tolerate it. </div> 
*/

#if 1
typedef struct NC_Dispatch NC_Dispatch; // This hangs too
#else
#include <netcdf.h>
#endif

int main() {
}
