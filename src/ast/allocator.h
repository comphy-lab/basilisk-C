typedef struct _Allocator Allocator;  
Allocator * new_allocator();
void * allocate (Allocator * a, long size);
void free_allocator (Allocator * a);
