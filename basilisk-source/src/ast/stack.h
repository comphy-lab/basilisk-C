typedef struct _Stack Stack;

Stack * stack_new (int size);

void   stack_push    (Stack * s, void * p);
void * stack_pop     (Stack * s);
void * stack_index   (Stack * s, int i);
void * stack_indexi  (Stack * s, int i);
void   stack_destroy (Stack * s);

void * stack_set_push (Stack * s, void * push (Stack *, void *));
void * stack_set_data (Stack * s, void * data);
void * stack_get_data (const Stack * s);

void * fast_stack_find (Stack * s, const char * key);
