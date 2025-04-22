%rename(_delete) delete;
%inline %{
  extern scalar new_scalar (const char * name);
  extern void delete (scalar * list);
%}
