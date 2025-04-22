%{
  #define SWIG_FILE_WITH_INIT
  #define face
  #define vertex

  typedef int scalar;
  typedef struct {
    scalar x, y;
  } vector;
  typedef struct {
    vector x, y;
  } tensor;

  extern void _init_solver();
%}

%rename(_scalar) scalar;
typedef int scalar;

%rename(_vector) vector;
typedef struct {
  scalar x, y;
} vector;

%rename(_tensor) tensor;
typedef struct {
  vector x, y;
} tensor;

// scalar lists
%typemap(in) scalar * {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (scalar *) malloc((size+1)*sizeof(scalar));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyInt_Check(o))
	$1[i] = PyInt_AsLong(PyList_GetItem($input,i));
      else {
	PyErr_SetString(PyExc_TypeError,"list must contain scalars");
	free($1);
	return NULL;
      }
    }
    $1[i] = -1;
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// This cleans up the scalar * array we malloc'd before the function call
%typemap(freearg) scalar * {
  free((char *) $1);
}

%exception py_scalar_init {
  $action
  if (result) {
     PyErr_SetString(PyExc_RuntimeError, "initialization failed");
     return NULL;
  }
}
%inline %{
  extern void init_grid (int n);
  extern void free_grid (void);
  struct _origin { double x, y, z; };
  extern void origin (struct _origin p);
  extern void size (double L);

  extern int py_scalar_init (scalar s, PyObject * f);
  extern int py_register_event (PyObject * action, PyObject * i, 
				PyObject * t);
%}

%include "numpy.i"

%init %{
  import_array();
  _init_solver();
%}

%pythoncode %{
from numpy import empty_like
class scalar(int):
    def __new__(cls,i=None):
        if i == None:
            i = new_scalar("python")
        return int.__new__(cls,i)
    def __del__(self):
        if callable(_delete):
            _delete([self])
    def __setattr__(self, name, value):
        if name == "f":
            py_scalar_init(self,value)
        else:
            self.__dict__[name] = value
    def f(self,x,y=0):
        try:
            ndim = x.ndim
        except AttributeError:
            return interpolate(self,x,y)
        if ndim == 1:
            return _interpolate1D(self,x,x.size)
        elif ndim == 2:
            z = empty_like(x)
            _interpolate2D(self,x,y,z)
            return z
    def norm(self):
        return normf(self)
    def stats(self):
        return statsf(self)

class vector:
    def __init__(self,v):
        self.x = scalar(v.x)
        self.y = scalar(v.y)

class tensor:
    def __init__(self,t):
        self.x = vector(t.x)
        self.y = vector(t.y)

def event (action, i = None, t = None):
    py_register_event(action, i, t)

from random import uniform
def noise():
    return uniform(-1.,1.)
%}
