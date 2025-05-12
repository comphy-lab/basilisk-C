%rename(_poisson) poisson;
%inline %{
  typedef struct {
    int i;              // number of iterations
    double resb, resa;  // maximum residual before and after the iterations
    double sum;         // sum of r.h.s.
  } mgstats;

  struct Poisson {
    scalar a, b;
    vector alpha;
    scalar _lambda;
    double tolerance;
    int nrelax;
    scalar * res;
  };
  
  extern mgstats poisson (struct Poisson p);
%}

%pythoncode %{
def poisson(a,b,alpha=None,lambda0=0,tolerance=1e-3):
    p = Poisson()
    p.a = a
    p.b = b
    if alpha != None: p.alpha = alpha
    p._lambda = lambda0
    p.tolerance = tolerance
    return _poisson(p)
%}

%{
  extern void mg_cycle (scalar * a, scalar * res, scalar * da,
			void (* relax) (scalar * da, scalar * res, 
				 int depth, void * data),
			void * data,
			int nrelax, int minlevel, int maxlevel);
  
  extern int NITERMAX;
  extern int NITERMIN;
  extern double TOLERANCE;

  extern mgstats mg_solve (scalar * a, scalar * b,
			   double (* residual) (scalar * a, scalar * b, 
						scalar * res,
						void * data),
			   void (* relax) (scalar * da, scalar * res, 
					   int depth, 
					   void * data),
			   void * data,
                           int nrelax,
                           scalar * res);
  
  extern mgstats project (face vector u, scalar p, 
			  face vector alpha, double dt, int nrelax);
%}
