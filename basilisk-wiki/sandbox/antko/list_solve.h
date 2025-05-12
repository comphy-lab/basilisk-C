/**
# Helper macro to invert a list of (linear) spatial operators

The macro below is a simple extension of the [solve macro](/src/solve.h) designed to invert linear equations. 
The macro also returns a (single) [multigrid statistics](poisson.h#mgstats). */

#include "poisson.h"

/**
## Implementation

The current implementation handles a maximum of five, possibly coupled, linear equations.
*Mandatory*: the size of the list must be provided through the `LIST_LEN` constant.
*/

const scalar NULL_SCALAR = {-1};

macro
mgstats list_solve 
 (scalar a, 
  double * func, 
  double * rhs,
  double * diags,
  scalar b = NULL_SCALAR,
  scalar c = NULL_SCALAR,
  scalar d = NULL_SCALAR,
  scalar e = NULL_SCALAR,
  int nrelax = 4,
  int minlevel = 0,
  double tolerance = TOLERANCE)
{{
    #ifndef LIST_LEN
        fprintf (stderr, "list_solve.h: CRITICAL: LIST_LEN not defined\n");
        abort();
    #endif
    if (LIST_LEN > 5) {
        fprintf (stderr, "list_solve.h:%d: error: list_solve() only handles up to 5 scalars\n", LINENO);
        return (mgstats){0};
    }
    scalar * list = NULL; 
    list = list_append(list, a);
    if (b.i != -1) list = list_append (list, b);
    if (c.i != -1) list = list_append (list, c);
    if (d.i != -1) list = list_append (list, d);
    if (e.i != -1) list = list_append (list, e);
    int len = list_len (list);
    if (len != LIST_LEN) {
        fprintf (stderr, "list_solve.h:%d: error: list length mismatch (%d -- should be %d)\n", LINENO, len, LIST_LEN);
        abort();
    }
    mgstats _s = (mgstats){0};
    scalar * _res_list, * _da_list;
    _da_list = list_clone (list);
    _res_list = list_clone (list);
    for (int b = 0; b < nboundary; b++)
        for (scalar _da in _da_list)
            _da.boundary[b] = _da.boundary_homogeneous[b];
    _s.nrelax = nrelax; 
    double _resb;
    {
        double maxres = 0.;
        foreach (reduction(max:maxres)) {
            for (int i = 0; i < len; i++) {
                scalar _res = _res_list[i];
                _res[] = rhs[i] - func[i];
                if (fabs (_res[]) > maxres)
                    maxres = fabs (_res[]);
            }
        }
        _resb = _s.resb = _s.resa = maxres;
    }
    for (_s.i = 0; _s.i < NITERMAX && (_s.i < NITERMIN || _s.resa > tolerance); _s.i++) {
        {
            for (scalar _res in _res_list) 
                restriction ({_res});
            int _maxlevel = grid->maxdepth;
            int _minlevel = min (minlevel, _maxlevel);
            for (int l = _minlevel; l <= _maxlevel; l++) {
                if (l == _minlevel)
                    foreach_level_or_leaf (l)
                        for (scalar _da in _da_list)
                            foreach_blockf (_da)
                                _da[] = 0.;
                else
                    foreach_level (l)
                        for (scalar _da in _da_list)
                            foreach_blockf (_da)
                                _da[] = bilinear (point, _da);
                for (scalar _da in _da_list)
                    boundary_level ({_da}, l);
                for (int i = 0; i < _s.nrelax; i++) {
                    foreach_level_or_leaf (l) {
                        scalar a = _da_list[0];
                        #if (LIST_LEN > 1)
                          scalar b = _da_list[1];
                        #endif
                        #if (LIST_LEN > 2) 
                          scalar c = _da_list[2];
                        #endif
                        #if (LIST_LEN > 3)
                          scalar d = _da_list[3];
                        #endif
                        #if (LIST_LEN > 4)
                          scalar e = _da_list[4];
                        #endif
                        for (int j = 0; j < len; j++) {
                            scalar _res = _res_list[j];
                            if (j == 0) a[] = 0;
                            #if (LIST_LEN > 1)
                                if (j == 1) b[] = 0;
                            #endif
                            #if (LIST_LEN > 2)
                                if (j == 2) c[] = 0;
                            #endif
                            #if (LIST_LEN > 3)
                                if (j == 3) d[] = 0;
                            #endif
                            #if (LIST_LEN > 4)
                                if (j == 4) e[] = 0;
                            #endif
                            double _n = _res[] - func[j], _d;
                            _d = diags[j];
                            if (j == 0) a[] = _n/_d;
                            #if (LIST_LEN > 1)
                                if (j == 1) b[] = _n/_d;
                            #endif
                            #if (LIST_LEN > 2)
                                if (j == 2) c[] = _n/_d;
                            #endif
                            #if (LIST_LEN > 3)
                                if (j == 3) d[] = _n/_d;
                            #endif
                            #if (LIST_LEN > 4)
                                if (j == 4) e[] = _n/_d;
                            #endif  
                        }
                    }
                    for (scalar _da in _da_list)
                        boundary_level ({_da}, l);
                }
            }
            foreach() {
                scalar a_iterator, _da_iterator;
                for (a_iterator, _da_iterator in list, _da_list)
                    foreach_blockf (a_iterator) 
                        a_iterator[] += _da_iterator[];
            }
        }
        {
            double maxres = 0.;
            foreach (reduction(max:maxres)) {
                for (int i = 0; i < len; i++) {
                    scalar _res = _res_list[i];
                    _res[] = rhs[i] - func[i];
                    if (fabs (_res[]) > maxres)
                        maxres = fabs (_res[]);
                }
            }
            _s.resa = maxres;
        }
        if (_s.resa > tolerance) {
            if (_resb/_s.resa < 1.2 && _s.nrelax < 100)
                _s.nrelax++;
            else if (_resb/_s.resa > 10 && _s.nrelax > 2)
                _s.nrelax--;
        }
        _resb = _s.resa;
    }
    _s.minlevel = minlevel;
    if (_s.resa > tolerance)
        fprintf (stderr,
            "list_solve.h:%d: warning: convergence for scalar list not reached after %d iterations\n"
            "  res: %g nrelax: %d\n", LINENO,
            _s.i, _s.resa, _s.nrelax),
            fflush (ferr);
    
    free (list);
    free (_res_list);
    free (_da_list);

    return _s;

}}