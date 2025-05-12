/**
# Treatment of periodicity and location shift
I report here just two handy functions so to avoid the (tedious) account
of periodicity and space shift of a given location.
*/

/**
### Periodic Wrap
- x is the original variable (location)
- c is the shift I want to apply, on the same axis
- L0 is the domain size. Here I consider a domain going from -L0/2 to +L0/2
*/

#define PERIODIC_WRAP(x, c, L0) (fmod((x)-(c)+0.5*(L0), (L0)) < 0 ? \
                                 fmod((x)-(c)+0.5*(L0), (L0)) + (L0) : \
                                 fmod((x)-(c)+0.5*(L0), (L0))) - 0.5*(L0)

#define xPW  PERIODIC_WRAP(x,+xShift,L0)
#define yPW  PERIODIC_WRAP(y,+yShift,L0)

#define xPWP  PERIODIC_WRAP(x,p().x,L0)
#define yPWP  PERIODIC_WRAP(y,p().y,L0)

/**
### Does it make easier to work on periodic boundaries?
'p' is the actual cell location (x,y,z), 'c' is the center wrt which I want to shift
PS stands for PeriodicShift.
*/
#define PS(p,c) ( fabs(p-c) > fabs(p-c+L0) ? p-c+L0 : \
                  fabs(p-c) > fabs(p-c-L0) ? p-c-L0 : \
                                             p-c)


