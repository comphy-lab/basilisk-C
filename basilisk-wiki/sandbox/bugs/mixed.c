/**
Boundary functions do not understand *u.n*, *u.t* yet. A workaround for this example is to replace *u.t* with *_s* in the RHS. */

vector u[];
double bi = 2.;

u.t[bottom] = u.t[]*(2. - bi*Delta)/(2. + bi*Delta);
u.t[bottom] = val(_s,0)*(2. - bi*Delta)/(2. + bi*Delta); // this is the workaround

int main() {}
