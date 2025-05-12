// in file src/examples/ginzburg-landau.c in line 44
// two variables are introduced:
mgstats mgd1, mgd2;

// while only mgd1 is used both in line 99:
mgd1 = diffusion (Ar, dt, r = r, beta = lambda);
// and later in line 104:
mgd1 = diffusion (Ai, dt, r = r, beta = lambda);
// while mgd2 is not used in that file
// (but I do not know whether it is changed somewhere else)

// I am not sure whether:
// 1. this is written as it should be;
// 2. mgd2 does not need to be introduced;
// 3. in line 104 mgd1 should be changed to mgd2;
