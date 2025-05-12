void streamfunction (vector u, scalar psi) {
	scalar omega[];
	vorticity (u, omega);
	boundary ({psi, omega});
	poisson (psi, omega);
}
