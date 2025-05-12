/**
# Contact angles on an embedded boundary

This header file implements contact angles for [VOF interfaces](/src/vof.h) on
[embedded boundaries](/src/embed.h) taken from the sandbox of M. Tavares.

Here, the contact angle is fixed to 90°
*/

extern scalar f;
const scalar contact_angle[] = pi / 2.;

static inline coord normal_contact(coord ns, coord nf, double angle)
{
	coord n;
	#if dimension == 2
		if (-ns.x * nf.y + ns.y * nf.x > 0){ // 2D cross product
			n.x = -ns.x * cos(angle) + ns.y * sin(angle);
			n.y = -ns.x * sin(angle) - ns.y * cos(angle);
		}
		else{
			n.x = -ns.x * cos(angle) - ns.y * sin(angle);
			n.y = ns.x * sin(angle) - ns.y * cos(angle);
		}
	#else // dimension == 3
		// 1. Calculate Rotation Axis (Cross Product)
		coord rot_axis;
		foreach_dimension()
		{
			rot_axis.x = ns.y * nf.z - ns.z * nf.y;
		}
		normalize(&rot_axis); // Normalize the rotation axis for unit length

		// 2. Store Original Position
		coord p1 = {ns.x, ns.y, ns.z};

		// 3. Construct Quaternion
		double re = cos(angle / 2.);						 // Real part (scalar)
		double x = rot_axis.x * sin(angle / 2.); // Imaginary i part
		double y = rot_axis.y * sin(angle / 2.); // Imaginary j part
		double z = rot_axis.z * sin(angle / 2.); // Imaginary k part

		// 4. Apply Rotation Matrix to Original Point
		n.x = re * re * p1.x + 2 * y * re * p1.z - 2 * z * re * p1.y + x * x * p1.x + 2 * y * x * p1.y + 2 * z * x * p1.z - z * z * p1.x - y * y * p1.x;
		n.y = 2 * x * y * p1.x + y * y * p1.y + 2 * z * y * p1.z + 2 * re * z * p1.x - z * z * p1.y + re * re * p1.y - 2 * x * re * p1.z - x * x * p1.y;
		n.z = 2 * x * z * p1.x + 2 * y * z * p1.y + z * z * p1.z - 2 * re * y * p1.x - y * y * p1.z + 2 * re * x * p1.y - x * x * p1.z + re * re * p1.z;
	#endif
	return n;
}

/**
This function is an adaptation of the
[reconstruction()](/src/fractions.h#reconstruction) function which
takes into account the contact angle on embedded boundaries. */

void reconstruction_contact(scalar f, vector n, scalar alpha)
{

	/**
	We first reconstruct the (n, alpha) fields everywhere, using the
	standard function. */

	reconstruction(f, n, alpha);

	/**
	In cells which contain an embedded boundary and an interface, we
	modify the reconstruction to take the contact angle into account. */

	foreach (){
		if (cs[] < 1. && cs[] > 0. && f[] < 1. && f[] > 0.){
			coord ns = facet_normal(point, cs, fs);
			normalize(&ns);

			coord nf;
			foreach_dimension()
				nf.x = n.x[];

			coord nc = normal_contact(ns, nf, contact_angle[]);

			foreach_dimension()
				n.x[] = nc.x;

			alpha[] = line_alpha(f[], nc);
		}
	}
	boundary({n, alpha});
}

event contact(i++){
	vector n[];
	scalar alpha[];

	reconstruction_contact(f, n, alpha);

	foreach (){
		if (cs[] == 0.){
			double fc = 0., sfc = 0.;
			coord o = {x, y, z};
			foreach_neighbor(){
				if (cs[] < 1. && cs[] > 0. && f[] < 1. && f[] > 0.)
				{ // standard version
					double coef = cs[] * (1. - cs[]) * f[] * (1. - f[]);
					if (coef == 0)
					{																	 // Case where the embedded boundary is strictly aligned with the mesh (coef=0)
						double mean = (f[] + cs[]) / 2.; // Approximation of the solid interface position to avoid division by zero
						coef = mean * (1. - mean) * f[] * (1. - f[]);
					}
					sfc += coef;

					coord nf;
					foreach_dimension()
						nf.x = n.x[];
					coord a = {x, y, z}, b;
					foreach_dimension()
					{
						a.x = (o.x - a.x) / Delta - 0.5,
						b.x = a.x + 1.;
					}
					fc += coef * rectangle_fraction(nf, alpha[], a, b);
				}
			}

			if (sfc > 0.)
			{
				f[] = fc / sfc;
			}
		}
	}
	boundary({f});
}

extern scalar *interfaces;

event properties(i = 0)
{
	for (scalar c in interfaces)
		if (c.height.x.i)
			heights(c, c.height);
}

event tracer_advection(i++)
{
	for (scalar c in interfaces)
		if (c.height.x.i)
			heights(c, c.height);
}