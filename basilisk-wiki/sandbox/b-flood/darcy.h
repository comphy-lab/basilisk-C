/**
# Darcy friction for Saint-Venant

When the Reynolds number is high (>2000), the stream becomes turbulent
and the friction can no longer be obtained analytically. Darcy and
Weisbach proposed an empirical relation describing the corresponding
friction. It can be written in its full form as:
$$ 
C_f = \frac{f}{8} \frac{q|q|}{h^2}, 
$$ 
where $f$, the "Darcy coefficient", is a free-parameter wich depends
on the nature of the soil. */

double f = 0.5; 

event darcy_friction (i++)
{
  foreach()
    if (h[] > dry)
      foreach_dimension()
	u.x[] /= 1 + dt*f*norm(u)/(8*h[]);
  boundary ((scalar *){u});
}

/**
## Link to the homepage

* [Homepage](Readme)
*/
