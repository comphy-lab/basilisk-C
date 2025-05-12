/**
# Manning friction for Saint-Venant

When the Reynolds number is high (>2000), the stream becomes turbulent
and the friction can no longer be obtained analytically. Manning
proposed an empirical relation describing the corresponding
friction. It can be written in its full form as: 
$$ C_f = - n^2 g \frac{q|q|}{h^{7/3}} $$ 
where $n$ is the "Manning coefficient", a
free-parameter wich depends on the nature of the soil. */

double n = 0; 
// fixme: nmanning[] is allocated even n is a constant
// this should use constant fields instead
scalar nmanning[];

event manning_friction (i++)
{
  foreach()
    if (h[] > dry) {      
      if (nmanning[] != 0)
	n = nmanning[];      
      double s = 1. + dt*G*sq(n)*norm(u)/pow(h[],4/3.);
      foreach_dimension()
	u.x[] /= s;
    }
  boundary ((scalar *){u});
}

/**
## Link to the homepage

* [Homepage](Readme)
*/
