/**
# Friction term : Manning in saint venant

When the Reynolds number is high (>2000), the stream becomes
turbulent. In this case, the friction term can't be easily solved
analytically. Manning proposed an empirical law describing this term,
it can be written in its full form as :

$$
Cf = - n^2 g \frac{q|q|}{h^{7/3}}
$$
where n is a free-parameter wich depends on the nature of the soil.
*/

/**
Edit (Leo Carriba Demange) : n = cst
*/

// Manning coefficient
double n = 0;
scalar nmanning[];

event manningsourceterm(i++){

	foreach(){
		if(h[] > dry){
			/* if( nmanning[] != 0 )
				//n = nmanning[];

			//double s = G*sq(n)*norm(u)/pow(h[],4/3.);  */  // initialize nmanning[] in "event init" instead
			double s = G*sq(nmanning[])*norm(u)/pow(h[],4/3.);
			foreach_dimension()
      			u.x[] /= 1+dt*s;
      	}
    }
}
