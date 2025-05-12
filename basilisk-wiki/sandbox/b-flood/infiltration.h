/**
# Infiltration in saint venant (mass equation)

## Green-ampt model
*/
scalar Inf[], Vol[], InfM[] ;
scalar K[], Ch[], Sat[];

event updateinf( i++ ){
  foreach(){      
    // Computing the maximum amount infiltrated in the soil
    if(Vol[] > 0)
      InfM[] = K[] * (1 + (Ch[]+h[]) * Sat[] / Vol[]);
    else
      InfM[] = HUGE ;
    if( dt > 0 ){ // fixing a little bug for the first step
      if( InfM[] > h[]/dt){
      	  Inf[] =  h[]/dt;
      	  foreach_dimension()
      	  	u.x[] = 0;
      	  h[] = 0;
      }
      else{
      	  Inf[] = InfM[];
      	  h[] -= Inf[]*dt;
      }
      if(dt < 1e9) // fixing a little bug when the whole area is dry at t=0
      	  Vol[] += Inf[]*dt;
    }
  }
}
/**
## Link to the homepage
* [Homepage](Readme)
*/

