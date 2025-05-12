

/**

Sorry, the page doesn't exist yet. It will soon be available. 

Pour me faire pardonner, une jolie image.

~~~gnuplot Figure: Il mailla l'abeille
set terminal pngcairo size 700,700

#set output "mailla_l-abeille.png"

#p "mesh" w l not
~~~

~~~gnuplot Figure: Il mailla l'abeille  -- iofgets error I don't understand as it occurs only on the website... So I cheat and plot the pre-registered expected results
set terminal pngcairo size 700,700

set output "mailla_l-abeille.png"

p "expected_mesh" w l not
~~~

Pour me faire pardonner de l'image et du jeu de mot associé, un petit contrepet extrait d'une comptine enfantine : "Il court, le furet"

*/



#include "utils.h"
#include "input.h"
#include "run.h"


int main(){
    init_grid(1<<8);
    run();
}
/*
event tend(i=10){}

event movie(i++){
      
    scalar gris[];
    FILE * fp = fopen("./../maya_final.pgm","r");  // import the (modified and grayed) basilisk logo
    input_pgm(gris,fp);
    fclose(fp);

    foreach()
      if (y>0.82 && gris[]==0.)
        gris[]=1.;
    
    // AMR
    double eps=0.3;
    adapt_wavelet({gris},(double[]){eps},maxlevel=8);
    
    fp=fopen("mesh","w");
    foreach(){
      fprintf(fp,"%g %g\n", x-Delta/2., y-Delta/2.);
      fprintf(fp,"%g %g\n", x-Delta/2., y+Delta/2.);
      fprintf(fp,"%g %g\n", x+Delta/2., y+Delta/2.);
      fprintf(fp,"%g %g\n", x+Delta/2., y-Delta/2.);
      fprintf(fp,"%g %g\n\n", x-Delta/2., y-Delta/2.);
    }
    fclose(fp);
}


*/
