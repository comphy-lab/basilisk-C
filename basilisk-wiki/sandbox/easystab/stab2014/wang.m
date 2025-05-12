%{
# Page de contribution : WANG jinshui
 %}
 
  %{
# Two ways of derivation
  Source : [two ways of derivation.m](sandbox/easystab/two ways of derivation)
  I use two ways to derive cos(x),one is [chebdif.m](http://basilisk.fr/sandbox/easystab/chebdif.m),the other is [finite difference](http://basilisk.fr/sandbox/easystab/differential_equation.m),and i find the same result. 

 ![method difference infinite](/diff.jpg)

 ![method chebychev](/chebychev.jpg)


 %}

 
  
 %{
# Derivation in 3D
I use [diffmat_3D](http://basilisk.fr/sandbox/easystab/diffmat_3D.m) to test an exponential function in [diffmat_3D_test_functions.m](http://basilisk.fr/sandbox/easystab/diffmat_3D_test_functions.m)
  
  ![test function](/3d.jpg)
 %} 
 
 %{
# I use the derivation 3D to solve equation of poisson 
I use derivation 3d to solve problem [poisson](http://basilisk.fr/sandbox/easystab/poisson3D.m) in [poissom.m](http://basilisk.fr/sandbox/easystab/stab2014/poisson.m)
![poisson](/poisson.jpg)


# Notes

De la part de Jérôme

~~~matlab
domaine	        valeur	note
connectivité 	2	2
recyclage	2	0
graphiques	2	1
théories	4	0
Originalité	4	0
note /14	14	3
~~~

Ta contribution principale est d'avoir testé de nouvelles fonction pour la dérivation en 2D et 3D. Ce genre de test correspond aux compétences des 2 premières séances de TP, alors que nous en avons fait 7! Le code de Poisson en 3D montre clairement que ça ne marche pas mais tu ne t'en rend même pas compte. Cela montre que tu n'as pas acquis les savoirs-faire de cette UE.