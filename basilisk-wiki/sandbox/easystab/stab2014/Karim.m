# My contributions

## Third order differentiation matrix

Based on the [diffmat](../diffmat.m) code, here the code that derivate three times a function :  [diffmat_thirdorder.m](diffmat_thirdorder.m)

## Differential equation with Chebychev

I used the finite differences and also the chebdif function to compute the differenciation matrices and to resolve a simple differential equation. Here the link : [differential_equation_chebychev.m]()

## Comparison between the second order and first order differentiation matrix 

Here we can compare the accuracy when we want to derivate a function twice by the fisrt derivative matrix or once by the second derivative matrix : [diffmat_comparison](diffmat_comparison.m)

# Diffusion equation

I used the [Vibrating String](../vibrating_string.m) code and modified it to model [the diffusion equation](Diffusion.m)

![](/Diffusion.gif)

# Advection global

In this section, there is :

- [Diffusion global](diffusion_global.m)
- [Advection diffusion global](advection_diffusion_global.m)

based on the [Advection global](../advection_global.m) code that solve the advection equation without using the march in time. The principle is to consider the time as a spatial variable.

# Notes

De la part de Jérôme

~~~matlab
domaine	        valeur	note
connectivité 	2	1
recyclage	2	1
graphiques	2	1
théories	4	0
Originalité	4	2
note /14	14	5
~~~

Ce que tu as fait pour la dérivée troisième est pas mal. Mais les poids de différentiation non centrés ne marchent pas. L'as tu vu? je ne crois pas. C'est pas mal aussi pour la comparaison entre les deux manière d'avoir la dérivée seconde, mais encore une fois, les résultats ne sont pas très bon. la manière dont tu présente tes résultats sans les discuter laisse le lecteur dans la confusion: la comparaison est-elle vraiment aussi mauvaise? Comment ça dépend du nombre de points de grille? Pas de discution. Pour la diffusion, tu fais en moins bien ce que d'autres ont déjà fait. La solution analytique que tu compare à ta simulation montre qu'il n'y a pas accord. Pourquoi? Pas de réponse. Il y a des liens orphelins sur ta page. Si c'était un rapport, on ne prendrais même pas le temps de le lire. 