%{
# Contributions

# Higer orde derivatives in diffmat_2D

[http://basilisk.fr/sandbox/easystab/stab2014/diffmat_2D_higher_order_derivatives.m]()
The limit of the method have been tested and an analyse have been done.

We have seen that we can't have a 16th orderderivative.

![Visualisation of the error according to derivative order](diff_2D_higer_order1.jpg)
![Visualisation of the error according to derivative order](diff_2D_higer_order2.jpg)
![Visualisation of the error according to derivative order](diff_2D_higer_order3.jpg)




# Flakner skan validation

We tried to validate the falkner skan code but found some unstable values for beta. We higlighted the problem we face to validate it with the Paulhosen profile and found values that are unstable for larage and small beta.

[http://basilisk.fr/sandbox/easystab/stab2014/falkner_skan_validation]()

The curves compared:

![The velocity profile](/sandbox/easystab/stab2014/falkner_skan_validation1.png)

The unstable values:

![The velocity profile](/sandbox/easystab/stab2014/falkner_skan_validation6.png)
![The velocity profile](/sandbox/easystab/stab2014/falkner_skan_validation7.png)


# Notes 

De la part de Jérôme

~~~matlab
domaine	          valeur	note
connectivité 	  2	        2
recyclage	  2	        2
graphiques	  2	        1
théories	  4	        0
Originalité	  4	        0
note /14	  14	        5
~~~

Dommage que tu n'as pas fait grand-chose. Et pourtant, tu avais l'air motivé au début. Ton anglais n'est pas non plus à la hauteur.


%}