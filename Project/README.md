


Pour compiler le programme, il faut exécuter les instructions suivantes dans le terminal:

1) mkdir build
2) cd build
3) cmake ..
4) make

Un exécutable myFem va être créé. Pour l'exécuter il faut entrer la ligne de commande suivante dans le terminal:
./myFem





Le programme lis le fichier mesh.txt qui contient le maillage et le fichier problem.txt qui contient les données physique du problème d'élasticité linéaire. Le programme écris le résultat du problème dans le fichier result.txt. Tous ces fichiers se trouve dans le dossier data.


Lorsqu'on lance le programme, au bout d'un certain moment, des données sur le problème s'affiche durant la résolution de celui-ci. Une fenêtre s'ouvre affichant la solution de ntre problème. En appuillant sur certaines touches, on peut modifier le contenu de la fenêtre pour visualiser d'autres aspects du problème:

1) La touche S permet d'afficher la spy matrice A du système linéaire.
2) La touche R permet d'exécuter le code sans renumérotation. 
3) La touche T permet d'exécuter le code avec une renumérotation selon l'axe x.
4) La touche T permet d'exécuter le code avec une renumérotation selon l'axe y.
5) La touche D permet d'afficher les domaines du maillages.
6) La touche D permet de parcourir les différents domaines du maillage.
7) La touche V permet d'afficher la solution du problème.
8) La touche X permet d'afficher les forces résultantes en x.
9) La touche Y permet d'afficher les forces résultantes en y.






Au début de la fonction main, on définie les conditions de Neumann et de Dirichlet.
Pour modifier les conditions frontières, on peut ajouter ou retirer les fonction suivantes:
1) femElasticityAddBoundaryCondition(theProblem, nom du domaine,NEUMANN_Y,value);
1) femElasticityAddBoundaryCondition(theProblem, nom du domaine,NEUMANN_X,value);
1) femElasticityAddBoundaryCondition(theProblem, nom du domaine,Dirichlet_Y,value);
1) femElasticityAddBoundaryCondition(theProblem, nom du domaine,Dirichlet_X,value);