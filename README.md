### Tas de sable abéliens

Pour lancer toutes les méthodes et consulter les temps / résultats / speedup :

	make test
	make run

Pour visualiser les algorithmes :

	make
	./bin/sandpiles-m<i>-c<j>-d<k>

où :

i représente la méthode utilisée et peut être 1 pour la version
	séquentielle de référence, 2 pour la version séquentielle
	euclienne optimisée, 3 pour la version parallèle;

j représente comment sont initialiser les cases et peut être 1 pour
	l'initialisation homogène des matrices et 2 pour l'initialisation
	d'un seul grand tas de sable;

k répresente la dimension de la matrice et peut être 128 ou 512
