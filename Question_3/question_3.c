#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#define X 80
#define Y 80

/*Étudier cette la transition de phase et rechercher numériquement 
la température critique Tc en fonction de la densité ρ = N/Ns. 
On pourra augmenter la taille du système à densité ρ constante pour vérifier 
que Tc n’en dépend que faiblement. */

/* Dans cette version de l'algorithme, beaucoup de modifications ont été effectuées pour minimiser le temps que prend le programme 
à tourner étant donné que nous allons simuler une grande grille. Principalement, nous allons retrouver un tableau des indices qui va
nous permettre de rendre le tableau principal aléatoire et aussi le tableau des 1 qui répertorie toutes les coordonnées des 1 dans 
le tableau principal. Même si ces tableaux et les foncitons ne sont pas évidentes, elles me permettent d'effectuer des opérations comme
rendre le tableau aléatoire ou bien choisir une partoicule aléatoirement sans faire une double boucle dans le tableau principal. 
Ces fonctions ne sont pour autant probablement pas la meilleure façon de faire mais c'est le mieux que j'ai pu trouver comme solution*/

const int N = 80;
const int Ns =1000;
int calc_nb_voisins( int tab[][Y], int x, int y);
int choisir_particule(int tab[][Y], int tab_1[][2], int *a, int *b);
void choisir_voisin(int tab[][Y], int *a, int*b, int*c, int*d);
double calculer_energie_case(int x, int y, double J, int tab[][Y]);
void effectuer_saut(int tab[][Y],int tab_1[][2],int x,int y,int vx,int vy,double proba,double p1, double dE, int particule);
void udpate_E_moy(double Enr);
void petit_saut_tab(int tab[][Y]);
void reinitialiser_fichier(const char *nom_fichier);
double calculer_energie_equilibre(const char *filename, int nb_energies);
void afficher_energie_equilibre_gnuplot();
double calculer_energie_configuration(int tab[][Y], double J);
void shuffle(int *array, int n);
double calculer_energie_etendue(int tab[][Y], int x, int y, double J);

int main(){
    printf("%d\n", Ns);
    clock_t debut, fin;
    double temps_ecoule;

    // Enregistrer le temps de début
    debut = clock();

    //Initialisation des variables
    int i, x, y,voisx, voisy, iteration_count, count_1;
    const double kb = 1; //.39 * pow(10.0, -23);
    int randx;
    int randy;
    double T = 0.1; /*Température*/
    double E1; /* nrj particule choisie*/
    double E2; /* nrj voisin */
    double dE; /*diff energie*/
    const double J = 1; /*constante J*/
    double p0; /*Proba de sauter au site voisin*/
    double p1; /*Pour calculer la proba*/
    int iteration = 100000; /*Nombre d'itérations*/
    const char *energie_equilibre_filename = "energie_equilibre.txt";
    double dT = 0.01; /* Pas de température */
    int compteur_voisin_nul;
    double energie_equilibre;
    double energie_config;

    /*Importer temps et régler la seed de drand*/
    time_t *tm = NULL;
    srand48(time(tm));


    /* Réinitaliser les fichiers */
    reinitialiser_fichier(energie_equilibre_filename);
    
    /* Initialisation des tableaux */
    int tab[X][Y]; /* Tableau principal */
    int indices[X * Y]; /* Tableau des indices pour rendre tab aléatoire à chaque nouvelle température */
    int tab_1[Ns][2]; /* Tableau qui tient à jour la positions de tous les 1 (toutes les particules) dans le tableau principal */


    /* On crée le tableau des indices, pour l'instant chaque case contient son indice */
    for (i = 0; i < X * Y; i++) {
        indices[i] = i;
    } 
    
    /* Passer différentes températures grâce à un while*/
    while(T<5){
        shuffle(indices, X * Y); /* Mélange aléatoirement le tableau des indices */

        /* Rendre le tableau principal aléatoire */
        count_1 = 0;
        for (i = 0; i < X * Y; i++) {
            x = indices[i] / Y;
            y = indices[i] % Y;
            if (i < Ns) {
                tab[x][y] = 1;
                tab_1[count_1][0] = x; /* Mettre à jour le tableau qui répertorie la position des 1 */
                tab_1[count_1][1] = y;
                count_1 += 1;
            } else {
                tab[x][y] = 0;
            }
        }

        reinitialiser_fichier("nrj.txt");
        compteur_voisin_nul = 0;

        /* Nous voulons créer un fichier qui va retenir l'énergie de la configuration donc nous avons besoin de la calculer au moins une fois, le reste du temps on va seulement ajouer dE*/
        energie_config = calculer_energie_configuration(tab, J);

        /*Boucle principale*/
        for(iteration_count = 0; iteration_count<iteration ; iteration_count++){

            /* Prendre une particule au hasard */ 
            int particule = choisir_particule(tab, tab_1, &randx, &randy); /*La fonction renvoie un entier qui correspond au numéro de la particule. Si c'est la première elle aura le numéro 0 etc...*/
            
            /* Prendre un voisin au hasard */
            choisir_voisin(tab, &voisx, &voisy, &randx, &randy); 
            
            if(tab[voisx][voisy] == 0){
                
                E1 = calculer_energie_etendue(tab, randx, randy, J); /* On n'a pas besoin de calculer l'énergie de toute la configuration mais seulemetn des voisins des voisins, on aurait pu encore améliorer le programme en calculant seulement l'énergie des voisins du voisin sélectionné au lieu de tous les voisins mais cela aurait compliqué l'algorithme */
                
                
                // On regarde si on échange voisin 
                tab[voisx][voisy] = 1;
                tab[randx][randy] = 0;

                E2 = calculer_energie_etendue(tab, randx, randy, J); /* On calcule la même énergie que E1 mais en ayant fait sauter la particule */

                dE = E2-E1;
        
                p0 = exp(-dE/(kb*T));
                p1 = drand48();
            
                // On remet les particules comme initialement
                tab[voisx][voisy] = 0;
                tab[randx][randy] = 1;

                effectuer_saut(tab, tab_1 ,randx, randy, voisx, voisy, p0, p1, dE, particule);

                if (dE < 0 || p1 < p0) {
                    energie_config += dE;
                }
                if(iteration_count > iteration-60){
                    udpate_E_moy(energie_config); // On met à jour le fichier seulement lorsque on arrive à la fin de notre nombre d'itération pour minimiser le temps que prend notre code à tourner
                    compteur_voisin_nul += 1;
                }            
            }
        }

        /* Mettre à jour le tableau des énergies d'équilibre avec la moyenne des 50 dernières valeurs si on a assez de valeurs */
        if(compteur_voisin_nul > 50){
            energie_equilibre = calculer_energie_equilibre("nrj.txt", 50);
        }else{
            energie_equilibre = calculer_energie_configuration(tab, J);
        }
        FILE *energie_equilibre_file = fopen(energie_equilibre_filename, "a");
        if (energie_equilibre_file == NULL) {
            printf("Erreur lors de l'ouverture du fichier %s.\n", energie_equilibre_filename);
            return -1;
        }
        fprintf(energie_equilibre_file, "%.2f %f\n", T, energie_equilibre);
        fclose(energie_equilibre_file);

        T = T + dT;
    }

    afficher_energie_equilibre_gnuplot(); 

    fin = clock();
    temps_ecoule = (double)(fin - debut) / CLOCKS_PER_SEC;
    printf("Temps écoulé : %.5f secondes\n", temps_ecoule);
}

void shuffle(int *array, int n) {
    srand(time(NULL));

    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

int calc_nb_voisins(int tab[][Y], int x, int y){
    int nb_vois = 0;
    int voisx1 = ( x + 1 ) % X;
    int voisy1 = ( y + 1 ) % Y;
    int voisx_1 = ( (x-1) + X ) % X;
    int voisy_1 = ( (y-1) + Y ) % Y;

    if(tab[x][voisy1] == 1){
        nb_vois+=1;
    }
    if(tab[x][voisy_1] == 1){
        nb_vois+=1;
    }
    if(tab[voisx1][y] == 1){
        nb_vois+=1;
    }
    if(tab[voisx_1][y] == 1){
        nb_vois+=1;
    }
    return nb_vois;
}

double calculer_energie_etendue(int tab[][Y], int x, int y, double J) {
    double nrj = 0;

    for (int i = x-2; i <= x+2; i++) {
        for (int j = y-2; j <= y+2; j++) {
            int nx = (i + X) % X;
            int ny = (j + Y) % Y;

            if (tab[nx][ny] == 1) {
                nrj += calculer_energie_case(nx, ny, J, tab);
            }
        }
    }

    return nrj;
}

int choisir_particule(int tab[][Y], int tab_1[][2], int *a, int *b){

    // Générer un nombre aléatoire entre 1 et le nombre total de particules
    int randnomb = rand() % Ns ;

    //Affecter les positions grâce à la liste tab_1
    *a = tab_1[randnomb][0];
    *b = tab_1[randnomb][1];

    return randnomb;

}

void choisir_voisin(int tab[][Y], int *a, int *b, int *c, int *d) {
    int voisins[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    int rand_voisin;

    do {
        rand_voisin = rand() % 4;
        *a = (*c + voisins[rand_voisin][0] + X) % X;
        *b = (*d + voisins[rand_voisin][1] + Y) % Y;
    } while (*a == *c && *b == *d);
}

double calculer_energie_case(int x, int y, double J, int tab[][Y]){
    double E;
    E = -J * calc_nb_voisins(tab, x, y);
    return E;
}

double calculer_energie_configuration(int tab[][Y], double J){
    double Enr = 0;
    int i,j;
    for(i = 0; i<N ; i++){
        for(j=0; j<N; j++){
            if(tab[i][j] == 1){
                Enr = Enr + calculer_energie_case(i, j, J, tab);
            } 
        }
    }
    return Enr;

}

void effectuer_saut(int tab[][Y],int tab_1[][2],int x,int y,int vx,int vy,double proba,double p1, double dE, int particule){
    if (tab[vx][vy] == 1) {
        // Le voisin est déjà occupé, on ne fait rien
        return;
    }
    if (tab[x][y] == 0){
        // On n'a pas de particule 
        return;
    }
     if (dE < 0 || p1 < proba) {
        tab[x][y] = 0;
        tab[vx][vy] = 1;
        tab_1[particule][0] = vx;
        tab_1[particule][1] = vy;
    } 
}

void udpate_E_moy(double Enr){
    FILE *fichier;
    fichier = fopen ("nrj.txt", "a");
    fprintf(fichier, "%f\n", Enr);
    fclose(fichier);
}

void petit_saut_tab(int tab[][Y]){
    FILE *fichier;
    fichier = fopen ("tableau.txt", "a");
    fprintf(fichier, "\n \n \n");
    fclose(fichier);
}

void reinitialiser_fichier(const char *nom_fichier) {
    FILE *fichier;

    // Ouvrir le fichier en mode écriture
    fichier = fopen(nom_fichier, "w");

    if (fichier == NULL) {
        printf("Erreur lors de l'ouverture du fichier.\n");
        return;
    }

    // Fermer le fichier, le contenu existant est effacé
    fclose(fichier);
}

double calculer_energie_equilibre(const char *filename, int nb_energies) {
    double somme_energies = 0.0;
    double energie;
    int count = 0;

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Erreur lors de l'ouverture du fichier %s.\n", filename);
        return -1;
    }

    // Parcourir le fichier jusqu'à la fin
    while (fscanf(file, "%lf", &energie) != EOF) {
        count++;
    }

    // Revenir en arrière de 'nb_energies' valeurs
    fseek(file, -nb_energies * (sizeof(double) + 1), SEEK_END);

    // Lire les 'nb_energies' dernières valeurs et les additionner
    for (int i = 0; i < nb_energies; i++) {
        fscanf(file, "%lf", &energie);
        somme_energies += energie;
    }

    fclose(file);

    // Calculer et retourner la moyenne
    return somme_energies / nb_energies;
}

void afficher_energie_equilibre_gnuplot() {
    system("gnuplot -e \"set terminal pngcairo size 800,600; set output 'energie_equilibre.png'; set xlabel \\\"Température\\\"; set ylabel \\\"Énergie d'équilibre\\\"; set grid; plot 'energie_equilibre.txt' using 1:2 with linespoints title \\\"Énergie d''équilibre\\\"\"");
}
