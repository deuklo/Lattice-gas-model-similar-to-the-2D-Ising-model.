#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#define X 20
#define Y 20

const int N = 20; // Dimension du tableau 
const int Ns =100; // Nombre de particules 
int calc_nb_voisins( int tab[][Y], int x, int y);
void enregistrer_tableau_tab(int tab[][Y]);
int choisir_particule(int tab[][Y], int tab_1[][2], int *a, int *b);
void choisir_voisin(int tab[][Y], int *a, int*b, int*c, int*d);
double calculer_energie_case(int x, int y, double J, int tab[][Y]);
void effectuer_saut(int tab[][Y],int tab_1[][2],int x,int y,int vx,int vy,double proba,double p1, double dE, int particule);
void afficher_gnuplot_tab();
void reinitialiser_fichier(const char *nom_fichier);
void shuffle(int *array, int n);
double calculer_energie_etendue(int tab[][Y], int x, int y, double J);

int main(){
    clock_t debut, fin;
    double temps_ecoule;

    // Enregistrer le temps de début
    debut = clock();

    //Initialisation des variables
    int i, x, y, voisx, voisy, iteration_count;
    const double kb = 1; //.39 * pow(10.0, -23);
    int randx;
    int randy;
    double T = 5; /* Température, pas forcément besoin de la définir dans cette question */
    double E1; /* nrj particule choisie*/
    double E2; /* nrj voisin */
    double dE; /*diff energie*/
    const double J = 1; /*constante J*/
    double p0; /*Proba de sauter au site voisin*/
    double p1; /*Pour calculer la proba*/
    int iteration = 100000; /*Nombre d'itérations*/
    int particule; // Numéro de particule choisie aléa

    /* Importer temps et régler la seed de drand */
    time_t *tm = NULL;
    srand48(time(tm));

    /* Réinitaliser les fichiers */
    reinitialiser_fichier("tableau.txt");
    
    /* Initialisation des tableaux */
    int tab[X][Y]; /* Tableau principal */
    int indices[X * Y]; /* Tableau des indices pour rendre tab aléatoire à chaque nouvelle température */
    int tab_1[Ns][2]; /* Tableau qui tient à jour la positions de tous les 1 (toutes les particules) dans le tableau principal */

    /* On crée le tableau des indices, pour l'instant chaque case contient son indice */
    for (int i = 0; i < X * Y; i++) {
        indices[i] = i;
    } 
    shuffle(indices, X * Y); /* Mélange aléatoirement le tableau des indices */
        
    /* Rendre le tableau principal aléatoire */
    int count_1 = 0;
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


    /*Boucle principale*/
    for(iteration_count = 0; iteration_count<iteration ; iteration_count++){

        /* Prendre une particule au hasard */ 
        particule = choisir_particule(tab, tab_1, &randx, &randy); /*La fonction renvoie un entier qui correspond au numéro de la particule. Si c'est la première elle aura le numéro 0 etc...*/        
        
        /* Prendre un voisin au hasard */
        choisir_voisin(tab, &voisx, &voisy, &randx, &randy); 
            if(tab[voisx][voisy] == 0){
            
            E1 = calculer_energie_etendue(tab, randx, randy, J); /* On n'a pas besoin de calculer l'énergie de toute la configuration mais seulemetn des voisins des voisins, on aurait pu encore améliorer le programme en calculant seulement l'énergie des voisins du voisin sélectionné au lieu de tous les voisins mais cela aurait compliqué l'algorithme */
                
                
            // On regarde si on échange voisin 
            tab[voisx][voisy] = 1;
            tab[randx][randy] = 0;

            E2 = calculer_energie_etendue(tab, randx, randy, J); /* On calcule la même énergie que E1 mais en ayant fait sauter la particule */
            dE = E2-E1;
    
            p0 = exp(-dE/(20*J));
            p1 = drand48();
        
            // On remet les particules comme initialement
            tab[voisx][voisy] = 0;
            tab[randx][randy] = 1;
            effectuer_saut(tab, tab_1 ,randx, randy, voisx, voisy, p0, p1, dE, particule);
        }
    }

    /* Afficher les tableaux sur gnuplot */
    enregistrer_tableau_tab(tab);
    afficher_gnuplot_tab();

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

void enregistrer_tableau_tab(int tab[][Y]) {
    int l, m;
    FILE *fichier;
    fichier = fopen("tableau.txt", "w");
    for (l = 0; l < N; l++) {
        for (m = 0; m < N; m++) {
            fprintf(fichier, "%d ", tab[l][m]);
        }
        fprintf(fichier, "\n");
    }
    fclose(fichier);
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

void afficher_gnuplot_tab() {
    system("gnuplot -e \"set terminal pngcairo size 800,800; set output 'tableau_visualisation.png'; set xrange [0:19]; set yrange [0:19]; set size square; set tics out nomirror; set key off; plot 'tableau.txt' matrix with image\"");
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
