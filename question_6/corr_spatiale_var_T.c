#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define X 80
#define Y 80

/*Ici on fait la corrélation spatiale en faisant varier la température*/

/*Pour faire tourner le programme, il faut inclure la bibliothèque de GSL*/

const int N = 80;
const int Ns =1000;
int calc_nb_voisins( int tab[][Y], int x, int y);
void enregistrer_tableau_tab(int tab[][Y]);
void enregistrer_tableau_g(double tab[][Y]);
int choisir_particule(int tab[][Y], int tab_1[][2], int *a, int *b);
void choisir_voisin(int tab[][Y], int *a, int*b, int*c, int*d);
double calculer_energie_case(int x, int y, double J, int tab[][Y]);
void effectuer_saut(int tab[][Y],int tab_1[][2],int x,int y,int vx,int vy,double proba,double p1, double dE, int particule);
void udpate_E_moy(double Enr);
void afficher_gnuplot_tab();
void afficher_gnuplot_g();
double calculer_densite_locale(int tab[][Y], int x, int y, int taille_fenetre);
double proportion_phase_condensee(int tab[][Y], int taille_fenetre, double seuil_densite);
void petit_saut_tab(int tab[][Y]);
void reinitialiser_fichier(const char *nom_fichier);
double calculer_energie_equilibre(const char *filename, int nb_energies);
void afficher_energie_equilibre_gnuplot();
double calculer_energie_configuration(int tab[][Y], double J);
void shuffle(int *array, int n);
double calculer_energie_etendue(int tab[][Y], int x, int y, double J);
void enregistrer_tableau_g_1D(double g_1D[][2]);

struct data
{
  double *t;
  double *y;
  size_t n;
};

double
expo(const double L, const double t)
{
  return exp(-t/L);
}

int
func_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double yi = d->y[i];
      double y = expo(a, ti);
      gsl_vector_set(f, i, yi - y);
    }

  return GSL_SUCCESS;
}

int
func_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double ei = exp(-ti/a);

      gsl_matrix_set(J, i, 0, -ei/a);

    }

  return GSL_SUCCESS;
}

int
func_fvv (const gsl_vector * x, const gsl_vector * v,
          void *params, gsl_vector * fvv)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);

  double va = gsl_vector_get(v, 0);

  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double ei = exp(-ti/a);
      double Daa = ei / (a*a);
    
      double sum;

      sum = va*va*Daa;

      gsl_vector_set(fvv, i, sum);
    }

  return GSL_SUCCESS;
}

void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double rcond;

  (void) params; /* not used */

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);
  fprintf(stderr, "iter %2zu: a = %.4f \n",
          iter,
          gsl_vector_get(x, 0)
        );
}

void
solve_system(gsl_vector *x, gsl_multifit_nlinear_fdf *fdf,
             gsl_multifit_nlinear_parameters *params)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;

  const size_t max_iter = 200;
  const double xtol = 1.0e-100;
  const double gtol = 1.0e-100;
  const double ftol = 1.0e-100;

  const size_t n = fdf->n;
  const size_t p = fdf->p;
  gsl_multifit_nlinear_workspace *work =
    gsl_multifit_nlinear_alloc(T, params, n, p);
  gsl_vector * f = gsl_multifit_nlinear_residual(work);
  gsl_vector * y = gsl_multifit_nlinear_position(work);
  int info;
  double chisq0, chisq, rcond;

  /* initialize solver */
  gsl_multifit_nlinear_init(x, fdf, work);

  /* store initial cost */
  gsl_blas_ddot(f, f, &chisq0);

  /* iterate until convergence */
  gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                              callback, NULL, &info, work);

  /* store final cost */
  gsl_blas_ddot(f, f, &chisq);

  /* store cond(J(x)) */
  gsl_multifit_nlinear_rcond(&rcond, work);

  gsl_vector_memcpy(x, y);

  /* print summary */

  fprintf(stderr, "NITER         = %zu\n", gsl_multifit_nlinear_niter(work));
  fprintf(stderr, "initial cost  = %.12e\n", chisq0);
  fprintf(stderr, "final cost    = %.12e\n", chisq);
  fprintf(stderr, "final x       = %.12e\n", gsl_vector_get(x, 0));

  gsl_multifit_nlinear_free(work);
}

int main(){
    printf("%d\n", Ns);
    clock_t debut, fin;
    double temps_ecoule;

    // Enregistrer le temps de début
    debut = clock();

    //Initialisation des variables
    int i, k, l, m, x, y ,voisx, voisy, iteration_count;
    const double kb = 1; //.39 * pow(10.0, -23);
    int randx;
    int randy;
    double Temperature = 1; /*Temperature*/
    double E1; /* nrj particule choisie*/
    double E2; /* nrj voisin */
    double dE; /*diff energie*/
    const double J = 1; /*constante J*/
    double p0; /*Proba de sauter au site voisin*/
    double p1; /*Pour calculer la proba*/
    int iteration = 100000; /*Nombre d'itérations*/

    /* Importer temps et régler la seed de drand */
    time_t *tm = NULL;
    srand48(time(tm));


    /* Réinitaliser les fichiers */
    reinitialiser_fichier("tableau.txt");
    reinitialiser_fichier("g_expo_T.txt");
    
    /* Initialisation des tableaux */
    int tab[X][Y]; /* Tableau principal */
    int indices[X * Y]; /* Tableau des indices pour rendre tab aléatoire à chaque nouvelle Temperature */
    int tab_1[Ns][2]; /* Tableau qui tient à jour la positions de tous les 1 (toutes les particules) dans le tableau principal */


    /* On crée le tableau des indices, pour l'instant chaque case contient son indice */
    for (i = 0; i < X * Y; i++) {
        indices[i] = i;
    } 
    while(Temperature<5){
        
        shuffle(indices, X * Y); /* Mélange aléatoirement le tableau des indices */
            
        /* Rendre le tableau principal aléatoire */
        for (i = 0; i < X * Y; i++) {
            x = indices[i] / Y;
            y = indices[i] % Y;
            if (i < Ns) {
                tab[x][y] = 1;
                tab_1[i][0] = x; /* Mettre à jour le tableau qui répertorie la position des 1 */
                tab_1[i][1] = y;
            } else {
                tab[x][y] = 0;
            }
        }

        reinitialiser_fichier("nrj.txt");

        /* Nous voulons créer un fichier qui va retenir l'énergie de la configuration donc nous avons besoin de la calculer au moins une fois, le reste du temps on va seulement ajouer dE*/
        double energie_config = calculer_energie_configuration(tab, J);

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
        
                p0 = exp(-dE/(kb*Temperature));
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
                }            
            }
        }

        /*Faire la corrélation spatiale*/
        double moy_g = 0;
        double g[X][Y];
        for( k=0; k<X; k++){
            for( l=0; l<Y; l++){
                g[k][l] = 0;
                for( m=0; m<Ns; m++){
                     x = tab_1[m][0];
                     y = tab_1[m][1];
                    g[k][l] += tab[x][y]*tab[(x+k)%X][(y+l)%Y];
                }g[k][l] = g[k][l] / Ns;
                moy_g += g[k][l];
            }
        }

        /* Enregistrer les différents tableaux */
        enregistrer_tableau_g(g);
        afficher_gnuplot_g();

        /* Centrer le tableau 2D */
        double g_temp[X][Y];

        for (int k = 0; k < X; k++) {
            for (int l = 0; l < Y; l++) {
                int new_k = (k + X / 2) % X;
                int new_l = (l + Y / 2) % Y;
                g_temp[new_k][new_l] = g[k][l];
            }
        }

        // Copier les valeurs centrées du tableau temporaire dans le tableau g
        for (int k = 0; k < X; k++) {
            for (int l = 0; l < Y; l++) {
                g[k][l] = g_temp[k][l];
            }
        }

        /* Transformer en tableau 1D */
        int max_dist = round(sqrt((X * X) + (Y * Y))) + 1;
        double g_1D[max_dist][2];

        for (int r = 0; r < max_dist; r++) {
            g_1D[r][0] = 0;
            g_1D[r][1] = 0;
        }

        int x_center = X / 2;
        int y_center = Y / 2;

        for (int k = 0; k < X; k++) {
            for (int l = 0; l < Y; l++) {
                int dx = k - x_center;
                int dy = l - y_center;
                int r = round(sqrt(dx * dx + dy * dy)); // Distance au centre (arrondie à l'entier le plus proche)

                g_1D[r][0] += g[k][l];
                g_1D[r][1] += 1;
            }
        }

        // Normaliser les valeurs dans g_1D
        for (int r = 0; r < max_dist; r++) {
            if (g_1D[r][1] != 0) {
                g_1D[r][0] /= g_1D[r][1];
            }
        }

        enregistrer_tableau_g_1D(g_1D);

        const size_t n = 6;  /* number of data points to fit */
        const size_t p = 1;    /* number of model parameters */

        const gsl_rng_type * T = gsl_rng_default;

        FILE *fichier;
        fichier = fopen("g_1D.txt", "r");

        gsl_vector *f = gsl_vector_alloc(n);
        gsl_vector *x = gsl_vector_alloc(p);
        gsl_multifit_nlinear_fdf fdf;
        gsl_multifit_nlinear_parameters fdf_params =
        gsl_multifit_nlinear_default_parameters();
        struct data fit_data;
        gsl_rng * r;
        size_t i;

        gsl_rng_env_setup ();
        r = gsl_rng_alloc (T);

        fit_data.t = malloc(n * sizeof(double));
        fit_data.y = malloc(n * sizeof(double));
        fit_data.n = n;

        /* generate synthetic data with noise */
        for (i = 0; i < n; ++i)
            {
            fscanf(fichier, "%lf %lf", &fit_data.t[i], &fit_data.y[i]);
            }   
        fclose(fichier);

        /* define function to be minimized */
        fdf.f = func_f;
        fdf.df = func_df;
        fdf.fvv = func_fvv;
        fdf.n = n;
        fdf.p = p;
        fdf.params = &fit_data;

        /* starting point */
        gsl_vector_set(x, 0, 1.0);

        fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
        solve_system(x, &fdf, &fdf_params);

        /* print data and model and put it in a file */
        FILE *file = fopen("g_expo_T.txt", "a");
        {
            double A = gsl_vector_get(x, 0);
            fprintf(file, "%lf %lf \n",Temperature, A);
        }
        fclose(file);
        gsl_vector_free(f);
        gsl_vector_free(x);
        gsl_rng_free(r);

        Temperature = Temperature + 0.1;

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

void enregistrer_tableau_g(double tab[][Y]) {
    int l, m;
    FILE *fichier;
    fichier = fopen("g.txt", "w");
    for (l = 0; l < N; l++) {
        for (m = 0; m < N; m++) {
            fprintf(fichier, "%f ", tab[(l+(X/2))%X][(m+(Y/2))%Y]);
        }
        fprintf(fichier, "\n");
    }
    fclose(fichier);
}

void enregistrer_tableau_g_1D(double g_1D[][2]){
    int l;
    FILE *fichier;
    fichier = fopen("g_1D.txt", "w");
    for (l = 0; l < N/2; l++) {
        fprintf(fichier, " %d %f\n", l, g_1D[l][0]);
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

void afficher_gnuplot_tab() {
    system("gnuplot -e \"set terminal pngcairo size 800,800; set output 'tableau_visualisation.png'; set xrange [0:79]; set yrange [0:79]; set size square; set tics out nomirror; set key off; plot 'tableau.txt' matrix with image\"");
}

void afficher_gnuplot_g() {
    system("gnuplot -e \"set terminal pngcairo size 800,800; set output 'g.png'; set xrange [0:79]; set yrange [0:79]; set size square; set tics out nomirror; set key off; set cbrange [0:1]; plot 'g.txt' matrix with image\"");
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
 
