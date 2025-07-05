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

/*Ici on va faire varier le nombre de particules Ns et analyser la fonction de corrélation temporelle*/

#define X 80
#define Y 80


struct data
{
  double *t;
  double *y;
  size_t n;
};


const int N = 80;
const int iteration = 100000; /*Nombre d'itérations*/
int calc_nb_voisins( int tab[][Y], int x, int y);
void enregistrer_tableau_tab(int tab[][Y]);
void enregistrer_tableau_g(double tab[][Y]);
int choisir_particule(int tab[][Y], int tab_1[][2], int *a, int *b, int Ns);
void choisir_voisin(int tab[][Y], int *a, int*b, int*c, int*d);
double calculer_energie_case(int x, int y, double J, int tab[][Y]);
void effectuer_saut(int tab[][Y],int tab_1[][2],int x,int y,int vx,int vy,double proba,double p1, double dE, int particule);
void udpate_E_moy(double Enr);
void afficher_gnuplot_tab();
void reinitialiser_fichier(const char *nom_fichier);
double calculer_energie_equilibre(const char *filename, int nb_energies);
void afficher_energie_equilibre_gnuplot();
double calculer_energie_configuration(int tab[][Y], double J);
void shuffle(int *array, int n);
double calculer_energie_etendue(int tab[][Y], int x, int y, double J);
void enregistrer_tableau_g_1D(double g_1D[][2]);
void enregistrer_tableau_g_moy(long double tab[]);
void afficher_g_fit();
double gaussian(const double a, const double b, const double c, const double t);
int func_f (const gsl_vector * x, void *params, gsl_vector * f);
int func_df (const gsl_vector * x, void *params, gsl_matrix * J);
void callback(const size_t iter, void *params,const gsl_multifit_nlinear_workspace *w);
void solve_system(gsl_vector *x, gsl_multifit_nlinear_fdf *fdf, gsl_multifit_nlinear_parameters *params);
int func_fvv (const gsl_vector * x, const gsl_vector * v,void *params, gsl_vector * fvv);

int main(){
    int Ns = 500;

    clock_t debut, fin;
    double temps_ecoule;

    // Enregistrer le temps de début
    debut = clock();

    //Initialisation des variables
    int l, j, particule_x, particule_y, voisx, voisy, iteration_count;
    const double kb = 1; //.39 * pow(10.0, -23);
    int randx;
    int randy;
    double Temperature = 1; /*Température*/
    double E1; /* nrj particule choisie*/
    double E2; /* nrj voisin */
    double dE; /*diff energie*/
    const double J = 1; /*constante J*/
    double p0; /*Proba de sauter au site voisin*/
    double p1; /*Pour calculer la proba*/
    int x_particule_corr;
    int y_partiule_corr;
    double nrj_t[iteration];


    /* Importer temps et régler la seed de drand */
    time_t *tm = NULL;
    srand48(time(tm));


    FILE *fichier2 = fopen("b_vs_Ns.txt", "w");
    while(Ns<6000){
    
    /* Initialisation des tableaux */
    int tab[X][Y]; /* Tableau principal */
    int indices[X * Y]; /* Tableau des indices pour rendre tab aléatoire à chaque nouvelle température */
    int tab_1[Ns][2]; /* Tableau qui tient à jour la positions de tous les 1 (toutes les particules) dans le tableau principal */
    int positions_initiales[Ns][2]; /* Tableau des positions initiales de quelques cases */
    int valeurs_initiales[Ns]; /* Valeurs initiales de ces quelques cases  */


    /* On crée le tableau des indices, pour l'instant chaque case contient son indice */
    for (l = 0; l < X * Y; l++) {
        indices[l] = l;
    } 
    shuffle(indices, X * Y); /* Mélange aléatoirement le tableau des indices */

        
    /* Rendre le tableau principal aléatoire */
    for (l = 0; l < X * Y; l++) {
        particule_x = indices[l] / X;
        particule_y = indices[l] % Y;
        if (l < Ns) {
            tab[particule_x][particule_y] = 1;
            tab_1[l][0] = particule_x; /* Mettre à jour le tableau qui répertorie la position des 1 */
            tab_1[l][1] = particule_y;
            positions_initiales[l][0] = particule_x;
            positions_initiales[l][1] = particule_y;
        } else {
            tab[particule_x][particule_y] = 0;
        }
    }


    /* Nous voulons créer un fichier qui va retenir l'énergie de la configuration donc nous avons besoin de la calculer au moins une fois, le reste du temps on va seulement ajouer dE*/
    double energie_config = calculer_energie_configuration(tab, J);

    int** g = malloc(Ns * sizeof(int*));
    for (l = 0; l < Ns; l++) {
        g[l] = malloc(iteration * sizeof(int));
    }
    
    /* Enregistrer ensuite la première valeur de ces cases */
    for(l=0; l<Ns; l++){
        particule_x = positions_initiales[l][0];
        particule_y = positions_initiales[l][1];
        if(tab[particule_x][particule_y] == 1){
            valeurs_initiales[l] = 1;
        }else{
            valeurs_initiales[l] = 0;
        }
    }

    /*Boucle principale*/
    for(iteration_count = 0; iteration_count<iteration ; iteration_count++){

        /* Prendre une particule au hasard */ 
        int particule = choisir_particule(tab, tab_1, &randx, &randy, Ns); /*La fonction renvoie un entier qui correspond au numéro de la particule. Si c'est la première elle aura le numéro 0 etc...*/        
        /* Prendre un voisin au hasard */
        choisir_voisin(tab, &voisx, &voisy, &randx, &randy); 
        /* Mettre à jour la corrélation */
            for(l = 0; l<Ns; l++){
                x_particule_corr = positions_initiales[l][0];
                y_partiule_corr = positions_initiales[l][1];
                if(tab[x_particule_corr][y_partiule_corr] == 1){
                    g[l][iteration_count] = valeurs_initiales[l]; // Ici on peut changer les 0 en -1 pour les corrélations si on le souhaite
                }else{
                    g[l][iteration_count] = 0;
                }
            }
            nrj_t[iteration_count] = energie_config;

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

    /*Maintenant que mon tableau g est rempli on peut moyenner sur touts les cases choisies (sur NS )*/
    long double g_moy[iteration];
    for(l=0; l<iteration; l++){
        g_moy[l] = 0;
        for(j=0; j<Ns; j++){
            g_moy[l] += g[j][l];
        }
       g_moy[l] /= Ns;
       
    }
    
    /*Enregistrer g_moy*/
    enregistrer_tableau_g_moy(g_moy);

    /* Afficher les tableaux sur gnuplot */
    enregistrer_tableau_tab(tab);
    afficher_gnuplot_tab();


    for (l = 0; l < Ns; l++) {
    free(g[l]);
    }   
    free(g);
    
    const size_t n = iteration;  /* number of data points to fit */
    const size_t p = 3;    /* number of model parameters */

    const gsl_rng_type * T = gsl_rng_default;

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

    /* Generate data  */
    for (l = 0; l < n; ++l)
    {
        fit_data.t[l] = l;
        fit_data.y[l] = g_moy[l];
    }   

    /* define function to be minimized */
    fdf.f = func_f;
    fdf.df = func_df;
    fdf.fvv = func_fvv;
    fdf.n = n;
    fdf.p = p;
    fdf.params = &fit_data;

    /* starting point */
    gsl_vector_set(x, 0, 1.0);
    gsl_vector_set(x, 1, 0.0);
    gsl_vector_set(x, 2, 1.0);

    fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
    solve_system(x, &fdf, &fdf_params);

    /* print data and model and put it in a file */
    FILE *file = fopen("g_gaussian.txt", "w");
    
    {
        double A = gsl_vector_get(x, 0);
        double B = gsl_vector_get(x, 1);
        double C = gsl_vector_get(x, 2);
        double densite = (double)Ns / (X * Y) ;
        double B_inv = 1/B;

        fprintf(fichier2, "%f %lf \n", densite, B_inv);

        for (l = 0; l < iteration; l++)
      {
        double yi = g_moy[l];
        double fi = gaussian(A, B, C, l);

        fprintf(file, "%d %lf %lf \n", l, yi, fi);
        
      }
    }
    fclose(file);
    
    gsl_vector_free(f);
    gsl_vector_free(x);
    gsl_rng_free(r);
    Ns += 500;
    }
    fclose(fichier2);
    afficher_g_fit();
    fin = clock();
    temps_ecoule = (double)(fin - debut) / CLOCKS_PER_SEC;
    printf("Temps écoulé : %.5f secondes\n", temps_ecoule);

    return 0;
}


/* model function: a * exp(-t*b)+c */
double
gaussian(const double a, const double b, const double c, const double t)
{
  const double z = t*b;
  return (a * exp(-z) + c);
} 

int
func_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  double b = gsl_vector_get(x, 1);
  double c = gsl_vector_get(x, 2);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double yi = d->y[i];
      double y = gaussian(a, b, c, ti);

      gsl_vector_set(f, i, yi - y);
    }

  return GSL_SUCCESS;
}

int
func_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  double b = gsl_vector_get(x, 1);
  double c = gsl_vector_get(x, 2);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double zi = b*ti;
      double ei = exp(-zi);

      gsl_matrix_set(J, i, 0, -ei);
      gsl_matrix_set(J, i, 1, a*ti*ei);
      gsl_matrix_set(J, i, 2, -1);
    }

  return GSL_SUCCESS;
}

int
func_fvv (const gsl_vector * x, const gsl_vector * v,
          void *params, gsl_vector * fvv)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  double b = gsl_vector_get(x, 1);
  double c = gsl_vector_get(x, 2);
  double va = gsl_vector_get(v, 0);
  double vb = gsl_vector_get(v, 1);
  double vc = gsl_vector_get(v, 2);
  size_t i;

  for (i = 0; i < d->n; ++i)
    {
      double ti = d->t[i];
      double zi = (ti - b) / c;
      double ei = exp(-0.5 * zi * zi);
      double Dab = -ti*ei;
      double Dac = -0;
      double Dbb = -a*ti*ti*ei;
      double Dbc = 0;
      double Dcc = 0;
      double sum;

      sum = 2.0 * va * vb * Dab +
            2.0 * va * vc * Dac +
                  vb * vb * Dbb +
            2.0 * vb * vc * Dbc +
                  vc * vc * Dcc;

      gsl_vector_set(fvv, i, sum);
    }

  return GSL_SUCCESS;
}

void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double avratio = gsl_multifit_nlinear_avratio(w);
  double rcond;

  (void) params; /* not used */

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stderr, "iter %2zu: a = %.4f, b = %.4f, c = %.4f, |a|/|v| = %.4f cond(J) = %8.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_vector_get(x, 1),
          gsl_vector_get(x, 2),
          avratio,
          1.0 / rcond,
          gsl_blas_dnrm2(f));
}

void
solve_system(gsl_vector *x, gsl_multifit_nlinear_fdf *fdf,
             gsl_multifit_nlinear_parameters *params)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  const size_t max_iter = 200;
  const double xtol = 1.0e-8;
  const double gtol = 1.0e-8;
  const double ftol = 1.0e-8;
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
  fprintf(stderr, "NFEV          = %zu\n", fdf->nevalf);
  fprintf(stderr, "NJEV          = %zu\n", fdf->nevaldf);
  fprintf(stderr, "NAEV          = %zu\n", fdf->nevalfvv);
  fprintf(stderr, "initial cost  = %.12e\n", chisq0);
  fprintf(stderr, "final cost    = %.12e\n", chisq);
  fprintf(stderr, "final x       = (%.12e, %.12e, %12e)\n",
          gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2));
  fprintf(stderr, "final cond(J) = %.12e\n", 1.0 / rcond);

  gsl_multifit_nlinear_free(work);
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

void enregistrer_tableau_g_moy(long double tab[]) {
    int l;
    FILE *fichier;
    fichier = fopen("g_moy.txt", "w");
    for (l = 0; l < iteration; l++) {
        fprintf(fichier, "%d %lf \n ", l, tab[l]);
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

int choisir_particule(int tab[][Y], int tab_1[][2], int *a, int *b, int Ns){

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

void reinitialiser_fichier(const char *nom_fichier) {
    FILE *fichier;

    // Ouvrir le fichier en mode écriture
    fichier = fopen(nom_fichier, "w");

    if (fichier == NULL) {
        printf("Erreur lors de l'ouverture du fichier.\n");
        return;
    }

    fclose(fichier);
}

void afficher_g_fit(){
    system("gnuplot -e \"set terminal pngcairo size 800,800; set output 'gaussian_visualisation.png'; plot 'g_gaussian.txt' using 1:2, 'g_gaussian.txt' using 1:3\"");
}