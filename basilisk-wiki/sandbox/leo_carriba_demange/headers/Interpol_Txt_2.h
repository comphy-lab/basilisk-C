// 20/04/2021
// Leo Carriba Demange
// Fonction pour interpoler f(x) a partir d'un jeu de donnees ( t, f(t) ) dans un fichier texte
// Fonction pour interpoler la derivee du jeu de donees ( t, f(t) ).

# define LONGUEUR_NOM_FICHIER 200 // longueur maximale du nom du fichier texte en entree

double interpolation_lineaire_textfile(char nom_fichier[LONGUEUR_NOM_FICHIER], double t)
// "nom_fichier" est le nom du fichier texte
// "t" est le temps ou on veut connaitre la hauteur d'eau ou le debit
{
    // declaration des variables
    double temps = -1.;
    double fonction = -1.;
    double temps_precedent = -1.;
    double fonction_precedent = -1.;
    double temps_suivant = -1.;
    double fonction_suivant = -1.;
    double fonction_retourne = -1.;

    fonction_retourne = -999.; // valeur retournee par defaut s'il y a un probleme

    // ouverture du fichier texte
    FILE* fichier = NULL;
    //char ligne[TAILLE_MAX_lIGNE] = "";
    fichier = fopen(nom_fichier, "r");

    if (fichier == NULL)
    {
        printf("erreur ! le fichier n'a pas pu etre ouvert a l'adresse %s\n", nom_fichier);
    }

    if (fichier != NULL)
    {
        if (fscanf(fichier, "%lg %lg", &temps, &fonction) != EOF) // On lit le fichier tant qu'on ne recoit pas d'erreur (NULL)
        {
            // on initialise l'instant precedent
            temps_precedent = temps;
            fonction_precedent = fonction;
        }

        while (fscanf(fichier, "%lg %lg", &temps, &fonction) != EOF) // On lit le fichier tant qu'on ne recoit pas d'erreur (NULL)
        {
            // on recupere l'instant suivant
            temps_suivant = temps;
            fonction_suivant = fonction;

            if (temps_precedent <= t && t <= temps_suivant)
            {
                // interpolation lineaire
                double theta = (t - temps_precedent) / (temps_suivant - temps_precedent);
                fonction_retourne = fonction_suivant * theta + fonction_precedent * (1 - theta);
                break;
            }

            // on prepare le tour de boucle suivant
            temps_precedent = temps_suivant;
            fonction_precedent = fonction_suivant;
        }

        // fermeture du fichier
        fclose(fichier);
    }


    return(fonction_retourne);
}

double derivee_textfile_forward(char nom_fichier[LONGUEUR_NOM_FICHIER], double t, double delta_t)
// "nom_fichier" est le nom du fichier texte
// "t" est le temps ou on veut connaitre la derivee de la hauteur d'eau ou du debit
// "delta_t" est l'intervalle sur lequel sera calculee la derivee
// La derivee est calculee par un schema forward
{
    // declaration des variables
    double derivee_retourne = -999.; // valeur retournee par defaut s'il y a un probleme
    double fonction_avant, fonction_apres;

    // evaluation de la fonction
    fonction_avant = interpolation_lineaire_textfile(nom_fichier, t);
    fonction_apres = interpolation_lineaire_textfile(nom_fichier, t + delta_t);

    // calcul de la derivee
    derivee_retourne = (fonction_apres - fonction_avant) / delta_t;

    return(derivee_retourne);
}




