#include <stdio.h>
#include <string.h>
#include "limace.h"
#include "erreurs.h"
#include "tai.h"


/* Codes de retour */
#define OK 0
#define HELP 1
#define ERR_NB_PARAM 2
#define ERR_MAT_IN 3
#define ERR_MAT_TYPE 4
#define ERR_MAT_DIM 5

void Syntaxe(void)
{
  Usage("Histogramme.mx\n"
        "-h\n");
}

void Description(void)
{
  Mesg("RÔLE\n");
  Mesg("\tCalcul du seuil de binarisation par la méthode d'Otsu\n");
  Mesg("ARGUMENTS\n");
  Mesg("\tHistogramme.mx : nom du fichier source (Matrix de Int 1 x 256)\n");
  Mesg("OPTION\n");
  Mesg("\t-h : affichage de l'aide\n");
  Mesg("CODES DE RETOUR\n");
  Code(OK,"succès");
  Code(HELP,"aide demandée");
  Code(ERR_NB_PARAM,"nombre d'arguments incorrect");
  Code(ERR_MAT_IN,"erreur durant la lecture du fichier en entrée");
  Code(ERR_MAT_TYPE,"mauvais type de matrice (Matrix de Int)");
  Code(ERR_MAT_DIM,"dimensions de matrice incorrectes (1x256)");
}

int main(int argc, char *argv[])
{
  Matrix Hist=NULL;
  int Seuil;

  /* Initialisation de la gestion des messages d'erreur */
  InitMesg(argv);
  /* Vérification du nombre d'arguments */
	if (argc!=2)
	{
	  Syntaxe();
	  return ERR_NB_PARAM;
	}
	if (argc==2)
	{
	  if (!strcmp(argv[1],"-h"))
	  {
	    Syntaxe();
	    Description();
	    return HELP;
	  }
	}
	/* Lecture de l'histogramme */
	Hist=MatReadAsc(argv[1]);
	if (Hist==NULL)
	{
		Erreur("Problème lors de la lecture de l'histogramme en entrée");
		return ERR_MAT_IN;
	}
	/* Vérification du type de la matrice */
	if (MatType(Hist)!=Int)
	{
    Erreur("L'histogramme doit être stocké dans une matrice de Int");
		MatFree(&Hist);
		return ERR_MAT_TYPE;
	}
	/* Vérification des dimensions de la matrice */
	if ((MatNbRow(Hist)!=1 || MatNbCol(Hist)!=256)
		&& (MatNbRow(Hist)!=256 || MatNbCol(Hist)!=1))
	{
		Erreur("L'histogramme doit être contenu dans une matrice 1x256 ou 256x1");
		MatFree(&Hist);
		return ERR_MAT_DIM;
	}
	/* Calcul du seuil d'Otsu */
	Seuil=Otsu(Hist);

	/* Écriture du seuil sur stdout */
	printf("%d\n",Seuil);
	/* Libération de la memoire */
	MatFree(&Hist);
	return OK;
}
