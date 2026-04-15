#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "limace.h"
#include "tai.h"
#include <limits.h>

#define DEBOGAGE

#ifdef DEBOGAGE
#define DEBUG fprintf(stderr,"Fichier %s, ligne %d\n",__FILE__,__LINE__);
#else
#define DEBUG
#endif

#define AFAIRE(ValeurRetour) \
  fprintf(stderr,"--> Fichier %s, ligne %d, corps de la fonction %s à écrire.\n",\
          __FILE__,__LINE__,__func__);\
  return ValeurRetour;  


/*
 * Conversion d'une image couleur en une image de niveaux de gris
 * Entrée : image initiale en couleur
 * Sortie : image de niveaux de gris résultat
 */
Image RGB2Gray(Image Im)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  Image Gray = ImAlloc(GrayLevel,NbRow,NbCol);
  unsigned char** gray = ImGetI(Gray);
  unsigned char** R = ImGetR(Im);
  unsigned char** G = ImGetG(Im);
  unsigned char** B = ImGetB(Im);
  for (int i =0;i<NbRow;i++){
    for (int j=0; j<NbCol;j++){
     gray[i][j]= round((0.299*R[i][j])+(0.587*G[i][j])+(0.114*B[i][j]));
    }
  }
  return Gray;
}


/*
 * Binarisation d'une image de niveaux de gris par seuillage global
 * Entrées : image de niveaux de gris initiale
             seuil (niveau de gris)
 * Sortie : image binaire
 */
Image Binarization(Image Im, unsigned char Threshold)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  unsigned char** gray = ImGetI(Im);
  Image Bin = ImAlloc(BitMap,NbRow,NbCol);
  unsigned char** bin = ImGetI(Bin);
  for (int i =0;i<NbRow;i++){
    for (int j=0; j<NbCol;j++){
      if (gray[i][j]<Threshold)
        bin[i][j]=0;
      else 
        bin[i][j]=1;
    }
  }
  return Bin;
}


/*
 * Inversion d'une image
 * Entrée : image initiale (binaire, niveaux de gris ou couleur)
 * Sortie : image résultat
 */
Image Inversion(Image Im)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  ImageType type = ImType(Im);
  if (type==Color){
    unsigned char** R = ImGetR(Im);
    unsigned char** G = ImGetG(Im);
    unsigned char** B = ImGetB(Im);
    for (int i =0;i<NbRow;i++){
      for (int j=0; j<NbCol;j++){
        R[i][j]=255-R[i][j];
        G[i][j]=255-G[i][j];
        B[i][j]=255-B[i][j];
      }
    }
  } else {
    unsigned char** I = ImGetI(Im);
    for (int i =0;i<NbRow;i++){
      for (int j=0; j<NbCol;j++){
        if (type==GrayLevel){
          I[i][j]=255-I[i][j];
        } else {
          I[i][j]=1-I[i][j];
        }
      }
    }
  }
  return Im;
}


/*
 * Calcul de l'histogramme d'une image de niveaux de gris
 * Entrée : image initiale (niveaux de gris)
 * Sortie : histogramme (matrice de int 1 x 256)
 */
Matrix Histogram(Image Im)
{
  Matrix mat =  MatAlloc(Int,1,256);
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  int ** M=MatGetInt(mat);
  for (int i=0; i < 256; i++){M[0][i]=0;}
  unsigned char** I = ImGetI(Im);
  for (int i =0; i<NbRow;i++){
    for (int j = 0; j<NbCol;j++){
      M[0][I[i][j]]+=1;
    }
  }
  return mat;
}


/*
 * Représentation d'un histogramme sous forme d'une image
 * Entrées : histogramme (matrice de int 1 x 256) et nombre de lignes de
 * l'image résultat (une échelle des niveaux de gris de 25 lignes est ajoutée
 * sous l'histogramme)
 * Sortie : image de niveaux de gris résultat
 */
Image Hist2Im(Matrix Hist, int NbLig)
{
	unsigned char **I;
	int *h,i,j,Max=0,NbCol=256,NbLig2=NbLig+25;
	Image Res;

	if (MatType(Hist)!=Int) return NULL;
  NbLig2=NbLig+25;
	Res=ImCAlloc(GrayLevel,NbLig2,NbCol);
  if (Res==NULL) return NULL;
	h=*MatGetInt(Hist);
	for (j=0;j<NbCol;j++)
		if (h[j]>Max) Max=h[j];
	I=ImGetI(Res);
	for (j=0;j<256;j++)
		for (i=NbLig-1;i>=(NbLig-NbLig*h[j]/Max);i--)
		    I[i][j]=255;
  for (j=0;j<256;j++)
    I[NbLig][j]=0;
  for (i=NbLig+1;i<NbLig2;i++)
    for (j=0;j<256;j++)
      I[i][j]=j;
	return Res;
}


/*
 * Calcul du seuil d'Otsu
 * Entrée : histogramme (matrice de int 1 x 256)
 * Sortie : seuil (niveau de gris) obtenu par la méthode d'Otsu
 */

 
int Cout(Matrix H, int S) {
  int **M = MatGetInt(H);
  int q1 = 0, q2 = 0, u1 = 0, u2 = 0;

  for (int k = 0; k < 256; k++) {
    if (k < S) {
      q1 += M[0][k];
      u1 += k * M[0][k];
    }else { 
      q2 += M[0][k];
      u2 += k * M[0][k];
    }
  }

  if (q1 == 0 || q2 == 0) return INT_MAX;
  u1 = u1 / q1;
  u2 = u2 / q2;

  int left = 0, right = 0;
  for (int k = 0; k < 256; k++) {
    if (k < S) 
      left  += M[0][k] * ((k - u1) * (k - u1));
    else 
      right += M[0][k] * ((k - u2) * (k - u2));
  }
  return left + right;
}


unsigned char Otsu(Matrix Hist){
  int seuil = 1;
  int min = Cout(Hist,1);
  for (int i = 2; i< 255; i++){
    int c = Cout(Hist,i);
    if (c<min){
      min = c;
      seuil =i;
    }
  }
  return seuil;
}



/*
 * Calcul de l'histogramme cumulé à partir de l'histogramme
 * Entrée : histogramme (matrice de int 1 x 256)
 * Sortie : histogramme cumulé (matrice de int 1 x 256)
 */
Matrix Hist2CumHist(Matrix Hist)
{
  Matrix mat =  MatAlloc(Int,1,256);
  int NbCol = MatNbCol(Hist);
  int ** M=MatGetInt(mat);
  int ** H=MatGetInt(Hist);
  M[0][0]=H[0][0];
  for (int i =1; i<NbCol;i++){
    M[0][i]= M[0][i-1]+H[0][i];
  }
  return mat;
}


/*
 * Application d'une transformation ponctuelle à une image de niveaux de gris
 * Entreés : image initiale (niveaux de gris) et
 * transformation ponctuelle (matrice de int 1 x 256)
 * Sortie : image de niveaux de gris transformée
 */
Image AppLUT(Image Im, Matrix LUT)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  int ** L=MatGetInt(LUT);
  Image App = ImAlloc(GrayLevel,NbRow,NbCol);
  unsigned char** A = ImGetI(App);
  unsigned char** S = ImGetI(Im);
  for (int i =0;i<NbRow;i++){
    for (int j=0; j<NbCol;j++){
      A[i][j]=L[0][S[i][j]];
    }
  }
  return App;
}


/*
 * Spécification d'histogramme
 * Entrées : histogramme cumulé de l'image et histogramme cumulé desiré
 * (on suppose que le dernier élément des deux histogrammes cumulés sont
 * les mêmes, c'est-à-dire qu'ils décrivent des images contenant le même nombre
 * de pixels)
 * Sortie : transformation ponctuelle (matrice 1 x 256)
 */
Matrix HistSpecif(Matrix CumHist, Matrix DesCumHist)
{
  Matrix mat =  MatAlloc(Int,1,256);
  int NbCol = MatNbCol(CumHist);
  int ** M=MatGetInt(mat);
  int ** I=MatGetInt(CumHist);
  int ** D=MatGetInt(DesCumHist);
  for (int i =0; i<NbCol;i++){
    int niv = 0;
    int best_diff = abs(I[0][i]-D[0][0]);
    for (int j=1; j<NbCol; j++){
      int diff = abs(I[0][i]-D[0][j]);
      if (diff <best_diff){
        niv = j;
        best_diff = diff;
      }
    }
    M[0][i]=niv; 
  }
  return mat;  
}


/*
 * Vérification de la validité d'une matrice représentant un élément
 * structurant binaire (pour l'érosion, la dilatation, etc.)
 * Entrée : matrice représentant un élément structurant
 * Sortie : 0 si la matrice est valide,
            SE_NOT_ODD si son nombre de lignes ou de colonnes n'est pas impair
            SE_NOT_INT si elle ne contient pas que des entiers
            SE_NOT_BIN si elle ne contient pas que des 0 et des 1
*/
int NotValidBinSE(Matrix StructuringElement)
{
  int **ES,NbLig,NbCol,i,j;

  if (MatType(StructuringElement)!=Int)
    return SE_NOT_INT;
  NbLig=MatNbRow(StructuringElement);
	if ((NbLig%2)!=1)
	  return SE_NOT_ODD;
  NbCol=MatNbCol(StructuringElement);
	if ((NbCol%2)!=1)
	  return SE_NOT_ODD;
  ES=MatGetInt(StructuringElement);
  for (i=0;i<NbLig;i++)
    for (j=0;j<NbCol;j++)
      if (ES[i][j]!=0 && ES[i][j]!=1)
        return SE_NOT_BIN;
  return 0;
}


/*
 * Vérification de la validité d'une matrice représentant un élément
 * structurant ternaire (pour la transformation "tout ou rien")
 * Entrée : matrice représentant un élément structurant
 * Sortie : 0 si la matrice est valide,
            SE_NOT_ODD si son nombre de lignes ou de colonnes n'est pas impair
            SE_NOT_INT si elle ne contient pas que des entiers
            SE_NOT_TERN si elle ne contient pas que des 0, des 1 et des 2
*/
int NotValidTernSE(Matrix StructuringElement)
{
  int **ES,NbLig,NbCol,i,j;

  if (MatType(StructuringElement)!=Int)
    return SE_NOT_INT;
  NbLig=MatNbRow(StructuringElement);
	if ((NbLig%2)!=1)
	  return SE_NOT_ODD;
  NbCol=MatNbCol(StructuringElement);
	if ((NbCol%2)!=1)
	  return SE_NOT_ODD;
  ES=MatGetInt(StructuringElement);
  for (i=0;i<NbLig;i++)
    for (j=0;j<NbCol;j++)
      if (ES[i][j]!=0 && ES[i][j]!=1 && ES[i][j]!=2)
        return SE_NOT_TERN;
  return 0;
}


/*
 * Amincissement d'une image binaire
 * Entreés : image binaire initiale et élément structurant (matrice de int
 * contenant uniquement des 0, des 1 et des 2 signifiant "peu importe")
 * Sortie : image binaire transformée
 */
Image Thinning(Image Im, Matrix StructuringElement)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  int MRow = (MatNbRow(StructuringElement)-1)/2;
  int MCol = (MatNbCol(StructuringElement)-1)/2;
  int ** E = MatGetInt(StructuringElement);

  unsigned char** I = ImGetI(Im);
  Image Copy = ImAlloc(BitMap, NbRow, NbCol);
  unsigned char **C = ImGetI(Copy);
  for (int i = 0; i < NbRow; i++)
    for (int j = 0; j < NbCol; j++)
      C[i][j] = I[i][j];

  for (int i = 1 ; i < NbRow-1; i++){
    for (int j = 1; j < NbCol-1 ; j++){
      int app = 1;
      for (int ie = -MRow; ie <= MRow; ie++){
        for (int je = -MCol; je <= MCol; je++){
          if(E[ie+MRow][je+MCol]!=2){
            app = app && I[i+ie][j+je]==E[ie+MRow][je+MCol];
          }
        }
      }
      C[i][j] = app ? 0 : I[i][j];
    }
  }
  return Copy;
}


/*
 * Erosion d'une image binaire
 * Entrées : image binaire et élément structurant (0, 1, 2="peu importe")
 * Sortie : image binaire érodée
 */
Image Erosion(Image Im, Matrix StructuringElement)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  int MRow = (MatNbRow(StructuringElement)-1)/2;
  int MCol = (MatNbCol(StructuringElement)-1)/2;
  int **E = MatGetInt(StructuringElement);
  unsigned char **I = ImGetI(Im);

  Image Out = ImAlloc(BitMap, NbRow, NbCol);
  unsigned char **O = ImGetI(Out);

  for (int i = 0; i < NbRow; i++)
    for (int j = 0; j < NbCol; j++)
      O[i][j] = I[i][j];

  for (int i = MRow; i < NbRow-MRow; i++) {
    for (int j = MCol; j < NbCol-MCol; j++) {
      int match = 1;
      for (int ie = -MRow; ie <= MRow; ie++) {
        for (int je = -MCol; je <= MCol; je++) {
          if (E[ie+MRow][je+MCol] == 1) {  // only check 1s, ignore 0s and 2s
            match = match && (I[i+ie][j+je] == 1);
          }
        }
      }
      O[i][j] = match ? 1 : 0;
    }
  }
  return Out;
}

/*
 * Dilatation d'une image binaire
 * Entrées : image binaire et élément structurant (0, 1, 2="peu importe")
 * Sortie : image binaire dilatée
 */
Image Dilatation(Image Im, Matrix StructuringElement)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  int MRow = (MatNbRow(StructuringElement)-1)/2;
  int MCol = (MatNbCol(StructuringElement)-1)/2;
  int **E = MatGetInt(StructuringElement);
  unsigned char **I = ImGetI(Im);

  Image Out = ImAlloc(BitMap, NbRow, NbCol);
  unsigned char **O = ImGetI(Out);

  for (int i = 0; i < NbRow; i++)
    for (int j = 0; j < NbCol; j++)
      O[i][j] = I[i][j];

  for (int i = MRow; i < NbRow-MRow; i++) {
    for (int j = MCol; j < NbCol-MCol; j++) {
      int match = 0;
      for (int ie = -MRow; ie <= MRow; ie++) {
        for (int je = -MCol; je <= MCol; je++) {
          if (E[ie+MRow][je+MCol] == 1) {
            // at least one pixel under SE is 1
            match = match || (I[i+ie][j+je] == 1);
          }
        }
      }
      O[i][j] = match ? 1 : 0;
    }
  }
  return Out;
}

/*
 * Ouverture : érosion suivie d'une dilatation
 */
Image Ouverture(Image Im, Matrix StructuringElement)
{
  Image eroded = Erosion(Im, StructuringElement);
  Image opened = Dilatation(eroded, StructuringElement);
  return opened;
}

/*
 * Fermeture : dilatation suivie d'une érosion
 */
Image Fermeture(Image Im, Matrix StructuringElement)
{
  Image dilated = Dilatation(Im, StructuringElement);
  Image closed  = Erosion(dilated, StructuringElement);
  return closed;
}

/*
 * Transformation "Tout ou Rien" (Hit-or-Miss)
 * Marque les pixels où le SE correspond exactement (1→1, 0→0, 2=peu importe)
 */
Image ToutOuRien(Image Im, Matrix StructuringElement)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  int MRow = (MatNbRow(StructuringElement)-1)/2;
  int MCol = (MatNbCol(StructuringElement)-1)/2;
  int **E = MatGetInt(StructuringElement);
  unsigned char **I = ImGetI(Im);

  Image Out = ImAlloc(BitMap, NbRow, NbCol);
  unsigned char **O = ImGetI(Out);

  for (int i = 0; i < NbRow; i++)
    for (int j = 0; j < NbCol; j++)
      O[i][j] = 0;

  for (int i = MRow; i < NbRow-MRow; i++) {
    for (int j = MCol; j < NbCol-MCol; j++) {
      int match = 1;
      for (int ie = -MRow; ie <= MRow; ie++) {
        for (int je = -MCol; je <= MCol; je++) {
          if (E[ie+MRow][je+MCol] != 2) {
            match = match && (I[i+ie][j+je] == E[ie+MRow][je+MCol]);
          }
        }
      }
      O[i][j] = match ? 1 : 0;
    }
  }
  return Out;
}

/*
 * Top-hat blanc : Im - Ouverture(Im, SE)
 * Détecte les petites régions blanches / saillances lumineuses
 */
Image TopHatBlanc(Image Im, Matrix StructuringElement)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  unsigned char **I = ImGetI(Im);

  Image opened = Ouverture(Im, StructuringElement);
  unsigned char **O = ImGetI(opened);

  Image Out = ImAlloc(BitMap, NbRow, NbCol);
  unsigned char **R = ImGetI(Out);

  for (int i = 0; i < NbRow; i++)
    for (int j = 0; j < NbCol; j++)
      R[i][j] = I[i][j] && !O[i][j]; // Im - Opening

  return Out;
}

/*
 * Top-hat noir : Fermeture(Im, SE) - Im
 * Détecte les petites régions noires / trous sombres
 */
Image TopHatNoir(Image Im, Matrix StructuringElement)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  unsigned char **I = ImGetI(Im);

  Image closed = Fermeture(Im, StructuringElement);
  unsigned char **C = ImGetI(closed);

  Image Out = ImAlloc(BitMap, NbRow, NbCol);
  unsigned char **R = ImGetI(Out);

  for (int i = 0; i < NbRow; i++)
    for (int j = 0; j < NbCol; j++)
      R[i][j] = C[i][j] && !I[i][j]; // Closing - Im

  return Out;
}

/*
 * Épaississement : inverse du thinning
 * Ajoute des pixels là où le SE correspond (Tout ou Rien → union avec Im)
 */
Image Epaississement(Image Im, Matrix StructuringElement)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  unsigned char **I = ImGetI(Im);

  Image tor = ToutOuRien(Im, StructuringElement);
  unsigned char **T = ImGetI(tor);

  Image Out = ImAlloc(BitMap, NbRow, NbCol);
  unsigned char **O = ImGetI(Out);

  for (int i = 0; i < NbRow; i++)
    for (int j = 0; j < NbCol; j++)
      O[i][j] = I[i][j] || T[i][j]; // union: add matched pixels

  return Out;
}

/*
 * Détection de contours
 * Un pixel est un contour s'il est à 1 et a au moins un voisin à 0 (4-connexité)
 * Equivalent à : Im - Erosion(Im, SE_croix)
 */
Image Contour(Image Im)
{
  int NbRow = ImNbRow(Im);
  int NbCol = ImNbCol(Im);
  unsigned char **I = ImGetI(Im);

  Image Out = ImAlloc(BitMap, NbRow, NbCol);
  unsigned char **O = ImGetI(Out);

  // 4-connectivity neighbors
  int di[] = {-1, 1,  0, 0};
  int dj[] = { 0, 0, -1, 1};

  for (int i = 0; i < NbRow; i++)
    for (int j = 0; j < NbCol; j++)
      O[i][j] = 0;

  for (int i = 1; i < NbRow-1; i++) {
    for (int j = 1; j < NbCol-1; j++) {
      if (I[i][j] == 1) {
        int border = 0;
        for (int d = 0; d < 4; d++) {
          if (I[i+di[d]][j+dj[d]] == 0)
            border = 1;
        }
        O[i][j] = border ? 1 : 0;
      }
    }
  }
  return Out;
}