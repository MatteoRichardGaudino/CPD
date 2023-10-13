/**
 *  Autore: Matteo Richard Gaudino
 *  Matricola: N86003226
 * 
 **/

#include <omp.h> // omp_get_wtime, omp parallel, ... 
#include <stdio.h> // printf
#include <stdlib.h> // malloc, rand, srand, atoi, ...
#include <string.h> // memset
#include <time.h> // time

const unsigned long long MAX = 1UL << 16;

// Genera una matrice con numeri casuali di dimensione row X col
long long** randomMatrix(unsigned int row, unsigned int col);
// Genera un vettore con numeri casuali di dimensione size
long long*  randomVector(unsigned int size);
// Esegue il prodotto Righe per colonne della matrice Matrix per il vettore colonna Vector
// Il vettore risultante viene scritto in Result
void matricePerVettore(long long **Matrix, unsigned int row, unsigned int col, long long *Vector, long long *Result);

int main(int argc, char *argv[]) { 
  long long **Matrix; // Matrice
  long long *Vector; // Vettore 
  long long *Result; // Vettore risultante

  unsigned int ROW = atoi(argv[1]); // = Result_size
  unsigned int COL = atoi(argv[2]); // Righe e colonne della matrice. COL = Vector_size, ROW = Result_size

  double t1, t2; // tempo di inizio e di fine
  
  srand(time(NULL)); // Genera numeri casuali

  // Creazione della matrice e del vettore
  Matrix = randomMatrix(ROW, COL);
  Vector = randomVector(COL);

  // Allocazione e inizializzazione del vettore risultante
  Result = malloc(ROW * sizeof(long long));
  Result = memset(Result, 0, ROW);

  
  t1 = omp_get_wtime(); // Tempo inizio

  matricePerVettore(Matrix, ROW, COL, Vector, Result); // Calcolo del prodotto

  t2 = omp_get_wtime(); // Tempo fine

  printf("%f\n", t2-t1); // Stampa il tempo impiegato
  
  // Stampa il Vettore Result
/*for (unsigned int i = 0; i < ROW; i++){ // Stampa il risultato
    printf("%lld ", Result[i]); 
  }
  printf("\n");*/
  return 0;
}

// ------- Funzioni -------

long long** randomMatrix(unsigned int row, unsigned int col){

  long long **Matrix;

  Matrix = malloc(row * (sizeof(long long*))); // Matrice allocata come vettore di puntatori
  for (unsigned int i = 0; i < row; i++){
    Matrix[i] = malloc(col * sizeof(long long));  
  }

  // Riempimento con elementi casuali < MAX
  for(unsigned int i = 0; i < row; i++)
    for(unsigned int j = 0; j < col; j++)    
      Matrix[i][j]= rand()%MAX;
    
  return Matrix;
}

long long*  randomVector(unsigned int size){
  long long *Vector;
  Vector = malloc(size * sizeof(long long));	
  for(unsigned int i = 0; i < size; i++){ // Riempimento con elementi casuali < MAX
    Vector[i]= rand()%MAX; 
  }   
  return Vector;     
}

void matricePerVettore(long long **Matrix, unsigned int row, unsigned int col, long long *Vector, long long *Result){
  unsigned int i, j;
  // Sezone parallela
  #pragma omp parallel for shared(Matrix, row, col, Vector, Result) private(i, j)
  for(i = 0; i < row; i++){ // Le righe di Matrix sono suddivise sui vari thread
    for(j = 0; j < col; j++)
      Result[i] += Matrix[i][j] * Vector[j];
  }   
}
