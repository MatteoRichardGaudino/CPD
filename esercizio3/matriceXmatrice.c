#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <time.h>

float* randomMatrix(int m, int n);
double matrXmatr(const float *A, int m, int n, const float *B, int k, float *C);
void crea_griglia(int dim, MPI_Comm *comm_grid, MPI_Comm *comm_row);

// ----------------------------------------

int main(int argc, char *argv[]){
   int worldRank, worldSize; // rank all'interno del COMM_WORLD e dimensione del comunicatore
   int gridRank, gridCoords[2]; // rank e coordinate all'interno del comunicatore griglia

   int m; // Righe e Colonne delle matrici
   int mSub; // Righe e colonne delle sottomatrici
   int matrixSize; // m*m
   int subMatrixSize; // m*m/worldSize

   float *A, *B; // Matrici di input
   float *C; // Matrice risultato
   float *subA, *subB, *subC; // Sottomatrici locali
   float *Temp; // Matrice di appoggio. Necessaria per non sovrascrivere subA durante la comunicazione in BMR

   double t = 0.F, start, end; // Tempo per impiegato per il calcolo e tempo totale con overhead

   MPI_Status status; // per send
   MPI_Request rqst; // per isend
   MPI_Comm comm_grid, comm_row; // Comunicatore di griglia e di riga



   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
   MPI_Comm_size(MPI_COMM_WORLD, &worldSize);


//------------------------ Controlli sull'input ---------------------------------
   int p = (int) sqrt(worldSize); // Lato della griglia quadrata di processori
   // Il numero di processori deve essere un quadrato perfetto per la BMR
   if (p*p != worldSize){
      perror("ERROR Processors must be perfect square.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
   }

   if(argc <= 1){
       perror("ERROR Argument required. usage: ./program <matrixDim>\n");
       MPI_Abort(MPI_COMM_WORLD, 1);
       return 1;
   }

   if(worldRank==0){      
      m = atoi(argv[1]);  // Righe e colonne delle matrici
      // Il numero di righe/colonne delle matrici deve essere multiplo di p
      if(m % p != 0){
         perror("ERROR rows/cols of the matrices must be multiple of sqrt nProc.\n");
         MPI_Abort(MPI_COMM_WORLD, 1);
         return 1;
      }
   }
//------------------------Generazione Matrici -----------------------------

    if (worldRank == 0){
      srand(time(NULL));   
      // Generazione delle matrici random
      A = randomMatrix(m, m);
      B = randomMatrix(m, m);
      // Alloca spazio per la matrice risultato
      C = malloc(m*m*sizeof(float));
   }

// ---------------------- Comunicazione input -----------------------------

   start = MPI_Wtime();

   MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&p, 1, MPI_INT, 0, MPI_COMM_WORLD);

    matrixSize = m*m;
    subMatrixSize = matrixSize/worldSize;
    mSub = m/p;


    crea_griglia(p, &comm_grid, &comm_row); // Crea la griglia di processori
   
   MPI_Comm_rank(comm_grid, &gridRank);
   MPI_Cart_coords(comm_grid, gridRank, 2, gridCoords);
   

   subA = malloc(subMatrixSize * sizeof(float)); // Sottomatrice di A
   Temp = malloc(subMatrixSize * sizeof(float)); // Blocco di appoggio
   subB = malloc(subMatrixSize * sizeof(float)); // Sottomatrice di B
   subC = malloc(subMatrixSize * sizeof(float)); // Sottomatrice di C
   memset(subC, 0, subMatrixSize); // Riempie C di zeri
    
   if(worldRank==0){
      int colOffset, rowOffset;
      for(int i = 1; i < worldSize; i++){
         rowOffset = (mSub) * (i/p); // Riga da cui partire
         colOffset = (mSub) * (i%p); // Colonna da cui partire
         for(int j = 0; j < mSub; j++){
            MPI_Send(A + (rowOffset * m) + colOffset, mSub, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            MPI_Send(B + (rowOffset * m) + colOffset, mSub, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            rowOffset++; // Prossima riga
         }
      }

      // Copia di subA
      for(int j = 0; j < mSub; j++){
         for(int z = 0; z < mSub; z++)
            *(subA + (mSub*j) + z) = *(A + (m*j) + z);
      }
      // Copia di subB
      for(int j = 0; j < mSub; j++){
         for(int z = 0; z < mSub; z++)
            *(subB + (mSub*j) + z) = *(B + (m*j) + z);
      }
   } else{
      // Ricezione delle proprie sottomatrici
      for(int j = 0; j < mSub; j++){
         MPI_Recv(subA + (mSub*j), mSub, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
         MPI_Recv(subB + (mSub*j), mSub, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
      }
   }
// ------------------------------- CALCOLO con BMR ------------------------------------------

   // SubB va inviato al processore sulla riga superiore e ricevuto dalla riga inferiore
   int upperRowRank, lowerRowRank;
   { // coords non è necessaro al di fuori di questo blocco
        int coords[2];
        coords[0] = (gridCoords[0] - 1); // Riga superiore
        coords[1] = gridCoords[1]; // Stessa colonna
        MPI_Cart_rank(comm_grid, coords, &upperRowRank);

        // SubB viene ricevuto dal processore sulla riga inferiore
        coords[0] = (gridCoords[0] + 1); // Riga inferiore, stessa colonna
        MPI_Cart_rank(comm_grid, coords, &lowerRowRank);
   }

   int bSenderCoords[2] = {gridCoords[0], gridCoords[0]}; // processore che effettua il broadcast sulla riga. Inizialmente sulla diagonale principale
   int bSenderRank; // Rank nel comunicator comm_grid
   MPI_Cart_rank(comm_row, bSenderCoords, &bSenderRank);


    if(gridCoords[0] == gridCoords[1]){ // Se si trova sulla diagonale si prepara ad inviare
      memcpy(Temp, subA, subMatrixSize * sizeof(float));
    }

    MPI_Bcast(Temp, subMatrixSize, MPI_FLOAT, bSenderRank, comm_row);

   t += matrXmatr(Temp, mSub, mSub, subB, mSub, subC);

    for(int i = 1; i < p; i++){
      bSenderCoords[1] += 1; // prossima diagonale
      MPI_Cart_rank(comm_row, bSenderCoords, &bSenderRank);

      if(gridCoords[1] == (bSenderCoords[1])%p){ // Se deve inviare
         memcpy(Temp, subA, mSub*mSub*sizeof(float));
      }

      MPI_Bcast(Temp, subMatrixSize, MPI_FLOAT, bSenderRank, comm_row); // Broadcast sulla stessa riga

        MPI_Isend(subB, subMatrixSize, MPI_FLOAT, upperRowRank, 0, comm_grid, &rqst); // Spedizione alla riga superiore

        MPI_Recv(subB, subMatrixSize, MPI_FLOAT, lowerRowRank, 0, comm_grid, &status); // Ricezione dalla riga inferiore


        t += matrXmatr(Temp, mSub, mSub, subB, mSub, subC); // prodotto
    }

// -------------------------------- FINE CALCOLO --------------------------------------

// -------------------------------- Comunicazione del risultato --------------------------------------

   if(worldRank==0){
       // Copia in C
       for (int i = 0; i < mSub; ++i) {
           memcpy(C + (i*m), subC + (i*mSub), mSub); // Copia ogni riga di subC in C
       }

       // Riceve dagli altri processori le loro parti di C
       int colOffset, rowOffset;
       for(int i = 1; i < worldSize; i++){
           rowOffset = (mSub) * (i/p); // Riga da cui partire
           colOffset = (mSub) * (i%p); // Colonna da cui partire
           for(int j = 0; j < mSub; j++){ // Ricezione per righe
               MPI_Recv(C + (rowOffset * m) + colOffset, mSub, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &status);
               rowOffset++; // passa alla prossima riga
           }
       }

   } else {
       for(int j = 0; j < mSub; j++){ // Invio per righe
           MPI_Send(subC + (mSub*j), mSub, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
       }
   }

   end = MPI_Wtime();

   double tMax; // Tempo massimo per il calcolo

   // Ricava il tempo massimo per il prodotto
   MPI_Reduce(&t, &tMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

   if (worldRank == 0) { // Stapa del risultato
       for (int i = 0; i < m; ++i) {
           for (int j = 0; j < m; ++j) {
               printf("%.2f\t", C[i*m + j]);
           }
           printf("\n");
       }

       printf("%f\t%f\n", tMax, end-start); // Tempo con e senza overhead
   }
   MPI_Finalize();
   return 0;
}

// -----------------------------------------------------------------------------------

// Alloca una matrice random di dimensione mxn
float* randomMatrix(int m, int n){
   float* M = malloc(m*n * sizeof(float));
   for(int i = 0; i < m*n; i++) {
       M[i] = ((float) rand())/((float) RAND_MAX/1000.0);
   }
   return M;
}

//crea una griglia bidimensionale periodica
void crea_griglia(int dim, MPI_Comm *comm_grid, MPI_Comm *comm_row){

   int dims[] = {dim, dim}; // Dimensione della griglia
   int period[] = {1, 1}; // La griglia deve essere periodica per ogni lato

   MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, 0, comm_grid); // Crrea una griglia di due dimensioni pxp periodica

   int remain[] = {0, 1};
   MPI_Cart_sub(*comm_grid, remain, comm_row); // Crea un comunicatore per broadcast su righe separate
}


// Prodotto righe per colonne
double matrXmatr(const float *A, int m, int n, const float *B, int k, float *C){
   double t_start, t_end;
    t_start = MPI_Wtime();

   for(int i = 0; i < m; i++){
       for(int z = 0; z < n; z++){
           for(int j = 0; j < k; j++){ // Prodotto scalare
               C[i*k+j] += A[i*n+z] * B[k*z+j];
           }
       }
   }

   t_end = MPI_Wtime();

   return t_end - t_start;
}
