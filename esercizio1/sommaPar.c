/**
 *  Autore: Matteo Richard Gaudino
 *  Matricola: N86003226
 * **/

#include <mpi.h> // Recv, Send, ...
#include <stdio.h> // std in/out/err
#include <stdlib.h> // malloc, atoi
#include <time.h> // Per la generazione di numeri random
#include <math.h> // log2 

/**
 * Utilizzo:
 *  mpiexec -machinefile <hostFile> -np <p> sommaPar <...>
 * 
 * Sistassi
 *  sommaPar strategy numToSum                          (numToSum > 20)
 *  sommaPar strategy numToSum num_1 ... num_n          (numToSum <= 20)
 * 
 * Semantica
 *  strategy stabilisce la strategia da utilizzare 1, 2 o 3.
 *     Nel caso il numero di processori non sia una potenza di 2 verrà utilizzata la strategia 1.
 *     Il valore di default è 1.
 *
 *
 *  Il programma prenderà in input <numToSum> numeri da sommare da riga di comando, se numToSum <=20.
 *  In alternativa, i numeri saranno generati casualmente.
 *     I numeri devono essere separati da uno spazio
 *
 * **/


int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    int worldSize, rank; // dimensione di MPI_COMM_WORLD e rank nel communicator
    MPI_Status status; // Per MPI_Recv

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    int commSteps = log2(worldSize); // per strategie 2 e 3
    unsigned long long pow2 = 1; // potenze di 2 per la strategia 2 e 3
// ----------------------- Lettura input ----------------------

    int nums_c = 0; // dimensione del vettore
    int* nums_v; // vettore da sommare

    int strategy = 1; // strategia per la comunicazione

    if(rank == 0){
        int opindex = 0; // indice per scorrere le opzioni
        if(argc < 3){
            perror("Too few options. Read documentation.\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        strategy = atoi(argv[1]);
        nums_c = atoi(argv[2]);

        if(strategy!=1 && strategy!=2 && strategy!=3){
            perror("Warning, strategy is a number other than 1, 2 or 3. Strategy is set to 1\n");
            strategy = 1;
        }

        nums_v = malloc(sizeof(int)*nums_c);

        if(nums_c <= 20){       // Prendi numeri da linea di comando
            if(argc - 3 != nums_c){
                perror("Too many options. Read documentation.\n");
                free(nums_v);
                MPI_Abort(MPI_COMM_WORLD, -1);
            }

            for(int i=0; i < nums_c; i++){
                nums_v[i] = atoi(argv[i+3]);
            }
        } else {                  // Genera numeri casuali
            srand(time(NULL)); // Genera numeri casuali
            unsigned long long max = 1UL << 32; // di massimo 32 bit

            for (int i = 0; i < nums_c; i++){
                nums_v[i] = rand()%max;
            }
        }

        if((worldSize & (worldSize - 1)) != 0 && strategy != 1) { // Verifica se il numero dei processori è una potenza di 2
            strategy = 1; // se non lo è si setta la strategia a 1
            perror("Warning strategy is not a power of 2, strategy is set to 1\n");
        }
    }

// ----------------------- Distribuzione input ----------------------

    MPI_Bcast(&nums_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&strategy, 1, MPI_INT, 0, MPI_COMM_WORLD);


    int rest = nums_c % worldSize;
    int numToSum = nums_c/worldSize + ((rank < rest)? 1: 0); // numeri da sommare per ogni processo

    if(rank == 0){
        int tmp = nums_c/worldSize; // Minimo numeri da sommare
        int cursor = numToSum;

        for (int i = 1; i < worldSize; ++i){
            int inumToSum = tmp + ((i < rest)? 1: 0); // aggiunge 1 se i < rest
            MPI_Send(nums_v + cursor, inumToSum, MPI_INT, i, 0, MPI_COMM_WORLD);

            cursor += inumToSum;
        }
    } else{
        nums_v = malloc(sizeof(int)*numToSum); // alloca spazio e riceve da p0
        MPI_Recv(nums_v, numToSum, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }

// ----------------------- Creazione dei timer ----------------------

    double t0, t1; // Per misurare i tempi
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime(); // è partita la misurazione

// ----------------------- Calcolo ----------------------

    long long sum = 0;
    for (int i = 0; i < numToSum; i++){
        sum += nums_v[i];
    }

/* ----------------------- Soluzione 1 ---------------------- */

    if(strategy == 1){
        if (rank != 0){ // se non è p0 invia
            MPI_Send(&sum, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
        } else{ // altrimenti ricevi e somma
            for (int i = 1; i < worldSize; i++){
                long long sum2;
                MPI_Recv(&sum2, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, &status);

                sum += sum2;
            }
        }
    } 
/* ----------------------- Soluzione 2 ---------------------- */
    else if(strategy == 2){ 
        for(int i = 0; i < commSteps; ++i){ // ripeti log2(p) volte
            if((rank % pow2) == 0){ // se partecipa alla comunicazione
                if((rank % (pow2 << 1)) == 0){ // se è un ricevente
                    long long sum2;
                    MPI_Recv(&sum2, 1, MPI_LONG_LONG, rank + pow2, 0, MPI_COMM_WORLD, &status);
                    sum += sum2;
                } else{ // altrimenti se è un mittente
                    MPI_Send(&sum, 1, MPI_LONG_LONG, rank - pow2, 0, MPI_COMM_WORLD);
                }
                pow2 <<= 1; // pow2 = pow2 * 2
            }
        }
    } 
/* ----------------------- Soluzione 3 ---------------------- */
    else if(strategy == 3){ 
        long long sum2;
        for(int i = 0; i < commSteps; ++i){ // come la strategia II, ma partecipano tutti
            if((rank % (pow2 << 1)) < pow2){
                MPI_Send(&sum, 1, MPI_LONG_LONG, rank + pow2, 0, MPI_COMM_WORLD);
                MPI_Recv(&sum2, 1, MPI_LONG_LONG, rank + pow2, 0, MPI_COMM_WORLD, &status);
            } else{
                MPI_Recv(&sum2, 1, MPI_LONG_LONG, rank - pow2, 0, MPI_COMM_WORLD, &status);
                MPI_Send(&sum, 1, MPI_LONG_LONG, rank - pow2, 0, MPI_COMM_WORLD);
                
            }
            sum += sum2;
            pow2 <<= 1;
        }
    }
    
/* ----------------------- Comunicazione dei risultati ---------------------- */
    t1 = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);

    double localTime = t1-t0;
    double maxTime;
    
    MPI_Reduce(&localTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); // calcola il tempo massimo impiegato

    if(rank == 0)
        printf("[%d] %lld %.10f\n", rank, sum, maxTime);

    free(nums_v);
    MPI_Finalize();
    return 0;
}