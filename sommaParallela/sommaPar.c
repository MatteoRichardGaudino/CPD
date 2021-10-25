#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    int worldSize, rank;
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    int commSteps = log2(worldSize); // per strategie 2 e 3
    unsigned long long pow2 = 1;

    int nums_c = 0; // dimensione del vettore
    int* nums_v; // vettore
    int strategy = 1; // strategia per la comunicazione

// ----------------------- Lettura input ----------------------
    if(argc < 3){
        if(rank == 0) printf("Error: usage mpiexec -np proc %s strategy filename\n", argv[0]);
        exit(-1);
    }

    if (rank == 0){
        strategy = atoi(argv[1]);

        FILE* test = fopen(argv[2], "r");
        
        fscanf(test, "%d" ,&nums_c);

        nums_v = malloc(sizeof(int)*nums_c);
        for (int i = 0; i < nums_c; i++){
            fscanf(test, "%d", &(nums_v[i]));
        }
        fclose(test);
    }

// ----------------------- Distribuzione input ----------------------

    MPI_Bcast(&nums_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&strategy, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int rest = nums_c % worldSize;
    int numToSum = nums_c/worldSize + ((rank < rest)? 1: 0);

    if(rank == 0){
        int tmp = nums_c/worldSize;
        int cursor = numToSum;
         
        for (int i = 1; i < worldSize; ++i){
            int inumToSum = tmp + ((i < rest)? 1: 0);
            MPI_Send(nums_v + cursor, inumToSum, MPI_INT, i, 0, MPI_COMM_WORLD);

            cursor += inumToSum;
        }
    } else{
        nums_v = malloc(sizeof(int)*numToSum);
        MPI_Recv(nums_v, numToSum, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }


    double t0, t1;
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();

    
// ----------------------- Calcolo ----------------------

    long long sum = 0;
    for (int i = 0; i < numToSum; i++){
        sum += nums_v[i];
    }

/* ----------------------- Soluzione 1 ---------------------- */
    if(strategy == 1){
        if (rank != 0){
            MPI_Send(&sum, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);
        } else{
            for (int i = 1; i < worldSize; i++){
                long long sum2;
                MPI_Status info;
                MPI_Recv(&sum2, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, &info);

                sum += sum2;
            }
        }
    } 
/* ----------------------- Soluzione 2 ---------------------- */
    else if(strategy == 2){ 
        for(int i = 0; i < commSteps; ++i){
            if((rank % pow2) == 0){
                if((rank % (pow2 << 1)) == 0){
                    long long sum2;
                    MPI_Recv(&sum2, 1, MPI_LONG_LONG, rank + pow2, 0, MPI_COMM_WORLD, &status);
                    sum += sum2;
                } else{
                    MPI_Send(&sum, 1, MPI_LONG_LONG, rank - pow2, 0, MPI_COMM_WORLD);
                }
                pow2 <<= 1;
            }
        }
    } 
/* ----------------------- Soluzione 3 ---------------------- */
    else if(strategy == 3){ 
        long long sum2;
        for(int i = 0; i < commSteps; ++i){
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
    

    MPI_Reduce(&localTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(strategy == 3){
        printf("[%d] %lld %.10f\n", rank, sum, localTime);
    }else if(rank == 0){
        printf("%lld %.10f\n", sum, maxTime);
    }

    MPI_Finalize();
    return 0;
}