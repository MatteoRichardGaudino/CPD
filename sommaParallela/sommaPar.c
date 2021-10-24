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

    int commSteps = log2(worldSize);
    unsigned long long pow2 = 1;

    int nums_c = 0;
    int* nums_v;

    int strategy = 1;

// ----------------------- Lettura input ----------------------
    if (rank == 0){
        strategy = atoi(argv[1]);
        nums_c = atoi(argv[2]);

        nums_v = malloc(sizeof(int)*nums_c);

        for (int i = 0; i < nums_c; i++){
            nums_v[i] = 1;
        }
    }
    

    /*if (rank == 0){
        FILE* test = fopen("test.txt", "r");
        
        fscanf(test, "%d" ,&nums_c);

        nums_v = malloc(sizeof(int)*nums_c);
        for (int i = 0; i < nums_c; i++){
            fscanf(test, "%d", &(nums_v[i]));
        }
    }*/

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

    int sum = 0;
    for (int i = 0; i < numToSum; i++){
        sum += nums_v[i];
    }

/* ----------------------- Soluzione 1 ---------------------- */
    if(strategy == 1){
        if (rank != 0){
            MPI_Send(&sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        } else{
            for (int i = 1; i < worldSize; i++){
                int sum2;
                MPI_Status info;
                MPI_Recv(&sum2, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &info);

                sum += sum2;
            }
        }
    } 
/* ----------------------- Soluzione 2 ---------------------- */
    else if(strategy == 2){ 
        for(int i = 0; i < commSteps; ++i){
            if((rank % pow2) == 0){
                if((rank % (pow2 << 1)) == 0){
                    int sum2;
                    MPI_Recv(&sum2, 1, MPI_INT, rank + pow2, 0, MPI_COMM_WORLD, &status);
                    sum += sum2;
                } else{
                    MPI_Send(&sum, 1, MPI_INT, rank - pow2, 0, MPI_COMM_WORLD);
                }
                pow2 <<= 1;
            }
        }
    } 
/* ----------------------- Soluzione 3 ---------------------- */
    else if(strategy == 3){ 
        int sum2;
        for(int i = 0; i < commSteps; ++i){
            if((rank % (pow2 << 1)) < pow2){
                MPI_Send(&sum, 1, MPI_INT, rank + pow2, 0, MPI_COMM_WORLD);
                MPI_Recv(&sum2, 1, MPI_INT, rank + pow2, 0, MPI_COMM_WORLD, &status);
            } else{
                MPI_Recv(&sum2, 1, MPI_INT, rank - pow2, 0, MPI_COMM_WORLD, &status);
                MPI_Send(&sum, 1, MPI_INT, rank - pow2, 0, MPI_COMM_WORLD);
                
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

    if(rank == 0){
        printf("%d =? %d %s %.10f\n", nums_c, sum, (sum == nums_c)? "OK": "ERROR", maxTime);
    }

    MPI_Finalize();
    return 0;
}