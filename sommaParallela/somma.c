#include <stdio.h>
#include <stdlib.h> // malloc()
#include <mpi.h> // MPI_Wtime()

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    unsigned int nums_c = 0;
    int* nums_v;
 
// ----------------------- Lettura input ----------------------

    FILE* test = fopen(argv[1], "r");
        
    fscanf(test, "%u" ,&nums_c);

    nums_v = malloc(sizeof(int)*nums_c);

    for (int i = 0; i < nums_c; i++){
        fscanf(test, "%d", &(nums_v[i]));
    }
    fclose(test);


    double start = MPI_Wtime();
// ----------------------- Somma ----------------------
    long long sum = 0;
    for (unsigned int i = 0; i < nums_c; i++){
        sum += nums_v[i];
    }
// ----------------------- Output ----------------------
    double end = MPI_Wtime();

    printf("%lld %.10f\n", sum, end-start);

    MPI_Finalize();
    return 0;
}