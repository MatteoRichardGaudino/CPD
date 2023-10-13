#include <stdio.h> 
#include <stdlib.h>
#include <time.h>


int main(int argc, char** argv){
    FILE* testFile = fopen(argv[1], "w");
    uint fileSize = atoi(argv[2]);
    unsigned long long max = 1UL << 32;

    srand(time(NULL));

    fprintf(testFile, "%d\n", fileSize);

    for (uint i = 0; i < fileSize; i++){
        fprintf(testFile, "%llu\n", rand()%max);
    }
    
    fclose(testFile);
    return 0;
}