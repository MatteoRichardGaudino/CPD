#include <stdio.h>

int main(){
    double med = 0;
    double in;
    for(int i = 0; i < 5; i++){
        scanf("%lf", &in);
        med += in;
    }
    printf("%lf\n", med/5);
}