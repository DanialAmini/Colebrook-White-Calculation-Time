



#include <iostream>
#include <chrono>
#include <sys/types.h>
#include <time.h>
//#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>


int main()
{
    const int N = 10000002;  //10mill logs takes 0.018 sec? like goddam :)) this is good, this is incredibly quick
    int dim = N;
    double *x = (double*) malloc(sizeof(double) * N);
    double *y = (double*)malloc(sizeof(double) * N);
    int i, j;

    //malloc uses heap memory allocation because otherwise you would get a stack overflow error. 
    // very large arrays do not fit in the stack, hence the error. 
    // heap is probably not "continuous" or something, I don't know. 

    clock_t start, finish;
    srand(time(NULL));


    for (i = 0; i < N; i++) { x[i] = 0.1 + ((double)rand() / (RAND_MAX)); }


    printf("size of x vector, N=%d\n", N);

    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = log(x[i]); }
        finish = clock();
        printf("time for log(x) is %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");
    
    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = x[i]*y[i]; }
        finish = clock();
        printf("time for x*y is %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");

    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = x[i] + y[i]; }
        finish = clock();
        printf("time for x+y is %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");

    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = x[i] / y[i]; }
        finish = clock();
        printf("time for x/y is %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");

    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = exp(x[i]); }
        finish = clock();
        printf("time for exp(x) is %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");

    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = pow(x[i],y[i]); }
        finish = clock();
        printf("time for x^y is %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");


    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = sqrt(x[i]); }
        finish = clock();
        printf("time for sqrt(x) is %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");

    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = x[i]*x[i]; }
        finish = clock();
        printf("time for x^2=x*x is %lf s\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");

    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = pow(x[i],2); }
        finish = clock();
        printf("time for x^2=pow(x,2) is %lf s (it recognizes integer power)\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");

    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = pow(x[i], 1/2); }
        finish = clock();
        printf("time for sqrt(x)=pow(x,1/2) is %lf s (incredibly quick)\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");

    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = pow(x[i], 1.201); }
        finish = clock();
        printf("time for x^1.201=pow(x,1.201) is %lf s (relatively slow)\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");

    for (j = 0; j < 5; j++) {
        start = clock();
        for (i = 0; i < N; i++) { y[i] = pow(1.201,x[i]); }
        finish = clock();
        printf("time for 1.201^x=pow(1.201,x) is %lf s (relatively slow)\n", (double)(finish - start) / CLOCKS_PER_SEC);
    }
    printf("\n");


    printf("\n");

    std::cout.precision(17);

    for (i = 0; i < N; i++) { y[i] = log(x[i]); }
    for (i = N - 10; i < N; i++) { std::cout << i << "\t" << x[i] << "\t" << y[i] << "\n"; }


    std::cout << "Hello World!\n";
}




//list of times
//
//size of x vector, N = 10000002
//time for log(x) is 0.017000 s
//time for x* y is 0.014000 s
//time for x + y is 0.015000 s
//time for x / y is 0.016000 s
//time for exp(x) is 0.019000 s
//time for x^ y is 0.049000 s
//time for sqrt(x) is 0.013000 s
//time for x ^ 2 = x * x is 0.011000 s
//time for x ^ 2 = pow(x, 2) is 0.011000 s(it recognizes integer power)
//time for sqrt(x) = pow(x, 1 / 2) is 0.006000 s(incredibly quick)
//time for x ^ 1.201 = pow(x, 1.201) is 0.056000 s(relatively slow)
//time for 1.201 ^ x = pow(1.201, x) is 0.051000 s(relatively slow)
