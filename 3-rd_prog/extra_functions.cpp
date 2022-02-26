#include "header.h"

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include <sys/time.h>
#include <sys/resource.h>

long double get_full_time()
{
    struct timeval buf;
    gettimeofday(&buf, 0);
    return (double)buf.tv_sec + (double)buf.tv_usec/1000000.0;
}

int Max(int i, int j){
    
    if (i>j) return i;
    else return j;
}

double choice(int k, int n, int i, int j){
    
    switch( k ){
            
        case (1):
            return (n - Max(i+1,j+1) + 1);
            break;
        case (2):
            return (Max(i+1, j+1));
            break;
        case (3):
            return (abs( i - j ));
            break;
        case (4):
            return (1.0/(i+j +1));
            break;
        
        default:
            return ERROR;
        }
};


int Matrix_Read(double *A, int n, char ** argv){
    
    FILE *in;
    in = fopen(*argv, "rt");
    
    if (in == NULL) {
        return -1;
    };
    
    for (int i = 0; i<n*n; i++){
        if ( fscanf(in, "%lf", &A[i]) != 1) return -1;
    }
    fclose(in);
    return SUCCESS;
};

void Matrix_Print(double *M, int m, int n){
    
    for (int i = 0; i< m; i++){
        
        for(int j = 0; j< n; j++)
        {
            
            printf("%10.3e ", M[i*m + j]);
        }
    printf("\n");
    }
};

double Norm(double *M, int row, int col){
    
    double current_sum;
    double current_max = 0;
    
    for (int j = 0; j<col; j++){
        current_sum = 0;
        
        for(int i = 0; i<row; i++)
            current_sum += abs(M[i*row + j]);
            
        if (current_sum > current_max)
            current_max = current_sum;
    }
    return (current_max);
};

double Norm_Nev(double *A, double *x, double *b, int n, double *C_nev){
    
    int i;
    for (i = 0; i< n; i++){
        
        int j;
        for (j = 0; j < n; j++) C_nev[i] += A[i*n + j] * x[j];
        
        C_nev[i] -= b[i];
            
    }
    return Norm(C_nev, 1, n) / Norm(b, 1, n);
};

double Norm_Razn(double *x, int n, double *X){
    
    int i;
    for (i = 0; i<n; i++)
        X[i] = (i+1)%2 - x[i];

    return Norm(X, 1, n);
};










