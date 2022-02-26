#include "header.h"

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<cmath>



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

double Norm(double *M, int n){
    
    double current_sum;
    double current_max = 0;
    
    for (int j = 0; j<n; j++){
        
        current_sum = 0;
        
        for(int i = 0; i<n; i++){
            
            current_sum += M[i*n + j];
        }
        if (current_sum > current_max) current_max = current_sum;
    }
    return (current_max);
};

double Norm_Nev(double *A, double *x, double *b, int n){
    
    double *C;
    C = new double[n];
    int i;
	double norma_upper;
	norma_upper = 0;
    for (i = 0; i< n; i++){
        
        int j;
        for (j = 0; j < n; j++) C[i] += A[i*n + j] * x[j];
        
        C[i] -= b[i];
	norma_upper += C[i];  
    }

	double nev;
	nev = norma_upper/Norm(b, n);
	
    
    return nev;
};

double Norm_Razn(double *x, int n, double *X){
 
    int i;

	double norma;
	norma = 0;

    for (i = 0; i<n; i++){

	X[i]= (i+1)%2 - x[i];
	norma+=X[i];
	}

	return norma;
};

//double Swap(int i, int j, int n, double *M, double eps){
//
//    int row;
//    for (row = i+1; row < n; row++){
//
//        if(M[row*n + j] > eps){
//
//            int t;
//            int column;
//            for (column = j; column < n; column++ ){
//
//                t = M[row*n + column];
//                M[row*n + column] = M[i*n + column];
//                M[i*j + column] = t;
//
//            }
//            return SUCCESS;
//        }
//    }
//    return SINGULAR_MATRIX;
//};
