#include "header.h"

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <cmath>

int Check_sym(double *A, int n){

    int i, j;
    
    for (i = 0; i < n ; i ++){
        
        for (j = i; j < n ; j++){
            
            if ( A[i*n + j] != A[j*n + i] ) return ERROR;
        }
    }
    return SUCCESS;
};

double Max(double i, double j){
    
    if (i>j) return i;
    else return j;
};

int Sgn (double a){
    
    if (a > 0)  return 1;
    if (a < 0)  return -1;
    if (a == 0) return 0;

    return SUCCESS;
};

void Matrix_Filling(double *A, int n, int k){
    
    int i, j;
    for (i = 0; i<n; i++){
        for (j = 0; j<n; j++) A[i*n + j] = choice(k, n, i, j);
    }
};

double choice(int k, int n, int i, int j){
    
    switch( k ){
            
        case (1):
            return (n - Max(i+1,j+1) + 1);
            break;
            
        case (2):
            if (i == j) {
                return (2);
                break;
            }
            
            if (abs(i-j) == 1){
                return (-1);
                break;
            }
            else{
                return 0;
                break;
            }
            
        case (3):
            if ( (j<n-1) && (i == j) ){
                return (1);
                break;
            }
            
            if (j == n-1){
                return (i);
                break;
            }
            
            if (i == n-1 ){
                return (j);
                break;
            }
            
            else {
                return (0);
                break;
            }
            
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
        return ERROR;
    };
    
    for (int i = 0; i<n*n; i++){
        
        if ( fscanf(in, "%lf", &A[i]) != 1 ){
            return ERROR;
        }
    }
    fclose(in);
    return SUCCESS;
};

void Matrix_Print(double *M, int m, int n){
    
    for(int i = 0; i< m; i++){
        
        for(int j = 0; j< n; j++) printf("%10.3e ", M[i*m + j]);
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

double Norm_inf(double *M, int n){

    double current_sum;
    double current_max = 0;

    for (int i = 0; i<n; i++){

        current_sum = 0;

        for(int j = 0; j<n; j++){

            current_sum += M[i*n + j];
        }
        if (current_sum > current_max) current_max = current_sum;
    }    
    return (current_max);
};
