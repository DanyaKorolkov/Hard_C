#include "header.h"

#include <cstdlib>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>

void find_B(double* B, int n, double* A){
    
    for (int i = 0; i<n; i++){
        
        for(int j = 0; j<(n+1)/2; j++) {
            
            B[i] += A[i*n + 2*j];
        }
    }
};

int find_LU_x_y (int n, double *A, double *B, double *x, double *L, double *U, double *y){

    double eps = (5e-16)*Norm(A, n);

    
    double sum;
    
    int i;
    int counter;
    int j;
    int k;
    
    // находим L
    for(counter = 0; counter < n; counter++){

        k = counter;
        for(i = counter; i<n; i++){
            
            sum = 0;
            
            for (j = 0; j<=k-1; j++) sum+=L[i*n + j]*U[j*n + k];

            L[i*n + k] = A[i*n + k] - sum;
            
            if ( (i == k) && (abs(L[i*n + k]) < eps) ){
                if (i == 0){
                    
                    return IMPOSSIBLE;
                }
                else{
                    
                    return FUNCTION_ERROR;
                }
            }
        }

        // находим U
        int i = counter;
        for(k = counter; k < n; k++){
            
            sum = 0;
            int j;
            
            for(j = 0; j<= i-1; j++) sum += L[i*n + j]*U[j*n + k];
            U[i*n + k] = (A[i*n + k] - sum) / L[i*n + i];
            
            if ( (i == k) && (U[i*n + k] != 1 )){

                return FUNCTION_ERROR;
            }
        }
    }
    

    // решаем Ly = B
    for (int i = 0; i < n; i++){
        
        sum = 0;
        
        for (int j = 0; j<i; j++) sum+= (L[i*n + j])*(y[j]);
        y[i] = (B[i] - sum)/(L[i*n + i]);

    }


    // решаем Ux = y
    for (int i = 1; i <= n; i++){
        
        sum = 0;
        
        for (int j = 1; j<i; j++) sum+= (U[(n-i)*n + n-j])*(x[n-j]);
        x[n-i] = (y[n-i] - sum)/(U[(n-i)*n + n-i]);
    }
    
    return SUCCESS;
};
