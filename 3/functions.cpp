#include "header.h"
#include <iostream>
#include <cstdlib>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include <pthread.h>

void synchronize(int total_threads){
    
    static pthread_mutex_t mutex        = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t condvar_in    = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t condvar_out   = PTHREAD_COND_INITIALIZER;
    static int threads_in = 0;
    static int threads_out = 0;
    
    pthread_mutex_lock(&mutex);
    
    threads_in++;
    
    if (threads_in >= total_threads){
        
        threads_out = 0;
        pthread_cond_broadcast(&condvar_in);
    }
    else
        while(threads_in < total_threads) pthread_cond_wait(&condvar_in, &mutex);
    
    threads_out++;
    if (threads_out >= total_threads)
        {
            threads_in = 0;
            pthread_cond_broadcast(&condvar_out);
        }
        else
            while (threads_out < total_threads)
                pthread_cond_wait(&condvar_out, &mutex);
        
        pthread_mutex_unlock(&mutex);
};

void find_B(double* B, int n, double* A){
    
    for (int i = 0; i<n; i++)
        
        for(int j = 0; j<(n+1)/2; j++)
            
            B[i] += A[i*n + 2*j];
};

int find_LU_x_y(int n, double *A, double *B, double *x, double *y, double *L, double *U, int thread_num, int total_threads, int *flag)
{
    
    double eps = (1e-15)*Norm(A, n, n);

    
    double sum;
    
    int i;
    int counter;
    int j;
    int k;
    
    int begin_col;
    int last_col;
    int begin_row;
    int last_row;
    
    // находим L
    for(counter = 0; counter < n; counter++){

        if (*flag == 1) return IMPOSSIBLE;
        
        synchronize(total_threads);
        
        k = counter;
        
        begin_col   = (n-counter)*thread_num/total_threads + counter;
    	last_col    = (n-counter)*(thread_num + 1)/total_threads + counter;
        
        for(i = begin_col; i<last_col; i++)
        {
            
            sum = 0;
            
            for (j = 0; j<=k-1; j++)
                sum+=L[i*n + j]*U[j*n + k];

            L[i*n + k] = A[i*n + k] - sum;
            
            if ( (i == k) && (abs(L[i*n + k]) < eps) )
                *flag = 1;
            
        }
        
        synchronize(total_threads);
        
        if (*flag == 1) return IMPOSSIBLE;
        
        synchronize(total_threads);

        // находим U
        begin_row   = (n-counter)*thread_num/total_threads + counter;
    	last_row    = (n-counter)*(thread_num + 1)/total_threads + counter;
        
        int i = counter;
        for(k = begin_row; k < last_row; k++){
            
            sum = 0;
            int j;
            
            for(j = 0; j<= i-1; j++)
                sum += L[i*n + j]*U[j*n + k];
            
            U[i*n + k] = (A[i*n + k] - sum) / L[i*n + i];
            
            if ( (i == k) && (U[i*n + k] != 1 ) )
                *flag = 1;
        }
        synchronize(total_threads);
    }
    
    if (*flag == 1)
        return IMPOSSIBLE;

    
    // решаем Ly = B
    if (thread_num == 0){
        
        for (int i = 0; i < n; i++){
            
            sum = 0;
            
            for (int j = 0; j<i; j++)
                sum+= (L[i*n + j])*(y[j]);
            
            y[i] = (B[i] - sum)/(L[i*n + i]);

        }

        // решаем Ux = y
        for (int i = 1; i <= n; i++){
            
            sum = 0;
            
            for (int j = 1; j<i; j++)
                sum+= (U[(n-i)*n + n-j])*(x[n-j]);
            
            x[n-i] = (y[n-i] - sum)/(U[(n-i)*n + n-i]);
        }
    }

    synchronize(total_threads);
    
    return SUCCESS;
}
