#include "header.h"

#include<iostream>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>

#include<pthread.h>
#include<thread>

using namespace std;

//стр 184
//номер 3

typedef struct _ARGS{

    double *A;
    double *B;
    double *L;
    double *U;
    double *x;
	double *y;
    int n;
    int thread_num;
    int total_threads;
    int *flag;
    long double thread_time;
} ARGS;

long int threads_total_time = 0;
pthread_mutex_t threads_total_time_mutex = PTHREAD_MUTEX_INITIALIZER;

void * find_LU_x_y_threaded(void *pa){
    
    ARGS *pargs = (ARGS*)pa;
    long double t;
    
    printf("\nThread %d started\n", pargs->thread_num + 1);
    
    synchronize(pargs->total_threads);
    
    t = get_full_time();
    
    find_LU_x_y(pargs->n, pargs->A, pargs->B, pargs->x, pargs->y, pargs ->L, pargs->U, pargs->thread_num, pargs ->total_threads, pargs->flag);
    
    t = get_full_time() - t;
    pargs -> thread_time = t;
    
    pthread_mutex_lock (&threads_total_time_mutex);
    threads_total_time += t;
    pthread_mutex_unlock (&threads_total_time_mutex);
    printf("\nThread %d finished, time = %Lf\n", pargs -> thread_num + 1, t);
    
    synchronize(pargs->total_threads);
    
    return NULL;
};



int main(int argc, char **argv) {
    
    if (!(argc == 5 || argc == 6)){
        printf ("Неверное кол-во входных аргументов.\n");
        return 0;
    }
    
    int n, m, k, nthreads;
    if (sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1 || sscanf(argv[3], "%d", &k) != 1 || sscanf(argv[4], "%d", &nthreads) != 1 || (k==0 ^ argc==6) == 1 || (k > 4) == 1 || (k!=0 ^ argc==5) == 1 || (nthreads <= 0) == 1)
    {
        printf("Неверный формат входных данных.\n");
        return 0;
    }
    
    double *A;
    double *L;
    double *U;
    double *x;
	double *y;
    pthread_t *threads;
    ARGS * args;
    long double t_full;
    int i;
    
    x = new double[n];
    
    try{
        
        args = new ARGS [nthreads];
        threads = new pthread_t [nthreads];
        
        A = new double[n*n];
        L = new double[n*n];
        U = new double[n*n];
		y = new double[n];
        
        if (k == 0){
            if (Matrix_Read(A, n, &argv[5]) == -1){
                printf("Входная матрица некорректна.\n");
                return 0;
            }
        }
        else{
            for (int i = 0; i<n; i++){
                for (int j = 0; j<n; j++) A[i*n + j] = choice(k, n, i, j);
            }
        }

        printf("Матрица A:\n");
        Matrix_Print(A, m, m);
        printf("\n");
        
        double* B;
        B = new double[n*n];
        find_B (B, n, A);
        

        int flag = 0;
            
        for (i = 0; i < nthreads; i++){
            args[i].A = A;
            args[i].B = B;
            args[i].L = L;
            args[i].U = U;
            args[i].x = x;
			args[i].y = y;
            args[i].n = n;
            args[i].thread_num = i;
            args[i].total_threads = nthreads;
            args[i].flag = &flag;
        }
        
        t_full = get_full_time();
        
        for(i = 0; i < nthreads; i++){
            if(pthread_create( threads + i, 0, find_LU_x_y_threaded, args + i) ){
                
            printf("Не удалось создать поток.\n");
            delete []A;
            delete []x;
			delete []y;
            delete []B;
            delete []L;
            delete []U;
            delete []args;
            delete []threads;
                
            return 10;
                
            }
        }
        
        for (i = 0; i < nthreads; i++){
            if (pthread_join( threads[i], 0)){
                
                printf("Не можем дождаться %d thread.\n", i);
                
                delete []A;
                delete []x;
				delete []y;
                delete []B;
                delete []L;
                delete []U;
                delete []args;
                delete []threads;
                
                return 10;
            }
        }
        
        if (flag == 1){
            
            printf("Невозможно применить LU разложение\n");
            
            delete []B;
            delete []x;
			delete []y;
            delete []A;
            delete []L;
            delete []U;
            delete []args;
            delete []threads;
            
            return 5;
        }
        
        t_full = get_full_time() - t_full;

        printf("\nВектор x:\n");
        Matrix_Print(x, 1, m);
        printf("\n");
        
        double* X;
        X = new double[n];
        printf("Норма разности: %10.3e\n",  Norm_Razn(x, n, X));
        delete []X;
        
        double *C_nev;
        C_nev = new double[n];
        printf("Норма невязки:%10.3e\n", Norm_Nev(A, x, B, n, C_nev));
        delete []C_nev;
        
        printf("\nОбщее время работы программы = %Lf\n", t_full);
        
        delete []B;
        delete []x;
		delete []y;
        delete []A;
        delete []L;
        delete []U;
        delete []args;
        delete []threads;
        
        return 0;
    }
    catch(const exception& e)
     {
         std::cerr << e.what() << endl;
         delete []A;
         delete []args;
         delete []threads;
         delete []x;
		 delete []y;
         delete []L;
         delete []U;
         return ERROR ;
     }
}




