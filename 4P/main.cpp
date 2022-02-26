#include "header.h"

#include<iostream>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include<unistd.h>

#include "mpi.h"

#define BUF_LEN 256

using namespace std;

//стр 232

int main(int argc, char **argv) {
    
    int n, m, k;
    
    double *A;
    double *B;
    double *L;
    double *U;
    double *x;
    double *y;
    double *X;
    int i, j;
    double time;
    
    double *Help;
    double *SuperHelp;
    
    int my_rank;
    int p;
    int tag =0;
    char message[BUF_LEN];
    int source;
    int dest;
    
    int max_rows;
    int flag = 1;
    
    int counter;
    int q;
    
    MPI_Status status;
    
    
    long double t_full;
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    if(my_rank != 0 ){
        
        sprintf(message, "Hello from process %d\n", my_rank);
        dest = 0;
        MPI_Send(message, strlen(message) + 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
    }
    else{
        for (source = 1; source < p; source++){
            MPI_Recv(message, BUF_LEN, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
            printf("%s\n", message);
        }
    }
    if ( my_rank == 0 ){
        
        if (!(argc == 4 || argc == 5)){
	    flag = 0;
            printf ("Неверное кол-во входных аргументов.\n");
        }
        
        if (sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1 || sscanf(argv[3], "%d", &k) != 1 || (k==0 ^ argc== 5) == 1 || (k!=0 ^ argc== 4) == 1 )
        {

	    flag = 0;
            printf("Неверный формат входных данных.\n");

        }
    }
    

	MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	
	if (flag == 0){
	
	printf("flag = 0\n");
        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return -1;
	}
	
    
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    max_rows = n / p + 1;
    
   
int mem_error = 0;
int er = 0;
try
	{
		x = new double[n];
	}
	catch (bad_alloc& e)
	{
		cout << e.what() << endl;
		mem_error = -1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &er, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	if (er<0)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}

	try
	{
		y = new double[n];
	}
	catch (bad_alloc& e)
	{
		cout << e.what() << endl;
		mem_error=-1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &er, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	if (er<0)
	{
		delete[] x;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}

	try
	{
		A = new double[n*max_rows];
	}
	catch (bad_alloc& e)
	{
		cout << e.what() << endl;
		mem_error=-1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &er, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	if (er<0)
	{
		delete[] x;
		delete[] y;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}

	try
	{
		B = new double[max_rows];
	}
	catch (bad_alloc& e)
	{
		cout << e.what() << endl;
		mem_error=-1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &er, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	if (er<0)
	{
		delete[] A;
		delete[] x;
		delete[] y;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}

	try
	{
		L = new double[max_rows *n];
	}
	catch (bad_alloc& e)
	{
		cout << e.what() << endl;
		mem_error=-1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &er, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	if (er<0)
	{
		delete[] A;
		delete[] x;
		delete[] y;
		delete[] B;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}

	try
	{
		U = new double[max_rows * n];
	}
	catch (bad_alloc& e)
	{
		cout << e.what() << endl;
		mem_error=-1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &er, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	if (er<0)
	{
		delete[] A;
		delete[] B;
		delete[] L;
		delete[] x;
		delete[] y;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}

	try
	{
		Help = new double[n];
	}
	catch (bad_alloc& e)
	{
		cout << e.what() << endl;
		mem_error=-1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &er, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	if (er<0)
	{
		delete[] A;
		delete[] L;
		delete[] B;
		delete[] x;
		delete[] y;
		delete[] U;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}

	try
	{
		SuperHelp = new double[n];
	}
	catch (bad_alloc& e)
	{
		cout << e.what() << endl;
		mem_error=-1;
		MPI_Barrier(MPI_COMM_WORLD);
		return -1;
	}
	MPI_Allreduce(&mem_error, &er, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	if (er<0)
	{
		delete[] A;
		delete[] L;
		delete[] B;
		delete[] x;
		delete[] y;
		delete[] U;
		delete[] Help;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return -1;
	}





        
   /*     if (! (x = new double[n]) ){
            fprintf(stderr, "Not enough memory!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
	
       if (! (y = new double[n]) ){
            fprintf(stderr, "Not enough memory!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (! (A = new double[max_rows * n])){
            fprintf(stderr, "Not enough memory!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        if (! (B = new double[max_rows])){
            fprintf(stderr, "Not enough memory!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        if (! (L = new double[max_rows * n])){
            fprintf(stderr, "Not enough memory!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        if (! (U = new double[max_rows * n])){
            fprintf(stderr, "Not enough memory!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (! (Help = new double[n])){
            fprintf(stderr, "Not enough memory!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }       
	if (! (SuperHelp = new double[n])){
            fprintf(stderr, "Not enough memory!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
      */  

        if (k == 0){
            



           if ( ( read_matrix_A(A, Help, n, argv[5], my_rank, p, B) ) == 0){
                
            if (my_rank == 0) printf("Error reading\n");
                
        	delete[]x;
			delete[]y;
	        delete[]A;
	        delete[]B;
	        delete[]L;
	        delete[]U;
			delete[]Help;
			delete[]SuperHelp;
                
		    MPI_Barrier(MPI_COMM_WORLD);
		    MPI_Finalize();
		    return -1;
            }
        }
        else{
            if (init_matrix_A (A, n, p, k, my_rank, B) != 1){
                
                if (my_rank == 0) printf("Error initialization\n");
                
              	delete [] x;
				delete [] y;
				delete [] A;
				delete [] B;
				delete [] L;
				delete [] U;
				delete [] Help;
				delete [] SuperHelp;
                
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Finalize();
                return -1;
            }
        }
        
        if (my_rank == 0) printf("Матрица A:\n");
        for(counter = 0; counter < m; counter++)
        {
            if(my_rank == counter % p)
                print_matrix(n, counter / p + 1, m, counter / p, A);
	

            MPI_Barrier(MPI_COMM_WORLD);
        }
        printf("\n");
	

	MPI_Barrier(MPI_COMM_WORLD);

	time = MPI_Wtime();

        if (Solution (my_rank, p, n, L, U, A, B, flag, Help, SuperHelp, x, y) != 1){
            if (my_rank == 0) printf("Невозможно применить метод\n");
            
            delete [] x;
			delete [] y;
			delete [] A;
			delete [] B;
			delete [] L;
			delete [] U;
			delete [] Help;
			delete [] SuperHelp;
            
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
            return 0;
        }

		time = MPI_Wtime() - time;

		if (my_rank == 0){
			printf("Вектор x =\n");
			for (i = 0; i < n; i ++) 
				printf("%10.3e ", x[i]);
		}

		 if (! (X = new double[n]) ){
				    fprintf(stderr, "Not enough memory!\n");
				    MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
		if (my_rank == 0)
			printf("\nНорма разности = %10.3e\n", Norm_Razn(x, n, X, max_rows));
	

		double norma_nev;
		norma_nev = Norm_Nev(A, x, B, n, my_rank, p);
	
		if(my_rank == 0)
		 	printf("Норма невязки = %10.3e\n", norma_nev);


		if(my_rank == 0) printf("Время работы = %.2f sec\n", time);
		


	delete [] x;
	delete [] y;
	delete [] A;
	delete [] B;
	delete [] X;
	delete [] L;
	delete [] U;
	delete [] Help;
	delete [] SuperHelp;

	MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 1;
}




