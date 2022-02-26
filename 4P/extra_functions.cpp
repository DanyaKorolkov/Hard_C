#include "header.h"

#include<iostream>

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include <sys/time.h>
#include <sys/resource.h>

#include"mpi.h"

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
            return 0;
        }
};
double Norm(double * A, int n, int size){
	int i, j;
	int sum = 0;
	double norm_total = 0;
	double mem = 0;
	for (int i = 0; i < size; i++)
	{
		sum = 0;
		for (int j = 0; j < n; j++)
			sum += abs(A[i*n+j]);

		if (sum > mem)
			mem = sum;
	}


	MPI_Allreduce(&mem, &norm_total, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
return norm_total;
}

double Norm_Nev(double *A, double *x, double *b, int n, int my_rank, int p){

    int max_rows = n/p;    
    if (my_rank < n % p) max_rows++;
    
    int i,j;
    double cur;
    double sum = 0;
    double total_sum = 0;
    for (i = 0; i < max_rows; i++){
       
		cur = 0;
		
		for (j = 0; j < n; j++)
			cur += A[i*n + j] * x[j];

		
		cur -= b[i];

		sum += abs(cur);
    }
    MPI_Allreduce(&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    double norm_b;
    sum = 0;
    for (i = 0; i < max_rows; ++i)
    {
        sum += abs(b[i]);
    }
    MPI_Allreduce(&sum, &norm_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return total_sum / norm_b;
};

double Norm_Razn(double *x, int n, double *X, int size){
    
    int i;
	double norm_razn = 0;

    for (i = 0; i < n; i++){

        norm_razn += abs( (i+1)%2 - x[i] );
	}

    return norm_razn;
};






int read_matrix_A (double *A, double *Help, int n, char * argv, int my_rank, int p, double *B){
    
    int flag = 1;
    int i, j, dest, row;
    
    MPI_Status status;
    
    FILE *in = nullptr;
    
    if(my_rank == 0){
        in = fopen(argv, "rt");
        if (!in || in == NULL){
            
            fprintf(stderr, "Error opening data %c\n", *argv);
            flag = 0;
        }
    }
    
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (flag == 0)
        return 0;
    
    for (i = 0; i < n; i++){
        
        if (my_rank == 0){
            
            for(j = 0; j < n; j++)
                if(! fscanf(in, "%lf", &Help[j]) ){
                    flag = 0;
                    fclose(in);
                    break;
            }
        }
        
        MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (flag == 0)
            return 0;
        
        dest = i % p;
        row = (i - my_rank)/p;
        
        if (my_rank == 0 && dest == 0)
            for(j = 0; j < n; ++j)
                A[row*n + j] = Help[j];
        
        else if (my_rank == 0)
            MPI_Send(Help, n, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
        
        else if (dest == my_rank)
        {
            MPI_Recv(Help, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            for(j = 0; j < n; ++j)
                A[row*n + j] = Help[j];
        }
    }
    

	int size;
	size = n / p;
	if(my_rank < n % p) size ++;
 	for(i = 0; i < size; i ++)
		for (j = 0; j< (n+1)/2; j++) B[i]=A[i*n + 2*j];
	



    if (my_rank == 0) fclose(in);
    return 1;
}


int init_matrix_A(double *A, int n, int p, int k, int my_rank, double *B){
    
    int i, j;
    int counter = 0;
    for (i = 0; i < n; ++i)
    {
        int dest;
        dest = i % p;
        if (dest == my_rank)
        {
            for (j = 0; j < n; ++j)
            {
                A[counter*n + j] = choice(k, n, i, j);
            }
            counter++;
        }
    }



	int size;
	size = n / p;
	if(my_rank < n % p) size ++;
 	for(i = 0; i < size; i ++){
		B[i] = 0;
		for (j = 0; j< (n+1) / 2; j++) B[i]+=A[i*n + j*2];
	
	}

    return 1;
}

void print_matrix(int n, int q, int m, int start, double* A){
    
    int colmin = (q > m ? m : q);
    int rowmin = (n > m ? m : n);
    int i, j;
    for (i = start; i < colmin; i++)
    {
        for (j = 0; j < rowmin; j++)
        {
            printf("%10.3e ", A[i*n+j]);
        }
        printf("\n");
    }
}
