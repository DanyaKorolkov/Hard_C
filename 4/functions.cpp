#include "header.h"

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include<unistd.h>

#include "mpi.h"

double Solution(int my_rank, int p, int n, double *L, double * U, double * A, double * B, int flag, double * Help, double * SuperHelp, double *x, double *y){
   
	
    int j, i, counter;
    double mem = 1;
    double sum;
    int max_rows = n/p;
    
    MPI_Status status;

    double eps;

	if (my_rank < n % p) max_rows++;

	MPI_Barrier(MPI_COMM_WORLD);

	eps = Norm(A, n, max_rows) * 1e-15;

	MPI_Barrier(MPI_COMM_WORLD);


    for (j = 0; j < n; j++)
    {
        if (j == 0)

        {

            if (my_rank == 0)
            {
                mem = A[0];
                if (abs(mem) < eps) flag = 0;
            }
            MPI_Bcast(&mem, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (flag == 0)
                return 0;

  

            if (my_rank == 0)
            {
                for (i = 0; i < n; ++i)
                    Help[i] = A[i];
            }

            MPI_Bcast(Help, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);


            for (i = 0; i < max_rows; i++)
            {	
                L[i * n] = A[i * n];
                U[i * n] = Help[i*p + my_rank] / mem;
            }

        }

        else
        {
            if (my_rank == j % p)
            {
                for (counter = 0; counter < n; counter++)
                {
                    Help[counter] = U[(j / p) * n + counter];
                }
            }

            MPI_Bcast(Help, n, MPI_DOUBLE, j % p, MPI_COMM_WORLD);

            for (i = 0; i < max_rows; i++)
            {
                if (p * i + my_rank >= j)
                { 
                    sum = 0;
                    for (counter = 0; counter < j; counter++)
                        sum += L[i * n + counter] * Help[counter];

                    L[i * n + j] = A[i * n + j] - sum;
                }
            }
            if (my_rank == j % p)
                mem = L[(j / p) * n + j];

            MPI_Bcast(&mem, 1, MPI_DOUBLE, j % p, MPI_COMM_WORLD);

            if (my_rank == j % p)
            {
                for (counter = 0; counter < n; counter++)
                    Help[counter] = L[(j / p) * n + counter];
            }
            MPI_Bcast(Help, n, MPI_DOUBLE, j % p, MPI_COMM_WORLD);

            if (my_rank == j % p)
            {
                for (counter = 0; counter < n; counter++)
                    SuperHelp[counter] = A[(j / p) * n + counter];
            }

            MPI_Bcast(SuperHelp, n, MPI_DOUBLE, j % p, MPI_COMM_WORLD);
            
            for (i = 0; i < max_rows; i++)
            {

                if (p * i + my_rank >= j)
                {
                    sum = 0;
                    for (counter = 0; counter < j; counter++){

                        sum += Help[counter] * U[i * n + counter];
						
					}

                    U[i * n + j] = (SuperHelp[i*p + my_rank] - sum) / mem;
                }
            }

        }


    }


MPI_Barrier(MPI_COMM_WORLD);


    // Ly = B - находим y
        
        
    for (i = 0; i < n; i ++)
    {   
        if (my_rank == i % p)
        {
                
            sum = 0;
                
            for (j = 0; j < i; j++)
            {
                sum += L[(i/p) * n + j] * y[j];
            }

		if(abs( L[(i/p) * n + i]) < eps) flag = 0;
              
		 
            y[i] = (B[i / p] - sum ) / L[(i/p) * n + i];
        }
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(y, n, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
    }

	MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (flag == 0) return 0;    


	MPI_Bcast(y, n, MPI_DOUBLE, (n - 1) % p, MPI_COMM_WORLD);

	//MPI_Barrier(MPI_COMM_WORLD);


    //Ux = y - находим x
    for (i = 0; i < n; i ++) Help[i] = 0;
        
    for(i = n - 1; i >= 0; i --)
    {
            
        //копируем "строку" U
        for (j = n - 1; j >= i; j--)
        {
            if(my_rank == j % p)
            {
                Help[j] = U[ (j / p) * n + i];
            }
        }
        MPI_Allreduce(Help, SuperHelp, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
/*if (my_rank == 0) printf("Матрица L:\n");
        for(counter = 0; counter < n; counter++)
        {
            if(my_rank == counter % p)
                print_matrix(n, counter / p + 1, n, counter / p, L);

			sleep(0.1);

            MPI_Barrier(MPI_COMM_WORLD);
        }
        printf("\n");
	sleep(0.5);

if (my_rank == 0) printf("Матрица U:\n");
        for(counter = 0; counter < n; counter++)
        {
            if(my_rank == counter % p)
                print_matrix(n, counter / p + 1, n, counter / p, U);

			sleep(0.1);

            MPI_Barrier(MPI_COMM_WORLD);
        }
        printf("\n");
	sleep(0.5);

    */        
        if (my_rank == i % p)
        {

            sum = 0;
            for (j = n-1; j > i; j--)
            	sum+=SuperHelp[j] * x[j];
                
            if(abs(SuperHelp[i]) < eps)
           	flag = 0;
       

            if (flag != 0)
                x[i] = (y[i] - sum) / SuperHelp[i];
        }
            
        MPI_Bcast(&flag, 1, MPI_INT, i % p, MPI_COMM_WORLD);

        if (flag == 0) return 0;
	
	MPI_Bcast(x, n, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
    }
    return 1;
}
