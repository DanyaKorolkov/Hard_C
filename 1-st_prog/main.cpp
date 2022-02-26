#include "header.h"

#include<iostream>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>

using namespace std;
int main1(int argc, char **argv);

int main(int argc, char **argv){
return main1(argc, argv);
}
int main1(int argc, char **argv) {
    if (!(argc == 4 || argc == 5)){
        printf ("Неверное кол-во входных аргументов.\n");
        return 0;
    }
    
    int n, m, k;
    if (sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1 || sscanf(argv[3], "%d", &k) != 1 || (k==0 ^ argc==5) == 1 || (k < 0) == 1|| (k > 4) == 1 || (k!=0 ^ argc==4) == 1)
    {
        printf("Неверный формат входных данных.\n");
        return 0;
    }
    
    try{
          
	double *A; 
        A = new double[n*n];
        
        if (k == 0){
            if (Matrix_Read(A, n, &argv[4]) == -1){
                printf("Входная матрица некорректна.\n");
		delete []A;
                return 0;
            }
        }
        else{
            for (int i = 0; i<n; i++){
                for (int j = 0; j<n; j++) A[i*n + j] = choice(k, n, i, j);
            }
        }
        printf("Матрица A:\n");
        Matrix_Print(A, m, n);
            
        double* B;
        B = new double[n*n];
        find_B (B, n, A);
        
    //    printf("\nМатрица B\n");
    //    Matrix_Print(B, 1, n);
    //    printf("\n");
        
        
        double* x;
        x = new double[n];

	double *L, *U, *y;
	L = new double[n*n];
	U = new double[n*n];
	y = new double[n*n];
	
        clock_t begin = clock();

	int res=find_LU_x_y(n, A, B, x, L, U, y);
	clock_t end = clock();

        if (res != SUCCESS ){
            if(res == IMPOSSIBLE){
                printf("Невозможно применить методо LU разложения, первый угловой минор = 0\n");
                delete []A;
                delete []B;
                delete []x;
		delete []L;
		delete []U;
		delete []y;
                return 0;
            }
            else{
                printf("Матрица вырождена.\n");
                delete []A;
                delete []B;
                delete []x;
		delete []L;
		delete []U;
		delete []y;
                return 0;
            }
        };
        
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        
        printf("\nВектор x:\n");
        Matrix_Print(x, 1, n);       

        double *X;
        X = new double[n];
        printf("Норма разности:\n");
        printf("%10.3e\n", Norm_Razn(x, n, X));
        
        double *C;
        C = new double[n];
        printf("Норма невязки:\n");
        printf("%10.3e\n", Norm_Nev(A, x, B, n));

	printf("Время работы: %lf\n", time_spent);
       
	delete []A;
	delete []B;
	delete []x;
	delete []L;
	delete []U;
	delete []y;
	delete []X;
	delete []C;
        
        return 0;
    }
    catch(const exception& e)
     {
         std::cerr << e.what() << endl;
         return ERROR ;
     }
}
    
