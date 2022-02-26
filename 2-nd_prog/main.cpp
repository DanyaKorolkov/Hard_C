#include "header.h"

#include<iostream>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>

using namespace std;

//6 номер
//пример ввода ./a 3 3 1 0 1 input1.txt

int main(int argc, char **argv){
    
    if (!(argc == 6 || argc == 7)){
        printf ("Неверное кол-во входных аргументов.\n");
        return 0;
    }
    
    int n, m, k, number;
    double epsilon;
    if (sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1 || sscanf(argv[3], "%lf", &epsilon) != 1 || sscanf(argv[4], "%d", &k) != 1 || (k < 0) == 1 || (k > 4) == 1|| (k==0 ^ argc==7) == 1 || (k!=0 ^ argc==6) == 1 || sscanf(argv[5], "%d", &number) != 1 || (number <= 0 ) == 1 || (number > n) == 1 )
    {
        printf("Неверный формат входных данных.\n");
        return 0;
    }
    
    double* A;
    
    try{

        A = new double[n*n];

        if (k == 0){
            if (Matrix_Read(A, n, &argv[6]) == ERROR){
                
                printf("Входная матрица некорректна.\n");
                delete []A;
                return 0;
            }
            if (Check_sym(A, n) == ERROR ){
                
                printf("Входная матрица не симметрична.\n");
                delete []A;
                return 0;
            }
        }
        else Matrix_Filling(A, n, k);
        
        printf("\nМатрица A:\n");
        Matrix_Print(A, m, n);
        printf("\n");

        clock_t time = clock();
        
        printf("%d-е собственное значение: %10.3e\n", number, Solution(n, A, epsilon, number));
        
        time = ( clock() - time ) / CLOCKS_PER_SEC;
        
        printf("\nВремя работы: %lu\n", time);
        
        delete []A;
    }
    
    catch(const exception& e)
     {
         std::cerr << e.what() << endl;
         delete []A;
         return ERROR ;
     }
}
