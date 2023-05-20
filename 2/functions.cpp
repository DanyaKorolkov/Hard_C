#include "header.h"

#include<cstdlib>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<cmath>

int Make_three_diag(int n, double *A, double eps){
    
    if (n < 2) return SUCCESS;
    
    double b_ij, b_ii, b_jj;
    int i, j, q;
    
    for ( i=1; i< n-1; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            if (abs(A[j*n + i-1]) >= eps)
            {
                
                b_ii= sqrt(A[i*n+i-1]*A[i*n+i-1] +A[j*n + i-1]*A[j*n + i-1]);
                
                if (b_ii>=eps)
                {
                    
                    double c= A[i*n+i-1]/b_ii;
                    double s = -A[j*n+i-1]/b_ii;
                    
                    A[i*n+i-1] = b_ii;
                    A[(i-1)*n+i]= b_ii;
                    A[j*n+i-1]= 0;
                    A[(i-1)*n+j]= 0;
                    
                    for (q = i+1; q<n; q++)
                    {
                        if (q!=j)
                        {
                            
                            b_ii= A[i*n+q]*c-A[j*n+q]*s;
                            b_ij= A[i*n+q]*s+A[j*n+q]*c;
                            A[q*n+i]= b_ii;
                            A[q*n+j]= b_ij;
                            A[i*n+q]= b_ii;
                            A[j*n+q]= b_ij;
                        }
                    }
                            
                    b_ii=(A[i*n+ i]*c - A[j*n+ i]*s)*c - (A[i*n+ j]*c - A[j*n+ j]*s)*s;
                    b_ij=(A[i*n+ i]*c - A[j*n+ i]*s)*s + (A[i*n+ j]*c - A[j*n+ j]*s)*c;
                    b_jj=(A[i*n+ i]*s + A[j*n+ i]*c)*s + (A[i*n+ j]*s + A[j*n+ j]*c)*c;
                    
                    A[i*n+ i]= b_ii;
                    A[i*n+ j]= b_ij;
                    A[j*n+ j]= b_jj;
                    A[j*n+ i]= b_ij;
                    
                    
                }
            }
        }
    }
    
    return SUCCESS;
};


double Solution (int n, double *A, double epsilon, int number){
    
    double eps = (1e-15)*Norm(A, n);
    
    double b;
    b = Norm_inf(A, n);
    if(b == 0){
        return 0;
    }
    
    return Rec(A, -b, b, epsilon, n, number, eps); // вмест a подставил -b
};

double Rec(double *A, double a, double b, double epsilon, int n, int number, double eps){
    
    if ( (b - a)  > epsilon ){
        
        double my_n = n_(A, n, (a+b) / 2, eps) ;
        
        if (my_n < number ) a = (a+b) / 2;
        else b = (a+b) / 2;
        
        return Rec(A, a, b, epsilon, n, number, eps);
    }
    return ( (a + b) / 2 );
};


int n_(double *A, int n, double c, double eps){
    
    int i;
    double a_max, b_max;
    a_max = 0;
    b_max = 0;
    for (i = 0; i < n; i++){
        
        if (abs(A[i*n + i] - c) > a_max) a_max = abs(A[i*n + i]);
        if (i!= n-1){
            if (abs(A[i*n + i + 1]) > b_max) b_max = abs(A[i*n + i + 1]);
        }
    }
    double alpha;
    alpha = 4*Max(a_max, b_max);
    
    int m;
    double x, y, a_help, b_help, gamma, u, v;

    x = (A[0] - c) / alpha;
    y = 1;
    
    if (Sgn(x) < 0) m = 1;
    else m = 0;
    
    int k;
    for(k = 1; k < n; k++){
        
        a_help = (A[k*n + k] - c) / alpha;
        b_help = A[k*n + k -  1] / alpha;
        
        if (abs(b_help) < eps) b_help = 0;
        
        gamma = (1/eps) /( Max( abs(x), abs( b_help*(b_help*y) ) ) );
        
        u = gamma * (a_help*x - b_help*b_help*y );
        v = gamma*x;
        
        if (Sgn(u*x) < 0) m++;
        
        x = u;
        y = v;
    }
    
    return m;
};


//for (i=0; i<n-1; ++i)
//    {
//    for (j=i+2; j<n; ++j)
//        {
//            if (abs(A[i*n + j]) > eps){
//
//                double y = A[i*n + i] - A[j*n + j] ;
//                double x = -2*A[i*n + j];
//
//                printf("x = %lf,y = %lf\n", x, y);
//
//                if (abs(y) < eps){
//                    c = 1 / sqrt(2);
//                    s = c;
//                }
//                else{
//
//                    t = x/y;
//                    if ( abs(t) > eps){
//
//                        c = sqrt(0.5 * (1 + abs(y) /sqrt(x*x + y*y)));
//                        s = Sgn(x*y)*abs(x)/(2*c*sqrt(x*x + y*y));
//
//                        printf("МЫ В IF c = %lf, s = %lf\n", c, s);
//                    }
//                    else{
//                        c = 1;
//                        s = 0;
//
//                        printf("МЫ В ELSE c = %lf, s = %lf\n", c, s);
//                    }
//                }
//
//                    b_ii = (c*c - s*s)*A[i*n+i];
//                    b_ij = s*c*(A[i*n + i] - A[j*n + j]) + (c*c - s*s)*A[i*n + j];
//                    b_jj = (c*c - s*s)*A[j*n+j];
//                    b_ji = b_ij;
//
//                    A[i*n +i] = b_ii;
//                    A[i*n +j] = b_ij;
//                    A[j*n +i] = b_ji;
//                    A[j*n +j] = b_jj;
//
//                }
//            }
//        }
