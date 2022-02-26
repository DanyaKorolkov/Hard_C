#ifndef Header_h
#define Header_h

#define SUCCESS 1
#define FUNCTION_ERROR 2
#define EXTRA_FUNCTION_ERROR 3
#define SINGULAR_MATRIX 4
#define IMPOSSIBLE 5
#define ERROR 6

// in extra_functions
int Check_sym(double *A, int n);
double Max (double i, double j);
int Sgn (double a);
int Matrix_Read(double *A, int n, char ** argv);
void Matrix_Print(double *M, int m, int n);
void Matrix_Filling(double *A, int n, int k);
double choice(int k, int n, int i, int j);
double Norm(double *M, int n);
double Norm_inf(double *M, int n);

//in functions
int Make_three_diag(int n, double *A, double eps);
int n_(double *A, int n, double c, double eps);
double Solution (int n, double *A, double epsilon, int number);
double Rec(double *A, double a, double b, double epsilon, int n, int number, double eps);

#endif
