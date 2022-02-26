#ifndef Header_h
#define Header_h

#define SUCCESS 1
#define FUNCTION_ERROR 2
#define EXTRA_FUNCTION_ERROR 3
#define SINGULAR_MATRIX 4
#define IMPOSSIBLE 5
#define ERROR 6

// in extra_functions
int Max (int i, int j);
int Matrix_Read(double *A, int n, char ** argv);
void Matrix_Print(double *M, int m, int n);
double Norm(double *M, int row, int col);
double choice(int k, int n, int i, int j);
double Swap(int i, int j, int n, double *M, double eps);
double Norm_Nev(double *A, double *x, double *b, int n, double *C_nev);
double Norm_Razn(double *x, int n, double *X);

long double get_full_time();
void * find_LU_x_y_threaded(void *pa);

// in functions
void synchronize(int total_threads);
void find_B(double* B, int n, double* A);
int find_LU_x_y(int n, double *A, double *B, double *x, double *y, double *L, double *U, int my_rank, int total_threads, int *flag);

#endif /* header_h */
