#ifndef Header_h
#define Header_h

// in extra_functions
int Max (int i, int j);
double Norm(double * A, int n, int size);
double choice(int k, int n, int i, int j);
double Norm_Nev(double *A, double *x, double *b, int n, int my_rank, int p);
double Norm_Razn(double *x, int n, double *X, int size);

//new
int read_matrix_A(double *A, double *Help, int n, char * argv, int my_rank, int p, double *B);
int init_matrix_A(double *A, int n, int p, int k, int my_rank, double *B);

void print_matrix(int n, int q, int m, int start, double* A);
void print_vector(double * vector, int first_row, int last_row, int p);

// in functions
double find_B(int q, int n, double* A);
double Solution(int my_rank, int p, int n, double *L, double * U, double * A, double * B, int flag, double * Help, double * SuperHelp, double *x, double *y);

#endif /* header_h */
