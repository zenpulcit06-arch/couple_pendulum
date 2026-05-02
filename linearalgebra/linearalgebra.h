#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H
typedef struct matrix{
    int rows;
    int cols;
    double** data;
} matrix;
matrix create_matrix(int rows, int cols);
void freemat(matrix* m);
matrix add(matrix* a, int n);
matrix scalarmult(matrix a, double scalar);
matrix mult(matrix* a, int n);
double determinant(matrix* m);
matrix transpose(matrix* m);
matrix cofactor(matrix* m);
matrix adj(matrix* m);
matrix inverse(matrix* m);
matrix rref(matrix* m);
matrix identity(int n);
matrix augment(matrix a, matrix b);
void print_matrix(matrix m);
matrix basis(matrix* m);
void entry_in_matrix(matrix* m, int row,int col);
#endif // LINEARALGEBRA_H