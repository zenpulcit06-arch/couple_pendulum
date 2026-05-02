#include <math.h>
#include "linearalgebra.h"
#include <stdlib.h>
#include <stdio.h>

matrix create_matrix(int rows, int cols){
    if(rows <= 0 || cols <= 0){
        fprintf(stderr, "Matrix dimensions must be positive\n");
        exit(EXIT_FAILURE);
    }
    matrix m;
    m.rows = rows;
    m.cols = cols;
    m.data = (double**) malloc(rows*sizeof(double*));
    if(m.data == NULL){
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < rows; i++){
        m.data[i] = (double*) malloc(cols*sizeof(double));
        if(m.data[i] == NULL){
            fprintf(stderr, "Memory allocation failed\n");
            // Free already allocated memory to prevent leaks
            for(int j = 0; j < i; j++){
                free(m.data[j]);
            }
            free(m.data);
            exit(EXIT_FAILURE);
        }
    }
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            m.data[i][j] = 0.0;
        }
    }
    return m;
}

void freemat(matrix* m){
    for(int i = 0; i < m->rows; ++i){
        free(m->data[i]);
    }
    free(m->data);
    m ->data = NULL;
}

matrix add(matrix* a, int n){
    if(n <= 0){
        fprintf(stderr, "Number of matrices must be greater than 0\n");
        exit(EXIT_FAILURE);
    }
    matrix result = create_matrix(a[0].rows, a[0].cols);

    for(int i = 0; i < n; ++i){
        if(a[i].rows != result.rows || a[i].cols != result.cols){
            fprintf(stderr, "All matrices must have the same dimensions\n");
            exit(EXIT_FAILURE);
        }
        for(int j = 0; j < a[i].rows; ++j){
            for(int k=0; k < a[i].cols; ++k){
                result.data[j][k] += a[i].data[j][k];
            }
        }
    }
    return result;
}

matrix scalarmult(matrix a, double scalar){
    matrix result = create_matrix(a.rows, a.cols);
    if (result.data == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < a.rows; ++i){
        for (int j =0; j<a.cols; ++j){
            result.data[i][j] = a.data[i][j] * scalar;
        }
    }
    return result;
}

matrix mult(matrix* a, int n){
    if(n <= 0){
        fprintf(stderr, "Number of matrices must be greater than 0\n");
        exit(EXIT_FAILURE);
    }
    if(a == NULL){
        fprintf(stderr, "Matrix array is NULL\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < n; i++){
        if(a[i].data == NULL || a[i].rows <= 0 || a[i].cols <= 0){
            fprintf(stderr, "Invalid matrix at index %d\n", i);
            exit(EXIT_FAILURE);
        }
    }
    matrix result = create_matrix(a[0].rows, a[0].cols);
    // Copy a[0] to result
    for(int i = 0; i < a[0].rows; ++i){
        for(int j = 0; j < a[0].cols; ++j){
            result.data[i][j] = a[0].data[i][j];
        }
    }
    if(n == 1){
        return result;
    }
    for(int i = 1; i < n; ++i){
        if(result.cols != a[i].rows){
            fprintf(stderr, "Incompatible dimensions for multiplication at matrix %d\n", i);
            freemat(&result);
            exit(EXIT_FAILURE);
        }
        matrix temp = create_matrix(result.rows, a[i].cols);
        for(int j = 0; j < temp.rows; ++j){
            for(int k = 0; k < temp.cols; ++k){
                for(int l = 0; l < result.cols; ++l){
                    temp.data[j][k] += result.data[j][l] * a[i].data[l][k];
                }
            }
        }
        freemat(&result);
        result = temp;
    }
    return result;
}

double determinant(matrix* m){
    if(m->rows != m->cols){
        fprintf(stderr, "Determinant is only defined for square matrices\n");
        exit(EXIT_FAILURE);
    
    }

    for(int i= 0; i < m->rows; ++i){
        if(m->data[i] == NULL){
            fprintf(stderr, "Matrix data is NULL at row %d\n", i);
            exit(EXIT_FAILURE);
        }      
    }

    for(int i = 0; i < m->rows; ++i){
        for(int j = 0; j < m->cols; ++j){
            if(isnan(m->data[i][j]) || isinf(m->data[i][j])){
                fprintf(stderr, "Matrix contains invalid number at (%d, %d)\n", i, j);
                exit(EXIT_FAILURE);
            }
        }
    }
    

	int sign = 1;
	double prev_pivot = 1.0;

	matrix result = create_matrix(m->rows,m->cols);
    int row = 0;

    for(int i = 0; i < m->rows; ++i){
        for(int j = 0; j < m->cols; ++j){
            result.data[i][j] = m->data[i][j];
        }
    }

    for (int i = 0; i < result.cols ; i++)
    {
        int max_row = -1;
        for(int k = row; k < result.rows; ++k){
            if(result.data[k][i] != 0){
                max_row = k;
                break;
            }

        }

        if(max_row == -1) continue;

        if(max_row != row){
    		double* temp = result.data[max_row];
    		result.data[max_row] = result.data[row];
    		result.data[row] = temp;
    		sign *= -1;
		}   

        double pivot = result.data[row][i];

        for (int j = row+1; j < result.rows; ++j)
        {
        	for (int l = 0; l < result.cols; ++l)
        		result.data[j][l] = (result.data[j][l] * pivot - result.data[j][i] * result.data[row][l]) / prev_pivot;
        }

        row++;
        prev_pivot = pivot;
    }

    double det =  sign*result.data[result.rows-1][result.cols-1];

    freemat(&result);

    return det;
}

matrix transpose(matrix* m){
    matrix result = create_matrix(m->cols, m->rows);
    if (result.data == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < m->cols; i++)
    {
        for (int j = 0; j < m->rows; j++)
        {
            result.data[i][j] = m->data[j][i];
        }
        
    }

    return result;
    
}
matrix cofactor(matrix* m){
    if(m->rows != m->cols){
        fprintf(stderr, "Cofactor is only defined for square matrices\n");
        exit(EXIT_FAILURE);
    }
    matrix result = create_matrix(m->rows, m->cols);

    if(result.data == NULL){
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for(int i= 0; i < m-> rows; ++i)
    {
        for(int j = 0; j < m->cols; ++j)
        {
            matrix sub = create_matrix(m->rows-1, m->cols-1);
            for(int k = 0; k < m->rows; ++k)
            {
                if(k==i) continue;
                for(int l = 0; l < m->cols; ++l)
                {
                    if(l ==j) continue;
                    sub.data[k < i ? k : k-1][l < j ? l : l-1] = m->data[k][l];

                }

            }

            double det = determinant(&sub);
            result.data[i][j] = ((i+j) % 2 == 0 ? 1 : -1) * det;
            freemat(&sub);
        }
    }
    return result;
}

matrix adj(matrix* m){
    matrix cofac = cofactor(m);
    matrix result = transpose(&cofac);
    freemat(&cofac);
    return result;
}

matrix inverse(matrix* m){
    double det = determinant(m);
    if(det == 0){
        fprintf(stderr, "Matrix is singular and cannot be inverted\n");
        exit(EXIT_FAILURE);
    }
    matrix adjoint = adj(m);
    matrix result = scalarmult(adjoint, 1.0/det);
    freemat(&adjoint);
    return result;
}

matrix rref(matrix* m){
    matrix result = create_matrix(m->rows,m->cols);
    int row = 0;

    for(int i = 0; i < m->rows; ++i){
        for(int j = 0; j < m->cols; ++j){
            result.data[i][j] = m->data[i][j];
        }
    }

    for (int i = 0; i < result.cols ; i++)
    {
        int max_row = -1;
        for(int k = row; k < result.rows; ++k){
            if(result.data[k][i] != 0){
                max_row = k;
                break;
            }

        }

        if(max_row == -1) continue;
        if(max_row != row){
    		double* temp = result.data[max_row];
    		result.data[max_row] = result.data[row];
    		result.data[row] = temp;
		} 

        double pivot = result.data[row][i];
        for(int j = 0; j < result.cols;++j)
        {
            result.data[row][j] /= pivot;
        }

        for (int j = 0; j < result.rows; j++)
        {
            if (j == row) continue;
            
            double factor = result.data[j][i];

            for(int l = 0; l < result.cols; ++l)
            {
                result.data[j][l] -= factor * result.data[row][l];
            }
        }

        row++;
        
    }

    return result;
    
}

matrix identity(int n){
    matrix result = create_matrix(n,n);
    if(result.data == NULL){
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < n; ++i){
        result.data[i][i] = 1.0;
    }

    return result;
}

matrix augment(matrix a, matrix b){
    if(a.rows != b.rows){
        fprintf(stderr, "Matrices must have the same number of rows to be augmented\n");
        exit(EXIT_FAILURE);
    }

    matrix result = create_matrix(a.rows, a.cols + b.cols);

    if(result.data == NULL){
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for(int i = 0 ; i < a.rows; ++i)
    {
        for(int j = 0; j < a.cols; ++j){
            result.data[i][j] = a.data[i][j];
        }

        for(int j = 0; j < b.cols; ++j)
        {
            result.data[i][a.cols + j] = b.data[i][j];
        }
    }

    return result;
}

void print_matrix(matrix m){
    for(int i = 0; i < m.rows; ++i){
        for(int j = 0; j < m.cols; ++j){
            printf("%lf ",m.data[i][j]);
        }
        printf("\n");
    }
}

matrix basis(matrix* m){
    matrix rref_m = rref(m);
        
    int* pvt_cols = (int*) malloc(m->rows * sizeof(int));
    if(pvt_cols == NULL){
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    int rank = 0;

    for(int i = 0; i < rref_m.rows;++i)
    {
        for(int j = 0; j < rref_m.cols; ++j)
        {
            if(rref_m.data[i][j] != 0){
                pvt_cols[rank++] = j;
                break;
            }
        }
    }
    freemat(&rref_m);

    matrix basis_mat = create_matrix(m->rows , rank);

    for(int i = 0; i < m->rows; ++i)
    {
        for(int j = 0; j < rank ; ++j)
        {
            basis_mat.data[i][j] = m->data[i][pvt_cols[j]];
        }
    }
    free(pvt_cols);
    return basis_mat;
}

void entry_in_matrix(matrix* m, int row,int col){
    if(row < 0 ||  col < 0|| m->rows <= row || m->cols <= col ){
        fprintf(stderr, "Row and column indices must be non-negative\n");
        exit(EXIT_FAILURE);
    }

    printf("Enter element at (%d,%d): ",row,col);
    scanf("%lf", &m->data[row][col]);
}