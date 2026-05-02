#include <stdio.h>
#include "linearalgebra.h"
#include <stdlib.h>

int main()
{
    matrix m1 = create_matrix(3,3), I3 = identity(3);
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j<3;++j){
            entry_in_matrix(&m1,i,j);
        }
    }
    matrix a = augment(m1, I3);
    matrix rref_a = rref(&a);
    print_matrix(a);
    printf("\n");
    print_matrix(rref_a);
    freemat(&m1);
    freemat(&I3);
    freemat(&a);
    freemat(&rref_a);
    return 0;
}