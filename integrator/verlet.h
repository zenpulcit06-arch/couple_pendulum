#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../linearalgebra/linearalgebra.h"

void inte_first_half(matrix* x , matrix* v , matrix* a ,double dt);

void inte_second_half(matrix*v , matrix* a, double dt);