#include "verlet.h"

void inte_first_half(matrix* x , matrix* v , matrix* a ,double dt)
{
    matrix arr1[] = {*x, scalarmult(*v,dt), scalarmult(*a, 0.5*dt*dt)};
    matrix temp1 = add(arr1,3);
    freemat(x);
    freemat(&arr1[1]);
    freemat(&arr1[2]);
    *x = temp1;

    matrix arr2[] =  {*v, scalarmult(*a, 0.5*dt)};
    matrix temp2 = add(arr2,2);
    freemat(v);
    freemat(&arr2[1]);
    *v = temp2;
}

void inte_second_half(matrix*v , matrix* a, double dt)
{
    matrix arr1[] = {*v, scalarmult(*a,0.5*dt)};
    matrix temp = add(arr1,2);
    freemat(v);
    freemat(&arr1[1]);
    *v = temp;
}