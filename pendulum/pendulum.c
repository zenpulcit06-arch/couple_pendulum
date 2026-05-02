#include "pendulum.h"
#include "math.h"

void compute_acceleration(pendulum* p1, pendulum* p2, double k, double g)
{
    p1->alpha = -(g/p1->length)*p1->theta - (k*p1->sop*(p1->theta - p2->theta))/(p1->mass*p1->length);
    p2->alpha = -(g/p2->length)*p2->theta - (k*p2->sop*(p2->theta - p1->theta))/(p2->mass*p2->length);

}