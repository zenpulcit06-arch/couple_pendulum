#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "raylib.h"
#include "integrator/verlet.h"
#include "linearalgebra/linearalgebra.h"
#include "pendulum/pendulum.h"

#define g 9.81

int main()
{
    pendulum p1,p2;
    double k, dt, total_time;
    printf("Enter mass of pendulums m1 and m2: ");
    scanf("%lf %lf", &p1.mass, &p2.mass);

    printf("Enter length of pendulum l1 and l2: ");
    scanf("%lf %lf", &p1.length, &p2.length);

    printf("Enter value of spring constant: ");
    scanf("%lf",&k);

    printf("Enter time step: ");
    scanf("%lf",&dt);

    printf("Enter Total time: ");
    scanf("%lf",&total_time);

    printf("Enter intial theta of pendulum in si theta1 and theta2: ");
    scanf("%lf %lf", &p1.theta, &p2.theta);

    printf("Enter intial omega of pendulum in si omega1 and omega2: ");
    scanf("%lf %lf", &p1.omega, &p2.omega);

    printf("Enter spring contact point of pendulum sop1 and sop2: ");
    scanf("%lf %lf", &p1.sop,&p2.sop);

    matrix x = create_matrix(1,2), v = create_matrix(1,2), a = create_matrix(1,2);

    x.data[0][0] = p1.theta;
    x.data[0][1] = p2.theta;

    v.data[0][0] = p1.omega;
    v.data[0][1] = p2.omega;

    compute_acceleration(&p1,&p2,k,g);

    a.data[0][0] =  p1.alpha;
    a.data[0][1] =  p2.alpha;

    long long steps = (long long) (total_time/dt);

    int scale = 200;  // pixels per meter
    int pivot1_x = 300, pivot1_y = 100;
    int pivot2_x = 500, pivot2_y = 100;

    InitWindow(800, 600, "Coupled Pendulum");

    double KE0 = 0.5*p1.mass*pow(p1.omega*p1.length, 2) + 
             0.5*p2.mass*pow(p2.omega*p2.length, 2);
    double PE0 = 0.5*k*pow(p1.sop*(p1.theta - p2.theta), 2);
    double E_total = KE0 + PE0;
    double PE_max = E_total;  // maximum PE can't exceed total energy
    if(PE_max == 0) PE_max = 1.0;  // avoid division by zero


    for (long long i = 0; i < steps && !WindowShouldClose(); i++)
    {
        BeginDrawing();
        ClearBackground(BLACK);
        
        inte_first_half(&x,&v,&a,dt);

        p1.theta = x.data[0][0];
        p2.theta = x.data[0][1];

        compute_acceleration(&p1,&p2,k,g);

        a.data[0][0] =  p1.alpha;
        a.data[0][1] =  p2.alpha;

        inte_second_half(&v,&a,dt);

        p1.omega = v.data[0][0];
        p2.omega = v.data[0][1];

        double PE = 0.5*k*pow(p1.sop*(p1.theta - p2.theta), 2);
        double normalized = PE / PE_max;  // 0 to 1
        float hue = 240.0f - normalized * 240.0f;  // 240 (blue) to 0 (red)
        Color mycolor = ColorFromHSV(hue, 1.0f, 1.0f);

        int bob1_x = pivot1_x + scale * sin(p1.theta);
        int bob1_y = pivot1_y + scale * cos(p1.theta);

        int bob2_x = pivot2_x + scale * sin(p2.theta);
        int bob2_y = pivot2_y + scale * cos(p2.theta);

        DrawLine(pivot1_x, pivot1_y, bob1_x, bob1_y, WHITE);
        DrawCircle(bob1_x, bob1_y, 10, RED);

        DrawLine(pivot2_x, pivot2_y, bob2_x, bob2_y, WHITE);
        DrawCircle(bob2_x, bob2_y, 10, BLUE);

        DrawLine(bob1_x, bob1_y, bob2_x, bob2_y, mycolor); 

        EndDrawing();


    }
    
    CloseWindow();
    freemat(&x);
    freemat(&v);
    freemat(&a);

    
    return 0;
}
