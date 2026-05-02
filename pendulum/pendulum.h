typedef struct pendulum{
    double theta;
    double omega;
    double alpha;
    double length;
    double mass;
    double sop; //spring_contact_point
} pendulum;

void compute_acceleration(pendulum* p1, pendulum* p2, double k, double g);