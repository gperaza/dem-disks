typedef struct disk
{
    int type; //0 mobile, 1 bottom, 2 wall, etc
    double x0; double y0; double w0; /* position x y and angular orientation*/
    double x1; double y1; double w1; /*higher derivatives*/
    double x2; double y2; double w2;
    double x3; double y3; double w3;
    double x4; double y4; double w4;
    double x5; double y5; double w5;
    double fx; double fy; double fw;
    double radius;
    double mass;
    double iMoment;
    int posFixed;
    int rotFixed;
    double xi; double yi; /*Initial positions of disks (useful for vibrating)*/
} disk;

typedef struct particleParameters
{
    double meanR;
    double kn; double kt; double pr; //elastic properties of the material
    double mu;
    double vGamma;
    double density;
} particleParameters;

typedef struct systemParameters
{
    long  nParticles;
    double box_w; double box_h;
    double gravity; double gravityAngle;
    double bGamma;
    double freq; double epsilon;
    double dimensionlessAc;
    long seed;
    double time;
    double timestep;
    double relaxTime; double runTime;
    double timeForGraph;
    double timeForWrite;
    double a[5];
    double c[6];
} systemParameters;

typedef struct neighbor_stats
{
    unsigned int touching;
    double s0;

} neighbor_stats;
