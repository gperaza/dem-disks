#include "common.h"

void predictor_positions(disk *particle) {
    double *a = global.a;

    particle->x0 +=
        particle->x1 * a[0] +
        particle->x2 * a[1] +
        particle->x3 * a[2] +
        particle->x4 * a[3] +
        particle->x5 * a[4];
    particle->x1 +=
        particle->x2 * a[0] +
        particle->x3 * a[1] +
        particle->x4 * a[2] +
        particle->x5 * a[3];
    particle->x2 +=
        particle->x3 * a[0] +
        particle->x4 * a[1] +
        particle->x5 * a[2];
    particle->x3 +=
        particle->x4 * a[0] +
        particle->x5 * a[1];
    particle->x4 +=
        particle->x5 * a[0];

    particle->y0 +=
        particle->y1 * a[0] +
        particle->y2 * a[1] +
        particle->y3 * a[2] +
        particle->y4 * a[3] +
        particle->y5 * a[4];
    particle->y1 +=
        particle->y2 * a[0] +
        particle->y3 * a[1] +
        particle->y4 * a[2] +
        particle->y5 * a[3];
    particle->y2 +=
        particle->y3 * a[0] +
        particle->y4 * a[1] +
        particle->y5 * a[2];
    particle->y3 +=
        particle->y4 * a[0] +
        particle->y5 * a[1];
    particle->y4 +=
        particle->y5 * a[0];

    return;
}

void predictor_rotations(disk *particle){
    double *a = global.a;

    particle->w0 +=
        particle->w1 * a[0] +
        particle->w2 * a[1] +
        particle->w3 * a[2] +
        particle->w4 * a[3] +
        particle->w5 * a[4];
    particle->w1 +=
        particle->w2 * a[0] +
        particle->w3 * a[1] +
        particle->w4 * a[2] +
        particle->w5 * a[3];
    particle->w2 +=
        particle->w3 * a[0] +
        particle->w4 * a[1] +
        particle->w5 * a[2];
    particle->w3 +=
        particle->w4 * a[0] +
        particle->w5 * a[1];
    particle->w4 +=
        particle->w5 * a[0];

    return;
}

void corrector_positions(disk *particle) {
    double *c = global.c;
    double correctorX, correctorY;

    correctorX = (particle->fx/particle->mass - particle->x2);
    particle->x0 += (c[0] * correctorX);
    particle->x1 += (c[1] * correctorX);
    particle->x2 += (c[2] * correctorX);
    particle->x3 += (c[3] * correctorX);
    particle->x4 += (c[4] * correctorX);
    particle->x5 += (c[5] * correctorX);

    correctorY = (particle->fy/particle->mass - particle->y2);
    particle->y0 += (c[0] * correctorY);
    particle->y1 += (c[1] * correctorY);
    particle->y2 += (c[2] * correctorY);
    particle->y3 += (c[3] * correctorY);
    particle->y4 += (c[4] * correctorY);
    particle->y5 += (c[5] * correctorY);

    return;
}

void corrector_rotations(disk *particle) {
    double *c = global.c;
    double correctorW;

    correctorW = (particle->fw/particle->iMoment - particle->w2);
    particle->w0 += (c[0] * correctorW);
    particle->w1 += (c[1] * correctorW);
    particle->w2 += (c[2] * correctorW);
    particle->w3 += (c[3] * correctorW);
    particle->w4 += (c[4] * correctorW);
    particle->w5 += (c[5] * correctorW);

    return;
}

void set_constants(double timestep) {
    double *a = global.a;
    double *c = global.c;

    a[0] = timestep;
    a[1]=a[0]*timestep/2.0;
    a[2]=a[1]*timestep/3.0;
    a[3]=a[2]*timestep/4.0;
    a[4]=a[3]*timestep/5.0;

    /*C[0] must be changed to 3/16 for forces dependent on velocity.*/
    c[0]=3.0/20.0*a[1]; //Should we try 3/16?
    c[1]=251.0/360.0*a[1]/a[0];
    c[2]=1.0;
    c[3]=11.0/18.0*a[1]/a[2];
    c[4]=1.0/6.0*a[1]/a[3];
    c[5]=1.0/60.0*a[1]/a[4];

    return;
}

extern gsl_rng * rgen;
void boundary_conditions(disk *particle, int b_cond, long kstep) {
    double epsilon = global.epsilon;
    double freq = global.freq;
    double time = global.time;
    double phase = global.phase;
    long stepsForVib = (long)(1/(freq*global.timestep));

    switch (b_cond) {
    case 1:
        /*Sinusoidal displacement*/
        if (global.vibrating) {
        particle->y0 = particle->yi + epsilon*sin(2*M_PI*freq*time - phase);
        particle->y1 = epsilon*2*M_PI*freq*cos(2*M_PI*freq*time - phase);
        }
        break;
    case 2:
        /*Random vibration*/
        if (kstep % stepsForVib == 0) {
            particle->y0 = particle->yi + epsilon*gsl_ran_flat(rgen, -1, 1);
            particle->x0 = particle->xi + epsilon*gsl_ran_flat(rgen, -1, 1);
        }
        break;
    default:
        /*Free relaxation*/
        printf("No excitation implemented\n");
        exit(0);
        break;
    }
    return;
}

void boundary_conditions_walls(long i) {
    double freq = global.freq;
    double time = global.time;
    double phase = global.phase;
    double Dy = global.epsilon*sin(2*M_PI*freq*time - phase);
    double Dx = 0;

    switch (global.bCondType) {
    case 1:
        /*Sinusoidal displacement*/
        if (global.vibrating) {
            particle[i].p1x = particle[i].p1x0 + Dx;
            particle[i].p1y = particle[i].p1y0 + Dy;
            particle[i].p2x = particle[i].p2x0 + Dx;
            particle[i].p2y = particle[i].p2y0 + Dy;

            particle[i].y1 = global.epsilon*
                2*M_PI*freq*cos(2*M_PI*freq*time - phase);
            particle[i].x1 = 0;
        }
        break;
    default:
        /*Free relaxation*/
        printf("No excitation implemented\n");
        exit(0);
        break;
    }
    return;
}
