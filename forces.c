#include "common.h"

void bulk_force(disk*, double, double);
void pair_force(disk*, disk*);
double normal_force_disk_disk(disk*, disk*, double, double, double, double);
double tangential_force_disk_disk(double, disk*, disk*, double, double, double);

void make_forces() {
    unsigned long i, j;
    unsigned long nParticles = global.nParticles;
    double gravityAngle = global.gravityAngle;
    double cga = cos(gravityAngle);
    double sga = sin(gravityAngle);

    for (i = 0; i < nParticles; i++) {
        particle[i].fx = 0; particle[i].fy = 0; particle[i].fw = 0;
    }

    for (i = 0; i < nParticles; i++) {

        if (particle[i].type == 0) {
            bulk_force(&particle[i], cga, sga);
        }

        for (j = i + 1; j < nParticles; j++) {
            if (particle[i].type == 0 || particle[j].type == 0) {
                pair_force(&particle[i], &particle[j]);
            }
        }
    }

    return;
}

void bulk_force(disk* particle, double cga, double sga) {
    double gravity = global.gravity;
    double bGamma = global.bGamma;

    particle->fx += -particle->mass*gravity*sga;
    particle->fy += -particle->mass*gravity*cga;

    particle->fx += - bGamma * particle->x1*particle->mass;
    particle->fy += - bGamma * particle->y1*particle->mass;

    return;
}

void pair_force(disk* p1, disk* p2) {

    double box_w = global.box_w;
    double box_h = global.box_h;

    double rx12, ry12, r12sq, r12;
    double normalForce = 0, tangentialForce = 0;
    double radsum = p1->radius + p2->radius;

    rx12 = p1->x0 - p2->x0;
    ry12 = p1->y0 - p2->y0;

    // Periodic boundary conditions on the RELATIVE VECTOR
    rx12  = rx12 - box_w*round( rx12/box_w );
    ry12  = ry12 - box_h*round( ry12/box_h );
    r12sq = rx12*rx12 + ry12*ry12;

    if (r12sq < radsum*radsum) {

        r12 = sqrt(r12sq);
        /*Normalized versor*/
        double rx12n = rx12/r12;
        double ry12n = ry12/r12;

        normalForce = normal_force_disk_disk(p1, p2, r12, rx12n, ry12n, radsum);
        tangentialForce = tangential_force_disk_disk(normalForce, p1, p2,
                                                     rx12, ry12, radsum);

        /*Add normal forces.*/
        p1->fx += rx12n*normalForce;
        p1->fy += ry12n*normalForce;
        p2->fx -= rx12n*normalForce;
        p2->fy -= ry12n*normalForce;
        /*Add tangential forces*/
        p1->fx += -ry12n*tangentialForce;
        p1->fy += rx12n*tangentialForce;
        p2->fx -= -ry12n*tangentialForce;
        p2->fy -= rx12n*tangentialForce;
        /*Add torques*/
        p1->fw -= tangentialForce*r12*p1->radius/(radsum);
        p2->fw -= tangentialForce*r12*p2->radius/(radsum);
    }

    return;
}

double normal_force_disk_disk(disk* p1, disk* p2, double r12,
                              double rx12n, double ry12n, double radsum) {
    double kn = diskParameters.kn;
    double vGamma = diskParameters.vGamma;

    double normalForce = 0;

    /*Relative velocities*/
    double vx12 = p1->x1 - p2->x1;
    double vy12 = p1->y1 - p2->y1;

    double effMass = p1->mass*p2->mass/(p1->mass + p2->mass);

    normalForce = kn*(radsum - r12);
    normalForce -= vGamma*sqrt(effMass)*(rx12n*vx12 + ry12n*vy12);
    normalForce = fmax(0, normalForce);

    return normalForce;
}

double tangential_force_disk_disk(double normalForce, disk* p1, disk* p2,
                                  double rx12, double ry12, double radsum){
    double kt = diskParameters.kt;
    double mu = diskParameters.mu;

    double tangentialForce = 0;
    double coulombLimit = mu*normalForce;
    double ss, d_ss, s0 = 0;

    /*Calculate tangential displacement.*/
    ss = p1->w0*p1->radius + p2->w0*p2->radius - atan2(ry12, rx12)*radsum;

    /*Check if they where touching before.*/
    if (1) {
        d_ss = ss - s0;
        d_ss = d_ss -
            2*M_PI*radsum*copysign((fabs(d_ss) > M_PI*radsum), d_ss);
    } else {
        s0 = ss;
        d_ss = 0;
    }

    tangentialForce = kt*d_ss;
    /*Check if contact is sliding*/
    if (fabs(tangentialForce) >= coulombLimit) {
        tangentialForce = copysign(coulombLimit, d_ss);
        d_ss = copysign(mu*normalForce/kt, d_ss);
        s0 = ss - d_ss;
    }


    return tangentialForce;
}
