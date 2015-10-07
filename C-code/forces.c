#include "common.h"
extern neighbor_stats *nStats;

double normal_force_disk_disk(long, long, double, double,
                              double, double);
double tangential_force_disk_disk(double, long, long, double,
                                  double, double, double, double);
void pair_force_disks(long, long);
void pair_force_disk_line(long, long);

void make_forces() {
    long i, j;
    long nParticles = global.nParticles;
    double cga = cos(global.gravityAngle);
    double sga = sin(global.gravityAngle);

    for (i = 0; i < nParticles; i++) {
        particle[i].fx = 0; particle[i].fy = 0; particle[i].fw = 0;
    }

    for (i = 0; i < nParticles; i++) {

        if (particle[i].type == 0) {
            bulk_force(i, cga, sga);
        }

        for (j = i + 1; j < nParticles; j++) {
            pair_force(i, j);
        }
    }

    return;
}

void bulk_force(long i, double cga, double sga) {
    double gravity = global.gravity;
    double bGamma = global.bGamma;

    /*A positive gravity angle tilts the system counter clockwise,
      i.e.  in this reference frame the gravity vector rotates
      clockwise.*/
    particle[i].fx += -particle[i].mass*gravity*sga;
    particle[i].fy += -particle[i].mass*gravity*cga;

    particle[i].fx += - bGamma * particle[i].x1*particle[i].mass;
    particle[i].fy += - bGamma * particle[i].y1*particle[i].mass;

    return;
}

void pair_force(long i, long j) {

    double box_w = global.box_w;
    double box_h = global.box_h;

    double rx12, ry12, r12sq, r12, radsum;
    double normalForce = 0, tangentialForce = 0;
    long l = i*global.nParticles + (j-1) - i*(i+3)/2;

#ifdef COLLISIONS
    nStats[l].pretouching = nStats[l].touching;
#endif

    if (particle[i].type + particle[j].type <= 1) {
        rx12 = particle[i].x0 - particle[j].x0;
        ry12 = particle[i].y0 - particle[j].y0;
        radsum = particle[i].radius + particle[j].radius;
    } else if (particle[i].type + particle[j].type == 2) {
        /* Wall index > disk index. Find the closest point on the
         line (p1x,p1y)-(p2x,p2y) using the projected length and subtract.*/
        double lDx = particle[j].p2x - particle[j].p1x;
        double lDy = particle[j].p2y - particle[j].p1y;
        double lenLineSqrd = lDx*lDx + lDy*lDy;
        double proy = ((particle[i].x0 - particle[j].p1x)*lDx +
                       (particle[i].y0 - particle[j].p1y)*lDy) / (lenLineSqrd);
        rx12 = particle[i].x0 - (particle[j].p1x + proy*(lDx));
        ry12 = particle[i].y0 - (particle[j].p1y + proy*(lDy));
        radsum = particle[i].radius;
    } else {
        return;
    }

    // Periodic boundary conditions on the RELATIVE VECTOR
    rx12  = rx12 - box_w*round( rx12/box_w );
    ry12  = ry12 - box_h*round( ry12/box_h );
    r12sq = rx12*rx12 + ry12*ry12;

    if (r12sq < radsum*radsum) {

        /*Link statistics*/
        global.linkCount++;
        if (nStats[l].touching == 0) {
            global.changingLinks++;
        }

        r12 = sqrt(r12sq);
        /*Normalized versor*/
        double rx12n = rx12/r12;
        double ry12n = ry12/r12;

        normalForce = normal_force_disk_disk(i, j, r12, rx12n, ry12n, radsum);
        tangentialForce = tangential_force_disk_disk(normalForce, i, j,
                                                     rx12, ry12, rx12n, ry12n,
                                                     radsum);

        /*Add normal forces.*/
        particle[i].fx += rx12n*normalForce;
        particle[i].fy += ry12n*normalForce;
        particle[j].fx -= rx12n*normalForce;
        particle[j].fy -= ry12n*normalForce;
        /*Add tangential forces (tangential versor (ry12n, -rx12n))*/
        particle[i].fx += ry12n*tangentialForce;
        particle[i].fy += -rx12n*tangentialForce;
        particle[j].fx -= ry12n*tangentialForce;
        particle[j].fy -= -rx12n*tangentialForce;
        /*Add torques*/
        particle[i].fw += tangentialForce*r12*particle[i].radius/(radsum);
        particle[j].fw += tangentialForce*r12*particle[j].radius/(radsum);
        //particle[i].fw = tangentialForce*particle[i].radius;
        //particle[j].fw = tangentialForce*particle[j].radius;

    } else {
        if (nStats[l].touching == 1) {
            global.changingLinks++;
        }
        nStats[l].touching = 0;
        nStats[l].nForce = nStats[l].tForce = 0;
        nStats[l].stretch = 0;
        nStats[l].sliding = 0;
    }

    return;
}

double normal_force_disk_disk(long i, long j, double r12,
                              double rx12n, double ry12n, double radsum) {
    double kn = diskParameters.kn;
    double vGamma = diskParameters.vGamma;
    double effMass;

    double normalForce = 0;

    /*Relative velocities*/
    double vx12 = particle[i].x1 - particle[j].x1;
    double vy12 = particle[i].y1 - particle[j].y1;

    if (particle[i].type == particle[j].type){
        effMass = particle[i].mass*particle[j].mass/(particle[i].mass
                                                     + particle[j].mass);
    } else if (particle[i].type != 0) {
        effMass = particle[j].mass;
    } else {
        effMass = particle[i].mass;
    }

    normalForce = kn*(radsum - r12);
    normalForce -= vGamma*sqrt(effMass)*(rx12n*vx12 + ry12n*vy12);
    normalForce = fmax(0, normalForce);

    global.potEnergyElasNorm += 0.5*kn*(radsum - r12)*(radsum - r12);

    return normalForce;
}

double tangential_force_disk_disk(double normalForce, long i,
                                  long j, double rx12, double ry12,
                                  double rx12n, double ry12n,
                                  double radsum){
    double kt = diskParameters.kt;
    double mu = diskParameters.mu;
    long l = i*global.nParticles + (j-1) - i*(i+3)/2;

    double tangentialForce;
    double coulombLimit = mu*normalForce;
    double ss, d_ss;

    /*Calculate tangential displacement.*/
    if (particle[i].type + particle[j].type <= 1) {
        ss = particle[i].w0*particle[i].radius
            + particle[j].w0*particle[j].radius
            - atan2(ry12, rx12)*radsum;
    } else {
        assert(particle[i].type + particle[j].type == 2);
        /*The tangent versor (tv) is (ry12n, -rx12n).*/
        /*Tangential displacement is calculated as
          R1w01 + (vec(ri) - vec(p1)) dot tversor.*/
        ss = particle[i].w0*particle[i].radius
            + (particle[i].x0 - particle[j].p1x)*ry12n -
            (particle[i].y0 - particle[j].p1y)*rx12n;
    }

    /*Check if they where touching before.*/
    if (nStats[l].touching == 1) {
        d_ss = ss - nStats[l].s0;
        d_ss = d_ss -
            2*M_PI*radsum*copysign((fabs(d_ss) > M_PI*radsum), d_ss);
    } else {
        nStats[l].s0 = ss;
        d_ss = 0;
        nStats[l].touching = 1;
    }

    tangentialForce = -kt*d_ss;
    /*Check if contact is sliding*/
    if (fabs(tangentialForce) >= coulombLimit) {
        tangentialForce = copysign(coulombLimit, -d_ss);
        d_ss = copysign(mu*normalForce/kt, d_ss);
        nStats[l].s0 = ss - d_ss;
        nStats[l].sliding = 1;
        global.slidingLinks++;
    } else {
        nStats[l].sliding = 0;
    }

    /*Link statistics*/
    nStats[l].stretch = d_ss;
    nStats[l].nForce = normalForce;
    nStats[l].tForce = tangentialForce;
    global.meanLinkSat += fabs(tangentialForce)/mu/normalForce;
    global.meanSqLinkSat +=
        tangentialForce/mu/normalForce*tangentialForce/mu/normalForce;

    global.potEnergyElasTg += 0.5*kt*d_ss*d_ss;

    return tangentialForce;
}
