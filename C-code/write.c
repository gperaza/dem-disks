#include "common.h"

extern FILE *fPhase, *fFirst, *fLast, *fLinkStat, *fLinks, *fEnergy;
extern FILE *fp;
extern links_3disks links3;
extern FILE *fCollisions;

void phase_plot(FILE*);
void write_3disk_linkstat(FILE*);
void write_3disk_avglinkstat();
void search_collisions();
void write_collision(long, int, unsigned long);
void save_collision_data(long, int);

void write_results() {

    phase_plot(fPhase);

    if (global.linkCount) {
        global.meanLinkSat = global.meanLinkSat/global.linkCount;
        global.meanSqLinkSat = global.meanSqLinkSat/global.linkCount;
    }
    fprintf(fLinks, "%e %lu %lu %lu %e %e\n",
            global.time,
            global.linkCount,
            global.slidingLinks,
            global.changingLinks,
            global.meanLinkSat,
            global.meanSqLinkSat
            );

    fprintf(fEnergy, "%e %e %e\n",
            global.time, global.potEnergyElasNorm, global.potEnergyElasTg);

    /*Calculate and write data for the 3Disk setup.*/
    /*0-time
      1-touching  2-sliding  3-nForce   4-tForce   5-stretch   6-s0
      7-touching  8-sliding  9-nForce  10-tForce  11-stretch  12-s0
    */
    if (global.nParticles == 3) {
        fprintf(fLinkStat,
                "%18.12e %u %u %18.12e %18.12e %18.12e %18.12e"
                " %u %u %18.12e %18.12e %18.12e %18.12e\n",
                global.time,
                nStats[0].touching, nStats[0].sliding,
                nStats[0].nForce, nStats[0].tForce,
                nStats[0].stretch, nStats[0].s0,
                nStats[1].touching, nStats[1].sliding,
                nStats[1].nForce, nStats[1].tForce,
                nStats[1].stretch, nStats[1].s0
                );
    }

    fflush(fLinks);
    fflush(fEnergy);
    fflush(fLinkStat);

#ifdef COLLISIONS
    search_collisions();
#endif
    return;
}

void phase_plot(FILE* fp) {
    long i;

    /*Print the header of each frame with relevant information.*/
    fprintf(fp, "#NewFrame\n"
            "#type:Disks3RadiiSinusoidalBottom\n"
            "#nParticles:%ld\n"
            "#time:%e\n"
            "#timestep:%e\n"
            "#box_w:%e\n"
            "#box_h:%e\n"
            "#freq:%e\n"
            "#dimensionlessAc:%e\n"
            "#epsilon:%e\n"
            "#gravity:%e\n"
            "#gravityAngle:%e\n"
            "#bGamma:%e\n"
            "#vGamma:%e\n"
            "#density:%e\n"
            "#kn:%e\n"
            "#kt:%e\n"
            "#mu:%e\n"
            "#bCondType:%d\n"
            "#EndOfHeader\n",
            global.nParticles,
            global.time,
            global.timestep,
            global.box_w,
            global.box_h,
            global.freq,
            global.dimensionlessAc,
            global.epsilon,
            global.gravity,
            global.gravityAngle,
            global.bGamma,
            diskParameters.vGamma,
            diskParameters.density,
            diskParameters.kn,
            diskParameters.kt,
            diskParameters.mu,
            global.bCondType
            );

    /*Print the state of the system*/
    /*
      0-id 1-type 2-radius 3-mass 4-iMoment
      5-x0 6-y0 7-w0 8-x1 9-y1 10-w1 11-x2 12-y2 13-w2
      14-x3 15-y3 16-w3 17-x4 18-y4 19-w4 20-x5 21-y5 22-w5
      23-fx 24-fy 25-fw
      26-xi 27-yi 28-wi 29-posFixed 30-rotFixed
      TODO: add fn and ft to this output to make it easier to analyze collisions.
    */
    for (i = 0; i < global.nParticles; i++) {
        fprintf(fp,"%ld %d %e %e %e"
                " %e %e %e %e %e %e %e %e %e"
                " %e %e %e %e %e %e %e %e %e"
                " %e %e %e"
                " %e %e %e %d %d"
                "\n",
                i, particle[i].type,
                particle[i].radius, particle[i].mass, particle[i].iMoment,
                particle[i].x0, particle[i].y0, particle[i].w0,
                particle[i].x1, particle[i].y1, particle[i].w1,
                particle[i].x2, particle[i].y2, particle[i].w2,
                particle[i].x3, particle[i].y3, particle[i].w3,
                particle[i].x4, particle[i].y4, particle[i].w4,
                particle[i].x5, particle[i].y5, particle[i].w5,
                particle[i].fx, particle[i].fy, particle[i].fw,
                particle[i].xi, particle[i].yi, particle[i].wi,
                particle[i].posFixed, particle[i].rotFixed
                );
    }

    fflush(fp);
    return;
}

void write_3disk_avglinkstat() {
    FILE *fp;

    assert(global.nParticles == 3);

    fp = fopen("3disk_avglinkstat.out", "w");
    fprintf(fp,  "%18.12e %18.12e %18.12e %18.12e %18.12e"
            " %18.12e %18.12e %18.12e %18.12e\n",
            links3.op_op/links3.lstkount, links3.op_sl/links3.lstkount,
            links3.op_cl/links3.lstkount, links3.sl_op/links3.lstkount,
            links3.sl_sl/links3.lstkount, links3.sl_cl/links3.lstkount,
            links3.cl_op/links3.lstkount, links3.cl_sl/links3.lstkount,
            links3.cl_cl/links3.lstkount
            );
    fflush(fp);
    fclose(fp);

    return;
}

void search_collisions() {
    int touch1 = nStats[0].touching;
    int touch2 = nStats[1].touching;
    int pretouch1 = nStats[0].pretouching;
    int pretouch2 = nStats[1].pretouching;
    static long i = 0;
    static int colliding = 0;
    static unsigned long collisionN = 0;

    if (!colliding) {  /*Search for a collision.*/
        if(touch1 + touch2 == 1 && pretouch1 + pretouch2 == 0) {
            colliding = (touch1 == 1) ? 1 : 2;
            i = 0;
            save_collision_data(i, colliding);
        }
    } else { /*Already colliding.*/
        if (touch1 == pretouch1 && touch2 == pretouch2) {
            /*Collision continues.*/
            i++;
            save_collision_data(i, colliding);
        } else if (!touch1 && !touch2) {
            /*Collision ends successfully.*/
            collisionN ++;
            i++;
            write_collision(i, colliding, collisionN);
            colliding = 0;
        } else {
            /*All other possibilities invalidate the collision.*/
            colliding = 0;
        }
    }

}

collision_temp_array coll_data[100000];
void save_collision_data(long i, int colliding) {
    i = i % 100000; /*To prevent segmentation fault.*/

    /*Calculate normal and tangential vectors*/
    double rx12 = particle[0].x0 - particle[colliding].x0;
    double ry12 = particle[0].y0 - particle[colliding].y0;
    /* Periodic boundary conditions on the RELATIVE VECTOR*/
    rx12  = rx12 - global.box_w*round( rx12/global.box_w );
    ry12  = ry12 - global.box_h*round( ry12/global.box_h );
    double r12 = sqrt(rx12*rx12 + ry12*ry12);
    /*Normal versor*/
    double rx12n = rx12/r12;
    double ry12n = ry12/r12;
    /*Relative velocities*/
    double vx12 = particle[0].x1 - particle[colliding].x1;
    double vy12 = particle[0].y1 - particle[colliding].y1;
    double vn = rx12n*vx12 + ry12n*vy12;
    double vt = ry12n*vx12 - rx12n*vy12;
    double gt = vt + particle[0].radius*particle[0].w1
        - particle[colliding].radius*particle[colliding].w1;
    double vt0 = ry12n*particle[0].x1 - rx12n*particle[0].y1;

    /*Fill an array with valuable data.*/
    coll_data[i].time = global.time;
    coll_data[i].Fn = nStats[colliding-1].nForce;
    coll_data[i].Ft = nStats[colliding-1].tForce;
    coll_data[i].gn = vn;
    coll_data[i].gt = gt;
    coll_data[i].x = particle[0].x0;
    coll_data[i].y = particle[0].y0;
    coll_data[i].theta = particle[0].w0;
    coll_data[i].vx = particle[0].x1;
    coll_data[i].vy = particle[0].y1;
    coll_data[i].w = particle[0].w1;
    coll_data[i].vt0 = vt0;
}

void write_collision(long i, int colliding, unsigned long collisionN) {
    FILE *fCollision;
    char fname[50];
    int j = 0;

    save_collision_data(i, colliding);

    sprintf(fname, "Collisions/collision_%ld.out", collisionN);
    fCollision = fopen(fname, "w");
    for (j = 0; j <= i; j++) {
        fprintf(fCollision, "%e %e %e %e %e %e %e %e %e %e %e %e\n",
                coll_data[j].time - coll_data[0].time,
                coll_data[j].Fn, coll_data[j].Ft,
                coll_data[j].gn, coll_data[j].gt,
                coll_data[j].x, coll_data[j].y, coll_data[j].theta,
                coll_data[j].vx, coll_data[j].vy, coll_data[j].w,
                coll_data[j].vt0);
    }
    fclose(fCollision);

    /*0-colN   1-tc   2-dgn   3-dgt   4-dvt0   5-dvw0*/
    fprintf(fCollisions,"%ld %e %e %e %e %e\n",
            collisionN,
            coll_data[i].time - coll_data[0].time,
            coll_data[i].gn - coll_data[0].gn,
            coll_data[i].gt - coll_data[0].gt,
            coll_data[i].vt0 - coll_data[0].vt0,
            coll_data[i].w - coll_data[0].w);
    fflush(fCollisions);

}
