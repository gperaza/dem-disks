#include "common.h"

extern FILE *fPhase, *fFirst, *fLast, *fLinkStat, *fLinks, *fEnergy;
extern FILE *fp;
extern links_3disks links3;

void phase_plot(FILE*);
void write_3disk_linkstat(FILE*);
void write_3disk_avglinkstat();

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

    //Calculate and write data for the 3Disk setup.
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

    return;
}

void phase_plot(FILE* fp) {
    long i;

    /*Print the header of each frame with relevant information.*/
    fprintf(fp, "#New Frame\n"
            "#type:Disks3RadiiSinusoidalBottom\n"
            "#nParticles:%ld\n"
            "#time:%e\n"
            "#timestep:%e\n"
            "#box_w:%lf\n"
            "#box_h:%lf\n"
            "#freq:%lf\n"
            "#dimensionlessAc:%lf\n"
            "#gravity:%lf\n"
            "#gravityAngle:%lf\n"
            "#bGamma:%lf\n"
            "#vGamma:%lf\n"
            "#density:%lf\n"
            "#kn:%lf\n"
            "#kt:%lf\n"
            "#mu:%lf\n"
            "#EndOfHeader\n",
            global.nParticles,
            global.time,
            global.timestep,
            global.box_w,
            global.box_h,
            global.freq,
            global.epsilon,
            global.gravity,
            global.gravityAngle,
            global.bGamma,
            diskParameters.vGamma,
            diskParameters.density,
            diskParameters.kn,
            diskParameters.kt,
            diskParameters.mu
            );

    /*Print the state of the system*/
    /*
      id type radius mass iMoment
      x0 y0 w0 x1 y1 w1 x2 y2 w2 x3 y3 w3 x4 y4 w4 x5 y5 w5
      fx fy fw
      xi yi posFixed rotFixed
    */
    for (i = 0; i < global.nParticles; i++) {
        fprintf(fp,"%ld %d %lf %lf %lf"
                " %lf %lf %lf %lf %lf %lf %lf %lf %lf"
                " %lf %lf %lf %lf %lf %lf %lf %lf %lf"
                " %lf %lf %lf"
                " %lf %lf %d %d"
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
                particle[i].xi, particle[i].yi,
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
    fclose(fp);

    return;
}
