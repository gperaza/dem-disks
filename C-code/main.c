/*
  File name: main.c

  This program simulates the time evolution of a set of disks in 2D.
*/
//#define GRAPHICS

#include "common.h"
#include <string.h>
#include <time.h>

disk *particle;
neighbor_stats *nStats;
particleParameters diskParameters;
systemParameters global;
links_3disks links3;

FILE *fPhase, *fFirst, *fLast, *fLinkStat, *fLinks, *fEnergy;
FILE *fp;

gsl_rng * rgen;

void get_input();
void clock_time(int);
void write_results();
void phase_plot(FILE*);
void write_3disk_avglinkstat();

/******************************************************************************/

void get_input() {

/******************************************************************************/
    char input[200];
    char type[200];
    char value[200];
    int inputs = 23, countInputs = 0;

    if ( (fp = fopen("input_file", "r")) == NULL ) {
        printf("Error opening input file.\n");
        exit(1);
    }

    while (fgets(input, sizeof(input), fp) != NULL) {
        sscanf(input,"%s %s",type, value);
        countInputs++;
        if (strcmp(type,"#nParticles") == 0) {
            sscanf(value, "%ld", &global.nParticles);
            printf("Particles = %ld.\n", global.nParticles);
        } else if (strcmp(type,"#seed") == 0) {
            sscanf(value, "%ld", &global.seed);
            printf("Seed for rng = %ld.\n", global.seed);
        } else if (strcmp(type,"#box_w") == 0) {
            sscanf(value, "%lf", &global.box_w);
            printf("box_w = %e disks.\n", global.box_w);
        } else if (strcmp(type,"#box_h") == 0) {
            sscanf(value, "%lf", &global.box_h);
            printf("box_h = %e disks.\n", global.box_h);
        } else if (strcmp(type,"#freq") == 0) {
            sscanf(value, "%lf", &global.freq);
            printf("Bottom frequency = %e Hz.\n", global.freq);
        } else if (strcmp(type,"#dimensionlessAc") == 0) {
            sscanf(value, "%lf", &global.dimensionlessAc);
            printf("Bottom acceleration = %e g's.\n", global.dimensionlessAc);
        } else if (strcmp(type,"#gravity") == 0) {
            sscanf(value, "%lf", &global.gravity);
            printf("Gravity = %e m/s^2.\n", global.gravity);
        } else if (strcmp(type,"#gravityAngle") == 0) {
            sscanf(value, "%lf", &global.gravityAngle);
            printf("Gravity angle = %e PI radians.\n", global.gravityAngle);
        } else if (strcmp(type,"#bGamma") == 0) {
            sscanf(value, "%lf", &global.bGamma);
            printf("Bulk gamma = %e. \n", global.bGamma);
        } else if (strcmp(type,"#timestep") == 0) {
            sscanf(value, "%lf", &global.timestep);
            printf("Time step = %e s. \n", global.timestep);
        } else if (strcmp(type,"#relaxTime") == 0) {
            sscanf(value, "%lf", &global.relaxTime);
            printf("Relaxing time = %e s. \n", global.relaxTime);
        } else if (strcmp(type,"#runTime") == 0) {
            sscanf(value, "%lf", &global.runTime);
            printf("Simulation time = %e s. \n", global.runTime);
        } else if (strcmp(type,"#timeForGraph") == 0) {
            sscanf(value, "%lf", &global.timeForGraph);
            printf("Time between graphics = %e s. \n", global.timeForGraph);
        } else if (strcmp(type,"#timeForWrite") == 0) {
            sscanf(value, "%lf", &global.timeForWrite);
            printf("Time between writes = %e s. \n", global.timeForWrite);
        } else if (strcmp(type,"#meanR") == 0) {
            sscanf(value, "%lf", &diskParameters.meanR);
            printf("Mean disk radius = %e m. \n", diskParameters.meanR);
        } else if (strcmp(type,"#density") == 0) {
            sscanf(value, "%lf", &diskParameters.density);
            printf("Disk's density = %e kg/m^3. \n", diskParameters.density);
        } else if (strcmp(type,"#kn") == 0) {
            sscanf(value, "%lf", &diskParameters.kn);
            printf("Normal elastic constant = %e. \n", diskParameters.kn);
        } else if (strcmp(type,"#pr") == 0) {
            sscanf(value, "%lf", &diskParameters.pr);
            printf("Disk's Poisson's ratio = %e. \n", diskParameters.pr);
        } else if (strcmp(type,"#mu") == 0) {
            sscanf(value, "%lf", &diskParameters.mu);
            printf("Disk's friction coefficient = %e. \n", diskParameters.mu);
        } else if (strcmp(type,"#vGamma") == 0) {
            sscanf(value, "%lf", &diskParameters.vGamma);
            printf("Viscoelastic gamma = %e. \n", diskParameters.vGamma);
        } else if (strcmp(type,"#bCondType") == 0) {
            sscanf(value, "%d", &global.bCondType);
            printf("Boundary conditions for bottom = %d. \n", global.bCondType);
        } else if (strcmp(type,"#thermalTime") == 0) {
            sscanf(value, "%lf", &global.thermalTime);
            printf("Thermalization time = %lf. \n", global.thermalTime);
        } else if (strcmp(type,"#relInitDisp") == 0) {
            sscanf(value, "%lf", &global.relInitDisp);
            printf("Relative initial displacement for bottom (-1,1) = %lf. \n",
                   global.relInitDisp);
        } else {
            printf("Unknown parameter.\n");
            exit(1);
        }
    }
    assert(countInputs == inputs);
    fclose(fp);

    diskParameters.kt = (1 - diskParameters.pr)/(1 - 0.5*diskParameters.pr)*
        diskParameters.kn;
    global.box_w *= 2*diskParameters.meanR;
    global.box_h *= 2*diskParameters.meanR;
    global.gravityAngle *= M_PI;
    if (global.bCondType == 1) {
        global.epsilon = global.dimensionlessAc*global.gravity/
            (4*M_PI*M_PI*global.freq*global.freq);
    } else if (global.bCondType == 2) {
        global.epsilon = global.dimensionlessAc*diskParameters.meanR;
    } else {
        printf("No excitation chosen");
        exit(0);
    }
    return;
}

/******************************************************************************/

int main(/*int argc, char *argv[]*/)

/******************************************************************************/
{
    time_t iTime = time(NULL);
    get_input();

    double timestep = global.timestep;
    long nStepsRelax = (long)(global.relaxTime/timestep);
    long nStepsThermal = (long)(global.thermalTime/timestep);
    long nStepsRun = (long)(global.runTime/timestep);
    long stepsForGraph = (long)(global.timeForGraph/timestep);
    long stepsForWrite = (long)(global.timeForWrite/timestep);
    global.time = -timestep*nStepsRelax - timestep*nStepsThermal;

    /*Set up constants for the gear integrator.*/
    set_constants(timestep);

    fPhase = fopen("phase_space.out", "w");
    fFirst = fopen("initial_phase_space.out", "w");
    fLast = fopen("ending_phase_space.out", "w");
    fLinkStat = fopen("3disk_linkstat.out", "w");
    fLinks = fopen("links.out", "w");
    fEnergy = fopen("elastic_energy.out", "w");

    /*Initialize the packing and make nParticles include walls.*/
    printf("Simulating for %ld particles.\n", global.nParticles);
    global.nParticles = init_system();
    printf("Number of particles including walls: %ld.\n", global.nParticles);

    /*Initialize linked cells*/
    init_cell();
    /*------------------------------------------------------------------------*/
    int relaxing = 1;
    long i;
    double gravityAngleAux = global.gravityAngle;
    global.gravityAngle = 0;
    double epsilonAux = global.epsilon;
    global.epsilon = 0;
    for (i = 0; i < nStepsRelax; i++) {
#ifdef GRAPHICS
        if (!(i % stepsForGraph)) {
            graphics(relaxing);
        }
#endif
        step(i);
        if (!(i % stepsForWrite)) {
            write_results();
        }
        global.bGamma = (nStepsRelax-1-i)/(nStepsRelax-1)*100000;
    }
    /*Tilt system and turn on vibration. Thermalization.
     Define phase such that vibration starts smoothly from zero.*/
    global.phase = 2*M_PI*global.freq*global.time - asin(global.relInitDisp);
    global.gravityAngle = gravityAngleAux;
    global.epsilon = epsilonAux;
    for (i = 0; i < nStepsThermal; i++) {
#ifdef GRAPHICS
        if (!(i % stepsForGraph)) {
            graphics(relaxing);
        }
#endif
        if (!(i % stepsForWrite)) {
            write_results();
        }
        step(i);
    }
    /*------------------------------------------------------------------------*/
    phase_plot(fFirst);
    for (i = 0; i < global.nParticles; i++) {
        if (particle[i].type == 0) {
            particle[i].wi = particle[i].w0;
            particle[i].xi = particle[i].x0;
            particle[i].yi = particle[i].y0;
        }
    }
    //Restart 3disk statistics to forget relaxation.
    links3.op_op = links3.op_sl = links3.op_cl = links3.sl_op = links3.sl_sl =
        links3.sl_cl = links3.cl_op = links3.cl_sl = links3.cl_cl = 0;
    links3.lstkount = 0;
    /*------------------------------------------------------------------------*/
    relaxing = 0;
    for (i = 0; i < nStepsRun; i++) {
#ifdef GRAPHICS
        if (!(i % stepsForGraph)) {
            graphics(relaxing);
        }
#endif
        if (!(i % stepsForWrite)) {
            write_results();
        }
        step(i);
    }

    phase_plot(fLast);
    if (global.nParticles == 3) {
        write_3disk_avglinkstat();
    }

    free(particle); free_cell(); free(nStats);
    fclose(fFirst); fclose(fPhase); fclose(fLast);
    fclose(fLinkStat); fclose(fLinks); fclose(fEnergy);
    clock_time((int)iTime);
    return 0;
}

/******************************************************************************/

long init_system() {

/******************************************************************************/
    long nParticles = global.nParticles;
    double box_w = global.box_w;
    double box_h = global.box_h;
    double meanR = diskParameters.meanR;
    double density = diskParameters.density;
    long seed = global.seed;

    long i, j;
    double layerY = 0;
    int layerN = 1;
    double sumWidth = 0;
    double wallR = meanR;
    long nBottom = (long)ceil(box_w/(2*wallR));
    long nDisks = nParticles;
    double sigma = meanR*0.25;
    disk tempParticle;

    rgen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rgen, seed);

    assert(nDisks > 0);
    assert(nBottom > 0);

    particle = (disk*)calloc(nParticles + nBottom, sizeof(disk));

    /*Setup the walls. Type = 1 for fixed disks*/
    for (i = nParticles; i < nParticles + nBottom; i++) {
        particle[i].type = 1;
        particle[i].radius = wallR;
        particle[i].x0 = sumWidth + particle[i].radius;
        particle[i].y0 = particle[i].radius;
        sumWidth += 2*particle[i].radius;
        /*Store initial positions for boundary conditions*/
        particle[i].xi = particle[i].x0;
        particle[i].yi = particle[i].y0;
        /*Set initial displacement of bottom*/
        particle[i].yi += global.epsilon*global.relInitDisp;
        /*Set disk properties*/
        particle[i].mass = density*particle[i].radius*particle[i].radius*M_PI;
        particle[i].iMoment =
            particle[i].mass*particle[i].radius*particle[i].radius/2;
    }

    /*Setup mobile particles. Type = 0 for mobile disk.*/
    /*Considering 3 different radius around a mean.*/
    for (i = 0; i < nDisks; i++) {
        particle[i].type = 0;
        if (i < (long)(nParticles/3))
            particle[i].radius = meanR - sigma;
        else if (i < (long)(2*nParticles/3))
            particle[i].radius = meanR + sigma;
        else
            particle[i].radius = meanR;
        particle[i].mass = density*particle[i].radius*particle[i].radius*M_PI;
        particle[i].iMoment =
            particle[i].mass*particle[i].radius*particle[i].radius/2;
    }

    /*Now make a random permutation (Knuth permutation).*/
    for (i = 0; i < nParticles; i++) {
        j = gsl_rng_uniform_int(rgen, i + 1);
        tempParticle = particle[i];
        particle[i] = particle[j];
        particle[j] = tempParticle;
    }

    /*Now set the x position such that contiguous particles barely touch.*/
    sumWidth = 0 + (layerN % 2)*meanR;
    layerY = 2*wallR + (meanR + sigma);
    for (i = 0; i < nParticles; i++) {
        if (sumWidth + 2*particle[i].radius > box_w) {
            layerY += 2*(meanR + sigma);
            layerN++;
            sumWidth = 0 + (layerN % 2)*meanR;
        }
        particle[i].x0 = sumWidth + particle[i].radius;
        particle[i].y0 = layerY;
        sumWidth += 2*particle[i].radius;
    }
    assert(layerY + (meanR + sigma)  < box_h);

    gsl_rng_free(rgen);

    nParticles += nBottom;
    nStats = (neighbor_stats*)calloc((nParticles*nParticles-nParticles)/2 ,
                             sizeof(neighbor_stats));
    for (i = 0; i < (nParticles*nParticles-nParticles)/2; i++) {
        nStats[i].touching = 0;
    }

    /*Restart RNG for simulation.*/
    rgen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rgen, global.seed);

    return (nParticles);
}

/******************************************************************************/

void step(long kstep) {

/******************************************************************************/

    long i = 0;
    long nParticles = global.nParticles;

    global.time += global.timestep;

    global.meanLinkSat = global.meanSqLinkSat = 0;
    global.linkCount = global.changingLinks = global.slidingLinks = 0;
    global.potEnergyElasNorm = global.potEnergyElasTg = 0;

    for (i = 0; i < nParticles; i++) {
        if (particle[i].type == 0) {
            if (particle[i].posFixed == 0) predictor_positions(&particle[i]);
            if (particle[i].rotFixed == 0) predictor_rotations(&particle[i]);
        }
        if (particle[i].type == 1) boundary_conditions(&particle[i],
                                                       global.bCondType, kstep);
    }
    make_forces();
    //make_forces_linked_cell();
    for (i = 0; i < nParticles; i++) {
        if (particle[i].type == 0) {
            if (particle[i].posFixed == 0) corrector_positions(&particle[i]);
            if (particle[i].rotFixed == 0) corrector_rotations(&particle[i]);
        }
        //Implement periodic boundary conditions if necessary (High vel drifts).
    }

    /*Link statistics for the 3disk setup*/
    if (global.nParticles == 3) {
        links3.lstkount++;

        links3.op_op += (1 - nStats[0].touching)*(1 - nStats[1].touching);
        links3.op_sl += (1 - nStats[0].touching)*(nStats[1].sliding);
        links3.op_cl += (1 - nStats[0].touching)*(nStats[1].touching)
            *(1 - nStats[1].sliding);

        links3.sl_op += (nStats[0].sliding)*(1 - nStats[1].touching);
        links3.sl_sl += (nStats[0].sliding)*(nStats[1].sliding);
        links3.sl_cl += (nStats[0].sliding)*(nStats[1].touching)
            *(1 - nStats[1].sliding);

        links3.cl_op += (1 - nStats[0].sliding)*(nStats[0].touching)
            *(1 - nStats[1].touching);
        links3.cl_sl += (1 - nStats[0].sliding)*(nStats[0].touching)
            *(nStats[1].sliding);
        links3.cl_cl += (1 - nStats[0].sliding)
            *(nStats[0].touching*nStats[1].touching)*(1 - nStats[1].sliding);
    }

    return;
}

/******************************************************************************/

void clock_time(int iTime) {

/******************************************************************************/
/*
  Purpose:

    CLOCK_TIME prints CPU and wall time of execution.

  Discussion:

    Prints CPU time and wall time in human readable format since a starting time
    given as a parameter.

  Modified:

    15 May 2014 by GGPM.

  Author:

    Gonzalo G. Peraza Mues (GGPM)

  Parameters:

    int iTime - Initial time.
*/
    int days, hours, minutes, seconds;
    int tTime = (int)time(NULL) - iTime;
    int t;

    days = tTime/(86400);
    t = tTime%(86400);
    hours = t/3600;
    t = t%3600;
    minutes = t/60;
    t = t%60;
    seconds = t;

    printf("CPU time = %f seconds.\n", (double)clock()/CLOCKS_PER_SEC);
    printf("Wall time = %d days, %d hours, %d minutes, %d seconds.\n", days,
           hours, minutes, seconds);
    return;
}
