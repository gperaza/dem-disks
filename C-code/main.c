/*
  File name: main.c

  This program simulates the time evolution of a set of disks in 2D.
*/

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
FILE *fCollisions;
FILE *fPhase2;
FILE *fStatsDisks, *fStatsPacking;

gsl_rng * rgen;

void get_input();
void clock_time(int);
void write_results();
void phase_plot(FILE*);
void write_3disk_avglinkstat();
void search_collisions();
void write_collision_stats();
void init_system();
void init_system_wedge();
void init_system_file(FILE*);
void pack_disks();
void write_stats();

void get_input() {
    char input[200];
    char type[200];
    char value[200];
    int inputs = 24, countInputs = 0;

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
            printf("box_w = %g disks.\n", global.box_w);
        } else if (strcmp(type,"#box_h") == 0) {
            sscanf(value, "%lf", &global.box_h);
            printf("box_h = %g disks.\n", global.box_h);
        } else if (strcmp(type,"#freq") == 0) {
            sscanf(value, "%lf", &global.freq);
            printf("Bottom frequency = %g Hz.\n", global.freq);
        } else if (strcmp(type,"#dimensionlessAc") == 0) {
            sscanf(value, "%lf", &global.dimensionlessAc);
            printf("Bottom acceleration = %g g's.\n", global.dimensionlessAc);
        } else if (strcmp(type,"#gravity") == 0) {
            sscanf(value, "%lf", &global.gravity);
            printf("Gravity = %g m/s^2.\n", global.gravity);
        } else if (strcmp(type,"#gravityAngle") == 0) {
            sscanf(value, "%lf", &global.gravityAngle);
            printf("Gravity angle = %g PI radians.\n", global.gravityAngle);
        } else if (strcmp(type,"#bGamma") == 0) {
            sscanf(value, "%lf", &global.bGamma);
            printf("Bulk damping ratio = %g. \n", global.bGamma);
        } else if (strcmp(type,"#relaxTime") == 0) {
            sscanf(value, "%lf", &global.relaxTime);
            printf("Relaxing time = %g s. \n", global.relaxTime);
        } else if (strcmp(type,"#runTime") == 0) {
            sscanf(value, "%lf", &global.runTime);
            printf("Simulation time = %g s. \n", global.runTime);
        } else if (strcmp(type,"#timeForGraph") == 0) {
            sscanf(value, "%lf", &global.timeForGraph);
            printf("Time between graphics = %g s. \n", global.timeForGraph);
        } else if (strcmp(type,"#timeForWriteRun") == 0) {
            sscanf(value, "%lf", &global.timeForWriteRun);
            printf("Time between writes during main simulation= %g s. \n",
                   global.timeForWriteRun);
        } else if (strcmp(type,"#timeForWriteThermal") == 0) {
            sscanf(value, "%lf", &global.timeForWriteThermal);
            printf("Time between writes during thermalization= %g s. \n",
                   global.timeForWriteThermal);
        } else if (strcmp(type,"#meanR") == 0) {
            sscanf(value, "%lf", &diskParameters.meanR);
            printf("Mean disk radius = %g m. \n", diskParameters.meanR);
        } else if (strcmp(type,"#density") == 0) {
            sscanf(value, "%lf", &diskParameters.density);
            printf("Disk's density = %g kg/m^3. \n", diskParameters.density);
        } else if (strcmp(type,"#kn") == 0) {
            sscanf(value, "%lf", &diskParameters.kn);
            printf("Normal elastic constant = %g. \n", diskParameters.kn);
        } else if (strcmp(type,"#pr") == 0) {
            sscanf(value, "%lf", &diskParameters.pr);
            printf("Elastic constants ratio = %g. \n", diskParameters.pr);
        } else if (strcmp(type,"#mu") == 0) {
            sscanf(value, "%lf", &diskParameters.mu);
            printf("Disk's friction coefficient = %g. \n", diskParameters.mu);
        } else if (strcmp(type,"#vGamma") == 0) {
            sscanf(value, "%lf", &diskParameters.vGamma);
            printf("Viscoelastic damping = %g. \n", diskParameters.vGamma);
        } else if (strcmp(type,"#bCondType") == 0) {
            sscanf(value, "%d", &global.bCondType);
            printf("Boundary conditions for bottom = %d. \n", global.bCondType);
        } else if (strcmp(type,"#thermalTime") == 0) {
            sscanf(value, "%lf", &global.thermalTime);
            printf("Thermalization time = %g. \n", global.thermalTime);
        } else if (strcmp(type,"#relInitDisp") == 0) {
            sscanf(value, "%lf", &global.relInitDisp);
            printf("Relative initial displacement for bottom (-1,1) = %g. \n",
                   global.relInitDisp);
        } else if (strcmp(type,"#wedge") == 0) {
            sscanf(value, "%d", &global.wedge);
            printf("Simulating a wedge? = %d. \n",
                   global.wedge);
        } else {
            printf("Unknown parameter.\n");
            exit(1);
        }
    }
    assert(countInputs == inputs);
    fclose(fp);

    diskParameters.kt = diskParameters.pr*diskParameters.kn;
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

int main(/*int argc, char *argv[]*/) {
    /*First, print the version of the program used.*/
    printf("Version: %s\n", VERSION);

    time_t iTime = time(NULL);
    /*Clean up.*/
    // if (system("rm *.svg")) printf("No *.svg to delete.\n");
    // if (system("rm *.out")) printf("No *.out to delete.\n");

    get_input();

    /* Automatically adjust timestep to be one hundredth of the
       normal period of a disk oscillation. */
    global.timestep = 0.01*2*M_PI/
        sqrt(diskParameters.kn/(diskParameters.density*M_PI*
                                diskParameters.meanR*diskParameters.meanR));
    printf("Adjusted timestep = %g\n", global.timestep);

    double timestep = global.timestep;
    long nStepsRelax = (long)(global.relaxTime/timestep);
    long nStepsThermal = (long)(global.thermalTime/timestep);
    long nStepsRun = (long)(global.runTime/timestep);
#ifdef GRAPHICS
    long stepsForGraph = (long)(global.timeForGraph/timestep);
#endif
    long stepsForWriteRun = (long)(global.timeForWriteRun/timestep);
    long stepsForWriteThermal = (long)(global.timeForWriteThermal/timestep);
    long stepsForPrint = (long)(1.0/timestep);
    global.time = -timestep*nStepsRelax - timestep*nStepsThermal;

    /*Set up constants for the gear integrator.*/
    set_constants(timestep);

    fPhase = fopen("phase_space.out", "w");
    fPhase2 = fopen("phase_space_2.out", "w");
    fFirst = fopen("initial_phase_space.out", "w");
    fLast = fopen("ending_phase_space.out", "w");
    fLinkStat = fopen("3disk_linkstat.out", "w");
    fLinks = fopen("links.out", "w");
    fEnergy = fopen("elastic_energy.out", "w");
    fStatsDisks = fopen("disk_stats.out", "w");
    fStatsPacking = fopen("packing_stats.out", "w");
#ifdef COLLISIONS
    fCollisions = fopen("collisions.out", "w");
#endif

    /*Initialize the packing and make nParticles include walls, we
      now consider a wedge as a special case when wedge == 1.*/
    FILE *input_fptr;
    printf("Simulating for %ld particles.\n", global.nParticles);
    if (global.wedge == 1) {
        global.nParticles = 3;
        init_system_wedge();
        printf("Simulating a wedge.\n");
    } else if ((input_fptr = fopen("input_phase_space.out", "r")) != NULL) {
        init_system_file(input_fptr);
        fclose(input_fptr);
    } else {
        init_system();
    }
    printf("Number of particles including walls: %ld.\n", global.nParticles);

    /*------------------------------------------------------------------------*/
#ifdef GRAPHICS
    int relaxing = 1;
#endif
    long i;
    double gravityAngleAux = global.gravityAngle;
    global.gravityAngle = 0;
    global.vibrating = 0;
    for (i = 0; i < nStepsRelax; i++) {
        if (!(i % stepsForPrint)) {
            printf("\rRelaxing %.0lf%%", (double)i/nStepsRelax*100);
            fflush(stdout);}
#ifdef GRAPHICS
        if (!(i % stepsForGraph)) {
            graphics(relaxing);
        }
#endif
        step(i);
        //global.bGamma = (nStepsRelax-1-i)/(nStepsRelax-1)*1000;
    }
    printf("\rRelaxing 100%%\n");
    /*Tilt system and turn on vibration. Thermalization.
      Define phase such that vibration starts smoothly from zero.*/
    global.phase = 2*M_PI*global.freq*global.time - asin(global.relInitDisp);
    global.gravityAngle = gravityAngleAux;
    global.vibrating = 1;
    for (i = 0; i < nStepsThermal; i++) {
        if (!(i % stepsForPrint)) {
            printf("\rThermalizing %.0lf%%", (double)i/nStepsThermal*100);}
        fflush(stdout);
#ifdef GRAPHICS
        if (!(i % stepsForGraph)) {
            graphics(relaxing);
        }
#endif
        if (!(i % stepsForWriteThermal)) {
            //write_results();
        }
        step(i);
    }
    printf("\rThermalizing 100%%\n");
    /*------------------------------------------------------------------------*/
    phase_plot(fFirst);
    for (i = 0; i < global.nParticles; i++) {
        if (particle[i].type == 0) {
            particle[i].wi = particle[i].w0;
            particle[i].xi = particle[i].x0;
            particle[i].yi = particle[i].y0;
        }
        particle[i].meanX = 0;
        particle[i].meanY = 0;
        particle[i].meanSX = 0;
        particle[i].meanSY = 0;
    }
    global.powerBottom = 0;
    global.powerVisc = 0;
    global.powerMu = 0;
    global.powerSpring = 0;
    global.potEnergyElasNorm = global.potEnergyElasTg = 0;
    global.potEnergyG = 0;
    global.kinEnergyTrans = global.kinEnergyRot =0;

    //Restart 3disk statistics to forget relaxation.
    links3.op_op = links3.op_sl = links3.op_cl = links3.sl_op = links3.sl_sl =
        links3.sl_cl = links3.cl_op = links3.cl_sl = links3.cl_cl = 0;
    links3.lstkount = 0;
    /*------------------------------------------------------------------------*/

#ifdef GRAPHICS
    relaxing = 0;
#endif

    for (i = 0; i < nStepsRun; i++) {
        if (!(i % stepsForPrint)) {
            printf("\rSimulating %.0lf%%", (double)i/nStepsRun*100);
            fflush(stdout);}
#ifdef GRAPHICS
        if (!(i % stepsForGraph)) {
            graphics(relaxing);
        }
#endif

#ifdef COLLISIONS
        search_collisions();
#endif

        if (!(i % stepsForWriteRun)) {
            write_results();
            if (i > 0) write_stats(stepsForWriteRun);
        }
        step(i);
    }
    printf("\rSimulating 100%%\n");

    phase_plot(fLast);
    if (global.nParticles == 3) {
        write_3disk_avglinkstat();
    }

#ifdef COLLISIONS
    write_collision_stats();
#endif

    free(particle); free(nStats);
    fclose(fFirst); fclose(fPhase); fclose(fPhase2); fclose(fLast);
    fclose(fLinkStat); fclose(fLinks); fclose(fEnergy);

#ifdef COLLISIONS
    fclose(fCollisions);
#endif

    clock_time((int)iTime);
    return 0;
}

/* long init_system() { */
/*     long nParticles = global.nParticles; */
/*     double box_w = global.box_w; */
/*     double box_h = global.box_h; */
/*     double meanR = diskParameters.meanR; */
/*     double density = diskParameters.density; */
/*     long seed = global.seed; */

/*     long i, j; */
/*     double layerY = 0; */
/*     int layerN = 1; */
/*     double sumWidth = 0; */
/*     double wallR = meanR; */
/*     long nBottom = (long)ceil(box_w/(2*wallR)); */
/*     double sigma = meanR*0.25; */
/*     disk tempParticle; */

/*     rgen = gsl_rng_alloc(gsl_rng_mt19937); */
/*     gsl_rng_set(rgen, seed); */

/*     assert(nParticles > 0); */
/*     assert(nBottom > 0); */

/*     particle = (disk*)calloc(nParticles + nBottom, sizeof(disk)); */

/*     /\*Setup the walls. Type = 1 for fixed disks*\/ */
/*     for (i = nParticles; i < nParticles + nBottom; i++) { */
/*         particle[i].type = 1; */
/*         particle[i].radius = wallR; */
/*         particle[i].x0 = sumWidth + particle[i].radius; */
/*         particle[i].y0 = particle[i].radius; */
/*         sumWidth += 2*particle[i].radius; */
/*         /\*Store initial positions for boundary conditions*\/ */
/*         particle[i].xi = particle[i].x0; */
/*         particle[i].yi = particle[i].y0; */
/*         /\*Set initial displacement of bottom*\/ */
/*         particle[i].y0 += global.epsilon*global.relInitDisp; */
/*         /\*Set disk properties*\/ */
/*         particle[i].mass = density*particle[i].radius*particle[i].radius*M_PI; */
/*         particle[i].iMoment = */
/*             particle[i].mass*particle[i].radius*particle[i].radius/2; */
/*     } */

/*     /\*Setup mobile particles. Type = 0 for mobile disk.*\/ */
/*     /\*Considering 3 different radius around a mean.*\/ */
/*     for (i = 0; i < nParticles; i++) { */
/*         particle[i].type = 0; */
/*         if (i < (long)(nParticles/3)) */
/*             particle[i].radius = meanR - sigma; */
/*         else if (i < (long)(2*nParticles/3)) */
/*             particle[i].radius = meanR + sigma; */
/*         else */
/*             particle[i].radius = meanR; */
/*         particle[i].mass = density*particle[i].radius*particle[i].radius*M_PI; */
/*         particle[i].iMoment = */
/*             particle[i].mass*particle[i].radius*particle[i].radius*0.5; */
/*         } */

/*         /\*Now make a random permutation (Knuth permutation).*\/ */
/*         for (i = 0; i < nParticles; i++) { */
/*             j = gsl_rng_uniform_int(rgen, i + 1); */
/*             tempParticle = particle[i]; */
/*             particle[i] = particle[j]; */
/*             particle[j] = tempParticle; */
/*         } */

/*         /\*Now set the x position such that contiguous particles barely touch.*\/ */
/*         sumWidth = 0 + (layerN % 2)*meanR; */
/*         layerY = 2*wallR + (meanR + sigma); */
/*         for (i = 0; i < nParticles; i++) { */
/*             if (sumWidth + 2*particle[i].radius > box_w) { */
/*                 layerY += 2*(meanR + sigma); */
/*                 layerN++; */
/*                 sumWidth = 0 + (layerN % 2)*meanR; */
/*             } */
/*             particle[i].x0 = sumWidth + particle[i].radius; */
/*             particle[i].y0 = layerY; */
/*             sumWidth += 2*particle[i].radius; */
/*         } */
/*         assert(layerY + (meanR + sigma)  < box_h); */

/*         gsl_rng_free(rgen); */

/*         nParticles += nBottom; */
/*         nStats = (neighbor_stats*)calloc((nParticles*nParticles-nParticles)/2 , */
/*                                          sizeof(neighbor_stats)); */
/*         for (i = 0; i < (nParticles*nParticles-nParticles)/2; i++) { */
/*             nStats[i].touching = 0; */
/*         } */

/*         /\*Restart RNG for simulation.*\/ */
/*         rgen = gsl_rng_alloc(gsl_rng_mt19937); */
/*         gsl_rng_set(rgen, global.seed); */

/*         return (nParticles); */
/* } */

void init_system() {
    /* Three radii. */
    long nDisks = global.nParticles;
    long nWalls = 4;
    double box_w = global.box_w;
    double box_h = global.box_h;
    global.box_w *= 1.1; global.box_h *= 1.1;
    global.nParticles += nWalls;
    double meanR = diskParameters.meanR;
    double sigma = meanR*0.25;
    double R1 = meanR - sigma, R2 = meanR, R3 = meanR + sigma;

    rgen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rgen, global.seed);

    /* Allocate memory for 4 walls and nDisks disks.*/
    particle = (disk*)calloc(global.nParticles, sizeof(disk));

    /* Set up box. Other geometries are possible.*/
    particle[nDisks].p1x = box_w*0.5; particle[nDisks].p1y = 0;
    particle[nDisks].nx = 0; particle[nDisks].ny = 1;
    particle[nDisks].x0 = box_w*0.5; particle[nDisks].y0 = 0;
    particle[nDisks+1].p1x = box_w; particle[nDisks+1].p1y = box_h*0.5;
    particle[nDisks+1].nx = -1; particle[nDisks+1].ny = 0;
    particle[nDisks+1].x0 = box_w; particle[nDisks+1].y0 = box_h*0.5;
    particle[nDisks+2].p1x = box_w*0.5; particle[nDisks+2].p1y = box_h;
    particle[nDisks+2].nx = 0;  particle[nDisks+2].ny = -1;
    particle[nDisks+2].x0 = box_w*0.5; particle[nDisks+2].y0 = box_h;
    particle[nDisks+3].p1x = 0; particle[nDisks+3].p1y = box_h*0.5;
    particle[nDisks+3].x0 = 0; particle[nDisks+3].y0 = box_h*0.5;
    particle[nDisks+3].nx = 1;  particle[nDisks+3].ny = 0;
    for (long i = nDisks; i < global.nParticles; i++) {
        particle[i].type = 2;
        /* Wall data: p1x p1y nx ny. */
        particle[i].p2x = particle[i].p1x + particle[i].ny;
        particle[i].p2y = particle[i].p1y - particle[i].nx;
        /* Store initial positions for boundary conditions */
        particle[i].p1x0 = particle[i].p1x;
        particle[i].p1y0 = particle[i].p1y;
        particle[i].p2x0 = particle[i].p2x;
        particle[i].p2y0 = particle[i].p2y;
        /* Set initial displacement of bottom for bottom wall. */
        if (particle[i].ny == 1) {
            particle[i].p1y += global.epsilon*global.relInitDisp;
            particle[i].p2y += global.epsilon*global.relInitDisp;
        }
    }

    /* Set up disk properties. Consider 3 different radii.*/
    for (long i = 0; i < nDisks; i++) {
        if       (i < (long)(nDisks/3))   particle[i].radius = R1;
        else if  (i < (long)(2*nDisks/3)) particle[i].radius = R2;
        else                              particle[i].radius = R3;

        particle[i].mass = M_PI*(particle[i].radius*particle[i].radius)*
            diskParameters.density;
        particle[i].iMoment = 0.5*particle[i].mass*
            (particle[i].radius*particle[i].radius);
        particle[i].type = 0;
    }
    /*Now make a random permutation (Knuth permutation).*/
    disk tempDisk;
    for (long i = 0; i < nDisks; i++) {
        long j = gsl_rng_uniform_int(rgen, i + 1);
        tempDisk = particle[i];
        particle[i] = particle[j];
        particle[j] = tempDisk;
    }

    nStats = (neighbor_stats*)calloc((global.nParticles*global.nParticles
                                      - global.nParticles)/2 ,
                                     sizeof(neighbor_stats));


    /*Restart RNG for simulation.*/
    rgen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rgen, global.seed);

    /* Setup a dense packing using an advancing front algorithm. */
    pack_disks();
}

void init_system_wedge() {
    double meanR = diskParameters.meanR;
    double density = diskParameters.density;
    long i;
    double wedgeAngle = M_PI/6;

    global.box_w = 3*2*diskParameters.meanR;
    global.box_h = 3*2*diskParameters.meanR;
    particle = (disk*)calloc(3, sizeof(disk));

    /*Setup the walls. Type = 2 for straight walls.*/
    for (i = 1; i < 3; i++) {
        particle[i].type = 2;
        double m = pow(-1, i)*tan(wedgeAngle);
        double b = -m*global.box_w/2;
        particle[i].p1x = 0;
        particle[i].p1y = b;
        particle[i].p2x = 1;
        particle[i].p2y = m+b;
        /*Store initial positions for boundary conditions*/
        particle[i].p1x0 = particle[i].p1x;
        particle[i].p1y0 = particle[i].p1y;
        particle[i].p2x0 = particle[i].p2x;
        particle[i].p2y0 = particle[i].p2y;
        /*Set initial displacement of bottom*/
        particle[i].p1y += global.epsilon*global.relInitDisp;
        particle[i].p2y += global.epsilon*global.relInitDisp;
    }

    /*Setup mobile particle. Type = 0 for mobile disk.*/
    particle[0].type = 0;
    particle[0].radius = meanR;
    particle[0].mass = density*particle[0].radius*particle[0].radius*M_PI;
    particle[0].iMoment =
        particle[0].mass*particle[0].radius*particle[0].radius/2;
    particle[0].x0 = global.box_w/2;
    particle[0].y0 = particle[0].radius*2;

    nStats = (neighbor_stats*)calloc(3, sizeof(neighbor_stats));
    for (i = 0; i < 3; i++) {
        nStats[i].touching = 0;
    }

    /*Restart RNG for simulation.*/
    rgen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rgen, global.seed);

    return;
}

void init_system_file(FILE *input_fptr) {
    /* Get input from file. */
    /* First line is nDisks nWalls box_w box_h*/
    /* File data as
       p1x p1y nx ny
       x y R mass iMoment */
    char input[500];
    long nDisks, nWalls;

    /* Read nDisks and nWalls. */
    if (fgets(input, sizeof(input), input_fptr) == NULL) exit(1);
    sscanf(input,"%ld %ld %lf %lf",
           &nDisks, &nWalls, &global.box_w, &global.box_h);
    global.nParticles = nDisks + nWalls;
    global.box_w *= 1.1; global.box_h *= 1.1;

    /* Allocate memory. */
    particle = (disk*)calloc(global.nParticles, sizeof(disk));

    /* Setup particles and walls. */
    for (long i = nDisks; i < global.nParticles; i++) {
        double p1x, p1y, nx, ny;
        particle[i].type = 2;
        /* Wall data: p1x p1y nx ny. */
        if (fgets(input, sizeof(input), input_fptr) == NULL) exit(1);
        sscanf(input,"%lf %lf %lf %lf", &p1x, &p1y, &nx, &ny);
        particle[i].p1x = p1x;
        particle[i].p1y = p1y;
        particle[i].p2x = p1x + ny;
        particle[i].p2y = p1y - nx;
        /* Store initial positions for boundary conditions */
        particle[i].p1x0 = particle[i].p1x;
        particle[i].p1y0 = particle[i].p1y;
        particle[i].p2x0 = particle[i].p2x;
        particle[i].p2y0 = particle[i].p2y;
        /* Set initial displacement of bottom for bottom wall. */
        if (ny == 1) {
            particle[i].p1y += global.epsilon*global.relInitDisp;
            particle[i].p2y += global.epsilon*global.relInitDisp;
        }
    }
    for (long i = 0; i < nDisks; i++) {
        /* Disk data: x y R mass iMoment */
        double x, y, R, mass, iMoment;
        if (fgets(input, sizeof(input), input_fptr) == NULL) exit(1);
        sscanf(input,"%lf %lf %lf %lf %lf", &x, &y, &R, &mass, &iMoment);
        particle[i].type = 0;
        particle[i].x0 = x;
        particle[i].y0 = y;
        particle[i].radius = R;
        particle[i].mass = mass;
        particle[i].iMoment = iMoment;
    }

    nStats = (neighbor_stats*)calloc((global.nParticles*global.nParticles
                                      - global.nParticles)/2 ,
                                     sizeof(neighbor_stats));


    /*Restart RNG for simulation.*/
    rgen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rgen, global.seed);
}

void step(long kstep) {
    long i = 0;
    long nParticles = global.nParticles;

    global.time += global.timestep;

    global.meanLinkSat = global.meanSqLinkSat = 0;
    global.linkCount = global.changingLinks = global.slidingLinks = 0;

    for (i = 0; i < nParticles; i++) {
        if (particle[i].type == 0) {
            if (particle[i].posFixed == 0) predictor_positions(&particle[i]);
            if (particle[i].rotFixed == 0) predictor_rotations(&particle[i]);
        }
        if (particle[i].type == 1) boundary_conditions(&particle[i],
                                                       global.bCondType, kstep);
        if (particle[i].type == 2) boundary_conditions_walls(i);
    }
    make_forces();
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

    /* Statistics */
    for (long i = 0; i < nParticles; i++) {
        particle[i].meanX += particle[i].x0;
        particle[i].meanY += particle[i].y0;
        particle[i].meanSX += particle[i].x0*particle[i].x0;
        particle[i].meanSY += particle[i].y0*particle[i].y0;
    }

    return;
}

void clock_time(int iTime) {
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
