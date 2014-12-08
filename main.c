/*
  File name: main.c

  This program simulates the time evolution of a set of disks in 2D

*/

#include "common.h"
#include <string.h>
#include <time.h>

disk *particle;
neighbor_stats *nStats;
particleParameters diskParameters;
systemParameters global;
FILE *fPhase, *fFirst, *fLast;
FILE *fp;

void get_input();
void clock_time(int);

void get_input() {
    char input[200];
    char type[200];
    char value[200];
    int inputs = 20, countInputs = 0;

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
            printf("box_w = %lf disks.\n", global.box_w);
        } else if (strcmp(type,"#box_h") == 0) {
            sscanf(value, "%lf", &global.box_h);
            printf("box_h = %lf disks.\n", global.box_h);
        } else if (strcmp(type,"#freq") == 0) {
            sscanf(value, "%lf", &global.freq);
            printf("Bottom frequency = %lf Hz.\n", global.freq);
        } else if (strcmp(type,"#dimensionlessAc") == 0) {
            sscanf(value, "%lf", &global.dimensionlessAc);
            printf("Bottom acceleration = %lf g's.\n", global.dimensionlessAc);
        } else if (strcmp(type,"#gravity") == 0) {
            sscanf(value, "%lf", &global.gravity);
            printf("Gravity = %lf m/s^2.\n", global.gravity);
        } else if (strcmp(type,"#gravityAngle") == 0) {
            sscanf(value, "%lf", &global.gravityAngle);
            printf("Gravity angle = %lf PI radians.\n", global.gravityAngle);
        } else if (strcmp(type,"#bGamma") == 0) {
            sscanf(value, "%lf", &global.bGamma);
            printf("Bulk gamma = %lf. \n", global.bGamma);
        } else if (strcmp(type,"#timestep") == 0) {
            sscanf(value, "%lf", &global.timestep);
            printf("Time step = %lf s. \n", global.timestep);
        } else if (strcmp(type,"#relaxTime") == 0) {
            sscanf(value, "%lf", &global.relaxTime);
            printf("Relaxing time = %lf s. \n", global.relaxTime);
        } else if (strcmp(type,"#runTime") == 0) {
            sscanf(value, "%lf", &global.runTime);
            printf("Simulation step = %lf s. \n", global.runTime);
        } else if (strcmp(type,"#timeForGraph") == 0) {
            sscanf(value, "%lf", &global.timeForGraph);
            printf("Time between graphics = %lf s. \n", global.timeForGraph);
        } else if (strcmp(type,"#timeForWrite") == 0) {
            sscanf(value, "%lf", &global.timeForWrite);
            printf("Time between writes = %lf s. \n", global.timeForWrite);
        } else if (strcmp(type,"#meanR") == 0) {
            sscanf(value, "%lf", &diskParameters.meanR);
            printf("Mean disk radius = %lf m. \n", diskParameters.meanR);
        } else if (strcmp(type,"#density") == 0) {
            sscanf(value, "%lf", &diskParameters.density);
            printf("Disk's density = %lf kg/m^3. \n", diskParameters.density);
        } else if (strcmp(type,"#kn") == 0) {
            sscanf(value, "%lf", &diskParameters.kn);
            printf("Normal elastic constant = %lf. \n", diskParameters.kn);
        } else if (strcmp(type,"#pr") == 0) {
            sscanf(value, "%lf", &diskParameters.pr);
            printf("Disk's Poisson's ratio = %lf. \n", diskParameters.pr);
        } else if (strcmp(type,"#mu") == 0) {
            sscanf(value, "%lf", &diskParameters.mu);
            printf("Disk's friction coefficient = %lf. \n", diskParameters.mu);
        } else if (strcmp(type,"#vGamma") == 0) {
            sscanf(value, "%lf", &diskParameters.vGamma);
            printf("Viscoelastic gamma = %lf. \n", diskParameters.vGamma);
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
    global.epsilon = global.dimensionlessAc*global.gravity/
        (4*M_PI*M_PI*global.freq*global.freq);

    return;
}

int main(/*int argc, char *argv[]*/)
{
    time_t iTime = time(NULL);
    get_input();

    double timestep = global.timestep;
    long nStepsRelax = (long)(global.relaxTime/timestep);
    long nStepsRun = (long)(global.runTime/timestep);
    long stepsForGraph = (long)(global.timeForGraph/timestep);
    long stepsForWrite = (long)(global.timeForWrite/timestep);
    global.time = -timestep*nStepsRelax;

    /*Set up constants for the gear integrator.*/
    set_constants(timestep);

    fPhase = fopen("phase_space.out", "w");
    fFirst = fopen("initial_phase_space.out", "w");
    fLast = fopen("ending_phase_space.out", "w");

    /*Initialize the packing and make nParticles include walls.*/
    printf("Simulating for %ld particles.\n", global.nParticles);
    global.nParticles = init_system();
    printf("Number of particles including walls: %ld.\n", global.nParticles);

    /*Initialize linked cells*/
    init_cell();

    int relaxing = 1;
    long i;
    for (i = 0; i < nStepsRelax; i++) {
        if (!(i % stepsForGraph)) {
            graphics(relaxing);
        }
        step();
     }

    //phase_plot(ffirst);

    relaxing = 0;
    for (i = 0; i < nStepsRun; i++) {
        if (!(i % stepsForGraph)) {
            graphics(relaxing);
        }
        if (!(i % stepsForWrite)) {
            //phase_plot(fphase);
        }
        step();
    }

    //phase_plot(flast);
    free(particle);
    free_cell();
    free(nStats);
    fclose(fFirst);
    fclose(fPhase);
    fclose(fLast);
    clock_time((int)iTime);
    return 0;
}

long init_system() {
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
    gsl_rng * rgen = gsl_rng_alloc(gsl_rng_mt19937);
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
            particle[i].radius = meanR;
        else
            particle[i].radius = meanR + sigma;
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

    return (nParticles);
}

void step() {

    long i = 0;
    long nParticles = global.nParticles;

    global.time += global.timestep;

    for (i = 0; i < nParticles; i++) {
        if (particle[i].type == 0) {
            if (particle[i].posFixed == 0) predictor_positions(&particle[i]);
            if (particle[i].rotFixed == 0) predictor_rotations(&particle[i]);
        }
        if (particle[i].type == 1) boundary_conditions(&particle[i]);
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

    return;
}

void phase_plot(FILE* fp) {
    long i;

    for (i = 0; i < global.nParticles; i++) {
        fprintf(fp,"Hello World.\n");
    }
    fflush(fp);
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
