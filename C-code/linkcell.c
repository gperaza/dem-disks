#include "common.h"

long nx, ny;
long **C; //screening array
long *E; //element array
double d; //size of cell

void init_cell() {
    long nParticles = global.nParticles;
    double box_w = global.box_w, box_h = global.box_h;
    long i, j;

    d = diskParameters.meanR*2;
    nx = (long)ceil(box_w/d);
    ny = (long)ceil(box_h/d);
    C = (long**) calloc(nx, sizeof(long*));
    for (i = 0; i < nx; i++) {
        C[i] = (long*)calloc(ny, sizeof(long));
    }
    E = (long*)calloc(nParticles, sizeof(long));

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            C[i][j] = -1;
        }
    }
    for (i = 0; i < nParticles; i++) {
        E[i] = -1;
    }
    return;
}

void free_cell() {
    long i;

    for (i = 0; i < nx; i++) {
        free(C[i]);
    }
    free(C);
    free(E);
}

void make_forces_linked_cell() {
    long i, j, k;
    long nParticles = global.nParticles;
    double gravityAngle = global.gravityAngle;
    double cga = cos(gravityAngle);
    double sga = sin(gravityAngle);
    long ix, iy;

    /*Create cells*/
    for (i = 0; i < nParticles; i++) {
        ix = (long)(particle[i].x0/d);
        iy = (long)(particle[i].y0/d);
        E[i] = C[ix][iy];
        C[ix][iy] = (long)i;
    }

    /*Set forces to zero*/
    for (i = 0; i < nParticles; i++) {
        particle[i].fx = 0; particle[i].fy = 0; particle[i].fw = 0;
    }

    /*Calculate forces*/
    for (i = 0; i < nParticles; i++) {

        /*Bulk forces*/
        if (particle[i].type == 0) {
            bulk_force(i, cga, sga);
        }

        /*Pair forces using linked cell*/
        ix = (long)(particle[i].x0/d);
        iy = (long)(particle[i].y0/d);
        if (C[ix][iy] < nParticles && C[ix][iy] != -1) {
            C[ix][iy] += nParticles;
            j = C[ix][iy] - nParticles;
            while(j != -1) {
                /*loop over own cell*/
                k = E[j];
                while (k != -1) {
                    pair_force(j, k);
                    k = E[k];
                }
                /*loop over neighbor cells*/
                k = C[(ix-1+nx)%nx][iy]%nParticles;
                while (k != -1) {
                    pair_force(j, k);
                    k = E[k];
                }
                k = C[(ix-1+nx)%nx][(iy-1+ny)%ny]%nParticles;
                while (k != -1) {
                    pair_force(j, k);
                    k = E[k];
                }
                k = C[ix][(iy-1+ny)%ny]%nParticles;
                while (k != -1) {
                    pair_force(j, k);
                    k = E[k];
                }
                k = C[(ix+1)%nx][(iy-1+ny)%ny]%nParticles;
                while (k != -1) {
                    pair_force(j, k);
                    k = E[k];
                }
                j = E[j];
            }
        }
    }

    /*Empty cells*/
    for (i = 0; i < nParticles; i++) {
        ix = (long)(particle[i].x0/d);
        iy = (long)(particle[i].y0/d);
        E[i] = -1;
        C[ix][iy] = -1;
    }

    return;
}
