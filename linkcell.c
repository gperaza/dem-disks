#include "common.h"

int nx, ny;
long unsigned *Y;
int *B;

void make_link_cell();

void init_algorithm() {
    long unsigned nParticles = global.nParticles;
    double box_w = global.box_w, box_h = global.box_h;

    double d = diskParameters.meanR*2;
    nx = ceil(box_w/d);
    ny = ceil(box_h/d);
    Y =(long unsigned*)calloc(nParticles, sizeof(long unsigned));
    B = (int*)calloc(ny, sizeof(int));


}

void make_link_cell()
{
    unsigned long nParticles = global.nParticles;
    long ix, iy;
    unsigned long i;
    double d = diskParameters.meanR*2;

    for (i = 0; i < nParticles; i++) {
        iy = (int)(particle[i].y0/d);

    }
}

double* allocate_2d_array(double **arr, unsigned long  n, unsigned long m)
{
    unsigned long i = 0;
    arr = (double**)malloc(n * sizeof(double*));
    double *arr_data = malloc( n * m * sizeof(double));
    for (i = 0; i < n; i++)
        (arr)[i] = arr_data + i * m ;
    return arr_data; //free point
}
void deallocate_2d_array(double** arr, double* arr_data){
    free(arr_data);
    free(arr);
}
