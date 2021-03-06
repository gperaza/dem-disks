#ifndef _common_h
#define _common_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "particle_structures.h"

extern disk *particle;
extern particleParameters diskParameters;
extern systemParameters global;
extern neighbor_stats *nStats;

void step(long);
void predictor_positions(disk*);
void predictor_rotations(disk*);
void corrector_positions(disk*);
void corrector_rotations(disk*);
void set_constants(double);
void make_forces();
void graphics(int);
void boundary_conditions(disk*, int, long);
void boundary_conditions_walls(long);
void bulk_force(long, double, double);
void pair_force(long, long);
#endif
