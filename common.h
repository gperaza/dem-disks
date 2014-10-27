#ifndef _common_h
#define _common_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include "particle_structures.h"

extern disk *particle;
extern particleParameters diskParameters;
extern systemParameters global;

void step();
void integrate();
void predictor_positions(disk*);
void predictor_rotations(disk*);
void corrector_positions(disk*);
void corrector_rotations(disk*);
unsigned long init_system();
void set_constants(double);
void make_forces();
void graphics(int);
void boundary_conditions(disk*);
#endif
