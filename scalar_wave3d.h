/*
 *  scalar_wave.h
 *  
 *
 *  Created by Geoffrey Ryan on 13-05-29.
 *  Copyright 2013 None. All rights reserved.
 *
 */

#include "cow.h"

void (*initialize)(cow_dfield*, double, double*);
void (*timestep)(cow_dfield*, double, double);

double min3(double a, double b, double c);
void initial_sine(cow_dfield *phi, double v, double *params);
void initial_sine_stand(cow_dfield *phi, double v, double *params);
void initial_gauss(cow_dfield *phi, double v, double *params);
void initial_impulse(cow_dfield *phi, double v, double *params);
void initial_sphericalwave(cow_dfield *phi, double v, double *params);
void initial_spherical_pulse(cow_dfield *phi, double v, double *params);
void initial_test(cow_dfield *phi, double v, double *params);
int run(cow_dfield *phi, double v, double cfl, double T);
void forward_euler(cow_dfield *phi, double v, double dt);
void rk2(cow_dfield *phi, double v, double dt);
void rk4(cow_dfield *phi, double v, double dt);
void leap_frog(cow_dfield *phi, double v, double dt);