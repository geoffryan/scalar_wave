/*
 *  scalar_wave.h
 *  
 *
 *  Created by Geoffrey Ryan on 13-05-29.
 *  Copyright 2013 None. All rights reserved.
 *
 */

#include "cow.h"

void (*timestep)(cow_dfield*, double, double, double, double);

void initialize(cow_dfield *phi, double a, double b, double width);
int run(cow_dfield *phi, double a, double b, double v, double cfl, double T);
void forward_euler(cow_dfield *phi, double v, double dt, double a, double b);
void rk2(cow_dfield *phi, double v, double dt, double a, double b);
void leap_frog(cow_dfield *phi, double v, double dt, double a, double b);
void write_phi(cow_dfield *phi, double time, char *fname);