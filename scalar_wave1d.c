/*
 *  scalar_wave.c
 *  
 *
 *  Created by Geoffrey Ryan on 13-05-29.
 *  Copyright 2013 None. All rights reserved.
 *
 */

#include "scalar_wave1d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void initial_sine(cow_dfield *phi, double a, double b, double v, double wavelength, double phase)
{
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	
	int i;
	double x;
	
	for(i=0; i<nx; i++)
	{
		x = (b - a) * cow_domain_positionatindex(dom, 0, i+ng) + a;
		field[2*(i+ng)] = sin(2 * (x-a) * M_PI / wavelength + phase);
		field[2*(i+ng)+1] = -(2 * M_PI * v) / wavelength * cos(2 * (x-a) * M_PI / wavelength + phase);
	}
	
	cow_dfield_syncguard(phi);
}

void initial_sine_stand(cow_dfield *phi, double a, double b, double v, double wavelength, double phase)
{
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	
	int i;
	double x;
	
	for(i=0; i<nx; i++)
	{
		x = (b - a) * cow_domain_positionatindex(dom, 0, i+ng) + a;
		field[2*(i+ng)] = sin(2 * (x-a) * M_PI / wavelength + phase);
		field[2*(i+ng)+1] = 0;
	}
	
	cow_dfield_syncguard(phi);
}

void initial_gauss(cow_dfield *phi, double a, double b, double v, double width, double mean)
{
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	
	int i;
	double x;
	
	for(i=0; i<nx; i++)
	{
		x = (b - a) * cow_domain_positionatindex(dom, 0, i+ng) + a;
		field[2*(i+ng)] = exp(-(x-mean)*(x-mean) / (2.0*width*width));
		field[2*(i+ng)+1] = 0;
	}
	
	cow_dfield_syncguard(phi);
}

int run(cow_dfield *phi, double a, double b, double v, double cfl, double T)
{
	cow_domain *dom = cow_dfield_getdomain(phi);	
	
	double dx = (b-a) * cow_domain_getgridspacing(dom, 0);
	double dt = cfl * dx / v;
	
	double t = 0;
	int nt = 0;
	double *data;
	
	char basename[] = "out/output_1d";
	char name[50];
	
	printf("dx: %lg, dt: %lg\n", dx, dt);
	
	data = cow_dfield_getdatabuffer(phi);
	
	sprintf(name, "%s_%d.h5", basename, nt);
	cow_dfield_write(phi, name);
	
	while(t < T)
	{
		timestep(phi, v, dt, a, b);
		data = cow_dfield_getdatabuffer(phi);

		t += dt;
		nt++;
		
		sprintf(name, "%s_%d.h5", basename, nt);
		cow_dfield_write(phi, name);		
	}
	
	printf("nt: %d, T: %lg\n", nt, t);
	
	return nt;
}

void forward_euler(cow_dfield *phi, double v, double dt, double a, double b)
{
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	
	double dx = (b-a) * cow_domain_getgridspacing(dom, 0);
	
	int i;
	double *dphi1 = (double *) malloc(nx * sizeof(double));
	double *dphiv1 = (double *) malloc(nx * sizeof(double));
	
	cow_dfield_syncguard(phi);
	
	for(i=0; i<nx; i++)
	{
		dphi1[i] = dt * field[2*(i+ng)+1];
		dphiv1[i] = dt*v*v/(dx*dx) * (field[2*(i+ng+1)] - 2*field[2*(i+ng)] + field[2*(i+ng-1)]);
	}
	
	for(i=0; i<nx; i++)
	{
		field[2*(i+ng)] += dphi1[i];
		field[2*(i+ng)+1] += dphiv1[i];
	}
	
	free(dphi1);
	free(dphiv1);
}

void rk2(cow_dfield *phi, double v, double dt, double a, double b)
{
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	
	double dx = (b-a) * cow_domain_getgridspacing(dom, 0);
	
	int i;
	
	cow_dfield *dphi1 = cow_dfield_new();
	cow_dfield_setdomain(dphi1, dom);
	cow_dfield_setname(dphi1, "phi");
	cow_dfield_addmember(dphi1, "phi");
	cow_dfield_addmember(dphi1, "phiv");
	cow_dfield_commit(dphi1);
	double *field1 = (double *) cow_dfield_getdatabuffer(dphi1);
	
	cow_dfield *dphi2 = cow_dfield_new();
	cow_dfield_setdomain(dphi2, dom);
	cow_dfield_setname(dphi2, "phi");
	cow_dfield_addmember(dphi2, "phi");
	cow_dfield_addmember(dphi2, "phiv");
	cow_dfield_commit(dphi2);
	double *field2 = (double *) cow_dfield_getdatabuffer(dphi2);
	
	cow_dfield_syncguard(phi); // sets periodic BC's
	
	for(i=0; i<nx; i++)
	{
		field1[2*(i+ng)] = dt * field[2*(i+ng)+1];
		field1[2*(i+ng)+1] = dt*v*v/(dx*dx) * (field[2*(i+ng+1)] - 2*field[2*(i+ng)] + field[2*(i+ng-1)]);
	}
	
	cow_dfield_syncguard(dphi1);
	
	for(i=0; i<nx; i++)
	{
		field2[2*(i+ng)] = dt * (field[2*(i+ng)+1] + 0.5*field1[2*(i+ng)+1]);
		field2[2*(i+ng)+1] = dt*v*v/(dx*dx) * (field[2*(i+ng+1)] - 2*field[2*(i+ng)] + field[2*(i+ng-1)]
								+ 0.5*(field1[2*(i+ng+1)] - 2*field1[2*(i+ng)] + field1[2*(i+ng-1)]));
	}
	
	for(i=0; i<nx; i++)
	{
		field[2*(i+ng)] += field2[2*(i+ng)];
		field[2*(i+ng)+1] += field2[2*(i+ng)+1];
	}
	
	cow_dfield_del(dphi1);
	cow_dfield_del(dphi2);
}

void rk4(cow_dfield *phi, double v, double dt, double a, double b)
{
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	
	double dx = (b-a) * cow_domain_getgridspacing(dom, 0);
	
	int i;
	
	cow_dfield *dphi1 = cow_dfield_new();
	cow_dfield_setdomain(dphi1, dom);
	cow_dfield_setname(dphi1, "phi");
	cow_dfield_addmember(dphi1, "phi");
	cow_dfield_addmember(dphi1, "phiv");
	cow_dfield_commit(dphi1);
	double *field1 = (double *) cow_dfield_getdatabuffer(dphi1);
	
	cow_dfield *dphi2 = cow_dfield_new();
	cow_dfield_setdomain(dphi2, dom);
	cow_dfield_setname(dphi2, "phi");
	cow_dfield_addmember(dphi2, "phi");
	cow_dfield_addmember(dphi2, "phiv");
	cow_dfield_commit(dphi2);
	double *field2 = (double *) cow_dfield_getdatabuffer(dphi2);
	
	cow_dfield *dphi3 = cow_dfield_new();
	cow_dfield_setdomain(dphi3, dom);
	cow_dfield_setname(dphi3, "phi");
	cow_dfield_addmember(dphi3, "phi");
	cow_dfield_addmember(dphi3, "phiv");
	cow_dfield_commit(dphi3);
	double *field3 = (double *) cow_dfield_getdatabuffer(dphi3);
	
	cow_dfield *dphi4 = cow_dfield_new();
	cow_dfield_setdomain(dphi4, dom);
	cow_dfield_setname(dphi4, "phi");
	cow_dfield_addmember(dphi4, "phi");
	cow_dfield_addmember(dphi4, "phiv");
	cow_dfield_commit(dphi4);
	double *field4 = (double *) cow_dfield_getdatabuffer(dphi4);
	
	cow_dfield_syncguard(phi); // sets periodic BC's
	
	for(i=0; i<nx; i++)
	{
		field1[2*(i+ng)] = dt * field[2*(i+ng)+1];
		field1[2*(i+ng)+1] = dt*v*v/(dx*dx) * (field[2*(i+ng+1)] - 2*field[2*(i+ng)] + field[2*(i+ng-1)]);
	}
	cow_dfield_syncguard(dphi1);
	
	for(i=0; i<nx; i++)
	{
		field2[2*(i+ng)] = dt * (field[2*(i+ng)+1] + 0.5*field1[2*(i+ng)+1]);
		field2[2*(i+ng)+1] = dt*v*v/(dx*dx) * (field[2*(i+ng+1)] - 2*field[2*(i+ng)] + field[2*(i+ng-1)]
											   + 0.5*(field1[2*(i+ng+1)] - 2*field1[2*(i+ng)] + field1[2*(i+ng-1)]));
	}
	cow_dfield_syncguard(dphi2);
	
	for(i=0; i<nx; i++)
	{
		field3[2*(i+ng)] = dt * (field[2*(i+ng)+1] + 0.5*field2[2*(i+ng)+1]);
		field3[2*(i+ng)+1] = dt*v*v/(dx*dx) * (field[2*(i+ng+1)] - 2*field[2*(i+ng)] + field[2*(i+ng-1)]
											   + 0.5*(field2[2*(i+ng+1)] - 2*field2[2*(i+ng)] + field2[2*(i+ng-1)]));
	}
	cow_dfield_syncguard(dphi3);
	
	for(i=0; i<nx; i++)
	{
		field4[2*(i+ng)] = dt * (field[2*(i+ng)+1] + field3[2*(i+ng)+1]);
		field4[2*(i+ng)+1] = dt*v*v/(dx*dx) * (field[2*(i+ng+1)] - 2*field[2*(i+ng)] + field[2*(i+ng-1)]
											   + field3[2*(i+ng+1)] - 2*field3[2*(i+ng)] + field3[2*(i+ng-1)]);
	}
	
	for(i=0; i<nx; i++)
	{
		field[2*(i+ng)] += (field1[2*(i+ng)] + 2*field2[2*(i+ng)] + 2*field3[2*(i+ng)] + field4[2*(i+ng)]) / 6.0;
		field[2*(i+ng)+1] += (field1[2*(i+ng)+1] + 2*field2[2*(i+ng)+1] + 2*field3[2*(i+ng)+1] + field4[2*(i+ng)+1]) / 6.0;
	}
	
	cow_dfield_del(dphi1);
	cow_dfield_del(dphi2);
	cow_dfield_del(dphi3);
	cow_dfield_del(dphi4);
}

void leap_frog(cow_dfield *phi, double v, double dt, double a, double b)
{
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	
	double dx = (b-a) * cow_domain_getgridspacing(dom, 0);
	
	int i;
	
	cow_dfield_syncguard(phi); 
	
	for(i=0; i<nx; i++)
	{
		field[2*(i+ng)] += dt * field[2*(i+ng)+1];
	}

	cow_dfield_syncguard(phi);
	
	for(i=0; i<nx; i++)
	{
		field[2*(i+ng)+1] += dt*v*v/(dx*dx) * (field[2*(i+ng+1)] - 2*field[2*(i+ng)] + field[2*(i+ng-1)]);
	}
}

void write_phi(cow_dfield *phi, double time, char *fname)
{
	double *data = (double *) cow_dfield_getdatabuffer(phi);
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getsize(dom, 0);
	int i;
	
	FILE *out = fopen(fname, "a");
	fprintf(out, "%lg", time);
	for(i=0; i<nx; i++)
		fprintf(out, " %lg %lg", data[2*(ng+i)], data[2*(ng+i)+1]);
	fprintf(out, "\n");
	fclose(out);
}

int main(int argc, char *argv[])
{	
	cow_init(argc, argv, 0);
	
	cow_domain *dom = cow_domain_new();
	cow_domain_setndim(dom, 1);
	cow_domain_setsize(dom, 0, 200);
	cow_domain_setguard(dom, 2);
	cow_domain_commit(dom);
	
	cow_dfield *phi = cow_dfield_new();
	cow_dfield_setdomain(phi, dom);
	cow_dfield_setname(phi, "phi");
	cow_dfield_addmember(phi, "phi");
	cow_dfield_addmember(phi, "phiv");
	cow_dfield_commit(phi);
	
	MPI_Comm comm;
	cow_domain_getcomm(dom, &comm);
	int r;
	MPI_Comm_rank(comm, &r);
	printf("Hello from %d!\n", r);
	
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getsize(dom, 0);
	
	printf("ng: %d, nx: %d\n", ng, nx);
	
	double a = -1;
	double b = 1;
	double cfl = 0.99;
	double T = 20.1;
	double v = 1.0;
	
	double lambda = (b - a) / 5;
	double phase = 0;
	initialize = &initial_sine;
	timestep = &rk4;
	
//	double lambda = 2*(b-a) / 2.0;
//	double phase = 0;
//	initialize = &initial_sine_stand;
//	timestep = &forward_euler;
	
//	double lambda = (b-a) / 10.0;
//	double phase = 0.5*(a+b);
//	initialize = &initial_gauss;
//	timestep = &leap_frog;
	
	int nt;
	
	//Initialize
	initialize(phi, a, b, v, lambda, phase);
	//Run
	nt = run(phi, a, b, v, cfl, T);
	//Output
	
	cow_dfield_del(phi);
	cow_domain_del(dom);
	
	cow_finalize();
	
	return 0;
}

