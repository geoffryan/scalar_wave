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

void initialize(cow_dfield *phi, double a, double b, double width)
{
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getsize(dom, 0);
	
	int i;
	
	double *buff = (double *) cow_dfield_getdatabuffer(phi);
	
	double dx = (b-a)/(nx-1);
	
	for(i=0; i<nx; i++)
	{
		double x = i*dx + a;
		buff[2*(i+ng)] = sin(2 * (x-a) * M_PI / width);
		buff[2*(i+ng)+1] = 0;
	}
	
	for(i=0; i<ng; i++)
	{
		buff[2*i] = 0;
		buff[2*i+1] = 0;
		buff[2*(i+ng+nx)] = 0;
		buff[2*(i+ng+nx)+1] = 0;
	}
}

int run(cow_dfield *phi, double a, double b, double v, double cfl, double T)
{
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getsize(dom, 0);
	
	double dx = (b-a)/(nx-1);
	double dt = cfl * dx / v;
	
	double t = 0;
	int nt;
	double *data;
	
	printf("dx: %lg, dt: %lg\n", dx, dt);
	
	data = cow_dfield_getdatabuffer(phi);
	
	FILE *out = fopen("output_1d.txt", "w");
	fclose(out);
	write_phi(phi, t, "output_1d.txt");
	
	nt = 0;
	
	while(t < T)
	{
		timestep(phi, v, dt, a, b);
		data = cow_dfield_getdatabuffer(phi);
		
		t += dt;
		write_phi(phi, t, "output_1d.txt");
		nt++;
	}
	
	printf("nt: %d, T: %lg\n", nt, t);
	
	return nt;
}

void forward_euler(cow_dfield *phi, double v, double dt, double a, double b)
{
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getsize(dom, 0);
	
	double dx = (b-a)/(nx-1);
	
	int i;
	double *dphi1 = (double *) malloc(nx * sizeof(double));
	double *dphiv1 = (double *) malloc(nx * sizeof(double));

	dphi1[0] = 0;
	dphiv1[0] = 0;
	dphi1[nx-1] = 0;
	dphiv1[nx-1] = 0;
	
	for(i=1; i<nx-1; i++)
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
	int nx = cow_domain_getsize(dom, 0);
	
	double dx = (b-a)/(nx-1);
	
	int i;
	double *dphi1 = (double *) malloc(nx * sizeof(double));
	double *dphiv1 = (double *) malloc(nx * sizeof(double));
	double *dphi2 = (double *) malloc(nx * sizeof(double));
	double *dphiv2 = (double *) malloc(nx * sizeof(double));
		
	dphi1[0] = 0;
	dphiv1[0] = 0;
	dphi1[nx-1] = 0;
	dphiv1[nx-1] = 0;
	dphi2[0] = 0;
	dphiv2[0] = 0;
	dphi2[nx-1] = 0;
	dphiv2[nx-1] = 0;
	
	for(i=1; i<nx-1; i++)
	{
		dphi1[i] = dt * field[2*(i+ng)+1];
		dphiv1[i] = dt*v*v/(dx*dx) * (field[2*(i+ng+1)] - 2*field[2*(i+ng)] + field[2*(i+ng-1)]);
	}
	
	for(i=1; i<nx-1; i++)
	{
		dphi2[i] = dt * (field[2*(i+ng)+1] + 0.5*dphiv1[i]);
		dphiv2[i] = dt*v*v/(dx*dx) * (field[2*(i+ng+1)] - 2*field[2*(i+ng)] + field[2*(i+ng-1)]
									  + 0.5*(dphi1[i+1]-2*dphi1[i]+dphi1[i-1]));
	}
	
	for(i=0; i<nx; i++)
	{
		field[2*(i+ng)] += dphi2[i];
		field[2*(i+ng)+1] += dphiv2[i];
	}
	
	free(dphi1);
	free(dphiv1);
	free(dphi2);
	free(dphiv2);
}

void leap_frog(cow_dfield *phi, double v, double dt, double a, double b)
{
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getsize(dom, 0);
	
	double dx = (b-a)/(nx-1);
	
	int i;
	
	field[2*ng] = 0;
	field[2*ng+1] = 0;
	field[2*(ng+nx-1)] = 0;
	field[2*(ng+nx-1)+1] = 0;
	
	for(i=1; i<nx-1; i++)
	{
		field[2*(i+ng)] += dt * field[2*(i+ng)+1];
	}
	
	for(i=1; i<nx-1; i++)
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
	cow_domain_setsize(dom, 0, 100);
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
	double width = 2*(b-a);
	double cfl = 0.8;
	double T = 1.0;
	double v = 1.0;
	
	int nt;
	
	//Initialize
	initialize(phi, a, b, width);
	timestep = &leap_frog;
	//Run
	nt = run(phi, a, b, v, cfl, T);
	//Output
	
	cow_dfield_del(phi);
	cow_domain_del(dom);
	
	cow_finalize();
	
	return 0;
}

