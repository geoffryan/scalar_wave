/*
 *  scalar_wave.c
 *  
 *
 *  Created by Geoffrey Ryan on 13-05-29.
 *  Copyright 2013 None. All rights reserved.
 *
 */

#include "scalar_wave3d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

double min3(double a, double b, double c)
{
	if(a < b)
	{
		if(a < c)
			return a;
		else
			return c;
	}
	else
	{
		if(b < c)
			return b;
		else
			return c;
	}
}

void initial_sine(cow_dfield *phi, double v, double *params)
{
	double wavelength = params[0];
	double phase = params[1];
	
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	int sx = cow_dfield_getstride(phi, 0);
	int sy = cow_dfield_getstride(phi, 1);
	int sz = cow_dfield_getstride(phi, 2);
	
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	
	int i, j, k;
	double x, y, z;
	double kx, ky, kz, w;
	
	kx = 2 * M_PI / wavelength;
	ky = 0;
	kz = 0;
	w = v * sqrt(kx*kx+ky*ky+kz*kz);
	
	for(i = ng; i < nx+ng; i++)
	{	
		x = cow_domain_positionatindex(dom, 0, i);
		for(j = ng; j < ny+ng; j++)
		{
			y = cow_domain_positionatindex(dom, 1, j);
			for(k = ng; k < nz+ng; k++)
			{
				z = cow_domain_positionatindex(dom, 2, k);
				field[sx*i+sy*j+sz*k] = sin(kx*x+ky*y+kz*z + phase);
				field[sx*i+sy*j+sz*k+1] = -w * cos(kx*x+ky*y+kz*z + phase);
			}
		}
	}
	
	cow_dfield_syncguard(phi);
}

void initial_sine_stand(cow_dfield *phi, double v, double *params)
{
	double wavelength = params[0];
	double phase = params[1];
	
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	int sx = cow_dfield_getstride(phi, 0);
	int sy = cow_dfield_getstride(phi, 1);
	int sz = cow_dfield_getstride(phi, 2);
	
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	
	int i, j, k;
	double x, y, z;
	double kx, ky, kz;
	
	kx = 2 * M_PI / wavelength;
	ky = 2 * M_PI / wavelength;
	kz = 2 * M_PI / wavelength;
	
	for(i = ng; i < nx+ng; i++)
	{	
		x = cow_domain_positionatindex(dom, 0, i);
		for(j = ng; j < ny+ng; j++)
		{
			y = cow_domain_positionatindex(dom, 1, j);
			for(k = ng; k < nz+ng; k++)
			{
				z = cow_domain_positionatindex(dom, 2, k);
				field[sx*i+sy*j+sz*k] = sin(kx*x+phase)*sin(ky*y+phase)*sin(kz*z+phase);
				field[sx*i+sy*j+sz*k+1] = 0;
			}
		}
	}
	
	cow_dfield_syncguard(phi);
}

void initial_gauss(cow_dfield *phi, double v, double *params)
{
	double width = params[0];
	double mean = params[1];
	
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	int sx = cow_dfield_getstride(phi, 0);
	int sy = cow_dfield_getstride(phi, 1);
	int sz = cow_dfield_getstride(phi, 2);
	
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	
	int i, j, k;
	double x, y, z;
	double inv_var = 1.0/(2.0*width*width);
	
	for(i = ng; i < nx+ng; i++)
	{	
		x = cow_domain_positionatindex(dom, 0, i);
		for(j = ng; j < ny+ng; j++)
		{
			y = cow_domain_positionatindex(dom, 1, j);
			for(k = ng; k < nz+ng; k++)
			{
				z = cow_domain_positionatindex(dom, 2, k);
				field[sx*i+sy*j+sz*k] = exp(-((x-mean)*(x-mean)+y*y+z*z)*inv_var);
				field[sx*i+sy*j+sz*k+1] = 0;
			}
		}
	}
	
	cow_dfield_syncguard(phi);
}

void initial_impulse(cow_dfield *phi, double v, double *params)
{
	double radius = params[0];
	double height = params[1];
	
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	int sx = cow_dfield_getstride(phi, 0);
	int sy = cow_dfield_getstride(phi, 1);
	int sz = cow_dfield_getstride(phi, 2);
	
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	
	int i, j, k;
	double x, y, z;
	
	for(i = ng; i < nx+ng; i++)
	{	
		x = cow_domain_positionatindex(dom, 0, i);
		for(j = ng; j < ny+ng; j++)
		{
			y = cow_domain_positionatindex(dom, 1, j);
			for(k = ng; k < nz+ng; k++)
			{
				z = cow_domain_positionatindex(dom, 2, k);
				if(x*x+y*y+z*z <= radius*radius)
					field[sx*i+sy*j+sz*k] = height;
				else
					field[sx*i+sy*j+sz*k] = 0;
				field[sx*i+sy*j+sz*k+1] = 0;
			}
		}
	}
	
	cow_dfield_syncguard(phi);
}

void initial_sphericalwave(cow_dfield *phi, double v, double *params)
{
	double wavelength = params[0];
	double scale = params[1];
	
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	int sx = cow_dfield_getstride(phi, 0);
	int sy = cow_dfield_getstride(phi, 1);
	int sz = cow_dfield_getstride(phi, 2);
	
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	
	int i, j, k;
	double x, y, z, r;
	double wavenum = 2 * M_PI / wavelength;
	
	for(i = ng; i < nx+ng; i++)
	{	
		x = cow_domain_positionatindex(dom, 0, i);
		for(j = ng; j < ny+ng; j++)
		{
			y = cow_domain_positionatindex(dom, 1, j);
			for(k = ng; k < nz+ng; k++)
			{
				z = cow_domain_positionatindex(dom, 2, k);
				r = sqrt(x*x + y*y + z*z);
				field[sx*i+sy*j+sz*k] = scale * sin(wavenum * r) / r;
				field[sx*i+sy*j+sz*k+1] = -scale * v * wavenum * cos(wavenum * r) / r;
			}
		}
	}
	
	cow_dfield_syncguard(phi);
}

void initial_spherical_pulse(cow_dfield *phi, double v, double *params)
{
	double width = params[0];
	double scale = params[1];
	double r0 = params[2];
	
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	int sx = cow_dfield_getstride(phi, 0);
	int sy = cow_dfield_getstride(phi, 1);
	int sz = cow_dfield_getstride(phi, 2);
	
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	
	int i, j, k;
	double x, y, z, r;
	double inv_var = 1.0/(2.0*width*width);
	
	for(i = ng; i < nx+ng; i++)
	{	
		x = cow_domain_positionatindex(dom, 0, i);
		for(j = ng; j < ny+ng; j++)
		{
			y = cow_domain_positionatindex(dom, 1, j);
			for(k = ng; k < nz+ng; k++)
			{
				z = cow_domain_positionatindex(dom, 2, k);
				r = sqrt(x*x+y*y+z*z);
				
				field[sx*i+sy*j+sz*k] = scale/r * (exp(-inv_var*(r-r0)*(r-r0)) - exp(-inv_var*(-r-r0)*(-r-r0)));
				field[sx*i+sy*j+sz*k+1] = 2.0*scale*v*inv_var/r * ((r-r0)*exp(-inv_var*(r-r0)*(r-r0)) - (-r-r0)*exp(-inv_var*(-r-r0)*(-r-r0)));
			}
		}
	}
	
	cow_dfield_syncguard(phi);
}

void initial_test(cow_dfield *phi, double v, double *params)
{
	double width = params[0];
	
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	int sx = cow_dfield_getstride(phi, 0);
	int sy = cow_dfield_getstride(phi, 1);
	int sz = cow_dfield_getstride(phi, 2);
	
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	
	int i, j, k;
	double x, y, z;
	double inv_var = 1.0/(2.0*width*width);
	
	for(i = ng; i < nx+ng; i++)
	{	
		x = cow_domain_positionatindex(dom, 0, i);
		for(j = ng; j < ny+ng; j++)
		{
			y = cow_domain_positionatindex(dom, 1, j);
			for(k = ng; k < nz+ng; k++)
			{
				z = cow_domain_positionatindex(dom, 2, k);
				
				field[sx*i+sy*j+sz*k] = exp(-((x-1.0)*(x-1.0)+(y-1.0)*(y-1.0)+(z-1.0)*(z-1.0))*inv_var)*sin(2*M_PI*x*4.0)*z;
				field[sx*i+sy*j+sz*k+1] = 0;
			}
		}
	}
	
	cow_dfield_syncguard(phi);
}

int run(cow_dfield *phi, double v, double cfl, double T)
{
	cow_domain *dom = cow_dfield_getdomain(phi);	
	int size = cow_domain_getsize(dom, 0);
	int cartrank = cow_domain_getcartrank(dom);
	
	double ax = cow_domain_getlowercoord(dom, 0);
	double bx = cow_domain_getuppercoord(dom, 0);
	double ay = cow_domain_getlowercoord(dom, 1);
	double by = cow_domain_getuppercoord(dom, 1);
	double az = cow_domain_getlowercoord(dom, 2);
	double bz = cow_domain_getuppercoord(dom, 2);
	
	double dx = cow_domain_getgridspacing(dom, 0);
	double dy = cow_domain_getgridspacing(dom, 1);
	double dz = cow_domain_getgridspacing(dom, 2);
	double dt = cfl / v * min3(dx, dy, dz);
	
	double t = 0;
	int nt = 0;
	
	char basename[] = "out/output_3d";
	char name[50];
	char sumname[50];
	
	sprintf(name, "%s_%d_%lg_%d.h5", basename, size, t, nt);
	sprintf(sumname, "%s_%d_summary.txt", basename, size);

	if(cartrank == 0)
	{
		FILE *sumfile = fopen(sumname, "w");
		fprintf(sumfile, "size: %d\n", size);
		fprintf(sumfile, "a: %lg %lg %lg\n", ax, ay, az);
		fprintf(sumfile, "b: %lg %lg %lg\n", bx, by, bz);
		fprintf(sumfile, "v: %lg\n", v);
		fprintf(sumfile, "cfl: %lg\n", cfl);
		fprintf(sumfile, "dx: %lg %lg %lg\n", dx, dy, dz);
		fprintf(sumfile, "dt: %lg\n", dt);
		fprintf(sumfile, "T: %lg\n", T);
		fclose(sumfile);
	}
	
	printf("dx: %lg, dy: %lg, dz: %lg, dt: %lg\n", dx, dy, dz, dt);
	
	cow_dfield_write(phi, name);
	
	while(t < T)
	{
		timestep(phi, v, dt);

		t += dt;
		nt++;
		
		sprintf(name, "%s_%d_%lg_%d.h5", basename, size, t, nt);
		cow_dfield_write(phi, name);		
	}
	
	printf("nt: %d, T: %lg\n", nt, t);
	
	return nt;
}

void forward_euler(cow_dfield *phi, double v, double dt)
{
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	int sx = cow_dfield_getstride(phi, 0);
	int sy = cow_dfield_getstride(phi, 1);
	int sz = cow_dfield_getstride(phi, 2);
	
	double dx = cow_domain_getgridspacing(dom, 0);
	double dy = cow_domain_getgridspacing(dom, 1);
	double dz = cow_domain_getgridspacing(dom, 2);
	
	double idx = dt*v*v / dx*dx;
	double idy = dt*v*v / dy*dy;
	double idz = dt*v*v / dz*dz;
	
	int i,j,k;
	
	cow_dfield *dphi1 = cow_dfield_new();
	cow_dfield_setdomain(dphi1, dom);
	cow_dfield_setname(dphi1, "phi");
	cow_dfield_addmember(dphi1, "phi");
	cow_dfield_addmember(dphi1, "phiv");
	cow_dfield_commit(dphi1);
	double *field1 = (double *) cow_dfield_getdatabuffer(dphi1);
	
	cow_dfield_syncguard(phi);
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)
			{
				field1[sx*i+sy*j+sz*k] = dt * field[sx*i+sy*j+sz*k+1];
				field1[sx*i+sy*j+sz*k+1] = idx*(field[sx*(i+1)+sy*j+sz*k] + field[sx*(i-1)+sy*j+sz*k] - 2*field[sx*i+sy*j+sz*k]);
				field1[sx*i+sy*j+sz*k+1] += idy*(field[sx*i+sy*(j+1)+sz*k] + field[sx*i+sy*(j-1)+sz*k] - 2*field[sx*i+sy*j+sz*k]);
				field1[sx*i+sy*j+sz*k+1] += idz*(field[sx*i+sy*j+sz*(k+1)] + field[sx*i+sy*j+sz*(k-1)] - 2*field[sx*i+sy*j+sz*k]);
			}
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)
			{
				field[sx*i+sy*j+sz*k] += field1[sx*i+sy*j+sz*k];
				field[sx*i+sy*j+sz*k+1] += field1[sx*i+sy*j+sz*k+1];
			}
	
	cow_dfield_del(dphi1);
}

void rk2(cow_dfield *phi, double v, double dt)
{
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	int sx = cow_dfield_getstride(phi, 0);
	int sy = cow_dfield_getstride(phi, 1);
	int sz = cow_dfield_getstride(phi, 2);
	
	double dx = cow_domain_getgridspacing(dom, 0);
	double dy = cow_domain_getgridspacing(dom, 1);
	double dz = cow_domain_getgridspacing(dom, 2);
	
	double idx = dt*v*v / dx*dx;
	double idy = dt*v*v / dy*dy;
	double idz = dt*v*v / dz*dz;
	
	int i,j,k;
	
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
	
	cow_dfield_syncguard(phi);
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)
			{
				field1[sx*i+sy*j+sz*k] = dt * field[sx*i+sy*j+sz*k+1];
				field1[sx*i+sy*j+sz*k+1] = idx*(field[sx*(i+1)+sy*j+sz*k] + field[sx*(i-1)+sy*j+sz*k] - 2*field[sx*i+sy*j+sz*k]);
				field1[sx*i+sy*j+sz*k+1] += idy*(field[sx*i+sy*(j+1)+sz*k] + field[sx*i+sy*(j-1)+sz*k] - 2*field[sx*i+sy*j+sz*k]);
				field1[sx*i+sy*j+sz*k+1] += idz*(field[sx*i+sy*j+sz*(k+1)] + field[sx*i+sy*j+sz*(k-1)] - 2*field[sx*i+sy*j+sz*k]);
			}
	
	cow_dfield_syncguard(dphi1);
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)
			{
				field2[sx*i+sy*j+sz*k] = dt * (field[sx*i+sy*j+sz*k+1]+0.5*field1[sx*i+sy*j+sz*k+1]);
				field2[sx*i+sy*j+sz*k+1] = idx*(field[sx*(i+1)+sy*j+sz*k] + field[sx*(i-1)+sy*j+sz*k] - 2*field[sx*i+sy*j+sz*k]
												+ 0.5*(field1[sx*(i+1)+sy*j+sz*k] + field1[sx*(i-1)+sy*j+sz*k] - 2*field1[sx*i+sy*j+sz*k]));
				field2[sx*i+sy*j+sz*k+1] += idy*(field[sx*i+sy*(j+1)+sz*k] + field[sx*i+sy*(j-1)+sz*k] - 2*field[sx*i+sy*j+sz*k]
												+ 0.5*(field1[sx*i+sy*(j+1)+sz*k] + field1[sx*i+sy*(j-1)+sz*k] - 2*field1[sx*i+sy*j+sz*k]));
				field2[sx*i+sy*j+sz*k+1] += idz*(field[sx*i+sy*j+sz*(k+1)] + field[sx*i+sy*j+sz*(k-1)] - 2*field[sx*i+sy*j+sz*k]
												+ 0.5*(field1[sx*i+sy*j+sz*(k+1)] + field1[sx*i+sy*j+sz*(k-1)] - 2*field1[sx*i+sy*j+sz*k]));
			}
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)
			{
				field[sx*i+sy*j+sz*k] += field2[sx*i+sy*j+sz*k];
				field[sx*i+sy*j+sz*k+1] += field2[sx*i+sy*j+sz*k+1];
			}
	
	cow_dfield_del(dphi1);
	cow_dfield_del(dphi2);
}

void rk4(cow_dfield *phi, double v, double dt)
{
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	int sx = cow_dfield_getstride(phi, 0);
	int sy = cow_dfield_getstride(phi, 1);
	int sz = cow_dfield_getstride(phi, 2);
	
	double dx = cow_domain_getgridspacing(dom, 0);
	double dy = cow_domain_getgridspacing(dom, 1);
	double dz = cow_domain_getgridspacing(dom, 2);
	
	double idx = dt*v*v / dx*dx;
	double idy = dt*v*v / dy*dy;
	double idz = dt*v*v / dz*dz;
	
	int i,j,k;
	
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
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)
			{
				field1[sx*i+sy*j+sz*k] = dt * field[sx*i+sy*j+sz*k+1];
				field1[sx*i+sy*j+sz*k+1] = idx*(field[sx*(i+1)+sy*j+sz*k] + field[sx*(i-1)+sy*j+sz*k] - 2*field[sx*i+sy*j+sz*k]);
				field1[sx*i+sy*j+sz*k+1] += idy*(field[sx*i+sy*(j+1)+sz*k] + field[sx*i+sy*(j-1)+sz*k] - 2*field[sx*i+sy*j+sz*k]);
				field1[sx*i+sy*j+sz*k+1] += idz*(field[sx*i+sy*j+sz*(k+1)] + field[sx*i+sy*j+sz*(k-1)] - 2*field[sx*i+sy*j+sz*k]);
			}
	cow_dfield_syncguard(dphi1);
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)
			{
				field2[sx*i+sy*j+sz*k] = dt * (field[sx*i+sy*j+sz*k+1]+0.5*field1[sx*i+sy*j+sz*k+1]);
				field2[sx*i+sy*j+sz*k+1] = idx*(field[sx*(i+1)+sy*j+sz*k] + field[sx*(i-1)+sy*j+sz*k] - 2*field[sx*i+sy*j+sz*k]
												+ 0.5*(field1[sx*(i+1)+sy*j+sz*k] + field1[sx*(i-1)+sy*j+sz*k] - 2*field1[sx*i+sy*j+sz*k]));
				field2[sx*i+sy*j+sz*k+1] += idy*(field[sx*i+sy*(j+1)+sz*k] + field[sx*i+sy*(j-1)+sz*k] - 2*field[sx*i+sy*j+sz*k]
												 + 0.5*(field1[sx*i+sy*(j+1)+sz*k] + field1[sx*i+sy*(j-1)+sz*k] - 2*field1[sx*i+sy*j+sz*k]));
				field2[sx*i+sy*j+sz*k+1] += idz*(field[sx*i+sy*j+sz*(k+1)] + field[sx*i+sy*j+sz*(k-1)] - 2*field[sx*i+sy*j+sz*k]
												 + 0.5*(field1[sx*i+sy*j+sz*(k+1)] + field1[sx*i+sy*j+sz*(k-1)] - 2*field1[sx*i+sy*j+sz*k]));
			}
	cow_dfield_syncguard(dphi2);
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)
			{
				field3[sx*i+sy*j+sz*k] = dt * (field[sx*i+sy*j+sz*k+1]+0.5*field2[sx*i+sy*j+sz*k+1]);
				field3[sx*i+sy*j+sz*k+1] = idx*(field[sx*(i+1)+sy*j+sz*k] + field[sx*(i-1)+sy*j+sz*k] - 2*field[sx*i+sy*j+sz*k]
												+ 0.5*(field2[sx*(i+1)+sy*j+sz*k] + field2[sx*(i-1)+sy*j+sz*k] - 2*field2[sx*i+sy*j+sz*k]));
				field3[sx*i+sy*j+sz*k+1] += idy*(field[sx*i+sy*(j+1)+sz*k] + field[sx*i+sy*(j-1)+sz*k] - 2*field[sx*i+sy*j+sz*k]
												 + 0.5*(field2[sx*i+sy*(j+1)+sz*k] + field2[sx*i+sy*(j-1)+sz*k] - 2*field2[sx*i+sy*j+sz*k]));
				field3[sx*i+sy*j+sz*k+1] += idz*(field[sx*i+sy*j+sz*(k+1)] + field[sx*i+sy*j+sz*(k-1)] - 2*field[sx*i+sy*j+sz*k]
												 + 0.5*(field2[sx*i+sy*j+sz*(k+1)] + field2[sx*i+sy*j+sz*(k-1)] - 2*field2[sx*i+sy*j+sz*k]));
			}
	cow_dfield_syncguard(dphi3);
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)
			{
				field4[sx*i+sy*j+sz*k] = dt * (field[sx*i+sy*j+sz*k+1]+field2[sx*i+sy*j+sz*k+1]);
				field4[sx*i+sy*j+sz*k+1] = idx*(field[sx*(i+1)+sy*j+sz*k] + field[sx*(i-1)+sy*j+sz*k] - 2*field[sx*i+sy*j+sz*k]
												+ field2[sx*(i+1)+sy*j+sz*k] + field2[sx*(i-1)+sy*j+sz*k] - 2*field2[sx*i+sy*j+sz*k]);
				field4[sx*i+sy*j+sz*k+1] += idy*(field[sx*i+sy*(j+1)+sz*k] + field[sx*i+sy*(j-1)+sz*k] - 2*field[sx*i+sy*j+sz*k]
												 + field2[sx*i+sy*(j+1)+sz*k] + field2[sx*i+sy*(j-1)+sz*k] - 2*field2[sx*i+sy*j+sz*k]);
				field4[sx*i+sy*j+sz*k+1] += idz*(field[sx*i+sy*j+sz*(k+1)] + field[sx*i+sy*j+sz*(k-1)] - 2*field[sx*i+sy*j+sz*k]
												 + field2[sx*i+sy*j+sz*(k+1)] + field2[sx*i+sy*j+sz*(k-1)] - 2*field2[sx*i+sy*j+sz*k]);
			}
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)
			{
				field[sx*i+sy*j+sz*k] += (field1[sx*i+sy*j+sz*k]+2*field2[sx*i+sy*j+sz*k]+2*field3[sx*i+sy*j+sz*k]+field4[sx*i+sy*j+sz*k]) / 6.0;
				field[sx*i+sy*j+sz*k+1] += (field1[sx*i+sy*j+sz*k+1]+2*field2[sx*i+sy*j+sz*k+1]+2*field3[sx*i+sy*j+sz*k+1]+field4[sx*i+sy*j+sz*k+1]) / 6.0;
			}
	
	
	cow_dfield_del(dphi1);
	cow_dfield_del(dphi2);
	cow_dfield_del(dphi3);
	cow_dfield_del(dphi4);
}

void leap_frog(cow_dfield *phi, double v, double dt)
{
	double *field = (double *) cow_dfield_getdatabuffer(phi);
	cow_domain *dom = cow_dfield_getdomain(phi);
	int ng = cow_domain_getguard(dom);
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	int sx = cow_dfield_getstride(phi, 0);
	int sy = cow_dfield_getstride(phi, 1);
	int sz = cow_dfield_getstride(phi, 2);
	
	double dx = cow_domain_getgridspacing(dom, 0);
	double dy = cow_domain_getgridspacing(dom, 1);
	double dz = cow_domain_getgridspacing(dom, 2);
	
	double idx = dt*v*v / dx*dx;
	double idy = dt*v*v / dy*dy;
	double idz = dt*v*v / dz*dz;
	
	int i,j,k;
	
	cow_dfield_syncguard(phi); 
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)
			{
				field[sx*i+sy*j+sz*k] += dt * field[sx*i+sy*j+sz*k+1];
			}

	cow_dfield_syncguard(phi);
	
	for(i = ng; i < nx+ng; i++)
		for(j = ng; j < ny+ng; j++)
			for(k = ng; k < nz+ng; k++)				
			{
				field[sx*i+sy*j+sz*k+1] += idx * (field[sx*(i+1)+sy*j+sz*k] - 2*field[sx*i+sy*j+sz*k] + field[sx*(i-1)+sy*j+sz*k]);
				field[sx*i+sy*j+sz*k+1] += idy * (field[sx*i+sy*(j+1)+sz*k] - 2*field[sx*i+sy*j+sz*k] + field[sx*i+sy*(j-1)+sz*k]);
				field[sx*i+sy*j+sz*k+1] += idz * (field[sx*i+sy*j+sz*(k+1)] - 2*field[sx*i+sy*j+sz*k] + field[sx*i+sy*j+sz*(k-1)]);
			}
}

int main(int argc, char *argv[])
{	
	cow_init(argc, argv, 0);
	
	int size = 128;
	int ng = 2;
	double a[3] = {-1,-1,-1};
	double b[3] = {1, 1, 1};
	
	if(argc > 1)
	{
		size = (int) strtol(argv[1], NULL, 10);
		if(size == 0)
		{
			printf("Hey, ya gotta give me an positive integer size numbnuts!\n");
			return 0;
		}
	}
	
	cow_domain *dom = cow_domain_new();
	cow_domain_setndim(dom, 3);
	cow_domain_setsize(dom, 0, size);
	cow_domain_setsize(dom, 1, size);
	cow_domain_setsize(dom, 2, size);
	cow_domain_setextent(dom, 0, a[0], b[0]);
	cow_domain_setextent(dom, 1, a[1], b[1]);
	cow_domain_setextent(dom, 2, a[2], b[2]);
	cow_domain_setguard(dom, ng);
	cow_domain_commit(dom);
	
	cow_dfield *phi = cow_dfield_new();
	cow_dfield_setdomain(phi, dom);
	cow_dfield_setname(phi, "phi");
	cow_dfield_addmember(phi, "phi");
	cow_dfield_addmember(phi, "phiv");
	cow_dfield_commit(phi);
	
	int nx = cow_domain_getnumlocalzonesinterior(dom, 0);
	int ny = cow_domain_getnumlocalzonesinterior(dom, 1);
	int nz = cow_domain_getnumlocalzonesinterior(dom, 2);
	
	printf("size: %d, ng: %d, nx: %d, ny: %d, nz: %d\n", size, ng, nx, ny, nz);
	printf("Sx: %d\n", cow_dfield_getstride(phi, 0));
	printf("Sy: %d\n", cow_dfield_getstride(phi, 1));
	printf("Sz: %d\n", cow_dfield_getstride(phi, 2));
	
	double cfl = 0.1;
	double T = 1.0;
	double v = 1.0;
	double params[5];
	
	params[0] = 0.3;
	params[1] = 0;
	initialize = &initial_sine;
	timestep = &rk4;
	
//	params[0] = 0.05*(b[0]-a[0]);
//	params[1] = 1.0;
//	initialize = &initial_impulse;
//	timestep = &rk4;
	
//	params[0] = 0.1;
//	params[1] = 1.0;
//	params[2] = 0.3;
//	initialize = &initial_spherical_pulse;
//	timestep = &rk4;
	
	int nt;
	
	//Initialize
	initialize(phi, v, params);
	//Run
	nt = run(phi, v, cfl, T);
	//Output
	
	cow_dfield_del(phi);
	cow_domain_del(dom);
	
	cow_finalize();
	
	return 0;
}

