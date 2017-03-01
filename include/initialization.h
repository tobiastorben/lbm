#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int* mapObstacleCells(int* bbCellMat,int nBBcells, int nx, int ny);
double* initUx(int nx,int ny, double uIn);
double* initUy(int nx,int ny);
double* initRho(int nx, int ny, double rho);
double* initFin(int nf,int nx, int ny, double* ex, double* ey,
				double* ux, double* uy, double* w, double* rho);
double* initFout(int nx,int ny, int nf);
int* readObstacle(int* nBBCells_p, int* nx_p, int* ny_p, char* path);
void nonDimensionalize(int nx, int ny, double dtPhys, double dxPhys, double nuPhys, double u0Phys, double* u0_p, double* tau_p);
