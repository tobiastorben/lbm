#ifndef STRUCTS_H
#define STRUCTS_H

#include <pthread.h>

typedef struct {
	double w[9],ex[9],ey[9]; 
    int nx,ny,exI[9],eyI[9],opposite[9];
} LatticeConsts;

typedef struct {
	int nIter,cyclesPerWrite,startWrite,nBBcells,nThreads,*bbCellMat,*bbCells,*outputSelect;
	double u0,u0Phys,dxPhys,dtPhys,nuPhys,rhoPhys,tau;
	char *outDir, *obstaclePath;
} SimParams;

typedef struct {
double *rho,*ux,*uy,*fIn,*fOut;
} FlowData;

typedef struct {
	pthread_t thread;
	LatticeConsts* lc;
	SimParams* params;
	FlowData* flow;
	int startX,endX,count;
} ThreadData;

typedef struct {
	SimParams* params;
	double *uxCpy,*uyCpy,*rhoCpy;
	int nx,ny,iter;
} PrintData;

#endif