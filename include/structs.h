#ifndef STRUCTS_H
#define STRUCTS_H

#include <pthread.h>

typedef struct {
	double w[9],ex[9],ey[9]; 
    int nx,ny,exI[9],eyI[9],opposite[9],
	westShiftArray[6],northShiftArray[6],eastShiftArray[6],southShiftArray[6];
	
} LatticeConsts;

typedef struct {
	int nIter,cyclesPerWrite,startWrite,nBBcells,nThreads,*bbCellMat,*bbCells,*outputSelect;
	double startVelX,startVelY,uRef,dxPhys,dtPhys,nuPhys,rhoPhys,tau;
	char *outDir, *obstaclePath;
} SimParams;

typedef struct {
double *rho,*ux,*uy,*fIn,*fOut;
} FlowData;

typedef struct {
	void (*westFun)(FlowData*,LatticeConsts*,double,double);
	void (*northFun)(FlowData*,LatticeConsts*,int,int,double,double);
	void (*eastFun)(FlowData*,LatticeConsts*,double,double);
	void (*southFun)(FlowData*,LatticeConsts*,int,int,double,double);
	double westBC[2],northBC[2],eastBC[2],southBC[2];	
	int westBCType,northBCType,eastBCType,southBCType;
} BoundaryData;

typedef struct {
	pthread_t thread;
	LatticeConsts* lc;
	SimParams* params;
	FlowData* flow;
	BoundaryData* bcdata;
	int startX,endX;
} ThreadData;

typedef struct {
	SimParams* params;
	double *uxCpy,*uyCpy,*rhoCpy;
	int nx,ny,iter;
} PrintData;


#endif