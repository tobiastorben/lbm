#ifndef STRUCTS_H
#define STRUCTS_H

#include <pthread.h>

//Lattice constants
typedef struct {
	double w[9],ex[9],ey[9]; 
    int nx,ny,exI[9],eyI[9],opposite[9],
	westShiftArray[6],northShiftArray[6],eastShiftArray[6],southShiftArray[6];
	
} LatticeConsts;

//Simulation parameters
typedef struct {
	int nIter,cyclesPerWrite,startWrite,nBBcells,nThreads,blockSize,*bbCellMat,*bbCells,*outputSelect;
	double startVelX,startVelY,uRef,dxPhys,dtPhys,nuPhys,rhoPhys,tau;
	char *outDir, *obstaclePath;
} SimParams;

//Field variables of the flow
typedef struct {
double *rho,*ux,*uy,*fIn,*fOut;
} FlowData;

//Boundary conditions
typedef struct {
	void (*westFun)(FlowData*,LatticeConsts*,double,double);
	void (*northFun)(FlowData*,LatticeConsts*,int,int,double,double);
	void (*eastFun)(FlowData*,LatticeConsts*,double,double);
	void (*southFun)(FlowData*,LatticeConsts*,int,int,double,double);
	double westBC[2],northBC[2],eastBC[2],southBC[2];	
	int westBCType,northBCType,eastBCType,southBCType;
} BoundaryData;

//Data to pass to a solver thread
typedef struct {
	pthread_t thread;
	LatticeConsts* lc;
	SimParams* params;
	FlowData* flow;
	BoundaryData* bcdata;
	int startX,endX;
} ThreadData;

//Data to pass to a printer thread
typedef struct {
	SimParams* params;
	double *uxCpy,*uyCpy,*rhoCpy;
	int nx,ny,iter;
} PrintData;

#endif