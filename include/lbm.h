#ifndef LBM_H
#define LBM_H



typedef struct {
	double w[9],ex[9],ey[9]; 
    int nx,ny,exI[9],eyI[9],opposite[9];
} LatticeConsts;

typedef struct {
	int nIter,cyclesPerWrite,startWrite,nBBcells,*bbCellMat,*bbCells,*outputSelect;
	double u0,u0Phys,dxPhys,dtPhys,nuPhys,tau;
	char *outDir, *obstaclePath;
} SimParams;

typedef struct {
double *rho,*ux,*uy,*fIn,*fOut;
} FlowData;

#endif