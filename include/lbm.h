#ifndef LBM_H
#define LBM_H



typedef struct {
	int nx,ny;
	double w[9];
	double ex[9];
    double ey[9]; 
    int exI[9];
    int eyI[9];
    int opposite[9];
} LatticeConsts;

typedef struct {
	int nIter,cyclesPerWrite,startWrite,nBBcells,*bbCellMat,*bbCells,*outputSelect;
	double u0,u0Phys,dxPhys,dtPhys,nuPhys,tau;
	char* outDir, *obstaclePath;
} SimParams;

typedef struct {
double *rho,*ux,*uy,*fIn,*fOut;
} FlowData;

#endif