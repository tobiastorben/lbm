#include "postProcessing.h"

double* calcF(int ny, SimParams* params, double* rho) {
	double p,*F;
	int i,j,k,*bbCells,*bbCellMat,nBBcells;
	bbCells = params->bbCells;
	bbCellMat = params->bbCellMat;
	nBBcells = params->nBBcells;
	F = calloc(2,sizeof(double));

	for (k = 0; k < nBBcells; k++) {
		i = bbCells[k];
		j = bbCells[1*nBBcells +k];
		
		if (bbCellMat[(i+1)*ny + j] == 0) {
			p = rho[(i+1)*ny + j]-1;
			F[0] -= p;
		}
		
		if (bbCellMat[(i-1)*ny + j] == 0) {
			p = rho[(i-1)*ny + j]-1;
			F[0] += p;
		}
		
		if (bbCellMat[i*ny + j+1] == 0) {
			p = rho[i*ny + j+1]-1;
			F[1] -= p;
		}
		
		if (bbCellMat[i*ny + j-1] == 0) {
			p = rho[i*ny + j-1]-1;
			F[1] += p;
		}
	}
	return F;
}
