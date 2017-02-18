#include "postProcessing.h"

double* calcF(double* F, double* rho, int nx, int ny, int* bbCells, int nBBcells, int* bbCellMat) {
	double lift = 0;
	double drag = 0;
	double p;
	F[0] = 0;
	F[1] = 0;
	int i,j,k;
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
