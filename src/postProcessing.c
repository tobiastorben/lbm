#include "postProcessing.h"

double* calcF(int ny, SimParams* params, double* rho) {
	double p,*F,c,dx,dt,rhoPhys;
	int i,j,k,*bbCells,*bbCellMat,nBBcells;
	bbCells = params->bbCells;
	bbCellMat = params->bbCellMat;
	nBBcells = params->nBBcells;
	dx = params->dxPhys;
	dt = params->dtPhys;
	rhoPhys = params->rhoPhys;
	c = dx/dt;
	F = calloc(2,sizeof(double));

	for (k = 0; k < nBBcells; k++) {
		i = bbCells[k];
		j = bbCells[1*nBBcells +k];
		
		if (bbCellMat[(i+1)*ny + j] == 0) {
			p = c*c*rhoPhys*(1.0/3.0)*(rho[(i+1)*ny + j]-1);
			F[0] -= p*dx;
		}
		
		if (bbCellMat[(i-1)*ny + j] == 0) {
			p = c*c*rhoPhys*(1.0/3.0)*(rho[(i-1)*ny + j]-1);
			F[0] += p*dx;
		}
		
		if (bbCellMat[i*ny + j+1] == 0) {
			p = c*c*rhoPhys*(1.0/3.0)*(rho[i*ny + j+1]-1);
			F[1] -= p*dx;
		}
		
		if (bbCellMat[i*ny + j-1] == 0) {
			p = c*c*rhoPhys*(1.0/3.0)*(rho[i*ny + j-1]-1);
			F[1] += p*dx;
		}
	}
	return F;
}
