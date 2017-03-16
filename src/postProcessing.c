#include "postProcessing.h"

double* calcF(int ny, SimParams* params, double* rho, double* ux, double* uy) {
	double p,*F,c,dx,dt,rhoPhys,visc,fric;
	int i,j,k,*bbCells,*bbCellMat,nBBcells;
	bbCells = params->bbCells;
	bbCellMat = params->bbCellMat;
	nBBcells = params->nBBcells;
	dx = params->dxPhys;
	dt = params->dtPhys;
	rhoPhys = params->rhoPhys;
	visc = params->nuPhys;
	c = dx/dt;
	F = calloc(2,sizeof(double));

	for (k = 0; k < nBBcells; k++) {
		i = bbCells[k];
		j = bbCells[1*nBBcells +k];
		
		if (bbCellMat[(i+1)*ny + j] == 0) {
			p = c*c*rhoPhys*(1.0/3.0)*(rho[(i+1)*ny + j]-1);
			fric = c*rhoPhys*visc*(uy[(i+2)*ny + j] - uy[(i+1)*ny + j]);
			F[0] -= p*dx;
			F[1] += fric;
		}
		
		if (bbCellMat[(i-1)*ny + j] == 0) {
			p = c*c*rhoPhys*(1.0/3.0)*(rho[(i-1)*ny + j]-1);
			fric = c*rhoPhys*visc*(uy[(i-2)*ny + j] - uy[(i-1)*ny + j]);
			F[0] += p*dx;
			F[1] += fric;
		}
		
		if (bbCellMat[i*ny + j+1] == 0) {
			p = c*c*rhoPhys*(1.0/3.0)*(rho[i*ny + j+1]-1);
			fric = c*rhoPhys*visc*(ux[i*ny + j+2] - ux[i*ny + j+1]);
			F[1] -= p*dx;
			F[0] += fric;
		}
		
		if (bbCellMat[i*ny + j-1] == 0) {
			p = c*c*rhoPhys*(1.0/3.0)*(rho[i*ny + j-1]-1);
			fric = c*rhoPhys*visc*(ux[i*ny + j-2] - ux[i*ny + j-1]);
			F[1] += p*dx;
			F[0] += fric;
		}
	}
	return F;
}
