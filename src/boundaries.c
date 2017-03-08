#include "boundaries.h"

void inlet(FlowData* flow, LatticeConsts* lc, SimParams* params) {
	double y, sum1, sum2, rhoJ, coeff,H,*fIn,*ux,*uy;
	int nx,ny;
	
	nx = lc->nx;
	ny = lc->ny;	
	ux = flow->ux;
	uy = flow->uy;
	fIn = flow->fIn;
	
	H = (double) ny-2.0;
	coeff = 4.0*(params->u0)/(H*H);
	for (int j = 1; j < (ny-1);j++) {
		y = j-0.5;
		
		//Poiseuille and Zou/He BC
		ux[j] = coeff*y*(H-y);
		uy[j] = 0;
		sum1 = fIn[j] + fIn[2*nx*ny +j] + fIn[4*nx*ny+j];
		sum2 = fIn[3*nx*ny +j] + fIn[6*nx*ny + j] + fIn[7*nx*ny + j];
		rhoJ = (1.0/(1.0-ux[j]))*(sum1 + 2.0*sum2);
		fIn[1*nx*ny+j] = fIn[3*nx*ny+j] + (2.0/3)*rhoJ*ux[j];
		fIn[5*nx*ny+j] = fIn[7*nx*ny+j] + 0.5*(fIn[4*nx*ny+j]-fIn[2*nx*ny+j])
					   +(1.0/6)*(rhoJ*ux[j]);
		fIn[8*nx*ny+j] = fIn[6*nx*ny+j] + 0.5*(fIn[2*nx*ny+j]-fIn[4*nx*ny+j])
						+(1.0/6)*(rhoJ*ux[j]);
		flow->rho[j] = rhoJ;
	}
}

void outlet(FlowData* flow, LatticeConsts* lc) {
	double sum1,sum2,*fIn,*ux,*uy;
	int nx,ny;
	
	nx = lc->nx;
	ny = lc->ny;	
	fIn = flow->fIn;
	ux = flow->ux;
	uy = flow->uy;
	
	for (int j = 1; j < (ny-1);j++) {
		sum1 = fIn[(nx-1)*ny+j] + fIn[2*nx*ny+(nx-1)*ny+j] + fIn[4*nx*ny+(nx-1)*ny+j];
		sum2 = fIn[1*nx*ny+(nx-1)*ny+j] + fIn[5*nx*ny+(nx-1)*ny+j] + fIn[8*nx*ny+(nx-1)*ny+j];
		ux[(nx-1)*ny +j] = -1.0 +sum1 + 2.0*sum2;
		uy[(nx-1)*ny +j] = 0;
		fIn[3*nx*ny+(nx-1)*ny + j] = fIn[1*nx*ny+(nx-1)*ny + j] - (2.0/3)*ux[(nx-1)*ny +j];
		fIn[7*nx*ny+(nx-1)*ny + j] = fIn[5*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[2*nx*ny+(nx-1)*ny + j]-fIn[4*nx*ny+(nx-1)*ny + j])
					   -(1.0/6)*ux[(nx-1)*ny +j];
		fIn[6*nx*ny+(nx-1)*ny + j] = fIn[8*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[4*nx*ny+(nx-1)*ny + j]-fIn[2*nx*ny+(nx-1)*ny + j])
					  -(1.0/6)*ux[(nx-1)*ny +j];
		flow->rho[(nx-1)*ny +j] = 1.0;
	}
}

void bounce(FlowData* flow, LatticeConsts* lc, SimParams* params){
	int i,j,k,l,nx,ny,*opp;
	
	nx = lc->nx;
	ny = lc->ny;
	opp = lc->opposite;
	
	//Obstacle
	for (l = 0; l < params->nBBcells; l++) {
		for (k = 0; k < 9; k++){
				i = params->bbCells[l];
				j = params->bbCells[params->nBBcells+l];
				flow->fOut[nx*ny*k + ny*i + j] = flow->fIn[nx*ny*opp[k] + ny*i + j];
		}
	}
	
	//Top/Bottow wall
	for (i = 0; i < nx; i++) {
		for (k = 0; k < 9; k++){
		flow->fOut[k*nx*ny + ny*i] = flow->fIn[nx*ny*opp[k] + ny*i];
			flow->fOut[k*nx*ny +i*ny+ny-1] = flow->fIn[nx*ny*opp[k] + ny*i + ny-1];
		}
	}
}