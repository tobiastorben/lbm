#include "boundaries.h"


void southXY(FlowData* flow, LatticeConsts * lc, SimParams* params,int startX, int endX) {
	double *fIn,*rho,uyBC, uxBC;
	int nx,ny,nxny;
	
	uyBC = 0.0; uxBC = 0.05;
	fIn = flow->fIn;
	rho = flow->rho;
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	
	for (int i = startX; i <= endX;i++) {
		rho[ny*i] = (1.0/(1.0-uyBC))*(fIn[i*ny] + fIn[nxny + i*ny] + fIn[3*nxny + i*ny]
					+ 2.0*(fIn[4*nxny + i*ny] + fIn[7*nxny + i*ny] + fIn[8*nxny + i*ny]));
		fIn[2*nxny + i*ny] = fIn[4*nxny + i*ny] + (2.0/3.0)*(rho[i*ny])*uyBC;
		fIn[5*nxny + i*ny] = fIn[7*nxny + i*ny] - 0.5*(fIn[nxny + i*ny]-fIn[3*nxny + i*ny]) + 0.5*(rho[i*ny])*uxBC + (1.0/6.0)*(rho[i*ny])*uyBC;
		fIn[6*nxny + i*ny] = fIn[8*nxny + i*ny] + 0.5*(fIn[nxny + i*ny]-fIn[3*nxny + i*ny]) - 0.5*(rho[i*ny])*uxBC + (1.0/6.0)*(rho[i*ny])*uyBC;
		flow->ux[ny*i] = uxBC;
		flow->uy[ny*i] = uyBC;		
	}
}

void southPX(FlowData* flow, LatticeConsts * lc, SimParams* params,int startX, int endX) {
	double *fIn,*uy,rhoBC, uxBC;
	int nx,ny,nxny;
	
	rhoBC = 1.0; uxBC = 0.0;
	fIn = flow->fIn;
	uy = flow->uy;
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	
	for (int i = startX; i <= endX;i++) {
		flow->ux[ny*i] = uxBC;
		flow->rho[ny*i] = rhoBC;
		uy[i*ny] = 1 - (fIn[i*ny] + fIn[nxny + i*ny] + fIn[3*nxny + i*ny]
					+ 2.0*(fIn[4*nxny + i*ny] + fIn[7*nxny + i*ny] + fIn[8*nxny + i*ny]))/rhoBC;
		fIn[2*nxny + i*ny] = fIn[4*nxny + i*ny] + (2.0/3.0)*(rhoBC)*uy[ny*i];
		fIn[5*nxny + i*ny] = fIn[7*nxny + i*ny] - 0.5*(fIn[nxny + i*ny]-fIn[3*nxny + i*ny]) + 0.5*(rhoBC)*uxBC + (1.0/6.0)*(rhoBC)*uy[ny*i];
		fIn[6*nxny + i*ny] = fIn[8*nxny + i*ny] + 0.5*(fIn[nxny + i*ny]-fIn[3*nxny + i*ny]) - 0.5*(rhoBC)*uxBC + (1.0/6.0)*(rhoBC)*uy[ny*i];
		flow->ux[ny*i] = uxBC;
	}
}


void northXY(FlowData* flow, LatticeConsts * lc, SimParams* params,int startX, int endX) {
	double *fIn,*rho,uyBC, uxBC;
	int nx,ny,nxny;
	
	uyBC = 0; uxBC = 0.05;
	fIn = flow->fIn;
	rho = flow->rho;
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	
	for (int i = startX; i <= endX;i++) {
		rho[ny*i + ny -1] = (1.0/(1.0+uyBC))*(fIn[ny*i + ny -1] + fIn[nxny + ny*i + ny -1] + fIn[3*nxny + ny*i + ny -1]
					+ 2.0*(fIn[2*nxny + ny*i + ny -1] + fIn[5*nxny + ny*i + ny -1] + fIn[6*nxny + ny*i + ny -1]));
		fIn[4*nxny + ny*i + ny -1] = fIn[2*nxny + ny*i + ny -1] - (2.0/3.0)*(rho[ny*i + ny -1])*uyBC;
		fIn[8*nxny + ny*i + ny -1] = fIn[6*nxny + ny*i + ny -1] + 0.5*(fIn[3*nxny + ny*i + ny -1]-fIn[nxny + ny*i + ny -1]) + 0.5*(rho[ny*i + ny -1])*uxBC - (1.0/6.0)*(rho[ny*i + ny -1])*uyBC;
		fIn[7*nxny + ny*i + ny -1] = fIn[5*nxny + ny*i + ny -1] + 0.5*(fIn[nxny + ny*i + ny -1]-fIn[3*nxny + ny*i + ny -1]) - 0.5*(rho[ny*i + ny -1])*uxBC - (1.0/6.0)*(rho[ny*i + ny -1])*uyBC;
		flow->ux[ny*i + ny-1] = uxBC;
		flow->uy[ny*i + ny-1] = uyBC;		
	}
}

void northPX(FlowData* flow, LatticeConsts * lc, SimParams* params,int startX, int endX) {
	double *fIn,*uy,rhoBC,uxBC;
	int nx,ny,nxny;
	
	rhoBC = 1.01; uxBC = 0;
	fIn = flow->fIn;
	uy = flow->uy;
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	
	for (int i = startX; i <= endX;i++) {
		flow->ux[ny*i + ny-1] = uxBC;
		flow->rho[ny*i + ny-1] = rhoBC;
		uy[ny*i + ny-1] = -1.0 + (fIn[ny*i + ny -1] + fIn[nxny + ny*i + ny -1] + fIn[3*nxny + ny*i + ny -1]
					+ 2.0*(fIn[2*nxny + ny*i + ny -1] + fIn[5*nxny + ny*i + ny -1] + fIn[6*nxny + ny*i + ny -1]))/rhoBC;					
		fIn[4*nxny + ny*i + ny -1] = fIn[2*nxny + ny*i + ny -1] - (2.0/3.0)*(rhoBC)*uy[ny*i + ny-1];
		fIn[8*nxny + ny*i + ny -1] = fIn[6*nxny + ny*i + ny -1] + 0.5*(fIn[3*nxny + ny*i + ny -1]-fIn[nxny + ny*i + ny -1]) + 0.5*(rhoBC)*uxBC - (1.0/6.0)*(rhoBC)*uy[ny*i + ny-1];
		fIn[7*nxny + ny*i + ny -1] = fIn[5*nxny + ny*i + ny -1] + 0.5*(fIn[nxny + ny*i + ny -1]-fIn[3*nxny + ny*i + ny -1]) - 0.5*(rhoBC)*uxBC - (1.0/6.0)*(rhoBC)*uy[ny*i + ny-1];		
	}
}

void westXY(FlowData* flow, LatticeConsts* lc, SimParams* params) {
	double sum1, sum2, rhoJ,*fIn,uxBC, uyBC;
	int nx,ny;
	uxBC = 0.0;uyBC = 0.0;
	nx = lc->nx;
	ny = lc->ny;	
	fIn = flow->fIn;

	for (int j = 1; j < (ny-1);j++) {
		flow->ux[j] = uxBC;
		flow->uy[j] = uyBC;
		sum1 = fIn[j] + fIn[2*nx*ny +j] + fIn[4*nx*ny+j];
		sum2 = fIn[3*nx*ny +j] + fIn[6*nx*ny + j] + fIn[7*nx*ny + j];
		rhoJ = (1.0/(1.0-uxBC))*(sum1 + 2.0*sum2);
		fIn[1*nx*ny+j] = fIn[3*nx*ny+j] + (2.0/3.0)*rhoJ*uxBC;
		fIn[5*nx*ny+j] = fIn[7*nx*ny+j] + 0.5*(fIn[4*nx*ny+j]-fIn[2*nx*ny+j]) + 0.5*(rhoJ*uyBC) + (1.0/6)*(rhoJ*uxBC);
		fIn[8*nx*ny+j] = fIn[6*nx*ny+j] + 0.5*(fIn[2*nx*ny+j]-fIn[4*nx*ny+j]) - 0.5*(rhoJ*uyBC) + (1.0/6)*(rhoJ*uxBC);
		flow->rho[j] = rhoJ;
	}
}

void westPY(FlowData* flow, LatticeConsts* lc, SimParams* params) {
	double sum1,sum2,*fIn,*ux,rhoBC, uyBC;
	int nx,ny;
	rhoBC = 1;uyBC = 0.0;
	nx = lc->nx;
	ny = lc->ny;	
	fIn = flow->fIn;
	ux = flow->ux;

	for (int j = 1; j < (ny-1);j++) {
		flow->rho[j] = rhoBC;
		flow->uy[j] = uyBC;
		sum1 = fIn[j] + fIn[2*nx*ny +j] + fIn[4*nx*ny+j];
		sum2 = fIn[3*nx*ny +j] + fIn[6*nx*ny + j] + fIn[7*nx*ny + j];
		ux[j] = 1.0-(sum1 + 2.0*sum2)/rhoBC;
		fIn[1*nx*ny+j] = fIn[3*nx*ny+j] + (2.0/3.0)*rhoBC*ux[j];
		fIn[5*nx*ny+j] = fIn[7*nx*ny+j] + 0.5*(fIn[4*nx*ny+j]-fIn[2*nx*ny+j]) + 0.5*(rhoBC*uyBC) + (1.0/6)*(rhoBC*ux[j]);
		fIn[8*nx*ny+j] = fIn[6*nx*ny+j] + 0.5*(fIn[2*nx*ny+j]-fIn[4*nx*ny+j]) - 0.5*(rhoBC*uyBC) + (1.0/6)*(rhoBC*ux[j]);
	}
}

void eastPY(FlowData* flow, LatticeConsts* lc) {
	double sum1,sum2,*fIn,rhoBC,uyBC,*ux;
	int nx,ny;
	
	nx = lc->nx;
	ny = lc->ny;	
	fIn = flow->fIn;
	ux = flow->ux;
	rhoBC = 1.01;uyBC = 0.0;
	
	for (int j = 1; j < (ny-1);j++) {
		flow->rho[(nx-1)*ny +j] = rhoBC;
		flow->uy[(nx-1)*ny +j] = uyBC;
		sum1 = fIn[(nx-1)*ny+j] + fIn[2*nx*ny+(nx-1)*ny+j] + fIn[4*nx*ny+(nx-1)*ny+j];
		sum2 = fIn[1*nx*ny+(nx-1)*ny+j] + fIn[5*nx*ny+(nx-1)*ny+j] + fIn[8*nx*ny+(nx-1)*ny+j];
		ux[(nx-1)*ny +j] = -1.0 + (sum1 + 2.0*sum2)/rhoBC;
		fIn[3*nx*ny+(nx-1)*ny + j] = fIn[1*nx*ny+(nx-1)*ny + j] - (2.0/3)*ux[(nx-1)*ny +j];
		fIn[7*nx*ny+(nx-1)*ny + j] = fIn[5*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[2*nx*ny+(nx-1)*ny + j]-fIn[4*nx*ny+(nx-1)*ny + j]) - 0.5*rhoBC*uyBC - (1.0/6)*rhoBC*ux[(nx-1)*ny +j];
		fIn[6*nx*ny+(nx-1)*ny + j] = fIn[8*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[4*nx*ny+(nx-1)*ny + j]-fIn[2*nx*ny+(nx-1)*ny + j]) + 0.5*rhoBC*uyBC - (1.0/6)*rhoBC*ux[(nx-1)*ny +j];
	}
}

void eastXY(FlowData* flow, LatticeConsts* lc) {
	double sum1,sum2,*fIn,uxBC,uyBC,rhoJ;
	int nx,ny;
	
	nx = lc->nx;
	ny = lc->ny;
	fIn = flow->fIn;
	uxBC = 0.0;uyBC = 0.0;
	
	for (int j = 1; j < (ny-1);j++) {
		flow->ux[(nx-1)*ny +j] = uxBC;
		flow->uy[(nx-1)*ny +j] = uyBC;
		sum1 = fIn[(nx-1)*ny+j] + fIn[2*nx*ny+(nx-1)*ny+j] + fIn[4*nx*ny+(nx-1)*ny+j];
		sum2 = fIn[1*nx*ny+(nx-1)*ny+j] + fIn[5*nx*ny+(nx-1)*ny+j] + fIn[8*nx*ny+(nx-1)*ny+j];
		rhoJ = (sum1 + 2.0*sum2)/(1.0-uxBC);
		fIn[3*nx*ny+(nx-1)*ny + j] = fIn[1*nx*ny+(nx-1)*ny + j] - (2.0/3)*uxBC;
		fIn[7*nx*ny+(nx-1)*ny + j] = fIn[5*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[2*nx*ny+(nx-1)*ny + j]-fIn[4*nx*ny+(nx-1)*ny + j]) - 0.5*rhoJ*uyBC - (1.0/6)*rhoJ*uxBC;
		fIn[6*nx*ny+(nx-1)*ny + j] = fIn[8*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[4*nx*ny+(nx-1)*ny + j]-fIn[2*nx*ny+(nx-1)*ny + j]) + 0.5*rhoJ*uyBC - (1.0/6)*rhoJ*uxBC;
		flow->rho[(nx-1)*ny +j] = rhoJ;
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
}