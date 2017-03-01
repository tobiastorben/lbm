#include "core.h"

void updateRho(FlowData* flow, LatticeConsts* lc) {
	double sum;
	int nx,ny;
	
	nx =lc->nx;
	ny = lc->ny;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++){
			sum = 0;
			for (int k = 0; k < 9; k++) {
				sum += flow->fIn[nx*ny*k + ny*i + j];
			}
			flow->rho[ny*i + j] = sum;
		}
	}
}

void updateU(FlowData* flow, LatticeConsts* lc){
	double sumX, sumY;
	int nx,ny,nxny;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;	

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++){
			sumX = 0;
			sumY = 0;
			for (int k = 0; k < 9; k++) {
				sumX += flow->fIn[nxny*k + ny*i + j]*lc->ex[k];
				sumY += flow->fIn[nxny*k + ny*i + j]*lc->ey[k];
			}
			flow->ux[ny*i + j] = sumX/flow->rho[ny*i + j];
			flow->uy[ny*i + j] = sumY/flow->rho[ny*i + j];
		}
	} 				 
}



void collide(FlowData* flow, LatticeConsts* lc, SimParams* params){
	double u, fEq, uSq,rhoIJ, uxIJ, uyIJ;
	int i,j,k,nx,ny, nxny;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	
	for (i = 0; i < nx; i++){
		for (j = 0; j < ny; j++){
			uxIJ = flow->ux[ny*i + j];
			uyIJ = flow->uy[ny*i + j];
			rhoIJ = flow->rho[ny*i + j];
			uSq = 1.5*(uxIJ*uxIJ+uyIJ*uyIJ);
			for (k = 0; k < 9; k++){
		u = 3.0*(lc->ex[k]*uxIJ + lc->ey[k]*uyIJ);
		fEq = rhoIJ*(lc->w[k])*(1.0+u+0.5*u*u-uSq);
		flow->fOut[nxny*k + ny*i + j] = flow->fIn[nxny*k + ny*i + j]-(flow->fIn[nxny*k + ny*i + j]-fEq)/params->tau;
		}
	}
	}			 
}



void stream(FlowData* flow, LatticeConsts* lc) {
	//Ugly solution. Improve?
	int i,j,k,nx,ny,nxny,*ex,*ey;
	double *fIn,*fOut;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nxny;
	ex = lc->exI;
	ey = lc->eyI;
	fIn = flow->fIn;
	fOut = flow->fOut;
	
	//Interior
	for (i = 1; i < nx-1; i++) {
		for (j = 1; j < ny-1; j++) {
			for (k = 0; k < 9; k++) {
				fIn[nxny*k+ny*i+j] = fOut[k*nxny+(i-ex[k])*ny +j-ey[k]];
			}
		}
	}	
	
	//Outlet
	for (j = 1; j < ny-1; j++) {
		for (k = 0; k < 9; k++) {
			if (ex[k]==-1) fIn[k*nxny+(nx-1)*ny+j] = fOut[k*nxny+j-ey[k]];
			else fIn[k*nxny+(nx-1)*ny+j] = fOut[k*nxny + (nx-1-ex[k])*ny+j-ey[k]];
		}
	}
	//Inlet
	for (j = 1; j < ny-1; j++) {
		for (k = 0; k < 9; k++) {
			if (ex[k]==1) fIn[k*nxny +j] = fOut[k*nxny+(nx-1)*ny+j-ey[k]];
			else fIn[k*nxny +j] = fOut[k*nxny-ex[k]*ny+j-ey[k]];
		}
	}
	//Top wall
	for (i = 1; i < nx-1; i++) {
		for (k = 0; k < 9; k++) {
			if (ey[k]==-1) fIn[k*nxny +i*ny + ny-1] = fOut[k*nxny +(i-ex[k])*ny];
			else fIn[k*nxny +i*ny + ny-1] = fOut[k*nxny +(i-ex[k])*ny +ny-1-ey[k]];
		}
	}
	//Bottom wall
	for (i = 1; i < nx-1; i++) {
		for (k = 0; k < 9; k++) {
			if (ey[k]==1) fIn[k*nxny+ny*i] = fOut[k*nxny +(i-ex[k])*ny+ny-1];
			else fIn[k*nxny + ny*i] = fOut[k*nxny + (i-ex[k])*ny-ey[k]];
		}
	}
	//Top right corner
	for (k = 0; k < 9; k++) {
			if (ex[k]==-1){
				if (ey[k]==-1) fIn[k*nxny+(nx-1)*ny+ny-1] = fOut[k*nxny];
				else fIn[k*nxny+(nx-1)*ny+ny-1] = fOut[k*nxny +ny-1-ey[k]];
			}
			else if (ey[k]==-1) fIn[k*nxny+(nx-1)*ny+ny-1] = fOut[k*nxny + (nx-1-ex[k])*ny];
			else fIn[k*nxny+(nx-1)*ny+ny-1] = fOut[k*nxny + (nx-1-ex[k])*ny +ny-1-ey[k]];
		}		
	//Bottom right corner
	for (k = 0; k < 9; k++) {
			if (ex[k]==-1){
				if (ey[k]==1) fIn[k*nxny+(nx-1)*ny] = fOut[k*nxny + ny-1];
				else fIn[k*nxny+(nx-1)*ny] = fOut[k*nxny-ey[k]];
				}
			else if (ey[k]==1) fIn[k*nxny+(nx-1)*ny] = fOut[k*nxny + (nx-1-ex[k])*ny + ny-1];
			else fIn[k*nxny+(nx-1)*ny] = fOut[k*nxny + (nx-1-ex[k])*ny -ey[k]];
		}	
	//Top left corner
	for (k = 0; k < 9; k++) {
			if (ex[k]==1){
				if (ey[k]==-1) fIn[k*nxny+ ny-1] = fOut[k*nxny+(nx-1)*ny];
				else fIn[k*nxny+ ny-1] = fOut[k*nxny+(nx-1)*ny+ny-1-ey[k]];
			} 
			else if (ey[k]==-1) fIn[k*nxny+ ny-1] = fOut[k*nxny -ex[k]*ny];
			else fIn[k*nxny+ ny-1] = fOut[k*nxny -ex[k]*ny +ny-1-ey[k]];
		}		
	//Bottom left corner
	for (k = 0; k < 9; k++) {
			if (ex[k]==1){
				if (ey[k]==1) fIn[k*nxny] = fOut[k*nxny+(nx-1)*ny+ny-1];
				else fIn[k*nxny] = fOut[k*nxny+(nx-1)*ny+-ey[k]];
				}
			else if (ey[k]==1) fIn[k*nxny] = fOut[k*nxny-ex[k]*ny+ny-1];
			else fIn[k*nxny] = fOut[k*nxny-ex[k]*ny-ey[k]];
		}
}