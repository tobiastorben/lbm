#include "core.h"

void step(FlowData* flow, LatticeConsts* lc, SimParams* params) {
	//Update macroscopic variables	
	updateRho(flow,lc);		
	updateU(flow,lc);
	
	//Apply BCs
	inlet(flow,lc,params);//Poiseuille (Zou/He)		
	outlet(flow,lc);//Constant pressure (Zou/He)
	
	//Collide (Bhatnagar-Gross-Kroot model)
	collide(flow,lc,params);//Particle-Particle collisions		
	bounce(flow,lc,params);//Bounce-back collision with boundary
	
	//Streaming
	stream(flow,lc);
}

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
	double sumX, sumY,*fIn,*ex,*ey,*rho,*ux,*uy;
	int nx,ny,nxny;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	ex = lc->ex;
	ey = lc->ey;
	fIn = flow->fIn;
	rho = flow->rho;
	ux = flow->ux;
	uy = flow->uy;
	
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++){
			sumX = 0;
			sumY = 0;
			for (int k = 0; k < 9; k++) {
				sumX += fIn[nxny*k + ny*i + j]*ex[k];
				sumY += fIn[nxny*k + ny*i + j]*ey[k];
			}
			ux[ny*i + j] = sumX/rho[ny*i + j];
			uy[ny*i + j] = sumY/rho[ny*i + j];
		}
	} 				 
}



void collide(FlowData* flow, LatticeConsts* lc, SimParams* params){
	double u,fEq,uSq,rhoIJ,uxIJ,uyIJ,*ux,*uy,*ey,*ex,*fOut,*fIn,*rho,*w,tau;
	int i,j,k,nx,ny, nxny;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	ex = lc->ex;
	ey = lc->ey;
	fIn = flow->fIn;
	fOut = flow->fOut;
	rho = flow->rho;
	ux = flow->ux;
	uy = flow->uy;
	w = lc->w;
	tau = params->tau;
	
	for (i = 0; i < nx; i++){
		for (j = 0; j < ny; j++){
			uxIJ = ux[ny*i + j];
			uyIJ = uy[ny*i + j];
			rhoIJ = rho[ny*i + j];
			uSq = 1.5*(uxIJ*uxIJ+uyIJ*uyIJ);
			for (k = 0; k < 9; k++){
		u = 3.0*(ex[k]*uxIJ + ey[k]*uyIJ);
		fEq = rhoIJ*(w[k])*(1.0+u+0.5*u*u-uSq);
		fOut[nxny*k + ny*i + j] = fIn[nxny*k + ny*i + j]-(fIn[nxny*k + ny*i + j]-fEq)/tau;
		}
	}
	}			 
}



void stream(FlowData* flow, LatticeConsts* lc) {
	int i,j,k,nx,ny,nxny,*ex,*ey;
	double *fIn,*fOut;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
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