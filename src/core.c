#include "core.h"
#include <stdio.h>

void step(FlowData* flow, LatticeConsts* lc, SimParams* params, ThreadData* tdata) {
	int nThreads = params->nThreads;
	
	//streamBlockBoundaries(flow,lc,params);
	stream(flow,lc);	

	pthread_create(&(tdata[0].thread),NULL,updateFirstBlock,&tdata[0]);
	pthread_create(&(tdata[nThreads-1].thread),NULL,updateLastBlock,&tdata[nThreads-1]);
	
	for (int i = 1; i < nThreads-1; i++) {
		pthread_create(&(tdata[i].thread),NULL,updateBlock,&tdata[i]);
	}
		
	for (int i = 0; i < nThreads; i++) {
		pthread_join(tdata[i].thread,NULL);
	}

	bounce(flow,lc,params);//Bounce-back collision with boundary
}

void* updateBlock(void* tdata_void) {
	ThreadData* tdata = (ThreadData*) tdata_void;
	BoundaryData* bcdata = tdata->bcdata;
	//streamBlockInterior(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	updateRho(tdata->flow,tdata->lc,tdata->startX,tdata->endX);		
	updateU(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	(*(bcdata->southFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX,bcdata->southBC[0], bcdata->southBC[1]);
	(*(bcdata->northFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX,bcdata->northBC[0], bcdata->northBC[1]);
	collide(tdata->flow,tdata->lc,tdata->params,tdata->startX,tdata->endX);
	return NULL;
}

void* updateFirstBlock(void* tdata_void) {
	ThreadData* tdata = (ThreadData*) tdata_void;
	BoundaryData* bcdata = tdata->bcdata;
	//streamBlockInterior(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	updateRho(tdata->flow,tdata->lc,tdata->startX,tdata->endX);		
	updateU(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	(*(bcdata->southFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX,bcdata->southBC[0], bcdata->southBC[1]);
	(*(bcdata->northFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX,bcdata->northBC[0], bcdata->northBC[1]);
	(*(bcdata->westFun))(tdata->flow,tdata->lc,bcdata->westBC[0], bcdata->westBC[1]);	
	collide(tdata->flow,tdata->lc,tdata->params,tdata->startX,tdata->endX);
	return NULL;
}

void* updateLastBlock(void* tdata_void) {
	ThreadData* tdata = (ThreadData*) tdata_void;
	BoundaryData* bcdata = tdata->bcdata;
	//streamBlockInterior(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	updateRho(tdata->flow,tdata->lc,tdata->startX,tdata->endX);		
	updateU(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	(*(bcdata->southFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX,bcdata->southBC[0], bcdata->southBC[1]);
	(*(bcdata->northFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX,bcdata->northBC[0], bcdata->northBC[1]);
	(*(bcdata->eastFun))(tdata->flow,tdata->lc,bcdata->eastBC[0], bcdata->eastBC[1]);
	collide(tdata->flow,tdata->lc,tdata->params,tdata->startX,tdata->endX);
	return NULL;
}

void updateRho(FlowData* flow, LatticeConsts* lc, int startX, int endX) {
	double sum;
	int nx,ny;
	
	nx =lc->nx;
	ny = lc->ny;
	for (int i = startX; i <= endX; i++) {
		for (int j = 0; j < ny; j++){
			sum = 0;
			for (int k = 0; k < 9; k++) {
				sum += flow->fIn[nx*ny*k + ny*i + j];
			}
			flow->rho[ny*i + j] = sum;
		}
	}
}

void updateU(FlowData* flow, LatticeConsts* lc, int startX, int endX){
	double sumX,sumY,*fIn,*ex,*ey,*rho,*ux,*uy;
	int nx,ny,nxny,i,j,k;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	ex = lc->ex;
	ey = lc->ey;
	fIn = flow->fIn;
	rho = flow->rho;
	ux = flow->ux;
	uy = flow->uy;
	
	for (i = startX; i <= endX; i++) {
		if (i != 0 && i != 199) {
		for (j = 0; j < ny; j++){
			sumX = 0;
			sumY = 0;
			for (k = 0; k < 9; k++) {
				sumX += fIn[nxny*k + ny*i + j]*ex[k];
				sumY += fIn[nxny*k + ny*i + j]*ey[k];
			}
		
			ux[ny*i + j] = sumX/rho[ny*i + j];
			uy[ny*i + j] = sumY/rho[ny*i + j];
		}
	} 	
	}	
}



void collide(FlowData* flow, LatticeConsts* lc, SimParams* params, int startX, int endX){
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

	
	for (i = startX; i <= endX; i++){
		for (j = 0; j < ny; j++){
			uxIJ = ux[ny*i + j];
			uyIJ = uy[ny*i + j];
			rhoIJ = rho[ny*i + j];
			uSq = 1.5*(uxIJ*uxIJ+uyIJ*uyIJ);
			for (k = 0; k < 9; k++){
				u = 3.0*(ex[k]*uxIJ + ey[k]*uyIJ);
				fEq = rhoIJ*(w[k])*(1.0+u+0.5*u*u-uSq);
				fOut[nxny*k + ny*i + j] = fIn[nxny*k + ny*i + j]-(fIn[nxny*k + ny*i + j]-fEq)/tau;//Change to *omega for speed? 
			}
		}
	}			 
}

void streamBlockInterior(FlowData* flow,LatticeConsts* lc,int startX, int endX) {
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
	for (i = startX+2; i <= endX; i++) {
		for (j = 1; j < ny-1; j++) {
			for (k = 0; k < 9; k++) {
				fIn[nxny*k+ny*i+j] = fOut[k*nxny+(i-ex[k])*ny +j-ey[k]];
			}
		}
	}	
}

void streamBlockBoundaries(FlowData* flow,LatticeConsts* lc,SimParams* params){
	int i,j,k,nx,ny,nxny,*ex,*ey,blockSize,nThreads;
	double *fIn,*fOut;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	ex = lc->exI;
	ey = lc->eyI;
	fIn = flow->fIn;
	fOut = flow->fOut;
	nThreads = params->nThreads;
	blockSize = nx/nThreads;
	
	for (i = blockSize; i < nx-1; i+=blockSize) {
		for (j = 1; j < ny-1; j++) {
			for (k = 0; k < 9; k++) {
				fIn[nxny*k+ny*i+j] = fOut[k*nxny+(i-ex[k])*ny +j-ey[k]];
				fIn[nxny*k+ny*(i+1)+j] = fOut[k*nxny+(i+1-ex[k])*ny +j-ey[k]];
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
 			if (ex[k]==-1) 1;//fIn[k*nxny+(nx-1)*ny+j] = fOut[k*nxny+j-ey[k]];
 			else fIn[k*nxny+(nx-1)*ny+j] = fOut[k*nxny + (nx-1-ex[k])*ny+j-ey[k]];
 		}
 	}
 	//Inlet
 	for (j = 1; j < ny-1; j++) {
 		for (k = 0; k < 9; k++) {
 			if (ex[k]==1) 1;//fIn[k*nxny +j] = fOut[k*nxny+(nx-1)*ny+j-ey[k]];
 			else fIn[k*nxny +j] = fOut[k*nxny-ex[k]*ny+j-ey[k]];
 		}
 	}
 	//Top wall
 	for (i = 1; i < nx-1; i++) {
 		for (k = 0; k < 9; k++) {
 			if (ey[k]==-1) 1;//fIn[k*nxny +i*ny + ny-1] = fOut[k*nxny +(i-ex[k])*ny];
 			else fIn[k*nxny +i*ny + ny-1] = fOut[k*nxny +(i-ex[k])*ny +ny-1-ey[k]];
 		}
 	}
 	//Bottom wall
 	for (i = 1; i < nx-1; i++) {
 		for (k = 0; k < 9; k++) {
 			if (ey[k]==1) 1;//fIn[k*nxny+ny*i] = fOut[k*nxny +(i-ex[k])*ny+ny-1];
 			else fIn[k*nxny + ny*i] = fOut[k*nxny + (i-ex[k])*ny-ey[k]];
 		}
 	}
 	///Top right corner
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
 