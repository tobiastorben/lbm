#include "core.h"

void step(FlowData* flow, LatticeConsts* lc, SimParams* params, ThreadData* tdata) {
	int nThreads = params->nThreads;
	
	streamBlockBoundaries(flow,lc,params);

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
	streamBlockInterior(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	(*(bcdata->southFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX,bcdata->southBC[0], bcdata->southBC[1]);
	(*(bcdata->northFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX,bcdata->northBC[0], bcdata->northBC[1]);
	updateRho(tdata->flow,tdata->lc,tdata->startX,tdata->endX);		
	updateU(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	collide(tdata->flow,tdata->lc,tdata->params,tdata->startX,tdata->endX);
	return NULL;
}

void* updateFirstBlock(void* tdata_void) {
	ThreadData* tdata = (ThreadData*) tdata_void;
	BoundaryData* bcdata = tdata->bcdata;
	streamFirstBlockInterior(tdata->flow,tdata->lc,tdata->endX);
	(*(bcdata->southFun))(tdata->flow,tdata->lc,tdata->startX+1,tdata->endX,bcdata->southBC[0], bcdata->southBC[1]);
	(*(bcdata->northFun))(tdata->flow,tdata->lc,tdata->startX+1,tdata->endX,bcdata->northBC[0], bcdata->northBC[1]);
	(*(bcdata->westFun))(tdata->flow,tdata->lc,bcdata->westBC[0], bcdata->westBC[1]);	
	updateRho(tdata->flow,tdata->lc,tdata->startX+1,tdata->endX);		
	updateU(tdata->flow,tdata->lc,tdata->startX+1,tdata->endX);
	collide(tdata->flow,tdata->lc,tdata->params,tdata->startX,tdata->endX);
	return NULL;
}

void* updateLastBlock(void* tdata_void) {
	ThreadData* tdata = (ThreadData*) tdata_void;
	BoundaryData* bcdata = tdata->bcdata;
	streamLastBlockInterior(tdata->flow,tdata->lc,tdata->startX);
	(*(bcdata->southFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX-1,bcdata->southBC[0], bcdata->southBC[1]);
	(*(bcdata->northFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX-1,bcdata->northBC[0], bcdata->northBC[1]);
	(*(bcdata->eastFun))(tdata->flow,tdata->lc,bcdata->eastBC[0], bcdata->eastBC[1]);
	updateRho(tdata->flow,tdata->lc,tdata->startX,tdata->endX-1);		
	updateU(tdata->flow,tdata->lc,tdata->startX,tdata->endX-1);
	collide(tdata->flow,tdata->lc,tdata->params,tdata->startX,tdata->endX);
	return NULL;
}

void updateRho(FlowData* flow, LatticeConsts* lc, int startX, int endX) {
	double sum;
	int nx,ny;
	
	nx =lc->nx;
	ny = lc->ny;
	for (int i = startX; i <= endX; i++) {
		for (int j = 1; j < ny-1; j++){
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
		for (j = 1; j < ny-1; j++){
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
	int i,j,k,l,nx,ny,nxny,*ex,*ey,*northShiftArray,*southShiftArray;
	double *fIn,*fOut;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	ex = lc->exI;
	ey = lc->eyI;
	fIn = flow->fIn;
	fOut = flow->fOut;
	northShiftArray = lc->northShiftArray;
	southShiftArray = lc->southShiftArray;
	
	//Interior
	for (i = startX+2; i <= endX; i++) {
		for (j = 1; j < ny-1; j++) {
			for (k = 0; k < 9; k++) {
				fIn[nxny*k+ny*i+j] = fOut[k*nxny+(i-ex[k])*ny +j-ey[k]];
			}
		}
		//North/South boundary
		for (l = 0; l < 6; l++) {
			k = northShiftArray[l];
			fIn[nxny*k+ny*i+ny-1] = fOut[k*nxny+(i-ex[k])*ny +ny-1-ey[k]];
			k = southShiftArray[l];
			fIn[nxny*k+ny*i] = fOut[k*nxny+(i-ex[k])*ny-ey[k]];
		}
	}		
}

void streamFirstBlockInterior(FlowData* flow,LatticeConsts* lc,int endX) {
	int i,j,k,l,nx,ny,nxny,*ex,*ey,*northShiftArray,*southShiftArray,*westShiftArray;
	double *fIn,*fOut;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	ex = lc->exI;
	ey = lc->eyI;
	fIn = flow->fIn;
	fOut = flow->fOut;
	northShiftArray = lc->northShiftArray;
	southShiftArray = lc->southShiftArray;
	westShiftArray = lc->westShiftArray;
	
	//Interior
	for (i = 1; i <= endX; i++) {
		for (j = 1; j < ny-1; j++) {
			for (k = 0; k < 9; k++) {
				fIn[nxny*k+ny*i+j] = fOut[k*nxny+(i-ex[k])*ny +j-ey[k]];
			}
		}
		//North/South nodes
		for (l = 0; l < 6; l++) {
			k = northShiftArray[l];
			fIn[nxny*k+ny*i+ny-1] = fOut[k*nxny+(i-ex[k])*ny +ny-1-ey[k]];
			k = southShiftArray[l];
			fIn[nxny*k+ny*i] = fOut[k*nxny+(i-ex[k])*ny-ey[k]];
		}
	}
	
	//West nodes
	for (j = 1; j < ny-1; j++) {
		for (l = 0; l < 6; l++) {
			k = westShiftArray[l];
			fIn[nxny*k+j] = fOut[k*nxny-ex[k]*ny+j-ey[k]];
		}
	}
}

void streamLastBlockInterior(FlowData* flow,LatticeConsts* lc, int startX) {
	int i,j,k,l,nx,ny,nxny,*ex,*ey,*northShiftArray,*southShiftArray,*eastShiftArray;
	double *fIn,*fOut;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	ex = lc->exI;
	ey = lc->eyI;
	fIn = flow->fIn;
	fOut = flow->fOut;
	northShiftArray = lc->northShiftArray;
	southShiftArray = lc->southShiftArray;
	eastShiftArray = lc->eastShiftArray;
	
	//Interior
	for (i = startX; i < nx-1; i++) {
		for (j = 1; j < ny-1; j++) {
			for (k = 0; k < 9; k++) {
				fIn[nxny*k+ny*i+j] = fOut[k*nxny+(i-ex[k])*ny +j-ey[k]];
			}
		}
		//North/South nodes
		for (l = 0; l < 6; l++) {
			k = northShiftArray[l];
			fIn[nxny*k+ny*i+ny-1] = fOut[k*nxny+(i-ex[k])*ny +ny-1-ey[k]];
			k = southShiftArray[l];
			fIn[nxny*k+ny*i] = fOut[k*nxny+(i-ex[k])*ny-ey[k]];
		}
	}
	
	//East nodes
	for (j = 1; j < ny-1; j++) {
		for (l = 0; l < 6; l++) {
			k = eastShiftArray[l];
			fIn[nxny*k+ny*(nx-1)+j] = fOut[k*nxny+(nx-1-ex[k])*ny +j-ey[k]];
		}
	}
}

void streamBlockBoundaries(FlowData* flow,LatticeConsts* lc,SimParams* params){
	int i,j,k,l,nx,ny,nxny,*ex,*ey,blockSize,nThreads,*northShiftArray,*southShiftArray;
	double *fIn,*fOut;
	
	nx = lc->nx;
	ny = lc->ny;
	nxny = nx*ny;
	ex = lc->exI;
	ey = lc->eyI;
	northShiftArray = lc->northShiftArray;
	southShiftArray = lc->southShiftArray;
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
		//North/South boundary
		for (l = 0; l < 6; l++) {
			k = northShiftArray[l];
			fIn[nxny*k+ny*i+ny-1] = fOut[k*nxny+(i-ex[k])*ny +ny-1-ey[k]];
			fIn[nxny*k+ny*(i+1)+ny-1] = fOut[k*nxny+(i+1-ex[k])*ny +ny-1-ey[k]];
			k = southShiftArray[l];
			fIn[nxny*k+ny*i] = fOut[k*nxny+(i-ex[k])*ny-ey[k]];
			fIn[nxny*k+ny*(i+1)] = fOut[k*nxny+(i+1-ex[k])*ny-ey[k]];
		}
	}
}