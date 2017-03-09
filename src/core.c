#include "core.h"
#include <stdio.h>

void step(FlowData* flow, LatticeConsts* lc, SimParams* params, ThreadData* tdata) {//TODO : Fix for 1 thread
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
	streamBlockInterior(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	updateRho(tdata->flow,tdata->lc,tdata->startX,tdata->endX);		
	updateU(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	collide(tdata->flow,tdata->lc,tdata->params,tdata->startX,tdata->endX);
	return NULL;
}

void* updateFirstBlock(void* tdata_void) {
	ThreadData* tdata = (ThreadData*) tdata_void;
	streamBlockInterior(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	updateRho(tdata->flow,tdata->lc,tdata->startX,tdata->endX);		
	updateU(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	inlet(tdata->flow,tdata->lc,tdata->params);//Poiseuille (Zou/He)		
	collide(tdata->flow,tdata->lc,tdata->params,tdata->startX,tdata->endX);
	return NULL;
}

void* updateLastBlock(void* tdata_void) {
	ThreadData* tdata = (ThreadData*) tdata_void;
	streamBlockInterior(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	updateRho(tdata->flow,tdata->lc,tdata->startX,tdata->endX);		
	updateU(tdata->flow,tdata->lc,tdata->startX,tdata->endX);
	outlet(tdata->flow,tdata->lc);//Constant pressure (Zou/He)
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
	
	for (int i = startX; i <= endX; i++) {
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
	for (i = startX+1; i <= endX; i++) {
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
			}
		}
	}
}