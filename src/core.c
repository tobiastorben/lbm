#include "core.h"

//------------------------------------------------------------------------------
//LBM    	Function: step
//------------------------------------------------------------------------------
//PURPOSE:	Progress the simulation one time step. The domain is divided into
//			vertical strips, called blocks. The calculations for each of these are
//			executed in parallell. It also alls streamBlockBoundaries and bounce,
//			which are not computed in paralell, but requires very little CPU time.
//USAGE:	step(flow,lc,params,tdata)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			params	SimParams* 		The parameters of the simulation		
//			tdata	ThreadData*		Array of structs. Each struct contains one
//									pthread and the data it needs to execute.
//.............................................................................
//CALLS:	streamBlockBoundaries	Stream two columns between blocks
//			bounce					Bounce-back collision with obstacle
//			updateFirstBlock		Does the calculations for the west block
//			updateLastBlock			Does the calculations for the east block
//			updateBlock				Does the calculations for remaining blocks
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************

void step(FlowData* flow, LatticeConsts* lc, SimParams* params, ThreadData* tdata) {
	int nThreads = params->nThreads;
	
	//Stream nodes between blocks
	streamBlockBoundaries(flow,lc,params);
	
	//Launch east and west solver threads
	pthread_create(&(tdata[0].thread),NULL,updateFirstBlock,&tdata[0]);
	pthread_create(&(tdata[nThreads-1].thread),NULL,updateLastBlock,&tdata[nThreads-1]);
	
	//Launch other solver threads
	for (int i = 1; i < nThreads-1; i++) {
		pthread_create(&(tdata[i].thread),NULL,updateBlock,&tdata[i]);
	}
	
	//Synchronize all threads
	for (int i = 0; i < nThreads; i++) {
		pthread_join(tdata[i].thread,NULL);
	}
	
	//Bounce-back collision with obstacle
	bounce(flow,lc,params);
}

//------------------------------------------------------------------------------
//LBM    	Function: updateBlock
//------------------------------------------------------------------------------
//PURPOSE:	Calls the necessary functions to do all the calculations for one block
//USAGE:	updateBlock(tdata_void)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			tdata_void	void*		Void pointer to tdata for current block
//.............................................................................
//CALLS:	streamBlockInterior		Propagate distribution function to adjecent
//									nodes
//			southFun				Applies BCs on south boundary
//			northFun				Applies BCs on north boundary
//			updateRho				Calculates density from new distribution
//			updateU					Calculates velocity from new distribution
//			collide					Perform the collision step, by the BGK method
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void* updateBlock(void* tdata_void) {
	ThreadData* tdata = (ThreadData*) tdata_void;
	BoundaryData* bcdata = tdata->bcdata;
	streamBlockInterior(tdata->flow,tdata->lc,tdata->startX,tdata->endX);//Streaming
	(*(bcdata->southFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX,bcdata->southBC[0], bcdata->southBC[1]);//South BC
	(*(bcdata->northFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX,bcdata->northBC[0], bcdata->northBC[1]);//North BC
	updateRho(tdata->flow,tdata->lc,tdata->startX,tdata->endX);//Update the density field	
	updateU(tdata->flow,tdata->lc,tdata->startX,tdata->endX);//Update the velocity field
	collide(tdata->flow,tdata->lc,tdata->params,tdata->startX,tdata->endX);//Collision step
	return NULL;
}

//------------------------------------------------------------------------------
//LBM    	Function: updateFirstBlock
//------------------------------------------------------------------------------
//PURPOSE:	Calls the necessary functions to to all the calculations for west block
//USAGE:	updateFirstBlock(tdata_void)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			tdata_void	void*			Void pointer to tdata for current block
//.............................................................................
//CALLS:	streamFirstBlockInterior	Propagate distribution function to adjecent
//										nodes for west block.
//			southFun					Applies BCs on south boundary
//			northFun					Applies BCs on north boundary
//			westFun						Applies BCs on west boundary
//			updateRho					Calculates density from new distribution
//			updateU						Calculates velocity from new distribution
//			collide						Perform the collision step, by the BGK method
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void* updateFirstBlock(void* tdata_void) {
	ThreadData* tdata = (ThreadData*) tdata_void;
	BoundaryData* bcdata = tdata->bcdata;
	streamFirstBlockInterior(tdata->flow,tdata->lc,tdata->endX);//Streaming west block
	(*(bcdata->southFun))(tdata->flow,tdata->lc,tdata->startX+1,tdata->endX,bcdata->southBC[0], bcdata->southBC[1]);
	(*(bcdata->northFun))(tdata->flow,tdata->lc,tdata->startX+1,tdata->endX,bcdata->northBC[0], bcdata->northBC[1]);
	(*(bcdata->westFun))(tdata->flow,tdata->lc,bcdata->westBC[0], bcdata->westBC[1]);//West BC
	updateRho(tdata->flow,tdata->lc,tdata->startX+1,tdata->endX);		
	updateU(tdata->flow,tdata->lc,tdata->startX+1,tdata->endX);
	collide(tdata->flow,tdata->lc,tdata->params,tdata->startX,tdata->endX);
	return NULL;
}

//------------------------------------------------------------------------------
//LBM    	Function: updateLastBlock
//------------------------------------------------------------------------------
//PURPOSE:	Calls the necessary functions to to all the calculations for east block
//USAGE:	updateLastBlock(tdata_void)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			tdata_void	void*			Void pointer to tdata for current block
//.............................................................................
//CALLS:	streamLastBlockInterior		Propagate distribution function to adjecent
//										nodes for west block.
//			southFun					Applies BCs on south boundary
//			northFun					Applies BCs on north boundary
//			eastFun						Applies BCs on east boundary
//			updateRho					Calculates density from new distribution
//			updateU						Calculates velocity from new distribution
//			collide						Perform the collision step, by the BGK method
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void* updateLastBlock(void* tdata_void) {
	ThreadData* tdata = (ThreadData*) tdata_void;
	BoundaryData* bcdata = tdata->bcdata;
	streamLastBlockInterior(tdata->flow,tdata->lc,tdata->startX);//Streaming east block
	(*(bcdata->southFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX-1,bcdata->southBC[0], bcdata->southBC[1]);
	(*(bcdata->northFun))(tdata->flow,tdata->lc,tdata->startX,tdata->endX-1,bcdata->northBC[0], bcdata->northBC[1]);
	(*(bcdata->eastFun))(tdata->flow,tdata->lc,bcdata->eastBC[0], bcdata->eastBC[1]);//East BC
	updateRho(tdata->flow,tdata->lc,tdata->startX,tdata->endX-1);		
	updateU(tdata->flow,tdata->lc,tdata->startX,tdata->endX-1);
	collide(tdata->flow,tdata->lc,tdata->params,tdata->startX,tdata->endX);
	return NULL;
}

//------------------------------------------------------------------------------
//LBM    	Function: updateRho
//------------------------------------------------------------------------------
//PURPOSE:	Calculates the density at each node, by summing the 9 distribution
//			functions.
//USAGE:	updateRho(flow,lc,startX,endX)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			startX	int				Index of first block column
//			endX	int				Index of last block column
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void updateRho(FlowData* flow, LatticeConsts* lc, int startX, int endX) {
	double sum;
	int nx,ny;
	
	nx =lc->nx;
	ny = lc->ny;
	
	//Update density field
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

//------------------------------------------------------------------------------
//LBM    	Function: updateU
//------------------------------------------------------------------------------
//PURPOSE:	Calculates the x and y velocity at each node, as a weighted sum of
//			the 9 distribution functions, with the lattice vector in the the
//			respective directions as weights. This sum is then divided by the
//			density
//USAGE:	updateU(flow,lc,startX,endX)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			startX	int				Index of first block column
//			endX	int				Index of last block column
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
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
	
	//Update velocity field
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

//------------------------------------------------------------------------------
//LBM    	Function: collide
//------------------------------------------------------------------------------
//PURPOSE:	Performs the collision step. The collision model used is the BGK model.
//			It relaxes the current distrubution function, towards the local
//			equilibrium distribution given by the Maxwell distribution.
//			This operation can be viewed as a Successive Overrelaxation (SOR),
//			where omega is equal to 1/tau. This therefore poses an absolute
//			limit for the relaxtion time, tau, to be larger than 0.5.
//USAGE:	collide(flow,lc,startX,endX)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			params	SimParams* 		The parameters of the simulation		
//			startX	int				Index of first block column
//			endX	int				Index of last block column
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
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
			uSq = 1.5*(uxIJ*uxIJ+uyIJ*uyIJ);//Factorize loop invariant calculations
			for (k = 0; k < 9; k++){
				u = 3.0*(ex[k]*uxIJ + ey[k]*uyIJ);
				fEq = rhoIJ*(w[k])*(1.0+u+0.5*u*u-uSq);//Equilibrium distribution
				fOut[nxny*k + ny*i + j] = fIn[nxny*k + ny*i + j]-(fIn[nxny*k + ny*i + j]-fEq)/tau;//Relax towards eq. distribution
			}
		}
	}			 
}

//------------------------------------------------------------------------------
//LBM    	Function: streamBlockInterior
//------------------------------------------------------------------------------
//PURPOSE:	Propagates the outgoing distributions to adjecent nodes. This reduces
//			to shifting fOut[i] in the direction of (ex[i],ey[i]). Only the
//			interior of the block is streamed here, as the boundaries requiers
//			special treatment since they have to access data outside the block.
//USAGE:	streamBlockInterior(flow,lc,startX,endX)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			startX	int				Index of first block column
//			endX	int				Index of last block column
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
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
//------------------------------------------------------------------------------
//LBM    	Function: streamFirstBlockInterior
//------------------------------------------------------------------------------
//PURPOSE:	Streams the west block. This is equivalent to streamBlockInterior,
//			exept for the west nodes. Here, only the distribution from the east
//			neighbours are streamed. The incoming distributions for the
//			others are treaded by the boundary conditions.
//USAGE:	streamFirstBlockInterior(flow,lc,endX)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			endX	int				Index of last block column
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
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

//------------------------------------------------------------------------------
//LBM    	Function: streamLastBlockInterior
//------------------------------------------------------------------------------
//PURPOSE:	Streams the east block. This is equivalent to streamBlockInterior,
//			exept for the east nodes. Here, only the distribution from the west
//			neighbours are streamed. The incoming distributions for the
//			others are treaded by the boundary conditions.
//USAGE:	streamFirstBlockInterior(flow,lc,endX)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			startX	int				Index of first block column
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
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
	for (i = startX+2; i < nx-1; i++) {
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

//------------------------------------------------------------------------------
//LBM    	Function: streamBlockBoundaries
//------------------------------------------------------------------------------
//PURPOSE:	All the threads are synchronized, and then a 2 node wide column
//			between every block is streamed. This is to prevent race conditions
//USAGE:	streamBlockBoundaries(flow,lc,params)
//ARGUMENTS:
//			Name 	 	Type     			Description
//.............................................................................
//			flow	FlowData*		The field variables of the flow
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			params	SimParams*		The parameters of the simulation		
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void streamBlockBoundaries(FlowData* flow,LatticeConsts* lc,SimParams* params){
	int i,j,k,l,nx,ny,nxny,*ex,*ey,blockSize,*northShiftArray,*southShiftArray;
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
	blockSize = params->blockSize;
	
	//Nodes between blocks, exept top and bottom
	for (i = blockSize; i < nx-1; i+=blockSize) {
		for (j = 1; j < ny-1; j++) {
			for (k = 0; k < 9; k++) {
				fIn[nxny*k+ny*i+j] = fOut[k*nxny+(i-ex[k])*ny +j-ey[k]];
				fIn[nxny*k+ny*(i+1)+j] = fOut[k*nxny+(i+1-ex[k])*ny +j-ey[k]];
			}
		}
		//Top and bottom
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