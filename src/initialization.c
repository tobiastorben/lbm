#include "initialization.h"

//This file contains all function to initialize the simulation, before it
//enters the main loop.

//------------------------------------------------------------------------------
//LBM    Function: initialize 
//------------------------------------------------------------------------------
//PURPOSE:	Initializes everything before it enters the main loop. This work
//			distributed to subfunctions, see call list.
//USAGE:	initialize(flow,params,lc,tdata,pdata,bcdata,inPath)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			params	SimParams* 		The parameters of the simulation
//			flow	FlowData*		The field variables of the flow
//			bcdata	BoundaryData*	The boundary conditions of the flow
//			pdata	PrintData*		All the data needed for the print thread
//			tdata	TreadData**		Pointer to array of structs. Each struct contains
//									one pthread and the data it needs to execute.
//			inPath	char*			The path to the input file		
//.............................................................................
//CALLS:	
//			setLatticeConstants		Set constants of D2Q9 lattice
//			parseInput				Parse input file
//			readObstacle			Read obstacle mask into matrix, and use it to
//									determine grid size
//			mapObstacleCells		map obstacle cells to their index in the grid
//			nonDimensionalize		Convert to lattice units, and set relaxation time
//			initRho					Allocate and initialize density matrix
//			initU					Allocate and initialize velocity field
//			initFOut				Allocate and initialize FOut 3D matrix

//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void initialize(FlowData* flow, SimParams* params, LatticeConsts* lc, ThreadData** tdata, PrintData* pdata, BoundaryData* bcdata, char* inPath) {
	int nx,ny,i,nThreads,rest;
	pthread_t *threads;
	
	printf("Initializing...\n\n");
	
	//D2Q9 lattice constants
	setLatticeConstants(lc);
	
	if (parseInput(inPath,params,bcdata)) {
		printf("Invalid simulation parameters. Exiting.");
		exit(1);
	}

	//Read geometry and set dimensions
	readObstacle(lc,params);//Read obstacle data from file
	mapObstacleCells(lc,params);//Array of indices to bounce back nodes
	
	//Cast to non-dimensional form
	nonDimensionalize(lc, params,bcdata);
			
	//ICs: Initialze flow as Poiseuille flow
	nx = lc->nx;
	ny = lc->ny;
	initRho(lc,flow);
	initU(lc,flow,params);
	flow->fIn = (double*) malloc(nx*ny*9* sizeof(double));
	initFOut(lc,flow);
	
	//Initialize thread data
	nThreads = params->nThreads;
	threads = (pthread_t*) malloc(nThreads*sizeof(pthread_t));	
	rest = nx % nThreads;
	*tdata = (ThreadData*) malloc(nThreads*sizeof(ThreadData));
	for (i = 0; i < nThreads; i++){
		(*tdata)[i].thread = threads[i];
		(*tdata)[i].params = params;
		(*tdata)[i].lc = lc;
		(*tdata)[i].flow = flow;
		(*tdata)[i].bcdata = bcdata;
		(*tdata)[i].startX = i*(nx-rest)/nThreads;
		(*tdata)[i].endX = (i+1)*((nx-rest)/nThreads)-1;
	}
	(*tdata)[nThreads-1].endX = nx-1;
	params->blockSize = (*tdata)[0].endX+1;

	//Initialize print data
	pdata->uxCpy = (double*) malloc(nx*ny*sizeof(double));
	pdata->uyCpy = (double*) malloc(nx*ny*sizeof(double));
	pdata->rhoCpy = (double*) malloc(nx*ny*sizeof(double));
	pdata->nx = nx;
	pdata->ny = ny;
	pdata->params = params;
	
}

//------------------------------------------------------------------------------
//LBM    Function: readObstacle 
//------------------------------------------------------------------------------
//PURPOSE:	Reads the obstacle file, and stores it in a boolean mask (bbCellMat).
//			 It also sets the mesh size (nx,ny) based on the number off rows and
//			columns of the obstacle.
//USAGE:	readObstacle(lc,params)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			params	SimParams* 		The parameters of the simulation		
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void readObstacle(LatticeConsts* lc, SimParams* params) {
	FILE* fp = fopen(params->obstaclePath,"r");
	char* line = malloc(10000);		
	char* token;
	int i,j,nBBCells,*bbCellMat;
	nBBCells=0;i=0;j = 0;
	
	//Traverse file once to determine the dimensions of the grid
	while (!feof(fp)){
		fgets(line,10000,fp);
		token = strtok(line, ",");		
		j = 0;
		while( token != NULL ){ 	    
		token = strtok(NULL, ",");
		j++;
	    }
		i++;
	}	
	lc->nx = i;
	lc->ny = j;
	
	//Read data into matrix
	bbCellMat = calloc((lc->nx)*(lc->ny),sizeof(int));
	fp = fopen(params->obstaclePath,"r");
	i = 0;
	while (!feof(fp)){
		fgets(line,10000,fp);
		token = strtok(line, ",");		
		j = 0;
		while( token != NULL ){ 	    
		if (atoi(token) == 1) {
			bbCellMat[i*(lc->ny) + j] = 1;
			nBBCells++;
		}
		token = strtok(NULL, ",");
		j++;
	    }
		i++;
	}
	params->nBBcells = nBBCells;
	params->bbCellMat = bbCellMat;
	free(line);

}

//------------------------------------------------------------------------------
//LBM    Function: mapObstacleCells
//------------------------------------------------------------------------------
//PURPOSE: 	Traverses the boolean mask bbCellMat, and creates a 2-row matrix,
//			bbCells, where each column contains the x and y index of a bounce back
//			cell. This allows for fast execution of the bounce method in the core
//			solver.
//USAGE:	mapObstacleCells(lc,params)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			params	SimParams* 		The parameters of the simulation		
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void mapObstacleCells(LatticeConsts* lc, SimParams* params) {
	int i,j,k;
	int* bbCells = (int*) malloc(2*(params->nBBcells)*sizeof(int));
	k = 0;
	for (i = 0; i < lc->nx; i++){
		for (j = 0; j < lc->ny; j++){    
			if (params->bbCellMat[i*(lc->ny) + j]) {
				bbCells[k] = i;
				bbCells[params->nBBcells+k] = j;
				k++;
			}
	    }
	}
	params->bbCells = bbCells;
}

//------------------------------------------------------------------------------
//LBM    Function: initU
//------------------------------------------------------------------------------
//PURPOSE:	Allocate memory for, and initialize the velocity field. The velocity
//			for all obtsacle cells are set to zero.
//USAGE:	initU(lc,flow,params)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			flow	FlowData*		The field variables of the flow
//			params	SimParams* 		The parameters of the simulation		
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void initU(LatticeConsts* lc, FlowData* flow, SimParams* params) {
	int i,j,nx,ny;
	double *uy,*ux;
	
	nx = lc->nx;
	ny = lc->ny;	
	uy = (double*) calloc(nx*ny,sizeof(double));
	ux = (double*) calloc(nx*ny,sizeof(double));
	
	for (i = 0; i < lc->nx; i++){
		for (j = 0; j < lc->ny; j++) {
			if (!(params->bbCellMat[i*(lc->ny) + j])){
			uy[(lc->ny)*i + j] = params->startVelY;
			ux[(lc->ny)*i + j] = params->startVelX;
			}
		}
	}
	flow->uy = uy;	
	flow->ux = ux;
}

//------------------------------------------------------------------------------
//LBM    Function: initFOut
//------------------------------------------------------------------------------
//PURPOSE:	Allocate memory for the fOut 3D matrix, and initialize it to the
//			equilibrium distribution for the initial velocity field.
//USAGE:	initFOut(lc,flow)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			flow	FlowData*		The field variables of the flow
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void initFOut(LatticeConsts* lc, FlowData* flow) {
	int i,j,k,nx,ny;
	double u, *fOut,*ux,*uy;
	
	nx = lc->nx;
	ny = lc->ny;
	ux = flow->ux;
	uy = flow->uy;
	fOut = (double*) malloc(nx*ny*9* sizeof(double));
	
	//Initialize particle distribution as equilibrium distrubution for
	//still fluid
	for (k = 0; k < 9; k++){
		for (i = 0; i < lc->nx; i++){
			for (j = 0; j < lc->ny; j++){
			u = 3.0*((lc->ex[k])*ux[ny*i + j] + (lc->ey[k])*uy[ny*i + j]);
			fOut[nx*ny*k + ny*i + j] = (flow->rho[ny*i + j])*(lc->w[k])*(1.0+u+0.5*u*u-1.5*(ux[ny*i + j]*ux[ny*i + j]+uy[ny*i + j]*uy[ny*i + j]));
			}
		}
	}
	flow->fOut = fOut;
}
//------------------------------------------------------------------------------
//LBM    Function: initRho
//------------------------------------------------------------------------------
//PURPOSE:	Allocate memory for the rho matrix, and initialize to unity.
//USAGE:	initUx(lc,flow)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			flow	FlowData*		The field variables of the flow
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void initRho(LatticeConsts* lc, FlowData* flow) {
	int i,j,nx,ny;
	double* rho;

	nx = lc->nx;
	ny = lc->ny;
	rho = (double*) malloc(nx*ny* sizeof(double));
	
	for (i = 0; i < nx; i++){
			for (j = 0; j < ny; j++){
			rho[ny*i + j] = 1.0;
			}
		}
	flow->rho = rho;
}

//------------------------------------------------------------------------------
//LBM    Function: nonDimensionalize
//------------------------------------------------------------------------------
//PURPOSE:	Convert the physical units of the input file, to lattice units. All
//			velocities and BC's and IC's are scaled to lattice units. The
//			relaxation time, tau, is calculated based on the Reynolds number,
//			the time step and the spatial step. It also prints central parameters
//			of the simulation, which are useful for choosing time step and grid size
//			for achieving good accuracy and stability.
//
//USAGE:	nonDimensionalize(lc,params,bcdata)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//			params	SimParams* 		The parameters of the simulation
//			bcdata	BoundaryData*	The boundary conditions of the flow	
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void nonDimensionalize(LatticeConsts* lc, SimParams* params, BoundaryData* bcdata) {
	double width,Re,t0,dx,dt,nu,u0,scale,rhoPhys;
	
	scale = (params->dtPhys)/(params->dxPhys);
	rhoPhys = params->rhoPhys;
	width = (params->dxPhys)*((lc->ny)-1);
	t0 = width/(params->uRef);
	Re = width*(params->uRef)/(params->nuPhys);
	dx = 1.0/(lc->ny-1);
	dt = (params->dtPhys)/t0;
	u0 = dt/dx;
	nu = dt/(dx*dx*Re);
	params->tau = 3.0*nu + 0.5;
	params->startVelX = params->startVelX*scale;
	params->startVelY = params->startVelY*scale;
	
	if (bcdata->westBCType) {
		bcdata->westBC[0] = 1.0 + 3.0*bcdata->westBC[0]*scale*scale/rhoPhys; 
		bcdata->westBC[1] = bcdata->westBC[1]*scale;
	}
	
	else {
		bcdata->westBC[0] = bcdata->westBC[0]*scale;
		bcdata->westBC[1] = bcdata->westBC[1]*scale;
	}
	
	if (bcdata->northBCType) {
		bcdata->northBC[0] = 1.0 + 3.0*bcdata->northBC[0]*scale*scale/rhoPhys; 
		bcdata->northBC[1] = bcdata->northBC[1]*scale;
	}
	
	else {
		bcdata->northBC[0] = bcdata->northBC[0]*scale;
		bcdata->northBC[1] = bcdata->northBC[1]*scale;
	}
	
	if (bcdata->eastBCType) {
		bcdata->eastBC[0] = 1.0 + 3.0*bcdata->eastBC[0]*scale*scale/rhoPhys; 
		bcdata->eastBC[1] = bcdata->eastBC[1]*scale;
	}
	
	else {
		bcdata->eastBC[0] = bcdata->eastBC[0]*scale;
		bcdata->eastBC[1] = bcdata->eastBC[1]*scale;
	}
	
	if (bcdata->southBCType) {
		bcdata->southBC[0] = 1.0 + 3.0*bcdata->southBC[0]*scale*scale/rhoPhys; 
		bcdata->southBC[1] = bcdata->southBC[1]*scale;
	}
	
	else {
		bcdata->southBC[0] = bcdata->southBC[0]*scale;
		bcdata->southBC[1] = bcdata->southBC[1]*scale;
	}
	
	printf("Numerical simulation parameters:\n");
	printf("Reynolds number (based on domain length in Y-direction): %7f\n", Re);
	printf("Lattice Mach number (based on refrance velocity): %.3f\n", (u0)*sqrt(3.0));
	printf("Speed/Accuracy ratio: %.2f\n", dt/(dx*dx));
	printf("Relaxation time: %.2f\n\n", params->tau);
	
	return;
}

//------------------------------------------------------------------------------
//LBM    Function: setLatticeConstants
//------------------------------------------------------------------------------
//PURPOSE:	Sets all the lattice constants for the D2Q9 BGK model
//
//USAGE:	setLatticeConstants(lc)
//ARGUMENTS:
//			Name 	 Type     		Description
//.............................................................................
//			lc		LatticeConsts*	The constants of the D2Q9 lattice
//.............................................................................
//Author: Tobias Valentin Rye Torben
//Date/Version: 29.03.2017
//******************************************************************************
void setLatticeConstants(LatticeConsts* lc) {
	double w[] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};//Weights
	double ex[] = {0,1.0,0,-1.0,0,1.0,-1.0,-1.0,1.0};//x component of lattice vectors 
	double ey[] = {0,0,1.0,0,-1.0,1.0,1.0,-1.0,-1.0};//y component of lattice vectors
	int exI[] = {0,1,0,-1,0,1,-1,-1,1};//int version for indexing
	int eyI[] = {0,0,1,0,-1,1,1,-1,-1};//int version for indexing
	int opposite[] = {0,3,4,1,2,7,8,5,6};//Index of opposite vector in lattice
	
	//Index of the lattice velocities that are pointing inwards at each boyndary
	int westShiftArray[] = {0,2,3,4,6,7};
	int northShiftArray[] = {0,1,2,3,5,6};
	int eastShiftArray[] = {0,1,2,4,5,8};
	int southShiftArray[] = {0,1,3,4,7,8};
	
	memcpy(lc->w,w,sizeof(lc->w));
	memcpy(lc->ex,ex,sizeof(lc->ex));
	memcpy(lc->ey,ey,sizeof(lc->ey));
	memcpy(lc->exI,exI,sizeof(lc->exI));
	memcpy(lc->eyI,eyI,sizeof(lc->eyI));
	memcpy(lc->opposite,opposite,sizeof(lc->opposite));
	memcpy(lc->westShiftArray,westShiftArray,sizeof(lc->westShiftArray));
	memcpy(lc->northShiftArray,northShiftArray,sizeof(lc->northShiftArray));
	memcpy(lc->eastShiftArray,eastShiftArray,sizeof(lc->eastShiftArray));
	memcpy(lc->southShiftArray,southShiftArray,sizeof(lc->southShiftArray));
}