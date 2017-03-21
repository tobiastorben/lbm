#include "initialization.h"

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

void initUx(LatticeConsts* lc, FlowData* flow, SimParams* params) {
	int i,j,nx,ny;
	double *ux;
	
	nx = lc->nx;
	ny = lc->ny;	
	ux = (double*) malloc(nx*ny * sizeof(double));
	for (i = 0; i < lc->nx; i++){
		for (j = 0; j < lc->ny; j++) {
			ux[(lc->ny)*i + j] = params->startVelX;
		}
	}
flow->ux = ux;	
}

void initUy(LatticeConsts* lc, FlowData* flow, SimParams* params) {
	int i,j,nx,ny;
	double *uy;
	
	nx = lc->nx;
	ny = lc->ny;	
	uy = (double*) malloc(nx*ny * sizeof(double));
	for (i = 0; i < lc->nx; i++){
		for (j = 0; j < lc->ny; j++) {
			uy[(lc->ny)*i + j] = params->startVelY;
		}
	}
flow->uy = uy;	
}


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

void nonDimensionalize(LatticeConsts* lc, SimParams* params, BoundaryData* bcdata) {
	double width,Re,t0,dx,dt,nu,u0,scale,rhoPhys;
	
	scale = (params->dtPhys)/(params->dxPhys);
	rhoPhys = params->rhoPhys;
	width = (params->dxPhys)*(lc->ny);
	t0 = width/(params->uRef);
	Re = width*(params->uRef)/(params->nuPhys);
	dx = 1.0/(lc->ny);
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
	printf("Reynolds number (based on domain length in Y-direction): %.1f\n", Re);
	printf("Lattice Mach number (based on refrance velocity): %.3f\n", (u0)*sqrt(3.0));
	printf("Speed/Accuracy ratio: %.2f\n", dt/(dx*dx));
	printf("LB dt: %.3e\n", dt);
	printf("Relaxation time: %.2f\n\n", params->tau);
	
	return;
}

void initialize(FlowData* flow, SimParams* params, LatticeConsts* lc, ThreadData** tdata, PrintData* pdata, BoundaryData* bcdata, char* inPath) {
	int nx,ny,i,nThreads;
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
	initUx(lc,flow,params);
	initUy(lc,flow,params);
	flow->fIn = (double*) malloc(nx*ny*9* sizeof(double));
	initFOut(lc,flow);
	
	//Initialize thread data
	nThreads = params->nThreads;
	threads = (pthread_t*) malloc(nThreads*sizeof(pthread_t));	
	*tdata = (ThreadData*) malloc(nThreads*sizeof(ThreadData));
	for (i = 0; i < nThreads; i++){
		(*tdata)[i].thread = threads[i];
		(*tdata)[i].params = params;
		(*tdata)[i].lc = lc;
		(*tdata)[i].flow = flow;
		(*tdata)[i].bcdata = bcdata;
		(*tdata)[i].startX = i*(((float) nx)/nThreads);
		(*tdata)[i].endX = (i+1)*(((float) nx)/nThreads)-1;
	}

	//Initialize print data
	pdata->uxCpy = (double*) malloc(nx*ny*sizeof(double));
	pdata->uyCpy = (double*) malloc(nx*ny*sizeof(double));
	pdata->rhoCpy = (double*) malloc(nx*ny*sizeof(double));
	pdata->nx = nx;
	pdata->ny = ny;
	pdata->params = params;
	
}

void setLatticeConstants(LatticeConsts* lc) {
	double w[] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};//Weights
	double ex[] = {0,1.0,0,-1.0,0,1.0,-1.0,-1.0,1.0};//x component of lattice vectors 
	double ey[] = {0,0,1.0,0,-1.0,1.0,1.0,-1.0,-1.0};//y component of lattice vectors
	int exI[] = {0,1,0,-1,0,1,-1,-1,1};//int version for indexing
	int eyI[] = {0,0,1,0,-1,1,1,-1,-1};//int version for indexing
	int opposite[] = {0,3,4,1,2,7,8,5,6};//Index of opposite vector in lattice
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