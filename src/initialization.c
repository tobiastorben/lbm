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
	char* line = malloc(1000);		
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
	double y,H,*ux;
	
	nx = lc->nx;
	ny = lc->ny;
	H = (double) lc->ny-2;
		
	//Initialize ux as Poiseuille flow
	ux = (double*) malloc(nx*ny * sizeof(double));
	for (i = 0; i < lc->nx; i++){
		for (j = 0; j < lc->ny; j++) {
			y = j-0.5;
			ux[(lc->ny)*i + j] = 0.0;//(4.0*(params->u0)/(H*H))*(y*H-y*y);
		}
	}
flow->ux = ux;	
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
	//for a Poiseuille flow
	for (k = 0; k < 9; k++){
		for (i = 0; i < lc->nx; i++){
			for (j = 0; j < lc->ny; j++){
			//u = 3.0*((lc->ex[k])*ux[ny*i + j] + (lc->ey[k])*uy[ny*i + j]);
			//fOut[nx*ny*k + ny*i + j] = (flow->rho[ny*i + j])*(lc->w[k])*(1.0+u+0.5*u*u-1.5*(ux[ny*i + j]*ux[ny*i + j]+uy[ny*i + j]*uy[ny*i + j]));
			fOut[nx*ny*k + ny*i + j] = lc->w[k];
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

void nonDimensionalize(LatticeConsts* lc, SimParams* params) {
	double width,Re,t0,dx,dt,nu;
	
	width = (params->dxPhys)*(lc->ny);
	t0 = width/(params->u0Phys);
	Re = width*(params->u0Phys)/(params->nuPhys);
	dx = 1.0/(lc->ny);
	dt = (params->dtPhys)/t0;
	params->u0 = 0.05;//dt/dx;
	nu = dt/(dx*dx*Re);
	params->tau = 0.692;//= 3.0*nu + 0.5;
	
	printf("Numerical simulation parameters:\n");
	printf("Reynolds number (based on channel width): %.1f\n", Re);
	printf("Lattice Mach number: %.3f\n", (params->u0)*sqrt(3.0));
	printf("Speed/Accuracy ratio: %.2f\n", dt/(dx*dx));
	printf("Relaxation time: %.2f\n\n", params->tau);
	
	return;
}

void initialize(FlowData* flow, SimParams* params, LatticeConsts* lc, ThreadData** tdata, PrintData* pdata, char* inPath) {
	int nx,ny,i,nThreads;
	pthread_t *threads;
	
	printf("Initializing...\n\n");
	
	//D2Q9 lattice constants
	setLatticeConstants(lc);
	
	if (parseInput(inPath,params)) {
		printf("Invalid simulation parameters. Exiting.");
		exit(1);
	}

	//Read geometry and set dimensions
	readObstacle(lc,params);//Read obstacle data from file
	mapObstacleCells(lc,params);//Array of indices to bounce back nodes
	
	//Cast to non-dimensional form
	nonDimensionalize(lc, params);
			
	//ICs: Initialze flow as Poiseuille flow
	nx = lc->nx;
	ny = lc->ny;
	initRho(lc,flow);
	flow->uy = (double*) calloc(nx*ny,sizeof(double));
	initUx(lc,flow,params);
	initFOut(lc,flow);
	flow->fIn = (double*) malloc(nx*ny*9* sizeof(double));
	memcpy(flow->fIn,flow->fOut,nx*ny*9*sizeof(double));
	
	//Initialize thread data
	nThreads = params->nThreads;
	threads = (pthread_t*) malloc(nThreads*sizeof(pthread_t));	
	*tdata = (ThreadData*) malloc(nThreads*sizeof(ThreadData));
	for (i = 0; i < nThreads; i++){
		(*tdata)[i].thread = threads[i];
		(*tdata)[i].params = params;
		(*tdata)[i].lc = lc;
		(*tdata)[i].flow = flow;
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
	
	memcpy(lc->w,w,sizeof(lc->w));
	memcpy(lc->ex,ex,sizeof(lc->ex));
	memcpy(lc->ey,ey,sizeof(lc->ey));
	memcpy(lc->exI,exI,sizeof(lc->exI));
	memcpy(lc->eyI,eyI,sizeof(lc->eyI));
	memcpy(lc->opposite,opposite,sizeof(lc->opposite));
}