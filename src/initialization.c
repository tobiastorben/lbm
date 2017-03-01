#include "initialization.h"

void mapObstacleCells(LatticeConsts* lc, SimParams* params) {
	int i,j,k;
	int* bbCells = malloc(2*(params->nBBcells)*sizeof(int));
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
	int i,j;
	int nBBCells = 0;
	
	//Traverse file one to determine the dimensions of the grid
	i = 0;
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
	int* bbCellMat = calloc((lc->nx)*(lc->ny),sizeof(int));
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
	double y,H;
	
	nx = lc->nx;
	ny = lc->ny;
	H = (double) lc->ny-2;
	
	
	//Initialize ux as Poiseuille flow
	double* ux = malloc(nx*ny * sizeof(double));
	for (i = 0; i < lc->nx; i++){
		for (j = 0; j < lc->ny; j++) {
			y = j-0.5;//Change??
			ux[(lc->ny)*i + j] = (4.0*(params->u0)/(H*H))*(y*H-y*y);
		}
	}
flow->ux = ux;	
}

void initFin(LatticeConsts* lc, FlowData* flow) {
	int i,j,k,nx,ny;
	double u, *fIn,*ux,*uy;
	
	nx = lc->nx;
	ny = lc->ny;
	ux = flow->ux;
	uy = flow->uy;
	//Allocate three dimensional array. Indexed as: [f*nx*ny+(nx-1)*ny+ny]
	fIn = malloc(nx*ny*9* sizeof(double));
	
	//Initialize particle distribution as equilibrium distrubution for
	//for a Poiseuille flow
	for (k = 0; k < 9; k++){
		for (i = 0; i < lc->nx; i++){
			for (j = 0; j < lc->ny; j++){
			u = 3.0*((lc->ex[k])*ux[ny*i + j] + (lc->ey[k])*uy[ny*i + j]);
			fIn[nx*ny*k + ny*i + j] = flow->rho[ny*i + j]*(lc->w[k])*(1.0+u+0.5*u*u-1.5*(ux[ny*i + j]*ux[ny*i + j]+uy[ny*i + j]*uy[ny*i + j]));
			}
		}
	}
	flow->fIn = fIn;
}

void initRho(LatticeConsts* lc, FlowData* flow) {
	int i,j,nx,ny;
	double* rho;

	nx = lc->nx;
	ny = lc->ny;
	rho = malloc(nx*ny* sizeof(double));
	
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
	params->u0 = dt/dx;
	nu = dt/(dx*dx*Re);
	params->tau = 3.0*nu + 0.5;
	
	printf("Numerical simulation parameters:\n");
	printf("Reynolds number (based on channel width): %.1f\n", Re);
	printf("Lattice Mach number: %.3f\n", (params->u0)*sqrt(3.0));
	printf("Speed/Accuracy ratio: %.2f\n", dt/(dx*dx));
	printf("Relaxation time: %.2f\n\n", params->tau);
	
	return;
}

void initialize(FlowData* flow, SimParams* params, LatticeConsts* lc, char* inPath) {
	int nx,ny;
	
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
			
	//ICs: Initialze flow as Poiseuille flow
	nx = lc->nx;
	ny = lc->ny;
	initRho(lc,flow);
	initUx(lc,flow,params);
	initFin(lc,flow);
	flow->uy = calloc(nx*ny,sizeof(double));
	flow->fOut = malloc(nx*ny*9* sizeof(double));
	
	//Cast to non-dimensional form
	nonDimensionalize(lc, params);
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