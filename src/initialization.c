#include "initialization.h"

int* mapObstacleCells(int* bbCellMat,int nBBcells, int nx, int ny) {
	int i = 0;
	int j = 0;
	int k = 0;
	int* bbCells = malloc(2*nBBcells*sizeof(int));
	for (i = 0; i < nx; i++){
		for (j = 0; j < ny; j++){    
			if (bbCellMat[i*ny + j]) {
				bbCells[k] = i;
				bbCells[1*nBBcells+k] = j;
				k++;
			}
	    }
	}
	return bbCells;
}

int* readObstacle(int* nBBCells_p, int* nx_p, int* ny_p, char* path) {
	FILE* fp = fopen(path,"r");
	char* line = malloc(10000);		
	char* token;
	int i = 0;
	int j;
	int nBBCells = 0;
	
	//Traverse file one to determine the dimensions of the grid
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
	*nx_p = i;
	*ny_p = j;
	fclose(fp);
	
	//Read data into matrix
	int* bbCellMat = calloc((*nx_p)*(*ny_p),sizeof(int));
	fp = fopen(path,"r");
	i = 0;
	while (!feof(fp)){
		fgets(line,10000,fp);
		token = strtok(line, ",");		
		j = 0;
		while( token != NULL ){ 	    
		if (atoi(token) == 1) {
			bbCellMat[i*(*ny_p) + j] = 1;
			nBBCells++;
		}
		token = strtok(NULL, ",");
		j++;
	    }
		i++;
	}
	*nBBCells_p = nBBCells;
	free(line);
	return bbCellMat;
}

double* initUx(int nx,int ny,double uIn) {
	double H = (double) ny-2;
	int i,j;
	double y;
	
	//Initialize ux as Poiseuille flow
	double* ux = malloc(nx*ny * sizeof(double));
	for (i = 0; i < nx; i++){
		for (j = 0; j < ny; j++) {
			y = j-0.5;//Change??
			ux[ny*i + j] = (4.0*uIn/(H*H))*(y*H-y*y);
		}
	}
return ux;	
}

double* initUy(int nx,int ny) {	
	//Zero initialize uy
	double* uy = calloc(nx*ny,sizeof(double));
	return uy;
}

double* initFin(int nf,int nx, int ny, double* ex, double* ey,
				double* ux, double* uy, double* w, double* rho) {
	int i,j,k;
	double u;
	//Allocate three dimensional array. Indexed as: [f*nx*ny+(nx-1)*ny+ny]
	double* fIn = malloc(nx*ny*nf* sizeof(double));
	
	//Initialize particle distribution as equilibrium distrubution for
	//for a Poiseuille flow
	for (k = 0; k < nf; k++){
		for (i = 0; i < nx; i++){
			for (j = 0; j < ny; j++){
			u = 3.0*(ex[k]*ux[ny*i + j] + ey[k]*uy[ny*i + j]);
			fIn[nx*ny*k + ny*i + j] = rho[ny*i + j]*w[k]*(1.0+u+0.5*u*u-1.5*(ux[ny*i + j]*ux[ny*i + j]+uy[ny*i + j]*uy[ny*i + j]));
			}
		}
	}
	return fIn;
}

double* initFout(int nx,int ny, int nf) {
	double* fOut = malloc(nx*ny*nf* sizeof(double));
	return fOut;
}

double* initRho(int nx, int ny, double initialRho) {
	int i,j;
	double* rho = malloc(nx*ny* sizeof(double));
	
	for (i = 0; i < nx; i++){
			for (j = 0; j < ny; j++){
			rho[ny*i + j] = initialRho;
			}
		}
	return rho;
}

void nonDimensionalize(int nx, int ny, double dtPhys, double dxPhys, double nuPhys, double u0Phys, double* u0_p, double* tau_p) {
	double width,Re,t0,dx,dt,nu;
	
	width = dxPhys*ny;
	t0 = width/u0Phys;
	Re = width*u0Phys/nuPhys;
	dx = 1.0/ny;
	dt = dtPhys/t0;
	*u0_p = dt/dx;
	nu = dt/(dx*dx*Re);
	*tau_p = 3.0*nu + 0.5;
	
	printf("Numerical simulation parameters:\n");
	printf("Reynolds number (based on channel width): %.1f\n", Re);
	printf("Lattice Mach number: %.3f\n", *u0_p*sqrt(3.0));
	printf("Speed/Accuracy ratio: %.2f\n", dt/(dx*dx));
	printf("Relaxation time: %.2f\n\n", *tau_p);
	
	return;
}