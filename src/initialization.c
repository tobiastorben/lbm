#include "initialization.h"

//Reads geometry of bounce back cells from a file, stores it in a vector
//of indices
int* readObstacleData(int nBBcells) {
	FILE* fp = fopen("obstacle.mask","r");
	char* line = malloc(nBBcells*100);		
	char* token;
	int* bbCells = malloc(2*nBBcells*sizeof(int));
	int i = 0;
	int j = 0;
	int k = 0;
	while (!feof(fp)){
		fgets(line,nBBcells*100,fp);
		token = strtok(line, ",");		
		int j = 0;
		while( token != NULL ){ 	    
		if (atoi(token)) {
			bbCells[k] = i;
			bbCells[1*nBBcells+k] = j;
			k++;
		}
		token = strtok(NULL, ",");
		j++;
	    }
	i++;
	}
	return bbCells;
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