#include "utilities.h"


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

void printVecD(double* vec, int n) {
	for (int i = 0; i < n; i++){
		printf("%8f, ", vec[i]);
	}
		
}
void printVecI(int* vec, int n) {
	for (int i = 0; i < n; i++){
		printf("%8d, ", vec[i]);
	}	
}

void printMatI(int* mat, int n, int m) {
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		printf("%8d ", mat[m*i + j]);
		}
	printf("\n");
	}
	printf("\n");
}
void printMatD(double* mat, int n, int m) {
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		printf("%f ", mat[m*i + j]);
		}
	printf("\n");
	}
printf("\n");	
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

void updateRho(double* rho, double* fIn, int nx, int ny, int nf) {
	double sum;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++){
			sum = 0;
			for (int k = 0; k < nf; k++) {
				sum += fIn[nx*ny*k + ny*i + j];
			}
			rho[ny*i + j] = sum;
		}
	}
}

void updateU(double* ux, double* uy, double* fIn,
			 double* rho, double* ex, double* ey, int nx, int ny, int nf){
	double sumX, sumY;
	int nxny = nx*ny;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++){
			sumX = 0;
			sumY = 0;
			for (int k = 0; k < nf; k++) {
				sumX += fIn[nxny*k + ny*i + j]*ex[k];
				sumY += fIn[nxny*k + ny*i + j]*ey[k];
			}
			ux[ny*i + j] = sumX/rho[ny*i + j];
			uy[ny*i + j] = sumY/rho[ny*i + j];
		}
	} 				 
}

void inlet(double* ux, double* uy, double* rho, double* fIn, double uIn, int nx, int ny) {
	double y, sum1, sum2, rhoJ;
	double H = (double) ny-2.0;
	double coeff = 4.0*uIn/(H*H);
	for (int j = 1; j < (ny-1);j++) {
		y = j-0.5;//Change?? 
		//Poiseuille and Zou/He BC
		//Consider to change indexing of fIn
		ux[j] = coeff*y*(H-y);
		uy[j] = 0;
		sum1 = fIn[j] + fIn[2*nx*ny +j] + fIn[4*nx*ny+j];
		sum2 = fIn[3*nx*ny +j] + fIn[6*nx*ny + j] + fIn[7*nx*ny + j];
		rhoJ = (1.0/(1.0-ux[j]))*(sum1 + 2.0*sum2);
		fIn[1*nx*ny+j] = fIn[3*nx*ny+j] + (2.0/3)*rhoJ*ux[j];
		fIn[5*nx*ny+j] = fIn[7*nx*ny+j] + 0.5*(fIn[4*nx*ny+j]-fIn[2*nx*ny+j])
					   +(1.0/6)*(rhoJ*ux[j]);
		fIn[8*nx*ny+j] = fIn[6*nx*ny+j] + 0.5*(fIn[2*nx*ny+j]-fIn[4*nx*ny+j])
						+(1.0/6)*(rhoJ*ux[j]);
		rho[j] = rhoJ;
	}
}

void outlet(double* ux, double* uy, double* rho, double* fIn, int nx, int ny) {
	double sum1, sum2;
	for (int j = 1; j < (ny-1);j++) {
		sum1 = fIn[(nx-1)*ny+j] + fIn[2*nx*ny+(nx-1)*ny+j] + fIn[4*nx*ny+(nx-1)*ny+j];
		sum2 = fIn[1*nx*ny+(nx-1)*ny+j] + fIn[5*nx*ny+(nx-1)*ny+j] + fIn[8*nx*ny+(nx-1)*ny+j];
		ux[(nx-1)*ny +j] = -1.0 +sum1 + 2.0*sum2;//UNCERTAIN!
		uy[(nx-1)*ny +j] = 0;
		fIn[3*nx*ny+(nx-1)*ny + j] = fIn[1*nx*ny+(nx-1)*ny + j] - (2.0/3)*ux[(nx-1)*ny +j];
		fIn[7*nx*ny+(nx-1)*ny + j] = fIn[5*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[2*nx*ny+(nx-1)*ny + j]-fIn[4*nx*ny+(nx-1)*ny + j])
					   -(1.0/6)*ux[(nx-1)*ny +j];
		fIn[6*nx*ny+(nx-1)*ny + j] = fIn[8*nx*ny+(nx-1)*ny + j] + 0.5*(fIn[4*nx*ny+(nx-1)*ny + j]-fIn[2*nx*ny+(nx-1)*ny + j])
					  -(1.0/6)*ux[(nx-1)*ny +j];
		rho[(nx-1)*ny +j] = 1.0;
	}
}

void collide(double* ux, double* uy, double*fIn, double* fOut, double* rho, double* ex,
			 double* ey, int nx,int ny, int nf, double tau, double* w){
	double u, fEq, uSq,rhoIJ, uxIJ, uyIJ;
	int i,j,k;
	int nxny = nx*ny;
	
		for (i = 0; i < nx; i++){
			for (j = 0; j < ny; j++){
				uxIJ = ux[ny*i + j];
				uyIJ = uy[ny*i + j];
				rhoIJ = rho[ny*i + j];
				uSq = 1.5*(uxIJ*uxIJ+uyIJ*uyIJ);
				for (k = 0; k < nf; k++){
			u = 3.0*(ex[k]*uxIJ + ey[k]*uyIJ);
			fEq = rhoIJ*w[k]*(1.0+u+0.5*u*u-uSq);
			fOut[nxny*k + ny*i + j] = fIn[nxny*k + ny*i + j]-tau*(fIn[nxny*k + ny*i + j]-fEq);
			}
		}
	}			 
}

void bounce(double* fIn, double* fOut,int nx,int ny,int nf, int* bbCells, int nBBcells, int* opposite){
	int i,j,k,l;
	
	//Obstacle
	for (l = 0; l < nBBcells; l++) {
		for (k = 0; k < nf; k++){
				i = bbCells[l];
				j = bbCells[1*nBBcells+l];
				fOut[nx*ny*k + ny*i + j] = fIn[nx*ny*opposite[k] + ny*i + j];
		}
	}
	
	//Top/Bottow wall
	for (i = 0; i < nx; i++) {
		for (k = 0; k < nf; k++){
			fOut[k*nx*ny + ny*i] = fIn[nx*ny*opposite[k] + ny*i];
			fOut[k*nx*ny +i*ny+ny-1] = fIn[nx*ny*opposite[k] + ny*i + ny-1];
		}
	}
}

void stream(double* fIn, double* fOut,int* ex,int* ey, int nx, int ny, int nf) {
	//Ugly solution. Improve?
	int i,j,k;
	int nxny = nx*ny;
	//Interior
	for (i = 1; i < nx-1; i++) {
		for (j = 1; j < ny-1; j++) {
			for (k = 0; k < nf; k++) {
				fIn[nxny*k+ny*i+j] = fOut[k*nxny+(i-ex[k])*ny +j-ey[k]];
			}
		}
	}	
	
	//Outlet
	for (j = 1; j < ny-1; j++) {
		for (k = 0; k < nf; k++) {
			if (ex[k]==-1) fIn[k*nxny+(nx-1)*ny+j] = fOut[k*nxny+j-ey[k]];
			else fIn[k*nxny+(nx-1)*ny+j] = fOut[k*nxny + (nx-1-ex[k])*ny+j-ey[k]];
		}
	}
	//Inlet
	for (j = 1; j < ny-1; j++) {
		for (k = 0; k < nf; k++) {
			if (ex[k]==1) fIn[k*nxny +j] = fOut[k*nxny+(nx-1)*ny+j-ey[k]];
			else fIn[k*nxny +j] = fOut[k*nxny-ex[k]*ny+j-ey[k]];
		}
	}
	//Top wall
	for (i = 1; i < nx-1; i++) {
		for (k = 0; k < nf; k++) {
			if (ey[k]==-1) fIn[k*nxny +i*ny + ny-1] = fOut[k*nxny +(i-ex[k])*ny];
			else fIn[k*nxny +i*ny + ny-1] = fOut[k*nxny +(i-ex[k])*ny +ny-1-ey[k]];
		}
	}
	//Bottom wall
	for (i = 1; i < nx-1; i++) {
		for (k = 0; k < nf; k++) {
			if (ey[k]==1) fIn[k*nxny+ny*i] = fOut[k*nxny +(i-ex[k])*ny+ny-1];
			else fIn[k*nxny + ny*i] = fOut[k*nxny + (i-ex[k])*ny-ey[k]];
		}
	}
	//Top right corner
	for (k = 0; k < nf; k++) {
			if (ex[k]==-1){
				if (ey[k]==-1) fIn[k*nxny+(nx-1)*ny+ny-1] = fOut[k*nxny];
				else fIn[k*nxny+(nx-1)*ny+ny-1] = fOut[k*nxny +ny-1-ey[k]];
			}
			else if (ey[k]==-1) fIn[k*nxny+(nx-1)*ny+ny-1] = fOut[k*nxny + (nx-1-ex[k])*ny];
			else fIn[k*nxny+(nx-1)*ny+ny-1] = fOut[k*nxny + (nx-1-ex[k])*ny +ny-1-ey[k]];
		}		
	//Bottom right corner
	for (k = 0; k < nf; k++) {
			if (ex[k]==-1){
				if (ey[k]==1) fIn[k*nxny+(nx-1)*ny] = fOut[k*nxny + ny-1];
				else fIn[k*nxny+(nx-1)*ny] = fOut[k*nxny-ey[k]];
				}
			else if (ey[k]==1) fIn[k*nxny+(nx-1)*ny] = fOut[k*nxny + (nx-1-ex[k])*ny + ny-1];
			else fIn[k*nxny+(nx-1)*ny] = fOut[k*nxny + (nx-1-ex[k])*ny -ey[k]];
		}	
	//Top left corner
	for (k = 0; k < nf; k++) {
			if (ex[k]==1){
				if (ey[k]==-1) fIn[k*nxny+ ny-1] = fOut[k*nxny+(nx-1)*ny];
				else fIn[k*nxny+ ny-1] = fOut[k*nxny+(nx-1)*ny+ny-1-ey[k]];
			} 
			else if (ey[k]==-1) fIn[k*nxny+ ny-1] = fOut[k*nxny -ex[k]*ny];
			else fIn[k*nxny+ ny-1] = fOut[k*nxny -ex[k]*ny +ny-1-ey[k]];
		}		
	//Bottom left corner
	for (k = 0; k < nf; k++) {
			if (ex[k]==1){
				if (ey[k]==1) fIn[k*nxny] = fOut[k*nxny+(nx-1)*ny+ny-1];
				else fIn[k*nxny] = fOut[k*nxny+(nx-1)*ny+-ey[k]];
				}
			else if (ey[k]==1) fIn[k*nxny] = fOut[k*nxny-ex[k]*ny+ny-1];
			else fIn[k*nxny] = fOut[k*nxny-ex[k]*ny-ey[k]];
		}
}

void csvWriteD(double* mat, int n, int m, char* path) {
	FILE* fp = fopen(path,"w");
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		fprintf(fp,"%20.20f", mat[m*i + j]);
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}

void csvWriteI(int* mat, int n, int m, char* path) {
	FILE* fp = fopen(path,"w");
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		fprintf(fp,"%d", mat[m*i + j]);
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}

void csvWriteLayer(double* mat, int n, int m, int k, char* path) {
		FILE* fp = fopen(path,"w");
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		fprintf(fp,"%20.20f", mat[k*n*m +m*i + j]);
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}