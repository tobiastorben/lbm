#include "utilities.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

//Reads geometry of bounce back cells from a file, stores it in a vector
//of indices
int** readObstacleData(int nBBcells) {
	FILE* fp = fopen("obstacle.csv","r");
	char* line = malloc(nBBcells*100);		
	char* token;
	int** bbCells = (int **)malloc(2 * sizeof(int*));
	for(int l = 0; l < 2; l++){
		bbCells[l] = (int *)malloc(nBBcells * sizeof(int));
	}
	int i = 0;
	int j = 0;
	int k = 0;
	while (!feof(fp)){
		fgets(line,nBBcells*100,fp);
		token = strtok(line, ",");		
		int j = 0;
		while( token != NULL ){ 	    
		if (atoi(token)) {
			bbCells[0][k] = i;
			bbCells[1][k] = j;
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

void printMatI(int** mat, int n, int m) {
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		printf("%8d ", mat[i][j]);
		}
	printf("\n");
	}
	printf("\n");
}
void printMatD(double** mat, int n, int m) {
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		printf("%f ", mat[i][j]);
		}
	printf("\n");
	}
printf("\n");	
}

double** initUx(int nx,int ny,double uIn) {
	double H = (double) ny-2;
	int i,j;
	double y;
	
	//Initialize ux as Poiseuille flow
	double** ux = (double **)malloc(nx * sizeof(double*));
	for(i = 0; i < nx; i++){
		ux[i] = (double *)malloc(ny * sizeof(double));
	}
	
	for (i = 0; i < nx; i++){
		for (j = 0; j < ny; j++) {
			y = j-0.5;//Change??
			ux[i][j] = (4.0*uIn/(H*H))*(y*H-y*y);
		}
	}
return ux;	
}

double** initUy(int nx,int ny) {	
	//Zero initialize uy
	double** uy = (double **)calloc(nx,sizeof(double*));
	for(int i = 0; i < nx; i++){
		uy[i] = (double *)calloc(ny,sizeof(double));
	}
	return uy;
}

double*** initFin(int nf,int nx, int ny, double* ex, double* ey,
				double** ux, double** uy, double* w, double** rho) {
	int i,j,k;
	double u;
	//Allocate three dimensional array. Indexed as: [f][nx-1][ny]
	double*** fIn = (double ***)malloc(nf * sizeof(double**));
	for(i = 0; i < nf; i++){
		fIn[i] = (double **)malloc(nx * sizeof(double*));
		for(j = 0; j < nx; j++){
		fIn[i][j] = (double*)malloc(ny * sizeof(double));
		}
	}
	//Initialize particle distribution as equilibrium distrubution for
	//for a Poiseuille flow
	for (k = 0; k < nf; k++){
		for (i = 0; i < nx; i++){
			for (j = 0; j < ny; j++){
			u = 3.0*(ex[k]*ux[i][j] + ey[k]*uy[i][j]);
			fIn[k][i][j] = rho[i][j]*w[k]*(1.0+u+0.5*u*u-1.5*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j]));
			}
		}
	}
	return fIn;
}

double*** initFout(int nx,int ny, int nf) {
	int i,j;
	double*** fOut = (double ***)malloc(nf * sizeof(double**));
	for(i = 0; i < nf; i++){
		fOut[i] = (double **)malloc(nx * sizeof(double*));
		for(j = 0; j < nx; j++){
		fOut[i][j] = (double*)malloc(ny * sizeof(double));
		}
	}
	return fOut;
}

double** initRho(int nx, int ny, double initialRho) {
	int i,j;
	double** rho = (double**) malloc(nx*sizeof(double*));
	for(int i = 0; i < nx; i++){
		rho[i] = (double *)malloc(ny*sizeof(double));
	}
	
	for (i = 0; i < nx; i++){
			for (j = 0; j < ny; j++){
			rho[i][j] = initialRho;
			}
		}
	return rho;
}

void updateRho(double** rho, double*** fIn, int nx, int ny, int nf) {
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++){
			rho[i][j] = 0;
			for (int k = 0; k < nf; k++) {
				rho[i][j] += fIn[k][i][j];
			}
		}
	}
}

void updateU(double** ux, double** uy, double*** fIn,
			 double** rho, double* ex, double* ey, int nx, int ny, int nf){
	for (int i = 0; i < nx; i++) {//ARE THE LIMITS CORRECT? (1:ny-1)??
		for (int j = 0; j < ny; j++){
			ux[i][j] = 0;
			uy[i][j] = 0;
			for (int k = 0; k < nf; k++) {
				ux[i][j] += fIn[k][i][j]*ex[k]/rho[i][j];
				uy[i][j] += fIn[k][i][j]*ey[k]/rho[i][j];
			}
		}
	} 				 
}

void inlet(double** ux, double** uy, double** rho, double*** fIn, double uIn, int nx, int ny) {
	double y, sum1, sum2;
	double H = (double) ny-2;
	for (int j = 1; j < (ny-1);j++) {
		y = j-0.5;//Change??
		//Poiseuille and Zou/He BC
		ux[0][j] = (4.0*uIn/(H*H))*(y*H-y*y);
		uy[0][j] = 0;
		sum1 = fIn[0][0][j] + fIn[2][0][j] + fIn[4][0][j];
		sum2 = fIn[3][0][j] + fIn[6][0][j] + fIn[7][0][j];
		rho[0][j] = (1.0/(1-ux[0][j]))*(sum1 + 2*sum2);
		fIn[1][0][j] = fIn[3][0][j] + (2.0/3)*rho[0][j]*ux[0][j];
		fIn[5][0][j] = fIn[7][0][j] + 0.5*(fIn[4][0][j]-fIn[2][0][j])
					   +0.5*(rho[0][j]*uy[0][j])//ALWAYS ZERO??
					   +(1.0/6)*(rho[0][j]*ux[0][j]);
		fIn[8][0][j] = fIn[6][0][j] + 0.5*(fIn[2][0][j]-fIn[4][0][j])
					  -0.5*(rho[0][j]*uy[0][j])//ALWAYS ZERO?
					   +(1.0/6)*(rho[0][j]*ux[0][j]);
	}
}

void outlet(double** ux, double** uy, double** rho, double*** fIn, int nx, int ny) {
	double sum1, sum2;
	for (int j = 1; j < (ny-1);j++) {
		rho[nx-1][j] = 1.0;
		sum1 = fIn[0][nx-1][j] + fIn[2][nx-1][j] + fIn[4][nx-1][j];
		sum2 = fIn[1][nx-1][j] + fIn[5][nx-1][j] + fIn[8][nx-1][j];
		ux[nx-1][j] = -1.0 + (1.0/(rho[nx-1][j]))*(sum1 + 2.0*sum2);//UNCERTAIN!
		uy[nx-1][j] = 0;
		fIn[3][nx-1][j] = fIn[1][nx-1][j] - (2.0/3)*rho[nx-1][j]*ux[nx-1][j];
		fIn[7][nx-1][j] = fIn[5][nx-1][j] + 0.5*(fIn[2][nx-1][j]-fIn[4][nx-1][j])
					   -0.5*(rho[nx-1][j]*uy[nx-1][j])//ALWAYS ZERO?
					   -(1.0/6)*(rho[nx-1][j]*ux[nx-1][j]);
		fIn[6][nx-1][j] = fIn[8][nx-1][j] + 0.5*(fIn[4][nx-1][j]-fIn[2][nx-1][j])
					  +0.5*(rho[nx-1][j]*uy[nx-1][j])//ALWAYS ZERO?
					  -(1.0/6)*(rho[nx-1][j]*ux[nx-1][j]);
	}
}

void collide(double** ux, double** uy, double***fIn, double*** fOut, double** rho, double* ex,
			 double* ey, int nx,int ny, int nf, double tau, double* w){
	double u, fEq;
	int i,j,k;
	for (k = 0; k < nf; k++){
		for (i = 0; i < nx; i++){
			for (j = 0; j < ny; j++){
			u = 3.0*(ex[k]*ux[i][j] + ey[k]*uy[i][j]);
			fEq = rho[i][j]*w[k]*(1+u+0.5*u*u
			-1.5*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j]));
			fOut[k][i][j] = fIn[k][i][j]-tau*(fIn[k][i][j]-fEq);
			}
		}
	}			 
}

void bounce(double*** fIn, double*** fOut,int nx,int ny,int nf, int** bbCells, int nBBcells, int* opposite){
	int i,j,k,l;
	
	//Obstacle
	for (l = 0; l < nBBcells; l++) {
		for (k = 0; k < nf; k++){
				i = bbCells[0][l];
				j = bbCells[1][l];
				fOut[k][i][j] = fIn[opposite[k]][i][j];
		}
	}
	
	//Top/Bottow wall
	for (i = 0; i < nx; i++) {
		for (k = 0; k < nf; k++){
			fOut[k][i][0] = fIn[opposite[k]][i][0];
			fOut[k][i][ny-1] = fIn[opposite[k]][i][ny-1];
		}
	}
}

void stream(double*** fIn, double*** fOut,int* ex,int* ey, int nx, int ny, int nf) {
	//Uglu solution. Improve?
	int i,j,k;

	//Interior
	for (i = 1; i < nx-1; i++) {
		for (j = 1; j < ny-1; j++) {
			for (k = 0; k < nf; k++) {
				fIn[k][i][j] = fOut[k][i-ex[k]][j-ey[k]];
			}
		}
	}	
	
	//Outlet
	for (j = 1; j < ny-1; j++) {
		for (k = 0; k < nf; k++) {
			if (ex[k]==-1) fIn[k][nx-1][j] = fOut[k][0][j-ey[k]];
			else fIn[k][nx-1][j] = fOut[k][nx-1-ex[k]][j-ey[k]];
		}
	}
	//Inlet
	for (j = 1; j < ny-1; j++) {
		for (k = 0; k < nf; k++) {
			if (ex[k]==1) fIn[k][0][j] = fOut[k][nx-1][j-ey[k]];
			else fIn[k][0][j] = fOut[k][-ex[k]][j-ey[k]];
		}
	}
	//Top wall
	for (i = 1; i < nx-1; i++) {
		for (k = 0; k < nf; k++) {
			if (ey[k]==-1) fIn[k][i][ny-1] = fOut[k][i-ex[k]][0];
			else fIn[k][i][ny-1] = fOut[k][i-ex[k]][ny-1-ey[k]];
		}
	}
	//Bottom wall
	for (i = 1; i < nx-1; i++) {
		for (k = 0; k < nf; k++) {
			if (ey[k]==1) fIn[k][i][0] = fOut[k][i-ex[k]][ny-1];
			else fIn[k][i][0] = fOut[k][i-ex[k]][-ey[k]];
		}
	}
	//Top right corner
	for (k = 0; k < nf; k++) {
			if (ex[k]==-1){
				if (ey[k]==-1) fIn[k][nx-1][ny-1] = fOut[k][0][0];
				else fIn[k][nx-1][ny-1] = fOut[k][0][ny-1-ey[k]];
			}
			else if (ey[k]==-1) fIn[k][nx-1][ny-1] = fOut[k][nx-1-ex[k]][0];
			else fIn[k][nx-1][ny-1] = fOut[k][nx-1-ex[k]][ny-1-ey[k]];
		}		
	//Bottom right corner
	for (k = 0; k < nf; k++) {
			if (ex[k]==-1){
				if (ey[k]==1) fIn[k][nx-1][0] = fOut[k][0][ny-1];
				else fIn[k][nx-1][0] = fOut[k][0][-ey[k]];
				}
			else if (ey[k]==1) fIn[k][nx-1][0] = fOut[k][nx-1-ex[k]][ny-1];
			else fIn[k][nx-1][0] = fOut[k][nx-1-ex[k]][-ey[k]];
		}	
	//Top left corner
	for (k = 0; k < nf; k++) {
			if (ex[k]==1){
				if (ey[k]==-1) fIn[k][0][ny-1] = fOut[k][nx-1][0];
				else fIn[k][0][ny-1] = fOut[k][nx-1][ny-1-ey[k]];
			} 
			else if (ey[k]==-1) fIn[k][0][ny-1] = fOut[k][-ex[k]][0];
			else fIn[k][0][ny-1] = fOut[k][-ex[k]][ny-1-ey[k]];
		}		
	//Bottom left corner
	for (k = 0; k < nf; k++) {
			if (ex[k]==1){
				if (ey[k]==1) fIn[k][0][0] = fOut[k][nx-1][ny-1];
				else fIn[k][0][0] = fOut[k][nx-1][-ey[k]];
				}
			else if (ey[k]==1) fIn[k][0][0] = fOut[k][-ex[k]][ny-1];
			else fIn[k][0][0] = fOut[k][-ex[k]][-ey[k]];
		}
}

void csvWriteD(double** mat, int n, int m, char* path) {
	FILE* fp = fopen(path,"w");
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		fprintf(fp,"%f", mat[i][j]);
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}

void csvWriteI(int** mat, int n, int m, char* path) {
	FILE* fp = fopen(path,"w");
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		fprintf(fp,"%d", mat[i][j]);
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}