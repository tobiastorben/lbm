#include "utilities.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

//Reads geometry of bounce back cells from a file, stores it in a vector
//of indices
int* readObstacleData(int nBBcells) {
	FILE* fp = fopen("obst.dat","r");
	char* line = malloc(nBBcells*100);
	fgets(line,nBBcells*100,fp);
	char* token = strtok(line, ",");
	int* ind = malloc(nBBcells*sizeof(int));
	int i = 0;
	 while( token != NULL ) 
   {
    ind[i] = atoi(token);
	token = strtok(NULL, ",");
	i++;
   }
	return ind;
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
}
void printMatD(double** mat, int n, int m) {
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		printf("%8f ", mat[i][j]);
		}
	printf("\n");
	}	
}

double** initUx(int nx,int ny,double uIn) {
	double H = (double) ny-1;
	int i,j;
	double y;
	
	//Initialize ux as Poiseuille flow
	double** ux = (double **)malloc(nx * sizeof(double*));
	for(i = 0; i < nx; i++){
		ux[i] = (double *)malloc(ny * sizeof(double));
	}
	
	for (i = 0; i < nx; i++){
		for (j = 0; j < ny; j++) {
			y = (double) j;
			ux[i][j] = (4*uIn/(H*H))*(y*H-y*y);
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

double*** initFin(int nf,int nx, int ny, int* ex, int* ey,
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
			u = 3*ex[k]*ux[i][j] + ey[k]*uy[i][j];
			fIn[k][i][j] = rho[i][j]*w[k]*(1+u+0.5*u*u-1.5*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j]));
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
			 double** rho, int* ex, int* ey, int nx, int ny, int nf){
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
	double H = (double) ny-1;
	for (int j = 1; j < (ny-1);j++) {
		y = (double) j;
		//Poiseuille and Zou/He BC
		ux[0][j] = (4*uIn/(H*H))*(y*H-y*y);
		uy[0][j] = 0;
		sum1 = fIn[0][0][j] + fIn[2][0][j] + fIn[4][0][j];
		sum2 = fIn[3][0][j] + fIn[6][0][j] + fIn[7][0][j];
		rho[0][j] = (1/(1-ux[0][j]))*(sum1 + 2*sum2);
		fIn[1][0][j] = fIn[3][0][j] + (2.0/3)*rho[0][j]*ux[0][j];
		fIn[5][0][j] = fIn[7][0][j] + 0.5*(fIn[4][0][j]-fIn[2][0][j])
					   +0.5*(rho[0][j]*uy[0][j])//ALWAYS ZERO??
					   +(1/6)*(rho[0][j]*ux[0][j]);
		fIn[8][0][j] = fIn[6][0][j] + 0.5*(fIn[2][0][j]-fIn[4][0][j])
					  -0.5*(rho[0][j]*uy[0][j])//ALWAYS ZERO?
					   +(1/6)*(rho[0][j]*ux[0][j]);
	}
}

void outlet(double** ux, double** uy, double** rho, double*** fIn, int nx, int ny) {
	double sum1, sum2;
	for (int j = 1; j < (ny-1);j++) {
		rho[nx-1][j] = 1;
		sum1 = fIn[0][nx-1][j] + fIn[2][nx-1][j] + fIn[4][nx-1][j];
		sum2 = fIn[1][nx-1][j] + fIn[5][nx-1][j] + fIn[8][nx-1][j];
		ux[nx-1][j] = -1 + sum1 + 2*sum2;//UNCERTAIN!
		uy[nx-1][j] = 0;
		fIn[3][nx-1][j] = fIn[1][nx-1][j] - (2.0/3)*rho[nx-1][j]*ux[nx-1][j];
		fIn[7][nx-1][j] = fIn[5][nx-1][j] + 0.5*(fIn[2][nx-1][j]-fIn[4][nx-1][j])
					   -0.5*(rho[nx-1][j]*uy[nx-1][j])//ALWAYS ZERO?
					   -(1/6)*(rho[nx-1][j]*ux[nx-1][j]);
		fIn[6][nx-1][j] = fIn[8][nx-1][j] + 0.5*(fIn[4][nx-1][j]-fIn[2][nx-1][j])
					  +0.5*(rho[nx-1][j]*uy[nx-1][j])//ALWAYS ZERO?
					  -(1/6)*(rho[nx-1][j]*ux[nx-1][j]);
	}
}

void collide(double** ux, double** uy, double***fIn, double*** fOut, double** rho, int* ex,
			 int* ey, int nx,int ny, int nf, double tau, double* w){
	double u;
	int i,j,k;
	for (k = 0; k < nf; k++){
		for (i = 0; i < nx; i++){
			for (j = 0; j < ny; j++){
			u = 3*ex[k]*ux[i][j] + ey[k]*uy[i][j];
			fOut[k][i][j] = fIn[k][i][j]-tau*rho[i][j]*w[k]*(1+u+0.5*u*u
			-1.5*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j]));
			}
		}
	}			 
}

void bounce(double*** fIn, double*** fOut,int nx,int ny,int nf, int** bbCells, int nBBcells, int* opposite){
	int i,j,k,l;
	for (l = 0; l < nBBcells; l++) {
		for (k = 0; k < nf; k++){
				i = bbCells[0][l];
				j = bbCells[1][l];
				fOut[k][i][j] = fIn[opposite[k]][i][j];
		}
	}
}		 