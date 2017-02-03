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

/*double*** initF(int nf,int nx, int ny, double uIn) {
	
	//Allocate three dimensional array. Indexed as: [f][nx][ny]
	double*** fIn = (double ***)malloc(nf * sizeof(double**));
	for(int i = 0; i < nf; i++){
		fIn[i] = (double **)malloc(nx * sizeof(double*));
		for(int j = 0; j < nx; j++){
		fIn[i][j] = (double*)malloc(ny * sizeof(double));
		}
	}
	
	
	for (int f = 0; f < nf; f++){
		for (int i = 0; i < nx; i++){
			for (int j = 0; j < ny; j++){
			
			}
		}
	}
	printf("The number is: %f\n", fIn[1][2][3]);
	return fIn;
}*/