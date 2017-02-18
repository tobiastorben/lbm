#include "utilities.h"
#include "postProcessing.h"

void printVecD(double* vec, int n) {
	for (int i = 0; i < n; i++){
		printf("%8f, ", vec[i]);
	}
	printf("\n");
		
}
void printVecI(int* vec, int n) {
	for (int i = 0; i < n; i++){
		printf("%8d, ", vec[i]);
	}
	printf("\n");	
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

void csvWriteD(double* mat, int n, int m, char* path) {
	FILE* fp = fopen(path,"w");
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		fprintf(fp,"%.3e", mat[m*i + j]);
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
		fprintf(fp,"%.3e", mat[k*n*m +m*i + j]);
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}

void writeTimeSeries(double* v, int n, char* path) {
	FILE* fp = fopen(path,"a");
	for (int i = 0; i < (n-1); i++){
		fprintf(fp,"%.5e, ", v[i]);
		}
	fprintf(fp,"%.5e\n", v[n-1]);
	fclose(fp);
}

void csvWritePres(double* rho, int n, int m, char* path) {
	FILE* fp = fopen(path,"w");
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		fprintf(fp,"%.3e", rho[m*i + j]-1);
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}

void csvWriteOmega(double* ux, double* uy, int n, int m, char* path) {
	FILE* fp = fopen(path,"w");
	double omega;
	for (int i = 0; i < (n-1); i++){
		for (int j = 0; j < (m-1); j++){
		omega = 0.5*(uy[(i+1)*m +j] - uy[(i-1)*m +j]	-(ux[i*m +j+1] - uy[i*m +j-1]));
		fprintf(fp,"%.3e", omega);
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}
void error(char* reason) {
	printf("ERROR! ");
	printf(reason);
	printf("\n");
	exit(-1);
}

void writeResults(double* ux, double* uy, double* rho,int nx,int ny, int* bbCells, int* bbCellMat, int nBBcells, int iter, char* outDir,  int* outputSelect) {
	char* fName = malloc(100);
	char* path = malloc(100);
	double F[2];
	if (outputSelect[0]) {
		strcpy(path,outDir);
		sprintf(fName,"\\ux%d.csv", iter);
		strcat(path,fName);
		csvWriteD(ux,nx,ny,path);
	}
	
	if (outputSelect[1]) {
		strcpy(path,outDir);
		sprintf(fName,"\\uy%d.csv", iter);
		strcat(path,fName);
		csvWriteD(uy,nx,ny,path);
	}
	
	if (outputSelect[2]) {
		strcpy(path,outDir);
		sprintf(fName,"\\p%d.csv", iter);
		strcat(path,fName);
		csvWritePres(rho,nx,ny,path);
	}
	
	if (outputSelect[3]) {
		strcpy(path,outDir);
		sprintf(fName,"\\omega%d.csv", iter);
		strcat(path,fName);
		csvWriteOmega(ux,uy,nx,ny,path);
	}
	
	if (outputSelect[4]) {
		strcpy(path,outDir);
		strcpy(fName, "\\F.csv");
		strcat(path,fName);
		calcF(F,rho,nx,ny,bbCells,nBBcells,bbCellMat);
		writeTimeSeries(F,2,path);
	}

	free(path);
	return;
}