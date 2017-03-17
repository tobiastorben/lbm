#include "utilities.h"

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

void writeU(double* u, int n, int m, double c, char* path) {
	FILE* fp;
	fp = fopen(path,"w");
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		fprintf(fp,"%.3e",c*u[m*i + j]);
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

void writePres(double* rho, int n, int m, double c, double rhoPhys, char* path) {
	double cSq;
	FILE* fp;
	
	fp = fopen(path,"w");
	cSq = c*c;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
		fprintf(fp,"%.3e", rhoPhys*cSq*(1.0/3.0)*(rho[m*i + j]-1));
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}

void writeVorticity(double* ux, double* uy, int n, int m, double dt, char* path) {
	FILE* fp;
	double vort,coeff;
	
	fp = fopen(path,"w");
	coeff = 1.0/(4.0*dt);
	
	for (int i = 0; i < (n-1); i++){
		for (int j = 0; j < (m-1); j++){
		vort = coeff*(uy[(i+1)*m +j] - uy[(i-1)*m+j] - (ux[i*m +j+1] - ux[i*m +j-1]));
		fprintf(fp,"%.3e", vort);
		if (j != m-1) fprintf(fp,",");
		}
	if (i != (n-1)) fprintf(fp,"\n");
	}
	fclose(fp);
}

void* launchWriteThread(void* pdata_void) {
	
	double *F,*ux,*uy,*rho,dx,dt,c,rhoPhys;
	int nx,ny,iter,*outputSelect;
	char *fName,*path,*outDir;
	PrintData* pdata;
	SimParams* params;
	
	pdata = (PrintData*) pdata_void;
	fName =  (char*) malloc(100);
	path = (char*) malloc(100);
	nx = pdata->nx;
	ny = pdata->ny;
	ux = pdata->uxCpy;
	uy = pdata->uyCpy;
	rho = pdata->rhoCpy;
	params = pdata->params;
	outputSelect = params->outputSelect;
	outDir = params->outDir;
	iter = pdata->iter;
	dt = params->dtPhys;
	dx = params->dxPhys;	
	rhoPhys = params->rhoPhys;
	c = dx/dt;
	
	if (outputSelect[0]) {
		strcpy(path,outDir);
		sprintf(fName,"/ux%d.csv", iter);
		strcat(path,fName);
		writeU(ux,nx,ny,c,path);
	}
	
	if (outputSelect[1]) {
		strcpy(path,outDir);
		sprintf(fName,"/uy%d.csv", iter);
		strcat(path,fName);
		writeU(uy,nx,ny,c,path);
	}
	
	if (outputSelect[2]) {
		strcpy(path,outDir);
		sprintf(fName,"/p%d.csv", iter);
		strcat(path,fName);
		writePres(rho,nx,ny,c,rhoPhys,path);
	}
	
	if (outputSelect[3]) {
		strcpy(path,outDir);
		sprintf(fName,"/vort%d.csv", iter);
		strcat(path,fName);
		writeVorticity(ux,uy,nx,ny,dt,path);
	}
	
	if (outputSelect[4]) {
		strcpy(path,outDir);
		strcpy(fName, "/F.csv");
		strcat(path,fName);
		F = calcF(ny,params,rho,ux,uy);
		writeTimeSeries(F,2,path);
	}
	free(path);
	return NULL;
}

void printProgression(int iter, int nIter) {
	double progression;
	progression = 100*iter/nIter;
	printf("%2.0f%%\b\b\b", progression);
	fflush(stdout);
}

void writeResults(FlowData* flow, LatticeConsts* lc, int iter, pthread_t* printThread, PrintData* pdata) {
	int nx,ny;
	pthread_join(*printThread,NULL);
	pdata->iter = iter;
	nx = lc->nx;
	ny = lc->ny;
	
	memcpy(pdata->uxCpy,flow->ux,nx*ny*sizeof(double));
	memcpy(pdata->uyCpy,flow->uy,nx*ny*sizeof(double));
	memcpy(pdata->rhoCpy,flow->rho,nx*ny*sizeof(double));
	pthread_create(printThread,NULL,launchWriteThread,pdata);
	
	return;
}