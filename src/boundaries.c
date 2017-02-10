#include "boundaries.h"

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