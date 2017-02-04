//Main program
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utilities.h"

//Indexing: ux(i,j) = ux(x = i*h,y=j*h);

int main(int argc, char** argv) {
	//Define grid
	int nx = 10;
	int ny = 5;
	//Define obstacle //Remove??
	int cylX  = nx/5+1;
	int cylY = ny/2 + 3;
	int cylR = ny/10+1;
	
	//Define flow parameters
	double uIn = 0.1;
	double Re = 100;
	double nu = uIn*cylR/Re;
	double tau = 1/(3*nu+0.5);//Relaxation parameter of BGK collision operator
	
	//Define simulations parameters
	int nIter = 4;
	int cyclesPerWrite = 5;
	
	
	//D2Q9 lattice constants
	int nf = 9;
	double w[] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};//Weights
	int ex[] = {0,1,0,-1,0,1,-1,-1,1};//x component of lattice vectors 
	int ey[] = {0,0,1,0,-1,1,1,1,-1};//y component of lattice vectors
	int opposite[] = {0,3,4,1,2,7,8,5,6};//Index of opposite vector in lattice
	
	//BC data
	int nBBcells = 2;
	//int** bbCells = readObstacleData(nBBcells);//Array of indices to bounce back nodes
	int** bbCells = (int **)malloc(2 * sizeof(int*));
	for(int  i = 0; i < 2; i++){
		bbCells[i] = (int *)malloc(2 * sizeof(int));
	};
	bbCells[0][0] = 1;
	bbCells[0][1] = 2;
	bbCells[1][0] = 3;
	bbCells[1][1] = 4;
	int indInlet = 1;
	int indOutlet = nx;
	
	//ICs: Initialze flow as Poiseuille flow
	double H = ny-2;
	double initialRho = 1;
	double** rho = initRho(nx,ny,initialRho);
	double** ux = initUx(nx,ny,uIn);
	double** uy = initUy(nx,ny);
	double*** fIn = initFin(nf,nx,ny,ex,ey,ux,uy,w,rho);
	double*** fOut = initFout(nx,ny,nf);

	//Main loop
	for (int i = 0; i < nIter; i++){		
		//Update macroscopic variables
		updateRho(rho,fIn,nf,nx,ny);
		updateU(ux,uy,fIn,rho,ex,ey,nx,ny,nf);
		
		//Apply BCs
		inlet(ux,uy,rho,fIn,uIn,nx,ny);//Poiseuille (Zou/He)
		outlet(ux,uy,rho,fIn,nx,ny);//Constant pressure (Zou/He)
		
		//Collide (Bhatnagar-Gross-Kroot model)
		collide(ux,uy,fIn,fOut,rho,ex,ey,nx,ny,nf,tau,w);//Particle-Particle collisions
		bounce(fIn,fOut,nx,ny,nf,bbCells,nBBcells,opposite);//Bounce-back collision with boundary
	}
	
	printMatD(rho,nx,ny);	
	return 0;
}




