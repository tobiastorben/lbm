//Main program
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utilities.h"

//Indexing: ux(i,j) = ux(x = i*h,y=j*h);

int main(int argc, char** argv) {
	//Define grid
	int nx = 400;
	int ny = 100;
	double tull;
	//Define flow parameters
	double uIn = 0.1;
	double Re = 100;
	double nu = uIn*2*(ny/10 +1)/Re;
	double tau = 1/(3*nu+0.5);//Relaxation parameter of BGK collision operator
	
	//Define simulations parameters
	int nIter = 10000;
	int cyclesPerWrite = 5;
	
	
	//D2Q9 lattice constants
	int nf = 9;
	double w[] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};//Weights
	double ex[] = {0,1.0,0,-1.0,0,1.0,-1.0,-1.0,1.0};//x component of lattice vectors 
	double ey[] = {0,0,1.0,0,-1.0,1.0,1.0,-1.0,-1.0};//y component of lattice vectors
	int cx[] = {0,1,0,-1,0,1,-1,-1,1};//x component of lattice vectors 
	int cy[] = {0,0,1,0,-1,1,1,-1,-1};//y component of lattice vectors
	int opposite[] = {0,3,4,1,2,7,8,5,6};//Index of opposite vector in lattice
	
	//BC data
	int nBBcells = 377;
	int** bbCells = readObstacleData(nBBcells);//Array of indices to bounce back nodes
	//ICs: Initialze flow as Poiseuille flow
	double initialRho = 1;
	double** rho = initRho(nx,ny,initialRho);
	double** ux = initUx(nx,ny,uIn);
	double** uy = initUy(nx,ny);
	double*** fIn = initFin(nf,nx,ny,ex,ey,ux,uy,w,rho);
	double*** fOut = initFout(nx,ny,nf);
	csvWriteD(fOut[7],nx,ny, "fout.csv");
		csvWriteD(ux,nx,ny, "ux.csv");
	csvWriteD(uy,nx,ny, "uy.csv");
	csvWriteD(rho,nx,ny, "rho.csv");
	csvWriteD(fIn[7],nx,ny, "fin.csv");
		printf("Postinit/Preupdate printed\n");
		getchar();
	
	//Main loop. Marches flow in time
	for (int i = 0; i < nIter; i++){	
		//Update macroscopic variables	
		updateRho(rho,fIn,nx,ny,nf);		
		updateU(ux,uy,fIn,rho,ex,ey,nx,ny,nf);
		csvWriteD(fOut[7],nx,ny, "fout.csv");
		csvWriteD(ux,nx,ny, "ux.csv");
	csvWriteD(uy,nx,ny, "uy.csv");
	csvWriteD(rho,nx,ny, "rho.csv");
	csvWriteD(fIn[7],nx,ny, "fin.csv");
		printf("Pre BC's/Post update printed\n");
		getchar();
	
		//Apply BCs
		inlet(ux,uy,rho,fIn,uIn,nx,ny);//Poiseuille (Zou/He)		
		outlet(ux,uy,rho,fIn,nx,ny);//Constant pressure (Zou/He)
		
		//Collide (Bhatnagar-Gross-Kroot model)
		csvWriteD(fOut[7],nx,ny, "fout.csv");
		csvWriteD(ux,nx,ny, "ux.csv");
	csvWriteD(uy,nx,ny, "uy.csv");
	csvWriteD(rho,nx,ny, "rho.csv");
	csvWriteD(fIn[7],nx,ny, "fin.csv");
		printf("Precollide/Post BC printed\n");
		getchar();
		collide(ux,uy,fIn,fOut,rho,ex,ey,nx,ny,nf,tau,w);//Particle-Particle collisions		
		csvWriteD(fOut[7],nx,ny, "fout.csv");
		csvWriteD(ux,nx,ny, "ux.csv");
	csvWriteD(uy,nx,ny, "uy.csv");
	csvWriteD(rho,nx,ny, "rho.csv");
	csvWriteD(fIn[7],nx,ny, "fin.csv");
		printf("Postcollide/Prebounce printed\n");
		getchar();		
		bounce(fIn,fOut,nx,ny,nf,bbCells,nBBcells,opposite);//Bounce-back collision with boundary
		csvWriteD(fOut[7],nx,ny, "fout.csv");
		csvWriteD(ux,nx,ny, "ux.csv");
	csvWriteD(uy,nx,ny, "uy.csv");
	csvWriteD(rho,nx,ny, "rho.csv");
	csvWriteD(fIn[7],nx,ny, "fin.csv");
		printf("Postbounce/Prestream printed\n");
		getchar();
		//Streaming
		stream(fIn,fOut,cx,cy,nx,ny,nf);
csvWriteD(ux,nx,ny, "ux.csv");
	csvWriteD(uy,nx,ny, "uy.csv");
	csvWriteD(rho,nx,ny, "rho.csv");
	csvWriteD(fIn[7],nx,ny, "fin.csv");		
		csvWriteD(fOut[7],nx,ny, "fout.csv");
		printf("Poststream printed\n");
		getchar();
	}
	csvWriteD(ux,nx,ny, "ux.csv");
	csvWriteD(uy,nx,ny, "uy.csv");
	csvWriteD(rho,nx,ny, "rho.csv");
	csvWriteD(fIn[0],nx,ny, "fin.csv");
	csvWriteD(fOut[0],nx,ny, "fout.csv");
	return 0;
}



