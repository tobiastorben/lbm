//Main program
#include "utilities.h"
#include "initialization.h"
#include "boundaries.h"
#include "core.h"
#include "postProcessing.h"
#include "inputParser.h"

int main(int argc, char* argv[]) {
	
	time_t t1,t2;
	double Re,u0,u0Phys,dxPhys,dtPhys,nuPhys,tau,initialRho,progression,*rho,*ux,*uy,*fIn,*fOut;
	int nIter,cyclesPerWrite,startWrite,nx,ny,nf,nBBcells,iter,*bbCellMat,*bbCells;
		
	//D2Q9 lattice constants
	double w[] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};//Weights
	double ex[] = {0,1.0,0,-1.0,0,1.0,-1.0,-1.0,1.0};//x component of lattice vectors 
	double ey[] = {0,0,1.0,0,-1.0,1.0,1.0,-1.0,-1.0};//y component of lattice vectors
	int exI[] = {0,1,0,-1,0,1,-1,-1,1};//int version for indexing
	int eyI[] = {0,0,1,0,-1,1,1,-1,-1};//int version for indexing
	int opposite[] = {0,3,4,1,2,7,8,5,6};//Index of opposite vector in lattice
	
	char* inPath = malloc(1000);char* obstaclePath = malloc(1000); char* outDir = malloc(1000);
	int outputSelect[] = {0,0,0,0,0};int outputSelectAll[] = {1,1,1,1,1};
		
	printf("Initializing...\n\n");

	if (argc == 2) inPath = argv[1];
	else strcpy(inPath,"..\\input\\input.in");
		
	if (parseInput(inPath,obstaclePath,&dxPhys,&dtPhys,&nuPhys,&u0Phys,&nIter,&cyclesPerWrite,&startWrite,outputSelect, outDir)) {
		error("Invalid simulation parameters. Exiting.");
	}

	//Read geometry and set dimensions
	bbCellMat = readObstacle(&nBBcells,&nx,&ny,obstaclePath);//Read obstacle data from file
	bbCells = mapObstacleCells(bbCellMat,nBBcells,nx,ny);//Array of indices to bounce back nodes
	
	//Cast to non-dimensional form
	nonDimensionalize(nx,ny,dtPhys,dxPhys,nuPhys,u0Phys,&u0,&tau);

	//ICs: Initialze flow as Poiseuille flow
	nf = 9;
	initialRho = 1;
	rho = initRho(nx,ny,initialRho);
	ux = initUx(nx,ny,u0);
	uy = initUy(nx,ny);
	fIn = initFin(nf,nx,ny,ex,ey,ux,uy,w,rho);
	fOut = initFout(nx,ny,nf);

	//Time execution
	t1 = clock();
	
	//Main loop. Marches flow in time
	printf("Solving flow:  ");
	for (iter = 0; iter < nIter; iter++){	
		//Update macroscopic variables	
		updateRho(rho,fIn,nx,ny,nf);		
		updateU(ux,uy,fIn,rho,ex,ey,nx,ny,nf);
		
		//Apply BCs
		inlet(ux,uy,rho,fIn,u0,nx,ny);//Poiseuille (Zou/He)		
		outlet(ux,uy,rho,fIn,nx,ny);//Constant pressure (Zou/He)
		
		//Collide (Bhatnagar-Gross-Kroot model)
		collide(ux,uy,fIn,fOut,rho,ex,ey,nx,ny,nf,tau,w);//Particle-Particle collisions		
		bounce(fIn,fOut,nx,ny,nf,bbCells,nBBcells,opposite);//Bounce-back collision with boundary
		
		//Streaming
		stream(fIn,fOut,exI,eyI,nx,ny,nf);
		
		//Write results and progression
		if (iter >= startWrite && iter % cyclesPerWrite == 0) {
			writeResults(ux,uy,rho,nx,ny,bbCells,bbCellMat,nBBcells,iter,outDir,outputSelect);
		}
		progression = 100*iter/nIter;
		printf("%2.0f%%\b\b\b", progression);
		fflush(stdout);

	}
	t2 = clock();
	printf("100%%\n\nElapsed time: %f s\n\n", ((float) t2-(float) t1)/1000.0);
	writeResults(ux,uy,rho,nx,ny,bbCells,bbCellMat,nBBcells,iter,outDir,outputSelectAll);
	printf("Results written to: %s\n", outDir);
	return 0;
}



