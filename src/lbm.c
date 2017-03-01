//Main program
#include "lbm.h"
#include "utilities.h"
#include "initialization.h"
#include "boundaries.h"
#include "core.h"
#include "postProcessing.h"
#include "inputParser.h"

int main(int argc, char* argv[]) {	
	SimParams params;
	FlowData flow;
	LatticeConsts lc;
	time_t t1,t2;
	double initialRho,progression;
	int iter;
	int outputSelectAll[] = {1,1,1,1,1};
	char* inPath = malloc(1000);

	if (argc == 2) inPath = argv[1];
	else strcpy(inPath,"..\\input\\input.in");
	
	initialize(&flow,&params,&lc,inPath);

	//Time execution
	t1 = clock();
	
	//Main loop. Marches flow in time
	printf("Solving flow:  ");
	for (iter = 0; iter < params.nIter; iter++){
		//Update macroscopic variables	
		updateRho(&flow,&lc);		
		updateU(&flow,&lc);
		
		//Apply BCs
		inlet(&flow,&lc,&params);//Poiseuille (Zou/He)		
		outlet(&flow,&lc);//Constant pressure (Zou/He)
		
		//Collide (Bhatnagar-Gross-Kroot model)
		collide(&flow,&lc,&params);//Particle-Particle collisions		
		bounce(&flow,&lc,&params);//Bounce-back collision with boundary
		
		//Streaming
		stream(&flow,&lc);
		
		//Write results and progression
		/*if (iter >= startWrite && iter % cyclesPerWrite == 0) {
			writeResults(ux,uy,rho,nx,ny,bbCells,bbCellMat,nBBcells,iter,outDir,outputSelect);
		}
		progression = 100*iter/nIter;
		printf("%2.0f%%\b\b\b", progression);
		fflush(stdout);*/
	}
	t2 = clock();
	//printf("100%%\n\nElapsed time: %f s\n\n", ((float) t2-(float) t1)/1000.0);
	//writeResults(ux,uy,rho,nx,ny,bbCells,bbCellMat,nBBcells,iter,outDir,outputSelectAll);
	//printf("Results written to: %s\n", outDir);
	return 0;
}



